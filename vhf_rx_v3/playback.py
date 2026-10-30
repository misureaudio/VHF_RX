"""Audio output stage for the NBFM receiver prototype (design doc §4.5).

Mixes per-channel demodulated audio (8 kHz float32, deviation in Hz units)
to the local computer's sound device (48 kHz mono float32) via
sounddevice/PortAudio (WASAPI shared mode on Windows).

Pipeline:
  worker thread (per chunk, ~11 ms):
    push(f, audio, squelch_open) per channel
      -> per-channel 5 ms linear gain ramp (click-free squelch gating)
      -> fixed gain 1/5000 x vol (3 kHz deviation -> 0.6 FS; no
         auto-normalization — it pumps with speech level changes)
      -> sum of open channels, clipped at +/-0.95
    flush()
      -> resample 8k -> 48k (resample_poly, 6x integer ratio)
      -> append to a bounded ring (2 s, drop oldest)
  PortAudio callback (10 ms): drains the ring, silence on starvation
    (starved frames are counted as underruns for the run report).

The mix/resample/ring core is pure numpy: it is unit-tested without a
sound device (tests/test_nbfm.py::test_audio_player) by calling the
callback directly with synthetic outdata.
"""
from __future__ import annotations

import collections
import threading

import numpy as np
from scipy.signal import resample_poly

FS_IN = 8000          # channel audio rate (Hz)
FS_OUT = 48000        # output rate: WASAPI-native, 6x integer ratio to FS_IN
GAIN_DEFAULT = 1.0 / 5000.0   # 3 kHz deviation -> 0.6 FS
CLIP = 0.95           # post-mix clip level
RAMP_S = 0.005        # 5 ms squelch open/close ramp (click-free gating)


class AudioPlayer:
    """Mixes squelch-gated channel audio into one local output stream.

    Worker side: push(f, audio, open_) per channel, flush() per chunk,
    drop(f) when a channel leaves the demod set. Output side: start() /
    stop() manage the sounddevice stream.
    """

    def __init__(self, fs_in: int = FS_IN, fs_out: int = FS_OUT,
                 vol: float = 1.0, ring_s: float = 2.0):
        self.fs_in = fs_in
        self.fs_out = fs_out
        self.gain = GAIN_DEFAULT * vol
        self._ramp_n = max(1, int(RAMP_S * fs_in))
        self._ring_max = int(ring_s * fs_out)
        self.lock = threading.Lock()
        self.ring: "collections.deque[np.ndarray]" = collections.deque()
        self.ring_n = 0
        self.underrun_frames = 0        # callback frames filled with silence
        self._open_chans: set[float] = set()
        self._ch_gain: dict[float, float] = {}
        self._acc = np.zeros(0, dtype=np.float64)
        self._up = fs_out // fs_in if fs_out % fs_in == 0 else None
        self.stream = None

    # -- worker side ------------------------------------------------------
    def push(self, f: float, audio: np.ndarray, open_: bool) -> None:
        """One channel's new 8 kHz block (deviation in Hz). Squelch-gated
        with a 5 ms linear ramp, scaled, accumulated for the mix."""
        if audio.size == 0:
            return
        if open_:
            self._open_chans.add(f)
        else:
            self._open_chans.discard(f)
        g0 = self._ch_gain.get(f, 0.0)
        target = 1.0 if open_ else 0.0
        if g0 < target:
            g1 = min(target, g0 + audio.size / self._ramp_n)
        elif g0 > target:
            g1 = max(target, g0 - audio.size / self._ramp_n)
        else:
            g1 = g0
        self._ch_gain[f] = g1
        if g0 <= 0.0 and g1 <= 0.0:
            return                      # closed and fully faded: nothing to mix
        x = audio.astype(np.float64)
        if g0 != g1:
            x = x * np.linspace(g0, g1, x.size)
        else:
            x = x * g0
        x *= self.gain
        if self._acc.size == 0:
            self._acc = x
        else:
            m = min(self._acc.size, x.size)
            self._acc = self._acc[:m] + x[:m]

    def flush(self) -> None:
        """End of a chunk: clip the accumulated mix, resample, ring-append."""
        if self._acc.size == 0:
            return
        mix = np.clip(self._acc, -CLIP, CLIP).astype(np.float32)
        self._acc = np.zeros(0, dtype=np.float64)
        if self._up is not None and self._up != 1:
            # padtype='symmetric' (old scipy 'odd') extends the block edge
            # instead of zero-padding: zero pads would dip the tone ~10%
            # at every 11 ms block edge
            out = resample_poly(mix, self._up, 1, padtype="symmetric").astype(np.float32)
        else:
            out = mix                   # integer-ratio unsupported (test hook)
        with self.lock:
            self.ring.append(out)
            self.ring_n += out.size
            while self.ring_n > self._ring_max:
                drop = self.ring.popleft()
                self.ring_n -= drop.size

    def drop(self, f: float) -> None:
        """Channel removed from the demod set: forget its state."""
        self._ch_gain.pop(f, None)
        self._open_chans.discard(f)

    def _drain(self, frames: int) -> np.ndarray:
        """Remove up to `frames` samples from the ring (caller holds lock)."""
        if self.ring_n == 0:
            return np.zeros(0, dtype=np.float32)
        take = min(frames, self.ring_n)
        out = np.empty(take, dtype=np.float32)
        pos = 0
        while pos < take:
            seg = self.ring.popleft()
            k = min(seg.size, take - pos)
            out[pos:pos + k] = seg[:k]
            pos += k
            if k < seg.size:
                self.ring.appendleft(seg[k:])
                break
        self.ring_n -= take
        return out

    # -- sounddevice side ---------------------------------------------------
    def _cb(self, outdata, frames, time_info, status) -> None:
        with self.lock:
            x = self._drain(frames)
        n = x.size
        if n < frames and self._open_chans:
            # starved while a squelch-open channel should be producing:
            # a real underrun. Silence with no open channel is by design
            # (nothing to play) and is NOT counted.
            self.underrun_frames += frames - n
        outdata[:n, 0] = x
        outdata[n:, 0] = 0.0

    def start(self, device=None) -> None:
        """Open the output stream on the local sound device. Raises on
        failure (no device / unsupported rate) — the caller degrades to
        record-only."""
        import sounddevice as sd
        self.stream = sd.OutputStream(samplerate=self.fs_out, channels=1,
                                      dtype="float32",
                                      blocksize=self.fs_out // 100,
                                      device=device, callback=self._cb)
        self.stream.start()

    def stop(self) -> None:
        if self.stream is not None:
            self.stream.stop()
            self.stream.close()
            self.stream = None
