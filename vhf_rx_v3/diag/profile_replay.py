"""Capture 6 s of live RSPdx IQ, then replay through the worker DSP path
with per-component timing. Diagnostics only.

Usage:  python diag/profile_replay.py   (needs the RSPdx streaming)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nbfm import BandSpectrum, ChannelDemod, Squelch
from sdrplay_api import AGC_50HZ, BW_5M, SDRPlay

FS = 6.048e6
CHUNK = 1 << 16
DUR = 6.0


def main():
    rx = SDRPlay()
    dev = rx.find_device()
    prof = rx.profile
    rx.get_params()
    rx.configure(fs_hz=FS, rf_hz=145.0e6, bw_khz=BW_5M,
                 lna=prof.lna_max_2m, gr_db=50,
                 agc=AGC_50HZ, agc_setpoint_dbfs=-20)
    cap = []
    lock = __import__("threading").Lock()

    def on_iq(xi, xq, n, reset):
        i = np.ctypeslib.as_array(xi, shape=(n,)).copy()
        q = np.ctypeslib.as_array(xq, shape=(n,)).copy()
        with lock:
            cap.append((i, q))

    def on_event(eid, tuner, gain):
        pass

    rx.start(on_iq, on_event)
    t0 = time.monotonic()
    while time.monotonic() - t0 < DUR:
        time.sleep(0.2)
    rx.stop()
    with lock:
        cap = list(cap)
    n_blocks = len(cap)
    n_samples = sum(i.size for i, q in cap)
    print(f"captured {n_blocks} blocks, {n_samples} samples "
          f"({n_samples / DUR / 1e6:.3f} MSPS delivered, avg block {n_samples / n_blocks:.0f})")

    # replay with timing
    spec = BandSpectrum(FS, 145.0e6, 144.0e6, 146.0e6)
    demod = ChannelDemod(145.5e6, FS)
    sq = Squelch()

    t_conv = t_concat = t_spec = t_demod = 0.0
    pending = []
    pending_n = 0
    t0 = time.monotonic()
    for i16, q16 in cap:
        tc = time.perf_counter()
        iq = ((i16.astype(np.float32) + 1j * q16.astype(np.float32))
              / prof.sample_scale).astype(np.complex64)
        t_conv += time.perf_counter() - tc
        pending.append(iq)
        pending_n += iq.size
        while pending_n >= CHUNK:
            tc = time.perf_counter()
            allc = np.concatenate(pending)
            t_concat += time.perf_counter() - tc
            tc = time.perf_counter()
            cp, chans, nf = spec.update(allc[:CHUNK])
            t_spec += time.perf_counter() - tc
            tc = time.perf_counter()
            audio, m = demod.process(allc[:CHUNK])
            sq.update(m["carrier_db"], audio.size * 2)
            t_demod += time.perf_counter() - tc
            rest = allc[CHUNK:]
            pending = [rest] if rest.size else []
            pending_n = rest.size
    t_replay = time.perf_counter() - t0
    stream_s = n_samples / FS
    print(f"replayed in {t_replay:.2f} s of CPU for {stream_s:.2f} s of stream "
          f"-> {t_replay / stream_s * 100:.0f}% of real time")
    print(f"  int16->c64 conversion : {t_conv / stream_s * 1000:7.1f} ms/s")
    print(f"  concatenate          : {t_concat / stream_s * 1000:7.1f} ms/s")
    print(f"  spectrum (FFT+grid)  : {t_spec / stream_s * 1000:7.1f} ms/s")
    print(f"  demod (1 channel)    : {t_demod / stream_s * 1000:7.1f} ms/s")


if __name__ == "__main__":
    main()
