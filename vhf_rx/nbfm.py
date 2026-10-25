"""NBFM DSP chain for the RSPdx VHF/UHF receiver prototype.

Chain per design doc (docs/RSPdx_VHF_NBFM_Receiver_Design.md):

  band-limit FIR (full-band)
    -> ActivityMonitor: 8192-pt FFT -> 12.5 kHz channel power grid
    -> per channel: NCO DDC -> /7 (32-tap) -> /9 (64-tap) -> /6 (128-tap BPF)
       -> 16 kHz complex
       -> quadrature discriminator (normalized) -> audio in Hz
       -> 20 Hz DC-block -> /2 (2nd-order Butterworth, 3.4 kHz) -> 8 kHz
       -> de-emphasis (bilinear H(s)=1/(1+s*75us))
    -> Squelch: narrowband carrier energy vs noise floor, hysteresis + timing
    -> CTCSS: Goertzel bank, 38 tones (ETSI TS 103 236)
"""
from __future__ import annotations

import numpy as np
from scipy.signal import firwin, lfilter, butter

# Numba-accelerated NCO (geometric recurrence, ~1 complex mult/sample).
# A numpy exp() per sample costs ~250 ms/s of stream; the recurrence costs
# ~1 ms/s. Phase drift over a 64 k-sample chunk is <0.007 rad (harmless:
# it shows up as a slow frequency offset removed by the DC blocker).
try:
    import numba as _numba

    @_numba.njit
    def _nco_fill(seq, w, ph0):
        s = np.exp(-1j * ph0)
        for k in range(seq.size):
            seq[k] = s
            s *= w

    # 2nd-order PLL FM demod (design doc §4.3). Per sample: PD (wrapped
    # phase error vs the NCO phase), loop filter (proportional + integral),
    # NCO phase update. The NCO frequency phi = Kp*e + integ IS the audio
    # (the loop's frequency estimate, rad/sample * scale = Hz); using the
    # integrator alone low-passes the audio and loses ~30% of the tone.
    # The arctan2 PD is amplitude-invariant, so AGC steps don't modulate
    # the loop gain.
    @_numba.njit
    def _pll_block(x, theta, integ, Kp, Ki, scale):
        n = x.size
        out = np.empty(n)
        two_pi = 2.0 * np.pi
        for k in range(n):
            e = np.arctan2(x[k].imag, x[k].real) - theta
            if e > np.pi:
                e -= two_pi
            elif e < -np.pi:
                e += two_pi
            integ += Ki * e
            phi = Kp * e + integ
            theta += phi
            if theta > np.pi:
                theta -= two_pi
            elif theta < -np.pi:
                theta += two_pi
            out[k] = phi * scale
        return out, theta, integ

    _HAVE_NUMBA = True
except ImportError:      # pragma: no cover - numba is in the project venv
    _HAVE_NUMBA = False

    def _pll_block(x, theta, integ, Kp, Ki, scale):
        out = np.empty(x.size)
        two_pi = 2.0 * np.pi
        re, im = x.real, x.imag
        for k in range(x.size):
            e = np.arctan2(im[k], re[k]) - theta
            if e > np.pi:
                e -= two_pi
            elif e < -np.pi:
                e += two_pi
            integ += Ki * e
            phi = Kp * e + integ
            theta += phi
            if theta > np.pi:
                theta -= two_pi
            elif theta < -np.pi:
                theta += two_pi
            out[k] = phi * scale
        return out, theta, integ


def pll_tracking_response(f: float, f_n: float, zeta: float,
                          fs: float = 16000.0) -> float:
    """|T_f(e^{jw})| of the 2nd-order PLL: input frequency -> NCO (output)
    frequency, EXACT discrete form. T_f(u) = [Kp(1-u)+Ki] /
    [(1-u)^2 + u(Kp(1-u)+Ki)] with u = e^{-jw}, w_n = 2pi*f_n/fs,
    Kp = 2*zeta*w_n, Ki = w_n^2. (The continuous-time (Kp s + Ki)/(s^2 +
    Kp s + Ki) matches at low f but is wrong near/above the loop BW —
    the z^-1 NCO delay adds a resonant peak.)"""
    w_n = 2.0 * np.pi * f_n / fs
    Kp = 2.0 * zeta * w_n
    Ki = w_n * w_n
    u = np.exp(-2j * np.pi * f / fs)
    num = Kp * (1.0 - u) + Ki
    den = (1.0 - u) ** 2 + u * (Kp * (1.0 - u) + Ki)
    return float(abs(num / den))

# ETSI TS 103 236 V1.1.1 — 38 standard CTCSS tones (Hz)
CTCSS_TONES = np.array([
    67.0, 69.3, 71.9, 74.4, 77.0, 79.7, 82.5, 85.4, 88.5, 91.5,
    94.8, 97.4, 100.0, 103.5, 107.2, 110.9, 114.8, 118.8, 123.0, 127.3,
    131.8, 136.5, 141.3, 146.2, 151.4, 156.7, 162.2, 167.9, 173.8, 179.9,
    186.2, 192.8, 203.5, 210.7, 218.1, 225.7, 233.6, 241.8, 250.3,
])


def _decimate(x: np.ndarray, h: np.ndarray, m: int) -> np.ndarray:
    return lfilter(h, 1.0, x)[::m]


# ---------------------------------------------------------------------------
# Measurement instrumentation (EN 300 086 test plan)
# ---------------------------------------------------------------------------
# ITU-T O.41 telephone psophometric weighting — APPROXIMATION.
# Zero-phase FIR fitted to the classic P.39/O.41 response anchors, normalized
# to 0 dB at 800 Hz (the O.41 calibration point). The ITU PDF sits behind a
# login wall, so the anchor values are from the widely reproduced classic
# curve; expect ~1-2 dB uncertainty on the weighting itself. For
# certification-grade numbers, verify with a calibrated psophometer.
_PSOPH_ANCHORS = [
    (100, -15.0), (200, -8.0), (300, -4.5), (400, -2.0), (500, -0.5),
    (600, 0.0), (800, 0.8), (1000, 1.2), (1200, 1.0), (1600, 0.2),
    (2000, -1.0), (2500, -2.5), (3000, -4.5), (4000, -8.5), (5000, -13.0),
    (6000, -18.0), (7000, -24.0),
]


def _psoph_fir(fs: float = 8000.0, ntaps: int = 257) -> np.ndarray:
    f_a = np.array([a for a, _ in _PSOPH_ANCHORS])
    d_a = np.array([b for _, b in _PSOPH_ANCHORS])
    half = ntaps // 2
    # DFT bin grid for the positive half: f_k = k*fs/ntaps, k = 0..ntaps/2
    f_grid = np.linspace(0.0, fs / 2, half + 1)
    # piecewise-linear in dB vs log10(f); flat at the first anchor below 100 Hz
    logf = np.log10(np.maximum(f_grid, 1e-9))
    db = np.interp(logf, np.log10(f_a), d_a)
    db -= db[int(np.argmin(np.abs(f_grid - 800.0)))]   # 0 dB at 800 Hz
    mag = 10.0 ** (db / 20.0)
    # odd-length linear-phase (type 1) FIR: H[k] = mag[k] for k <= N/2 and
    # H[N-k] = mag[k] (conjugate symmetry); H[129..256] = mag[128..1]
    H = np.concatenate([mag, mag[half:0:-1]])
    # ifft of a zero-phase response is a *wrapped* kernel (center at index 0);
    # ifftshift gives the linear-phase FIR (center tap at (N-1)/2, constant
    # group delay) that lfilter expects
    return np.fft.ifftshift(np.real(np.fft.ifft(H)))


_PSOPH_CACHE: dict[float, np.ndarray] = {}


def psoph_weight(x: np.ndarray, fs: float = 8000.0) -> np.ndarray:
    """Apply O.41 psophometric weighting.

    Single lfilter pass of the linear-phase FIR: magnitude = |H| exactly
    (a double pass would square it); the constant (N-1)/2-sample delay is
    irrelevant for the stationary SND/ND measurements.
    """
    h = _PSOPH_CACHE.get(fs)
    if h is None:
        h = _psoph_fir(fs)
        _PSOPH_CACHE[fs] = h
    return lfilter(h, 1.0, x.astype(np.float64))


def snd_nd(audio: np.ndarray, fs: float = 8000.0, f0: float = 1000.0) -> float:
    """Psophometrically weighted SND/ND per EN 300 086 (SINAD meter).

    SND = weighted RMS of the whole output; ND = weighted RMS with the f0
    fundamental removed (±1 bin notch; the demodulated audio of a single
    tone is a pure tone + harmonics, so this isolates noise+distortion).
    """
    w = psoph_weight(audio, fs)
    n = w.size
    # The psoph FIR (257 taps = 32 ms @ 8 kHz) starts from zero state on
    # every call: the first ~32 ms of `w` is a fill-up burst of the tone.
    # In the rectangular spectrum that burst smears into a 30+ Hz-wide
    # pedestal around f0 that the ±1-bin notch cannot remove, flooding the
    # ND. Discard the transient (only for windows long enough to matter).
    tr = 320
    if n > 4 * tr:
        w = w[tr:]
        n = w.size
    X = np.fft.rfft(w)
    k0 = int(round(f0 * n / fs))
    X[max(0, k0 - 1):k0 + 2] = 0.0
    nd = np.fft.irfft(X, n)
    snd = float(np.sqrt(np.mean(w ** 2)))
    ndv = float(np.sqrt(np.mean(nd ** 2)))
    if ndv <= 0:
        return 120.0
    return float(20.0 * np.log10(snd / ndv))


# lfilter denominator for pure-FIR use; float32 keeps the recursion in
# complex64 when the input is complex64
ONE32 = np.float32(1.0)


class ChannelDemod:
    """One 12.5 kHz channel: DDC + decimation + FM demod + de-emphasis + CTCSS."""

    def __init__(self, f0: float, fs: float = 6.048e6, fs_audio: int = 8000,
                 tau: float = 75e-6, dc_block_hz: float = 20.0,
                 demod: str = "pll", pll_f_n: float = 1750.0,
                 pll_zeta: float = 0.707):
        self.f0 = f0
        self.fs = fs
        self.fs_audio = fs_audio
        # FM demod: "pll" = 2nd-order digital PLL (design doc §4.3, loop
        # bandwidth 2-3 kHz) — linear, tracks carrier offset, rejects
        # wideband adjacent noise; "disc" = quadrature discriminator
        # (first-prototype fallback, sin(dphi) nonlinearity at high dev)
        self.demod = demod
        # PLL loop gains per sample: w_n = 2pi*f_n/fs_ch, Kp = 2*zeta*w_n,
        # Ki = w_n^2. f_n = 1750 Hz, zeta = 0.707 -> crossover ~2.4 kHz
        self.pll_f_n = pll_f_n
        self.pll_zeta = pll_zeta
        w_n = 2.0 * np.pi * pll_f_n / (fs / (7 * 9 * 6))
        self.pll_kp = 2.0 * pll_zeta * w_n
        self.pll_ki = w_n * w_n
        self.pll_theta = 0.0
        self.pll_integ = 0.0
        # NCO
        self.phase = 0.0
        self.step = 2.0 * np.pi * f0 / fs
        # decimation stages: 6.048M -> 864k -> 96k -> 16k  (7*9*6 = 378)
        # real (symmetric) FIRs applied to complex data: every stage is
        # DC-centred after the NCO, so a plain lowpass passes +/-fc correctly
        # float32 coeffs keep lfilter in complex64 (2x faster than c128)
        self.h1 = firwin(32, 0.40 * (fs / 7), window=("kaiser", 4.0), fs=fs).astype(np.float32)
        self.h2 = firwin(64, 0.40 * (fs / 63), window=("kaiser", 5.0),
                         fs=fs / 7).astype(np.float32)
        self.h3 = firwin(128, 7500.0, window=("kaiser", 8.6), fs=fs / 63).astype(np.float32)
        self.fs_mid = fs / 63          # 96 kHz
        self.fs_ch = fs / (7 * 9 * 6)  # 16 kHz
        # NCO recurrence state
        self._nco_w = np.complex64(np.exp(-1j * self.step))
        self._nco_seq = np.empty(1 << 16, dtype=np.complex64)
        # lfilter states carried across blocks (resetting them injects a
        # filter transient glitch at every block boundary -> 363 Hz comb)
        self.zi1 = np.zeros(len(self.h1) - 1, dtype=np.complex64)
        self.zi2 = np.zeros(len(self.h2) - 1, dtype=np.complex64)
        self.zi3 = np.zeros(len(self.h3) - 1, dtype=np.complex64)
        # global input sample counter: the decimation grids (mod 7 / 63 / 378)
        # are anchored to the STREAM, not to each block. 16384 is not
        # divisible by 378, so per-block slicing would drift the output rate
        # by 1.5% (a 1 kHz tone would land at 985 Hz)
        self.n_done = 0
        # DC blocker (1st order HP) at 16 kHz — IIR form y[n] = a(x[n]-x[n-1]+y[n-1])
        a_dc = 1.0 / (1.0 + 2.0 * np.pi * dc_block_hz / self.fs_ch)
        self.dc_a = a_dc
        self.dc_b = np.array([a_dc, -a_dc], dtype=np.float32)
        self.dc_aa = np.array([1.0, -a_dc], dtype=np.float32)
        self.dc_state = np.zeros(1, dtype=np.float64)
        # audio lowpass for /2 -> 8 kHz
        b, a = butter(2, 3400.0, btype="low", fs=self.fs_ch)
        self.lp_b, self.lp_a = b, a
        self.lp_state = np.zeros(2, dtype=np.float64)
        # de-emphasis: bilinear H(s) = 1/(1+s*tau)  (1st-order IIR, 1 state)
        c = 2.0 * tau * fs_audio
        self.de_b = np.array([1.0, 1.0]) / (1.0 + c)
        self.de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
        self.de_state = np.zeros(1)
        # discriminator state
        self.prev_x = None
        # carrier detection (512-pt FFT of 16 kHz complex, 31.25 Hz/bin)
        self.fft_n = 512
        self.buf = np.zeros(self.fft_n, dtype=np.complex64)
        self.buf_i = 0
        # CTCSS state
        self.ctcss = None

    # ------------------------------------------------------------------
    def set_ctcss(self, tone_hz: float | None):
        self.ctcss = CTCSS(fs=self.fs_audio, tone_hz=tone_hz) if tone_hz else None

    def carrier_db(self) -> float:
        """FM carrier detect: total power in the audio band (31-3 kHz) vs
        the upper band (4-8 kHz), where a wanted NBFM signal has nothing.
        An FM signal concentrates its power in the audio band (carrier
        sidebands); thermal noise is flat -> ~0 dB on noise, +dB on signal."""
        X = np.fft.fft(self.buf * np.hanning(self.fft_n))
        p = (X.real ** 2 + X.imag ** 2).astype(np.float32)
        audio_sum = float(np.sum(p[1:97]))     # 31 - 3031 Hz at 16 kHz
        ref_sum = float(np.sum(p[129:255]))    # 4.0 - 7.97 kHz
        if ref_sum <= 0:
            return -120.0
        return float(10.0 * np.log10((audio_sum / 96.0) / (ref_sum / 126.0)))

    # ------------------------------------------------------------------
    def process(self, iq: np.ndarray):
        """Feed a block of full-band complex IQ (fs = self.fs).
        Returns (audio float32 @8 kHz, metrics dict)."""
        n = iq.size
        # NCO DDC — self.phase is the phase of the FIRST sample of this block
        if _HAVE_NUMBA:
            if n > self._nco_seq.size:
                self._nco_seq = np.empty(n, dtype=np.complex64)
            _nco_fill(self._nco_seq[:n], self._nco_w, self.phase)
            x = iq * self._nco_seq[:n]
        else:
            ph = self.phase + self.step * np.arange(n, dtype=np.float64)
            x = iq * np.exp(-1j * ph)
        self.phase = (self.phase + self.step * n) % (2.0 * np.pi)
        # decimation: filters at full rate, then slice the globally
        # continuous grid (sample k kept iff global index divisible by 7/63/378)
        # a=np.float32(1.0) keeps lfilter in complex64 (a Python float would
        # promote to complex128: 1.7x slower)
        G = self.n_done
        y1, self.zi1 = lfilter(self.h1, ONE32, x, zi=self.zi1)
        i7 = (-G) % 7
        y7 = y1[i7::7]
        g7 = (G + i7) // 7
        y2, self.zi2 = lfilter(self.h2, ONE32, y7, zi=self.zi2)
        i9 = (-g7) % 9
        y9 = y2[i9::9]
        g63 = (G + i7 + 7 * i9) // 63
        y3, self.zi3 = lfilter(self.h3, ONE32, y9, zi=self.zi3)
        i6 = (-g63) % 6
        x = y3[i6::6]
        # G+i7+7*i9 is the block's first 63-grid index (=63*g63); stepping i6
        # positions on the 63-domain grid costs 63 input samples each, NOT 6:
        # with 6*i6, g16 was floor(g63/6) (off by 1 whenever i6 != 0, i.e. 5 of
        # 6 blocks) -> i8=(-g16)%2 picked the wrong 16k parity -> the audio
        # grid jumped 62.5 us (22.5 deg at 1 kHz) at most block boundaries,
        # smearing every tone with a broadband phase-noise skirt.
        g16 = (G + i7 + 7 * i9 + 63 * i6) // 378
        self.n_done = G + n
        self.x16_last = x          # debug/diagnostic hook (last 16k block)
        # carrier detector buffer
        k = min(self.fft_n - self.buf_i, x.size)
        if k > 0:
            self.buf[self.buf_i:self.buf_i + k] = x[:k]
            self.buf_i += k
            if self.buf_i == self.fft_n:
                self.buf_i = 0
        carr = self.carrier_db()
        if x.size == 0:
            return np.zeros(0, dtype=np.float32), {"carrier_db": carr,
                                                   "ctcss_tone": None, "audio_rms_hz": 0.0,
                                                   "max_step_rad": 0.0}
        if self.prev_x is None:
            self.prev_x = 0.0 + 0.0j
        xp = np.empty_like(x)
        xp[0] = self.prev_x
        xp[1:] = x[:-1]
        self.prev_x = x[-1]
        # max phase step of the 16k signal: for |step| > pi the phase
        # detector wraps (deviation aliasing, |dev| > 8 kHz) — a live
        # interference monitor and the EN 300 086 "no clipping" check
        step_max = float(np.abs(np.angle(x * np.conj(xp))).max())
        if self.demod == "pll":
            dev_hz, self.pll_theta, self.pll_integ = _pll_block(
                x, self.pll_theta, self.pll_integ, self.pll_kp, self.pll_ki,
                self.fs_ch / (2.0 * np.pi))
        else:
            # quadrature discriminator (normalized)
            a = np.imag(x * np.conj(xp))
            a /= (np.abs(x) * np.abs(xp) + 1e-12)
            dev_hz = a * self.fs_ch / (2.0 * np.pi)
        # DC block (1st-order IIR, state carried across blocks)
        y, self.dc_state = lfilter(self.dc_b, self.dc_aa, dev_hz, zi=self.dc_state)
        # lowpass + decimate /2 -> 8 kHz (globally continuous even-index grid)
        y, self.lp_state = lfilter(self.lp_b, self.lp_a, y, zi=self.lp_state)
        i8 = (-g16) % 2
        audio = y[i8::2]
        # de-emphasis
        audio, self.de_state = lfilter(self.de_b, self.de_a, audio, zi=self.de_state)
        # CTCSS
        ctcss_tone = None
        if self.ctcss is not None:
            ctcss_tone = self.ctcss.update(audio)
        metrics = {"carrier_db": carr, "ctcss_tone": ctcss_tone,
                   "audio_rms_hz": float(np.sqrt(np.mean(audio ** 2))) if audio.size else 0.0,
                   "max_step_rad": step_max}
        return audio.astype(np.float32), metrics


class Squelch:
    """Hysteresis + timing squelch on a carrier-detect ratio (dB)."""

    def __init__(self, open_db: float = 12.0, close_db: float = 9.0,
                 open_ms: float = 100.0, close_ms: float = 1000.0, fs: float = 16000.0):
        self.open_db = open_db
        self.close_db = close_db
        self.open_n = int(open_ms / 1000.0 * fs)
        self.close_n = int(close_ms / 1000.0 * fs)
        self.open = False
        self._cnt = 0

    def update(self, ratio_db: float, n_samples: int) -> bool:
        if not self.open:
            if ratio_db >= self.open_db:
                self._cnt += n_samples
                if self._cnt >= self.open_n:
                    self.open = True
                    self._cnt = 0
            else:
                self._cnt = 0
        else:
            if ratio_db <= self.close_db:
                self._cnt += n_samples
                if self._cnt >= self.close_n:
                    self.open = False
                    self._cnt = 0
            else:
                self._cnt = 0
        return self.open


class CTCSS:
    """Goertzel CTCSS detector on 8 kHz audio (ETSI TS 103 236 tones)."""

    def __init__(self, fs: float = 8000.0, tone_hz: float | None = None,
                 tones: np.ndarray = CTCSS_TONES, n: int = 2048,
                 snr_db: float = 20.0, n_windows: int = 3):
        self.fs = fs
        self.tones = np.atleast_1d(np.asarray(tones, dtype=float))
        self.want = tone_hz
        self.n = n
        self.snr_db = snr_db
        self.n_windows = n_windows
        self.buf = np.zeros(n, dtype=np.float32)
        self.buf_i = 0
        self.run = 0
        # Goertzel coefficients for tones and off-tone (f + 25 Hz)
        self.g_tone = self._goertzel(self.tones, fs)
        self.g_off = self._goertzel(self.tones + 25.0, fs)

    @staticmethod
    def _goertzel(freqs: np.ndarray, fs: float) -> np.ndarray:
        # (ntones, 2): [r, 0.5*r] per Goertzel (cos term, half cos term)
        w = 2.0 * np.pi * freqs / fs
        return np.stack([2.0 * np.cos(w), np.cos(w)], axis=1)

    @staticmethod
    def _goertzel_power(x: np.ndarray, g: np.ndarray) -> np.ndarray:
        out = np.empty(g.shape[0], dtype=np.float32)
        for i in range(g.shape[0]):
            r, hc = g[i]
            s1 = s2 = 0.0
            for s in x:
                s0 = s + r * s1 - s2
                s2 = s1
                s1 = s0
            out[i] = s1 * s1 + s2 * s2 - r * s1 * s2
        return out

    def update(self, audio: np.ndarray) -> float | None:
        """Feed 8 kHz audio. Returns the detected tone (Hz) or None."""
        need = self.n - self.buf_i
        if audio.size < need:
            self.buf[self.buf_i:self.buf_i + audio.size] = audio
            self.buf_i += audio.size
            return None
        self.buf[self.buf_i:self.buf_i + need] = audio[:need]
        self.buf_i += need
        if self.buf_i < self.n:
            return None
        x = self.buf.copy()
        self.buf[:] = 0.0
        self.buf_i = 0
        p_tone = self._goertzel_power(x, self.g_tone)
        p_off = self._goertzel_power(x, self.g_off)
        if self.want is not None:
            idx = int(np.argmin(np.abs(self.tones - self.want)))
            ratio_db = 10.0 * np.log10(p_tone[idx] / max(p_off[idx], 1e-20))
            if ratio_db >= self.snr_db:
                self.run += 1
            else:
                self.run = 0
            return float(self.tones[idx]) if self.run >= self.n_windows else None
        # scan all tones
        ratio = 10.0 * np.log10(p_tone / np.maximum(p_off, 1e-20))
        best = int(np.argmax(ratio))
        if ratio[best] >= self.snr_db:
            self.run += 1
        else:
            self.run = 0
        return float(self.tones[best]) if self.run >= self.n_windows else None


class BandSpectrum:
    """Full-band spectrum engine: ONE 8192-pt FFT per window, shared by the
    12.5 kHz channel-power grid (auto channel selection, noise floor) and
    the waterfall display rows.

    Channel bins are computed on an fftshifted frequency axis (each channel
    then maps to one contiguous bin block); the waterfall keeps the natural
    fft order (0..+fs/2, -fs/2..0) for the display."""

    def __init__(self, fs: float, f_center: float, f_min: float, f_max: float,
                 nfft: int = 8192, spacing: float = 12.5e3, half_bw: float = 6.25e3,
                 wf_decim: int = 4):
        self.fs = fs
        self.f_center = f_center
        self.nfft = nfft
        self.win = np.hanning(nfft).astype(np.float32)
        self._winsum = 0.5 * float(np.sum(self.win ** 2))
        self.freqs = f_center + np.fft.fftshift(np.fft.fftfreq(nfft, 1.0 / fs))
        chans = np.arange(np.ceil(f_min / spacing) * spacing,
                          f_max + spacing / 2, spacing)
        self.channels = chans
        self.ch_bins = []
        for fch in chans:
            idx = np.nonzero(np.abs(self.freqs - fch) <= half_bw)[0]
            self.ch_bins.append((int(idx[0]), int(idx[-1])))
        # vectorized channel-grid extraction (cumsum of the power spectrum)
        self._ch_lo = np.array([b[0] for b in self.ch_bins], dtype=np.int64)
        self._ch_hi = np.array([b[1] for b in self.ch_bins], dtype=np.int64)
        self._ch_n = (self._ch_hi - self._ch_lo + 1).astype(np.float64)
        self.buf = np.zeros(nfft, dtype=np.complex64)
        self.buf_i = 0
        self.noise_db = -120.0
        self.ch_power_last: np.ndarray | None = None
        # waterfall
        self.wf_decim = wf_decim
        self.wf_freqs = (f_center + np.fft.fftfreq(nfft, 1.0 / fs))[::wf_decim]
        self.rows: list[np.ndarray] = []
        self.floor = -120.0

    def update(self, iq: np.ndarray) -> tuple[np.ndarray | None, np.ndarray | None, float]:
        """Returns (channel powers linear, channel freqs, noise floor dB);
        None while the FFT buffer fills. Appends a waterfall row per window."""
        n = min(self.nfft - self.buf_i, iq.size)
        if n > 0:
            self.buf[self.buf_i:self.buf_i + n] = iq[:n]
            self.buf_i += n
        if self.buf_i < self.nfft:
            return None, None, self.noise_db
        X = np.fft.fft(self.buf * self.win)
        p = ((X.real ** 2 + X.imag ** 2) / self._winsum).astype(np.float32)
        self.buf[:] = 0.0
        self.buf_i = 0
        # channel grid (shifted axis) — vectorized via cumsum
        p_sh = np.fft.fftshift(p)
        cum = np.concatenate(([0.0], np.cumsum(p_sh, dtype=np.float64)))
        cp = ((cum[self._ch_hi + 1] - cum[self._ch_lo]) / self._ch_n).astype(np.float32)
        self.noise_db = float(10.0 * np.log10(np.median(cp) + 1e-30))
        self.ch_power_last = cp
        # waterfall row (fft order)
        db = 10.0 * np.log10(p + 1e-30)
        self.floor = 0.95 * self.floor + 0.05 * float(np.median(db))
        self.rows.append(db[::self.wf_decim])
        return cp, self.channels, self.noise_db

    def to_png(self, path: str) -> None:
        if not self.rows:
            return
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        img = np.array(self.rows)
        img = img - self.floor
        lo, hi = np.percentile(img, 1), np.percentile(img, 99.9)
        img = np.clip(img, lo, max(hi, lo + 10))
        fig, ax = plt.subplots(figsize=(12, 5))
        im = ax.imshow(img, aspect="auto", origin="lower",
                       extent=[self.wf_freqs[0] / 1e6, self.wf_freqs[-1] / 1e6,
                               0, len(self.rows) * self.nfft / self.fs],
                       cmap="magma", vmin=lo, vmax=max(hi, lo + 10))
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Time (s)")
        ax.set_title(f"Waterfall — {len(self.rows)} spectra, {self.nfft} pts")
        fig.colorbar(im, ax=ax, label="dB over floor")
        fig.tight_layout()
        fig.savefig(path, dpi=110)
        plt.close(fig)


# back-compat name (tests)
ActivityMonitor = BandSpectrum
