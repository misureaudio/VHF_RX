"""Offline tests for the NBFM DSP chain + ctypes struct layout.

Run:  D:/Source/hermes-dir/.venv/Scripts/python.exe tests/test_nbfm.py
No hardware required (struct test only loads the DLL).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from nbfm import (ActivityMonitor, ChannelDemod, CTCSS, Squelch,  # noqa: E402
                  pll_tracking_response, psoph_weight, snd_nd)

FS = 6.048e6
FAILS: list[str] = []

# must match the ChannelDemod defaults (design doc §4.3: 2nd-order PLL,
# loop bandwidth 2-3 kHz -> f_n = 1750 Hz, zeta = 0.707, crossover ~2.4 kHz)
PLL_F_N = 1750.0
PLL_ZETA = 0.707


def check(name: str, cond: bool, detail: str = ""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {detail}")
    if not cond:
        FAILS.append(name)


def fm_modulate(audio_hz: np.ndarray, f_carrier: float, fs: float,
                amp: float = 0.5) -> np.ndarray:
    """audio_hz: instantaneous frequency deviation in Hz. Returns complex baseband."""
    phi = 2.0 * np.pi * np.cumsum(audio_hz) / fs
    t = np.arange(audio_hz.size) / fs
    return amp * np.exp(1j * (phi + 2.0 * np.pi * f_carrier * t))


def feed(demod: ChannelDemod, iq: np.ndarray, block: int = 16384):
    audio = []
    for i in range(0, iq.size, block):
        a, _ = demod.process(iq[i:i + block])
        audio.append(a)
    return np.concatenate(audio) if audio else np.zeros(0)


def tone_freqs(x: np.ndarray, fs: float = 8000) -> dict:
    """Return {freq: amplitude} for the strongest 4 spectral lines.
    Hanning coherent gain 0.5 -> amplitude = 4|X[k]|/N."""
    X = np.fft.rfft(x * np.hanning(x.size))
    f = np.fft.rfftfreq(x.size, 1 / fs)
    p = np.abs(X)
    out = {}
    for _ in range(4):
        k = int(np.argmax(p))
        if p[k] <= 0:
            break
        out[float(f[k])] = float(4.0 * p[k] / x.size)
        p[max(0, k - 2):k + 3] = 0
    return out


def test_device_profiles():
    """v1: per-device profile table (no hardware — pure data checks).

    The RSP2pro has no device ID of its own in API v3: it enumerates as
    hwVer 2 (RSP2 family) and is configured through the RSP2-family
    structs. Selection must therefore classify by hwVer, and the sample
    scale / LNA limits must come from the profile, not be hard-coded."""
    import sdrplay_api as s
    # family classification (match by hwVer, not serial prefix)
    check("family: hwVer 4 -> rspdx", s.family_of(4) == "rspdx")
    check("family: hwVer 2 -> rsp2pro (RSP2 family)", s.family_of(2) == "rsp2pro")
    check("family: unknown hwVer -> None", s.family_of(99) is None)
    dx, p2 = s.PROFILES["rspdx"], s.PROFILES["rsp2pro"]
    # sample scale: RSPdx 14-bit in int16; RSP2/RSP2pro native 16-bit (measured)
    check("scale: RSPdx full scale 8192", dx.sample_scale == 8192.0)
    check("scale: RSP2/RSP2pro full scale 32768", p2.sample_scale == 32768.0)
    # LNA state counts per band (out of range fails Init with OutOfRange)
    check("LNA max: RSPdx 26 / 20", (dx.lna_max_2m, dx.lna_max_70cm) == (26, 20))
    check("LNA max: RSP2/RSP2pro 8 / 5", (p2.lna_max_2m, p2.lna_max_70cm) == (8, 5))
    # min API version per family
    check("min API: RSPdx 3.06, RSP2/RSP2pro 3.08",
          dx.min_api == 3.06 and p2.min_api == 3.08)
    # antenna routing values (RSP2 uses 5/6, not 0/1)
    check("antenna: RSP2 port A = 5, port B = 6",
          s.RSP2_ANT_A == 5 and s.RSP2_ANT_B == 6)
    check("antenna: RSP2 update reason 0x200 (main, not ext1)",
          s.UPD_RSP2_ANTENNA == 0x200)


def test_structs():
    import ctypes as ct
    import sdrplay_api as s
    # hand-computed from the 3.11 headers (MSVC x64 alignment); note
    # TunerParamsT is 72 bytes with rfFreq at offset 40 (double alignment)
    expected = {
        s.FsFreqT: 16, s.SyncUpdateT: 8, s.ResetFlagsT: 3,
        s.Rsp1aParamsT: 2, s.Rsp2ParamsT: 1, s.RspDuoParamsT: 4,
        s.RspDxParamsT: 12, s.DevParamsT: 64,
        s.GainValuesT: 12, s.GainT: 24, s.RfFreqT: 16,
        s.DcOffsetTunerT: 12, s.TunerParamsT: 72,
        s.DcOffsetT: 2, s.DecimationT: 3, s.AgcT: 20, s.ControlParamsT: 32,
        s.Rsp1aTunerParamsT: 1, s.Rsp2TunerParamsT: 16,
        s.RspDuoTunerParamsT: 12, s.RspDxTunerParamsT: 4,
        s.RxChannelParamsT: 144, s.DeviceT: 96, s.DeviceParamsT: 24,
        s.ErrorInfoT: 1540, s.StreamCbParamsT: 20, s.GainCbParamT: 16,
    }
    for cls, size in expected.items():
        check(f"sizeof({cls.__name__})=={size}", ct.sizeof(cls) == size,
              f"(ctypes: {ct.sizeof(cls)})")
    def _offsetof(cls, field):
        inst = cls()
        return ct.addressof(getattr(inst, field)) - ct.addressof(inst)

    check("TunerParamsT.rfFreq offset == 40",
          _offsetof(s.TunerParamsT, "rfFreq") == 40,
          f"(ctypes: {_offsetof(s.TunerParamsT, 'rfFreq')})")
    rx = s.SDRPlay()  # loads DLL, opens API
    check("api version >= 3.06", rx.api_version >= 3.06, f"(v{rx.api_version:.2f})")
    rx._lib.sdrplay_api_Close()


def test_tone_recovery():
    n = int(2 * FS)
    t = np.arange(n) / FS
    audio = 1800.0 * np.sin(2 * np.pi * 1000.0 * t)   # 1 kHz @ 1.8 kHz dev
    iq = fm_modulate(audio, 145.0e6, FS)
    d = ChannelDemod(145.0e6, FS)
    out = feed(d, iq)
    lines = tone_freqs(out)
    f1 = max(lines, key=lines.get)
    check("tone freq ~1 kHz", abs(f1 - 1000) < 20, f"(got {f1:.1f} Hz)")
    # theoretical amplitude: the PLL demod is LINEAR in the deviation, so
    # fundamental = D * |T_f(1 kHz)| (the exact discrete loop tracking
    # response, incl. its ~+1.7 dB bump near f_n); then the linear stages
    # (DC block, LP, de-emphasis) each apply their |H(1 kHz)|
    w = 2j * np.pi * 1000 / 16000
    h_lp = abs((d.lp_b[0] + d.lp_b[1]*np.exp(-w) + d.lp_b[2]*np.exp(-2*w)) /
               (d.lp_a[0] + d.lp_a[1]*np.exp(-w) + d.lp_a[2]*np.exp(-2*w)))
    h_dc = abs(d.dc_a * (1 - np.exp(-w)) / (1 - d.dc_a * np.exp(-w)))
    # digital de-emphasis (bilinear) evaluated at 1 kHz / 8 kHz
    c = 2 * 75e-6 * 8000
    de_b = np.array([1.0, 1.0]) / (1.0 + c)
    de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
    w8 = 2j * np.pi * 1000 / 8000
    h_de = abs((de_b[0] + de_b[1]*np.exp(-w8)) / (de_a[0] + de_a[1]*np.exp(-w8)))
    T = pll_tracking_response(1000.0, PLL_F_N, PLL_ZETA)
    expected = 1800.0 * T * h_de * h_lp * h_dc
    amp = lines[f1]
    check("deviation ~1.8 kHz (PLL linear + stages)",
          abs(amp - expected) / expected < 0.02,
          f"(got {amp:.0f} Hz, expect {expected:.0f})")


def test_de_emphasis():
    n = int(2 * FS)
    # 500 Hz and 2 kHz tones, EQUAL deviation, WITH 6 dB/oct pre-emphasis
    # applied at the AUDIO rate (8 kHz), as a real transmitter does.
    # Pre-emphasis cancels de-emphasis exactly. The PLL demod is linear in
    # the deviation, so each tone's net amplitude = 600*|Hpe(f)|*|Tf(f)|*
    # |Hlp(f)|*|Hdc(f)|*|Hde(f)| with |Hde| = 1/|Hpe| -> the 600 and the
    # pre-emphasis cancel in the ratio, leaving |Tf|*|Hlp|*|Hdc|.
    n8 = n // 756                      # 6.048 MHz = 756 * 8 kHz
    t8 = np.arange(n8) / 8000.0
    raw8 = 600.0 * (np.sin(2 * np.pi * 500 * t8) + np.sin(2 * np.pi * 2000 * t8))
    import scipy.signal as sps
    c = 2 * 75e-6 * 8000
    de_b = np.array([1.0, 1.0]) / (1.0 + c)
    de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
    pre8 = sps.lfilter(de_a, de_b, raw8)
    pre = np.repeat(pre8, 756)         # back to the RF rate for modulation
    iq = fm_modulate(pre, 434.6e6, FS)
    d = ChannelDemod(434.6e6, FS)
    out = feed(d, iq)[:16000]          # exactly 2.000 s: both tones on bins
    lines = tone_freqs(out)
    f5 = min(lines, key=lambda f: abs(f - 500))
    f2 = min(lines, key=lambda f: abs(f - 2000))
    ratio_db = 20 * np.log10(lines[f5] / lines[f2])
    # theoretical ratio: see the comment above — PLL linear, pre-emphasis
    # cancels de-emphasis exactly
    def h_lp_db(f):
        w = 2j * np.pi * f / 16000
        H = (d.lp_b[0] + d.lp_b[1]*np.exp(-w) + d.lp_b[2]*np.exp(-2*w)) / \
            (d.lp_a[0] + d.lp_a[1]*np.exp(-w) + d.lp_a[2]*np.exp(-2*w))
        return abs(H)
    def h_dc_db(f):
        w = 2j * np.pi * f / 16000
        H = d.dc_a * (1 - np.exp(-w)) / (1 - d.dc_a * np.exp(-w))
        return abs(H)
    num = (pll_tracking_response(500.0, PLL_F_N, PLL_ZETA)
           * h_lp_db(500) * h_dc_db(500))
    den = (pll_tracking_response(2000.0, PLL_F_N, PLL_ZETA)
           * h_lp_db(2000) * h_dc_db(2000))
    expected = 20 * np.log10(num / den)
    check("de-emphasis restores pre-emphasis (500 vs 2 kHz)",
          abs(ratio_db - expected) < 2.5,
          f"(500/2kHz = {ratio_db:+.1f} dB, expect {expected:+.1f} dB)")


def test_pll_offset_tracking():
    """PLL-specific (design doc §4.3): the loop locks to a carrier
    frequency offset and removes it from the audio. +150 Hz offset on a
    1 kHz tone: the tone is intact and the audio DC (after the 20 Hz
    DC blocker) is ~0."""
    n = int(1 * FS)
    t = np.arange(n) / FS
    d = ChannelDemod(145.0e6, FS)
    out = feed(d, fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t),
                              145.0e6 + 150.0, FS))
    lines = tone_freqs(out)
    f1 = max(lines, key=lines.get)
    check("PLL: 1 kHz tone with +150 Hz carrier offset present",
          abs(f1 - 1000) < 20, f"(got {f1:.1f} Hz)")
    dc = float(np.mean(out[1600:]))
    check("PLL: carrier offset removed (audio DC ~ 0 Hz)",
          abs(dc) < 5.0, f"(DC {dc:+.2f} Hz)")


def test_pll_low_freq_large_dev():
    """PLL startup regression (2026-10-25): 300 Hz tone at 1.8 kHz
    deviation (beta = 6, the hardest low-frequency voice case). The
    decimation FIRs start from zero state, so the first ~10 16k samples
    of a stream are near-zero partial sums whose phase is numerically
    garbage. Starting the loop there injects a one-shot phase kick (up to
    ~pi in the PD error); for this tone the kick pushed the loop out of
    the capture basin into a stable slip orbit — the tone then came out
    at 72.6 % of the tracking-response model with its harmonics at
    -6 dB, forever. The loop now amplitude-gates the PD (|x| < 0.05) and
    initializes the NCO to the first valid sample's phase (e0 = 0). The
    tone must match the model and no slip harmonics may appear."""
    n = int(2 * FS)
    t = np.arange(n) / FS
    d = ChannelDemod(145.0e6, FS)
    out = feed(d, fm_modulate(1800.0 * np.sin(2 * np.pi * 300 * t),
                              145.0e6, FS))[:16000]  # exactly 2.000 s
    N = out.size
    X = np.abs(np.fft.rfft(out * np.hanning(N)))

    def amp(f):
        return float(4.0 * X[int(round(f * N / 8000))] / N)

    w = 2j * np.pi * 300 / 16000
    h_lp = abs((d.lp_b[0] + d.lp_b[1]*np.exp(-w) + d.lp_b[2]*np.exp(-2*w)) /
               (d.lp_a[0] + d.lp_a[1]*np.exp(-w) + d.lp_a[2]*np.exp(-2*w)))
    h_dc = abs(d.dc_a * (1 - np.exp(-w)) / (1 - d.dc_a * np.exp(-w)))
    c = 2 * 75e-6 * 8000
    de_b = np.array([1.0, 1.0]) / (1.0 + c)
    de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
    w8 = 2j * np.pi * 300 / 8000
    h_de = abs((de_b[0] + de_b[1]*np.exp(-w8)) / (de_a[0] + de_a[1]*np.exp(-w8)))
    expected = 1800.0 * pll_tracking_response(300.0, PLL_F_N, PLL_ZETA) \
        * h_de * h_lp * h_dc
    check("PLL: 300 Hz @ 1.8 kHz dev (beta=6) matches tracking model",
          abs(amp(300.0) - expected) / expected < 0.02,
          f"(got {amp(300.0):.0f} Hz, expect {expected:.0f})")
    # slip-lock signature: the tone's harmonics at -6 dB (measured pre-fix)
    worst_h = max(20 * np.log10(amp(f) / amp(300.0))
                  for f in (600.0, 900.0, 1200.0, 1500.0, 1800.0))
    check("PLL: no slip harmonics (300 Hz beta=6, <= -30 dB)",
          worst_h <= -30.0, f"(worst harmonic {worst_h:+.1f} dB rel)")


def test_multitone_linearity():
    """Multitone (400/1000/2000 Hz @ 800 Hz deviation each, 2.4 kHz peak):
    the PLL demod is linear in the deviation, so the output must be the sum
    of the three single-tone responses — no intermodulation products in
    band. A single-tone test cannot see inter-tone interaction (loop
    nonlinearity, saturation, gain modulation): it shows up as new tones at
    the 2nd-order products of the three that fall in-band and don't
    coincide with a tone: 600 (1000-400), 800 (2*400), 1400 (400+1000),
    1600 (2000-400), 2400 (400+2000), 3000 (1000+2000)."""
    n = int(2 * FS)
    t = np.arange(n) / FS
    A = 800.0
    f_tones = (400.0, 1000.0, 2000.0)
    audio = A * (np.sin(2 * np.pi * 400 * t)
                 + np.sin(2 * np.pi * 1000 * t)
                 + np.sin(2 * np.pi * 2000 * t))
    d = ChannelDemod(145.0e6, FS)
    out = feed(d, fm_modulate(audio, 145.0e6, FS))[:16000]  # exactly 2.000 s
    N = out.size
    X = np.abs(np.fft.rfft(out * np.hanning(N)))

    def amp(f):
        return float(4.0 * X[int(round(f * N / 8000))] / N)

    def h_lp(f):
        w = 2j * np.pi * f / 16000
        return abs((d.lp_b[0] + d.lp_b[1]*np.exp(-w) + d.lp_b[2]*np.exp(-2*w)) /
                   (d.lp_a[0] + d.lp_a[1]*np.exp(-w) + d.lp_a[2]*np.exp(-2*w)))

    def h_dc(f):
        w = 2j * np.pi * f / 16000
        return abs(d.dc_a * (1 - np.exp(-w)) / (1 - d.dc_a * np.exp(-w)))

    c = 2 * 75e-6 * 8000
    de_b = np.array([1.0, 1.0]) / (1.0 + c)
    de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])

    def h_de(f):
        w8 = 2j * np.pi * f / 8000
        return abs((de_b[0] + de_b[1]*np.exp(-w8)) / (de_a[0] + de_a[1]*np.exp(-w8)))

    # 1) superposition: each tone in the mix must match that tone ALONE —
    #    the solo baseline isolates inter-tone interaction from the small
    #    full-chain model deviation near the loop crossover (the model is
    #    exact for the loop alone; the 6M NCO zero-order hold adds a small
    #    feedback delay, documented in diag/README.md: -0.6 % @ 1 kHz,
    #    -2.5 % @ 2 kHz, -5.8 % @ 3 kHz)
    for f in f_tones:
        solo = A * np.sin(2 * np.pi * f * t)
        d_s = ChannelDemod(145.0e6, FS)
        out_s = feed(d_s, fm_modulate(solo, 145.0e6, FS))[:16000]
        Xs = np.abs(np.fft.rfft(out_s * np.hanning(N)))
        amp_solo = float(4.0 * Xs[int(round(f * N / 8000))] / N)
        check(f"multitone: {f:.0f} Hz tone unchanged by the other two tones",
              abs(amp(f) - amp_solo) / amp_solo < 0.02,
              f"(mix {amp(f):.0f} vs solo {amp_solo:.0f})")
    # 2) intermodulation floor: the in-band 2nd-order products must sit at
    #    the numerical floor (the PLL is linear in deviation; the arctan2
    #    PD is amplitude-invariant, so no gain modulation either)
    tone_amp = float(np.mean([amp(f) for f in f_tones]))
    im_bins = (600.0, 800.0, 1400.0, 1600.0, 2400.0, 3000.0)
    worst_f, worst_db = 0.0, -1e9
    for f in im_bins:
        rel = 20 * np.log10(amp(f) / tone_amp)
        if rel > worst_db:
            worst_f, worst_db = f, rel
    check("multitone: IM products <= -40 dB below tones (PLL linear)",
          worst_db <= -40.0, f"(worst {worst_f:.0f} Hz at {worst_db:+.1f} dB rel)")
    # 3) the theoretical model (exact T_f x stages) for reference
    for f in f_tones:
        expected = (A * pll_tracking_response(f, PLL_F_N, PLL_ZETA)
                    * h_de(f) * h_lp(f) * h_dc(f))
        print(f"  [info] {f:.0f} Hz: mix {amp(f):.0f} vs model {expected:.0f} "
              f"({amp(f) / expected - 1:+.1%})")
    # 4) the discriminator fallback has real sin(delta-phi) folding at
    #    2.4 kHz peak deviation (max 16k phase step 0.94 rad) — report,
    #    don't gate (design doc §4.3: why the PLL is the default)
    d2 = ChannelDemod(145.0e6, FS, demod="disc")
    out2 = feed(d2, fm_modulate(audio, 145.0e6, FS))[:16000]
    X2 = np.abs(np.fft.rfft(out2 * np.hanning(N)))
    amp2 = lambda f: float(4.0 * X2[int(round(f * N / 8000))] / N)
    tone_amp2 = float(np.mean([amp2(f) for f in f_tones]))
    worst2 = max(20 * np.log10(amp2(f) / tone_amp2) for f in im_bins)
    print(f"  [info] discriminator IM floor at the same deviation: "
          f"{worst2:+.1f} dB rel (sin folding)")


def test_adjacent_channel():
    n = int(2 * FS)
    t = np.arange(n) / FS
    want = fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS, amp=0.5)
    intr = fm_modulate(1800.0 * np.sin(2 * np.pi * 2000 * t), 145.0125e6, FS, amp=0.5)
    d = ChannelDemod(145.0e6, FS)
    out = feed(d, want + intr)
    lines = tone_freqs(out)
    f1 = min(lines, key=lambda f: abs(f - 1000))
    f2 = min(lines, key=lambda f: abs(f - 2000))
    rej = 20 * np.log10(lines[f1] / lines[f2])
    check("adjacent (12.5 kHz) rejected >20 dB", rej > 20, f"(got {rej:.1f} dB)")


def test_squelch_and_ctcss():
    n = int(2 * FS)
    t = np.arange(n) / FS
    rng = np.random.default_rng(0)
    noise = 0.01 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    voice = 1500.0 * np.sin(2 * np.pi * 300 * t)
    ctcss_tone = np.sin(2 * np.pi * 88.5 * t)
    # build: 0.5 s noise, 1.5 s signal with CTCSS
    sig = np.zeros(n, dtype=np.complex128)
    sig[:n // 4] = noise[:n // 4]
    m = fm_modulate(voice + 200.0 * ctcss_tone, 145.5e6, FS, amp=0.5)
    sig[n // 4:] = m[:3 * n // 4] + noise[n // 4:]
    d = ChannelDemod(145.5e6, FS)
    d.set_ctcss(88.5)
    sq = Squelch(open_ms=100, close_ms=500, fs=16000)
    open_blocks = 0
    total_blocks = 0
    detected = None
    for i in range(0, n, 16384):
        a, met = d.process(sig[i:i + 16384])
        if sq.update(met["carrier_db"], a.size * 2):
            open_blocks += 1
        total_blocks += 1
        if met["ctcss_tone"]:
            detected = met["ctcss_tone"]
    check("squelch opens on signal, closed on noise",
          0 < open_blocks < total_blocks, f"({open_blocks}/{total_blocks} blocks open)")
    check("CTCSS 88.5 Hz detected", detected == 88.5, f"(got {detected})")
    # wrong tone (outside +/-1%) must NOT latch
    d2 = ChannelDemod(145.5e6, FS)
    d2.set_ctcss(85.4)
    wrong = None
    for i in range(0, n, 16384):
        _, met = d2.process(sig[i:i + 16384])
        if met["ctcss_tone"]:
            wrong = met["ctcss_tone"]
    check("CTCSS 85.4 Hz not latched on 88.5 Hz signal", wrong is None, f"(got {wrong})")


def test_activity():
    n = int(1 * FS)
    t = np.arange(n) / FS
    sig = fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t), 145.0125e6, FS, amp=0.3)
    am = ActivityMonitor(FS, 145.0e6, 144.0e6, 146.0e6)
    cp = None
    last_chans = None
    for i in range(0, n, 16384):
        c, chans, nf = am.update(sig[i:i + 16384])
        if c is not None:
            cp, last_chans = c, chans
    assert cp is not None and last_chans is not None
    k = int(np.argmin(np.abs(last_chans - 145.0125e6)))
    snr = 10 * np.log10(cp[k] / np.median(cp))
    check("activity monitor finds 145.0125 MHz", snr > 15, f"(SNR {snr:.1f} dB)")


# ---------------------------------------------------------------------------
# EN 300 086 V2.1.2 conformance (12.5 kHz channel separation)
# ---------------------------------------------------------------------------
MPFD_12P5 = 2500.0          # Table 3: max permissible deviation, 12.5 kHz
NORM_DEV = 0.6 * MPFD_12P5  # cl. 6.1: normal test modulation = 1 kHz @ 60% MPFD


def _demod_audio(want: np.ndarray, intr: np.ndarray | None):
    d = ChannelDemod(145.0e6, FS)
    return feed(d, want + intr if intr is not None else want)


def test_co_channel():
    """EN 300 086 cl. 8.3: co-channel rejection, 12.5 kHz spacing.

    Wanted: 1 kHz @ 60% MPFD (cl. 6.1). Unwanted: 400 Hz @ 60% MPFD, at the
    same nominal frequency, displaced by 0 / ±1.5 kHz / ±3 kHz. Raise the
    unwanted until the psophometric SND/ND drops to 14 dB; the level ratio
    (unwanted/wanted) at that point must lie between -12.0 and 0 dB (cl.
    8.3.3). Noise-free DSP chain (the standard's 6 dBµV wanted level is
    arbitrary here — only the ratio matters).
    """
    n = int(1 * FS)
    t = np.arange(n) / FS
    want = fm_modulate(NORM_DEV * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS, amp=0.5)

    def ratio_db(off: float) -> float:
        # cl. 8.3: "the amplitude of the unwanted test signal shall be
        # adjusted until the SND/ND ratio ... is reduced to 14 dB" — i.e.
        # start with the unwanted signal quiet (SND/ND high) and RAISE it.
        # Scanning down from 0 dB would mis-report a crossing above 0 dB as
        # a borderline pass (the SND/ND is already < 14 dB at equal level
        # for in-band offsets, because the discriminator demodulates both
        # signals).
        r_prev = db_prev = None
        for db10 in range(-400, 101, 20):        # -40, -38, ..., +10 dB
            db = db10 / 10.0
            intr = fm_modulate(NORM_DEV * np.sin(2 * np.pi * 400 * t),
                               145.0e6 + off, FS, amp=0.5 * 10 ** (db / 20))
            r = snd_nd(_demod_audio(want, intr), 8000, 1000.0)
            if r <= 14.0:
                if r_prev is None:
                    return db                    # already below 14 at -40 dB
                frac = (r_prev - 14.0) / (r_prev - r)
                return db_prev + (db - db_prev) * frac
            r_prev, db_prev = r, db
        return 10.0                              # not reached: worse than +10

    worst = 0.0
    for off in (0, 1500, -1500, 3000, -3000):
        cross = ratio_db(off)
        worst = min(worst, cross)
        check(f"co-channel crossing in [-12, 0] dB @ {off:+d} Hz",
              -12.0 <= cross <= 0.0, f"(crossing {cross:+.1f} dB)")
    print(f"  [info] co-channel rejection (lowest of five): {worst:+.1f} dB")


def test_adjacent_channel_selectivity():
    """EN 300 086 cl. 8.4: adjacent channel selectivity, 12.5 kHz spacing.

    Unwanted: 400 Hz @ 60% MPFD at one channel separation (±12.5 kHz). Raise
    it until the psophometric SND/ND drops to 14 dB; the level ratio must be
    >= 60.0 dB (Table 8, normal conditions).
    """
    n = int(1 * FS)
    t = np.arange(n) / FS
    want = fm_modulate(NORM_DEV * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS, amp=0.5)

    def ratio_db(off: float) -> float:
        r_prev = db_prev = None
        for db10 in range(-200, 801, 50):        # -20, -15, ..., +80 dB
            db = db10 / 10.0
            intr = fm_modulate(NORM_DEV * np.sin(2 * np.pi * 400 * t),
                               145.0e6 + off, FS, amp=0.5 * 10 ** (db / 20))
            r = snd_nd(_demod_audio(want, intr), 8000, 1000.0)
            if r <= 14.0:
                if r_prev is None:
                    return db
                frac = (r_prev - 14.0) / (r_prev - r)
                return db_prev + (db - db_prev) * frac
            r_prev, db_prev = r, db
        return 80.0                               # not reached: better than +80

    for off in (12500, -12500):
        cross = ratio_db(off)
        check(f"adjacent selectivity >= 60 dB @ {off/1000:+.1f} kHz",
              cross >= 60.0, f"(crossing {cross:+.1f} dB)")


def test_deviation_headroom():
    """Test plan 6.2: ±3 kHz (IARU standard) and ±5 kHz (EN 300 086 25 kHz
    MPFD headroom) 1 kHz tones — no clipping, correct fundamental, and the
    psophometric SND/ND (the standard's SINAD meter metric)."""
    n = int(1 * FS)
    t = np.arange(n) / FS
    for dev in (3000.0, 5000.0):
        d = ChannelDemod(145.0e6, FS)
        iq = fm_modulate(dev * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS)
        audio, step_max = [], 0.0
        for i in range(0, iq.size, 16384):
            a, met = d.process(iq[i:i + 16384])
            audio.append(a)
            if i > 0:                       # skip the decimation transient
                step_max = max(step_max, met["max_step_rad"])
        a = np.concatenate(audio)
        # "no clipping" for this chain: the PLL phase detector wraps at pi,
        # so the 16k phase step must stay < pi, i.e. |instantaneous
        # deviation| < 8 kHz (deviation aliasing)
        check(f"no deviation aliasing @ ±{dev/1000:.0f} kHz dev",
              step_max < np.pi, f"(max step {np.degrees(step_max):.0f} deg)")
        lines = tone_freqs(a)
        f1 = max(lines, key=lines.get)
        check(f"1 kHz fundamental present @ ±{dev/1000:.0f} kHz dev",
              abs(f1 - 1000) < 20, f"(got {f1:.1f} Hz)")
        # PLL theory: linear in deviation -> fundamental = D*|Tf(1 kHz)|*stages
        w = 2j * np.pi * 1000 / 16000
        h_lp = abs((d.lp_b[0] + d.lp_b[1]*np.exp(-w) + d.lp_b[2]*np.exp(-2*w)) /
                   (d.lp_a[0] + d.lp_a[1]*np.exp(-w) + d.lp_a[2]*np.exp(-2*w)))
        h_dc = abs(d.dc_a * (1 - np.exp(-w)) / (1 - d.dc_a * np.exp(-w)))
        c = 2 * 75e-6 * 8000
        de_b = np.array([1.0, 1.0]) / (1.0 + c)
        de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
        w8 = 2j * np.pi * 1000 / 8000
        h_de = abs((de_b[0] + de_b[1]*np.exp(-w8)) / (de_a[0] + de_a[1]*np.exp(-w8)))
        T = pll_tracking_response(1000.0, PLL_F_N, PLL_ZETA)
        expected = dev * T * h_de * h_lp * h_dc
        if dev == 3000.0:
            check(f"fundamental amplitude @ ±{dev/1000:.0f} kHz dev within 5%",
                  abs(lines[f1] - expected) / expected < 0.05,
                  f"(got {lines[f1]:.0f}, expect {expected:.0f})")
        r = snd_nd(a, 8000, 1000.0)
        if dev == 3000.0:
            check("SND/ND >= 20 dB @ ±3 kHz dev", r >= 20.0, f"(got {r:.1f} dB)")
        else:
            # ±5 kHz is the 25 kHz MPFD — OUT of spec for our 12.5 kHz
            # receiver (IARU FM.1 = ±3 kHz, EN 300 086 Table 3 = ±2.5 kHz).
            # At ±5 kHz the startup frequency error exceeds the loop
            # pull-in range, so the loop runs in a slip state that
            # preserves the deviation (audio = deviation + a DC offset
            # removed by the DC blocker); report, don't gate.
            print(f"  [info] ±5 kHz (out-of-spec headroom): fundamental "
                  f"{lines[f1]:.0f} Hz (linear expect {expected:.0f}), "
                  f"SND/ND {r:.1f} dB")


def test_audio_player():
    """Audio output stage (playback.py): squelch gating + 5 ms ramp,
    fixed gain, multi-channel mix + clip, 8k->48k resample, ring and
    callback. Pure numpy — the PortAudio callback is called directly,
    no sound device required."""
    from playback import AudioPlayer, FS_IN, FS_OUT, GAIN_DEFAULT
    fs = FS_IN
    n = 4000                                  # 500 ms
    t = np.arange(n) / fs
    dev = 1800.0 * np.sin(2 * np.pi * 1000.0 * t)   # 1 kHz @ 1.8 kHz dev

    def drain(p):
        with p.lock:
            out = np.concatenate(list(p.ring)) if p.ring else np.zeros(0)
            p.ring.clear()
            p.ring_n = 0
        return out

    # single open channel: gain, resample ratio, tone freq, ramp start
    p = AudioPlayer()
    for i in range(0, n, 200):
        p.push(145.5e6, dev[i:i + 200], True)
        p.flush()
    out = drain(p)
    check("player: 8k->48k resample ratio 6", out.size == n * 6, f"(got {out.size})")
    lines = tone_freqs(out, FS_OUT)
    f1 = max(lines, key=lines.get)
    check("player: tone at 1 kHz after resample", abs(f1 - 1000) < 8,
          f"(got {f1:.1f} Hz)")
    expected = 1800.0 * GAIN_DEFAULT          # 0.36 FS
    check("player: fixed gain 1/5000 (1.8 kHz dev -> 0.36 FS)",
          abs(lines[f1] - expected) / expected < 0.03,
          f"(got {lines[f1]:.3f}, expect {expected:.3f})")
    check("player: squelch ramp starts at 0", abs(out[0]) < 0.05 * expected,
          f"(first sample {out[0]:.4f})")

    # squelch-closed channel contributes nothing
    p2 = AudioPlayer()
    for i in range(0, n, 200):
        p2.push(145.5e6, dev[i:i + 200], False)
        p2.flush()
    check("player: squelch-closed channel contributes nothing",
          drain(p2).size == 0)

    # two identical open channels: linear sum, no clip (2 x 0.36 = 0.72)
    p3 = AudioPlayer()
    for i in range(0, n, 200):
        p3.push(145.5e6, dev[i:i + 200], True)
        p3.push(145.6e6, dev[i:i + 200], True)
        p3.flush()
    out3 = drain(p3)
    lines3 = tone_freqs(out3, FS_OUT)
    f3 = max(lines3, key=lines3.get)
    check("player: two-channel mix sums linearly (0.72 FS)",
          abs(lines3[f3] - 2.0 * expected) / (2.0 * expected) < 0.03,
          f"(got {lines3[f3]:.3f}, expect {2.0 * expected:.3f})")

    # overdrive: 2 x 3 kHz dev = 1.2 -> clipped at +/-0.95
    dev_big = 3000.0 * np.sin(2 * np.pi * 1000.0 * t)
    p4 = AudioPlayer()
    for i in range(0, n, 200):
        p4.push(145.5e6, dev_big[i:i + 200], True)
        p4.push(145.6e6, dev_big[i:i + 200], True)
        p4.flush()
    out4 = drain(p4)
    # the plateau sits at 0.95; the resample FIR rings ~2% at the hard
    # clip edge (negative sidelobes), so check plateau + no >FS peak,
    # not the exact peak
    peak4 = float(np.max(np.abs(out4)))
    p99 = float(np.percentile(np.abs(out4), 99))
    check("player: overdriven mix clipped at +/-0.95 (no >FS peak)",
          p99 > 0.93 and peak4 < 1.0,
          f"(99th pct {p99:.3f}, peak {peak4:.3f})")

    # callback: fills outdata from the ring; silence + underrun count when starved
    p5 = AudioPlayer()
    for i in range(0, n, 200):
        p5.push(145.5e6, dev[i:i + 200], True)
        p5.flush()
    frames = 480
    outdata = np.zeros((frames, 1), dtype=np.float32)
    p5._cb(outdata, frames, None, None)
    check("player: callback fills outdata from ring",
          float(np.abs(outdata).max()) > 0.1,
          f"(peak {float(np.abs(outdata).max()):.3f})")
    with p5.lock:
        p5.ring.clear()
        p5.ring_n = 0
    before = p5.underrun_frames
    p5._cb(outdata, frames, None, None)
    check("player: starved callback -> silence + underrun count",
          np.all(outdata == 0) and p5.underrun_frames == before + frames,
          f"(underruns {before} -> {p5.underrun_frames})")


if __name__ == "__main__":
    test_device_profiles()
    test_structs()
    test_tone_recovery()
    test_pll_offset_tracking()
    test_pll_low_freq_large_dev()
    test_multitone_linearity()
    test_de_emphasis()
    test_adjacent_channel()
    test_squelch_and_ctcss()
    test_activity()
    test_co_channel()
    test_adjacent_channel_selectivity()
    test_deviation_headroom()
    test_audio_player()
    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)} test(s): {FAILS}")
        sys.exit(1)
    print("ALL TESTS PASSED")
