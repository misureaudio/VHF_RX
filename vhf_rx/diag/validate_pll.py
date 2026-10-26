"""PLL demod validation (design §4.3 upgrade).

Checks, on the full 6.048 MSPS chain (DDC + decimation + demod + stages):
  1. tracking response: 300 / 1000 / 3000 Hz tones vs pll_tracking_response
  2. linearity at ±3/±5 kHz deviation: SND/ND, PLL vs discriminator
  3. carrier offset tracking: +150 Hz offset removed, tone intact
  4. loop settling: freq step 0 -> 1.8 kHz, settle time
  5. noise rejection: SND/ND on a locked carrier + noise, PLL vs disc

No hardware needed. Run from anywhere:
  D:/Source/hermes-dir/.venv/Scripts/python.exe diag/validate_pll.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from nbfm import (ChannelDemod, _pll_block, pll_tracking_response, snd_nd)  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "tests"))
from test_nbfm import fm_modulate, feed, FS  # noqa: E402

F_N = 1750.0
ZETA = 0.707


def tone_amp(audio, f0, fs=8000.0, win_s=1.0):
    a = audio[: int(win_s * fs)]
    w = np.hanning(a.size)
    Y = np.fft.rfft(a * w)
    f = np.fft.rfftfreq(a.size, 1 / fs)
    amp = 2 * np.abs(Y) / np.sum(w)
    k = int(round(f0 * a.size / fs))
    return float(amp[k - 1:k + 2].max())


print("=== 1. tracking response (1.8 kHz dev tone, PLL) ===")
d = ChannelDemod(145.0e6, FS)
for f0 in (300.0, 1000.0, 3000.0):
    n = int(1 * FS)
    t = np.arange(n) / FS
    a = feed(d, fm_modulate(1800.0 * np.sin(2 * np.pi * f0 * t), 145.0e6, FS))
    got = tone_amp(a, f0)
    # stages at f0 (16k domain: LP, DC; 8k domain: de-emphasis)
    w16 = 2j * np.pi * f0 / 16000
    h_lp = abs((d.lp_b[0] + d.lp_b[1] * np.exp(-w16) + d.lp_b[2] * np.exp(-2 * w16)) /
               (d.lp_a[0] + d.lp_a[1] * np.exp(-w16) + d.lp_a[2] * np.exp(-2 * w16)))
    h_dc = abs(d.dc_a * (1 - np.exp(-w16)) / (1 - d.dc_a * np.exp(-w16)))
    c = 2 * 75e-6 * 8000
    de_b = np.array([1.0, 1.0]) / (1.0 + c)
    de_a = np.array([1.0, (1.0 - c) / (1.0 + c)])
    w8 = 2j * np.pi * f0 / 8000
    h_de = abs((de_b[0] + de_b[1] * np.exp(-w8)) / (de_a[0] + de_a[1] * np.exp(-w8)))
    T = pll_tracking_response(f0, F_N, ZETA)
    exp = 1800.0 * T * h_lp * h_dc * h_de
    print(f"  {f0:5.0f} Hz: got {got:7.1f}  expect {exp:7.1f}  "
          f"({got / exp - 1:+.1%})  |T|={T:.3f}")
    d = ChannelDemod(145.0e6, FS)  # fresh PLL state per tone

print("=== 2. linearity: SND/ND at ±3/±5 kHz dev, PLL vs disc ===")
for dev in (3000.0, 5000.0):
    for dem in ("pll", "disc"):
        d = ChannelDemod(145.0e6, FS, demod=dem)
        n = int(2 * FS)
        t = np.arange(n) / FS
        a = feed(d, fm_modulate(dev * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS))
        r = snd_nd(a, 8000, 1000.0)
        print(f"  {dem:4s} @ ±{dev/1000:.0f} kHz: SND/ND = {r:5.1f} dB")

print("=== 3. carrier offset tracking (+150 Hz) ===")
# The PLL locks the offset into the NCO, so the tone must come through
# unchanged (vs the no-offset baseline) and the offset itself (a DC
# level in the deviation signal) must be rejected by the DC blocker.
d = ChannelDemod(145.0e6, FS)
n = int(1 * FS)
t = np.arange(n) / FS
base = feed(d, fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS))
amp_base = tone_amp(base, 1000.0)
d = ChannelDemod(145.0e6, FS)
a = feed(d, fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t), 145.0e6 + 150.0, FS))
got = tone_amp(a, 1000.0)
print(f"  1 kHz tone: no offset {amp_base:.1f} Hz, with +150 Hz offset "
      f"{got:.1f} Hz ({got / amp_base - 1:+.2%})")
print(f"  offset level in audio (DC blocker rejection): "
      f"{np.mean(a[4000:]):.2f} Hz (expect ~0)")

print("=== 4. loop settling (hard 500 Hz carrier step) ===")
# A constant carrier offset appears in the demodulated signal as a DC
# level (deviation in Hz), so the audio path shows: (4a) the loop's own
# transient — phi rings at f_n·sqrt(1-ζ²) ≈ 1.24 kHz while locking onto
# the step (measured on _pll_block alone, 16 kHz), then (4b) the DC
# blocker (20 Hz corner, τ ≈ 8 ms) rejecting the offset from the audio.
w_n = 2 * np.pi * F_N / 16000
Kp, Ki = 2 * ZETA * w_n, w_n ** 2
n16 = 2 * 16000
freq16 = np.zeros(n16)
freq16[n16 // 2:] = 500.0
x16 = np.exp(2j * np.pi * np.cumsum(freq16) / 16000).astype(np.complex64)
out, _, _, _ = _pll_block(x16, 0.0, 0.0, Kp, Ki, 16000.0 / (2 * np.pi), 0)
post = out[n16 // 2:]
final = float(np.mean(post[-1000:]))
over = np.where(np.abs(post - final) > 0.1 * abs(final))[0]
t4a = (over[-1] + 1) / 16000 if over.size else 0.0
print(f"  4a loop alone: phi settles to {final:.0f} Hz (expect ~500), "
      f"within ±10 % at {t4a * 1000:.2f} ms after step "
      f"(f_n = {F_N:.0f} Hz, ζ = {ZETA})")
d = ChannelDemod(145.0e6, FS)
n = int(2 * FS)
t = np.arange(n) / FS
audio_in = np.zeros(n)
audio_in[n // 2:] = 500.0
a = feed(d, fm_modulate(audio_in, 145.0e6, FS))
# 2 ms sliding DC estimate of the 8 kHz audio, from the nominal step
# (group delay ~0.9 ms); the offset step must decay to ~0
nominal = (n // 2) // 756
dc = np.array([np.mean(a[i:i + 16]) for i in range(nominal, a.size - 16)])
peak = float(np.max(np.abs(dc[:200])))
over = np.where(np.abs(dc) > 0.1 * peak)[0]
t4b = (over[-1] + 1) / 8000 if over.size else 0.0
print(f"  4b full chain: audio offset peak {peak:.0f} Hz, decays below "
      f"10 % at {t4b * 1000:.1f} ms (DC blocker τ ≈ 1/2π/20 Hz ≈ 8 ms)")

print("=== 5. PLL vs disc: SND/ND on a locked carrier + noise ===")
rng = np.random.default_rng(0)
n16 = 16000 * 10
t16 = np.arange(n16) / 16000
w_n = 2.0 * np.pi * F_N / 16000
Kp, Ki = 2 * ZETA * w_n, w_n ** 2
sig = np.exp(2j * np.pi * np.cumsum(1800.0 * np.sin(2 * np.pi * 1000 * t16)) / 16000)
sig_p = float(np.mean(np.abs(sig) ** 2))
for snr_db in (10.0, 20.0, 30.0):
    noise_p = sig_p / (10 ** (snr_db / 10))
    noise = (rng.standard_normal(n16) + 1j * rng.standard_normal(n16)) * np.sqrt(noise_p / 2)
    x = sig + noise
    out_pll, _, _, _ = _pll_block(x.astype(np.complex64), 0.0, 0.0, Kp, Ki, 16000.0 / (2 * np.pi), 0)
    xp = np.roll(x, 1)
    a_disc = np.imag(x * np.conj(xp)) / (np.abs(x) * np.abs(xp) + 1e-12) * 16000.0 / (2 * np.pi)
    r_pll = snd_nd(out_pll[16000:], 16000, 1000.0)
    r_disc = snd_nd(a_disc[16000:], 16000, 1000.0)
    print(f"  SNR {snr_db:3.0f} dB: PLL SND/ND = {r_pll:5.1f} dB   disc = {r_disc:5.1f} dB")
