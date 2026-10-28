"""Micro-probe: where does the FM phase go?"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from nbfm import ChannelDemod, _decimate  # noqa: E402

FS = 6.048e6


def fm_modulate(audio_hz, f_carrier, fs, amp=0.5):
    phi = 2.0 * np.pi * np.cumsum(audio_hz) / fs
    t = np.arange(audio_hz.size) / fs
    return amp * np.exp(1j * (phi + 2.0 * np.pi * f_carrier * t))


def top_peaks(x, fs, label, npeaks=8):
    # complex input: full (two-sided) FFT — rfft only accepts real input
    if np.iscomplexobj(x):
        X = np.fft.fft(x * np.hanning(x.size))
        f = np.fft.fftfreq(x.size, 1 / fs)
    else:
        X = np.fft.rfft(x * np.hanning(x.size))
        f = np.fft.rfftfreq(x.size, 1 / fs)
    p = np.abs(X).copy()
    print(f"  {label}: top peaks (Hz, rel dB):")
    for _ in range(npeaks):
        k = int(np.argmax(p))
        if p[k] <= 0:
            break
        print(f"    {f[k]:12.1f} Hz  {20*np.log10(p[k]/p.max()):7.1f} dB")
        p[max(0, k - 2):k + 3] = 0


n = 16384 * 8
t = np.arange(n) / FS
audio = 1800.0 * np.sin(2 * np.pi * 1000.0 * t)
iq = fm_modulate(audio, 145.0e6, FS)
print("iq: max|iq| =", np.max(np.abs(iq)), " rms =", np.sqrt(np.mean(np.abs(iq)**2)))

# DDC by direct multiplication with the *same* phase formula as fm_modulate
t2 = np.arange(n) / FS
x_direct = iq * np.exp(-1j * 2.0 * np.pi * 145.0e6 * t2)
print("direct DDC residual: phase step rms =",
      np.sqrt(np.mean(np.diff(np.angle(x_direct))**2)), "rad")
top_peaks(np.imag(x_direct * np.conj(np.roll(x_direct, 1))) /
          (np.abs(x_direct) * np.abs(np.roll(x_direct, 1)) + 1e-12), FS, "disc @6.048M")

# now via ChannelDemod's NCO
d = ChannelDemod(145.0e6, FS)
x_nco = iq * np.exp(-1j * (d.step * np.arange(n)))
print("NCO DDC residual: phase step rms =",
      np.sqrt(np.mean(np.diff(np.angle(x_nco))**2)), "rad")

# decimation
x = _decimate(x_nco, d.h1, 7)
print("after /7: rms =", np.sqrt(np.mean(np.abs(x)**2)))
x = _decimate(x, d.h2, 9)
print("after /9: rms =", np.sqrt(np.mean(np.abs(x)**2)))
x = _decimate(x, d.h3, 6)
print("after /6: rms =", np.sqrt(np.mean(np.abs(x)**2)), " size =", x.size)
top_peaks(x, FS / 378, "x @16k complex")

a = np.imag(x * np.conj(np.roll(x, 1))) / (np.abs(x) * np.abs(np.roll(x, 1)) + 1e-12)
print("disc @16k: rms =", np.sqrt(np.mean(a**2)), " max =", np.max(np.abs(a)))
top_peaks(a, FS / 378, "disc @16k")
