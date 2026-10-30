"""Per-second spectral lines of the wideband WAV (diag helper)."""
import wave
import numpy as np

w = wave.open("D:/Source/hermes-dir/output/audio_103.0000.wav")
fs = w.getframerate()
n = w.getnframes()
x = np.frombuffer(w.readframes(n), dtype=np.int16).astype(np.float64)
print("rate", fs, "dur", round(n / fs, 2), "s")
win = fs
for i in range(1, min(8, n // win - 1)):
    seg = x[i * win:(i + 1) * win]
    wgt = np.hanning(win)
    X = np.abs(np.fft.rfft(seg * wgt))
    f = np.fft.rfftfreq(win, 1 / fs)
    amp = 2 * X / np.sum(wgt)
    top = np.argsort(amp)[::-1][:6]
    print(f"t={i}s rms={np.sqrt(np.mean(seg**2)):.0f} top: "
          + ", ".join(f"{f[j]:.0f}Hz({amp[j]:.0f})" for j in top))
# full-file top 12 lines (Hann, whole file)
wgt = np.hanning(n)
X = np.abs(np.fft.rfft(x * wgt))
f = np.fft.rfftfreq(n, 1 / fs)
amp = 2 * X / np.sum(wgt)
top = np.argsort(amp)[::-1][:12]
print("full-file top lines:")
for j in top:
    print(f"  {f[j]:9.1f} Hz  amp {amp[j]:8.0f}  "
          f"{20*np.log10(amp[j]/amp[top[0]]):+5.1f} dB")
