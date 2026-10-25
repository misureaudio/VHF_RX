"""Stage-by-stage probe of the demod chain (diagnostic, not a test)."""
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


def fund_amp(x, f0, fs):
    """Amplitude of the component at f0 (Hanning, coherent-gain corrected)."""
    X = np.fft.rfft(x * np.hanning(x.size))
    f = np.fft.rfftfreq(x.size, 1 / fs)
    k = int(np.argmin(np.abs(f - f0)))
    return float(4.0 * np.abs(X[k]) / x.size)


def main():
    n = int(2 * FS)
    t = np.arange(n) / FS
    audio = 1800.0 * np.sin(2 * np.pi * 1000.0 * t)
    iq = fm_modulate(audio, 145.0e6, FS)
    d = ChannelDemod(145.0e6, FS)

    # replay process() but capture intermediate stages
    ph = 0.0
    step = d.step
    prev_x = 0.0 + 0.0j
    dev_all, xch_all = [], []
    for i in range(0, n, 16384):
        blk = iq[i:i + 16384]
        p = ph + step * np.arange(blk.size)
        x = blk * np.exp(-1j * p)
        ph = (ph + step * blk.size) % (2 * np.pi)
        x = _decimate(x, d.h1, 7)
        x = _decimate(x, d.h2, 9)
        x = _decimate(x, d.h3, 6)
        xch_all.append(x)
        xp = np.empty_like(x)
        xp[0] = prev_x
        xp[1:] = x[:-1]
        a = np.imag(x * np.conj(xp))
        a /= (np.abs(x) * np.abs(xp) + 1e-12)
        prev_x = x[-1]
        dev_all.append(a * d.fs_ch / (2 * np.pi))
    dev = np.concatenate(dev_all)
    print(f"discriminator fundamental @1kHz: {fund_amp(dev, 1000, 16000):.0f} Hz (expect ~1800, sin() nonlinearity -> ~1190)")

    # DC block
    yv = xv = 0.0
    a_dc = d.dc_a
    y = np.empty_like(dev)
    for i in range(dev.size):
        yv = a_dc * (dev[i] - xv + yv)
        y[i] = yv
        xv = dev[i]
    print(f"after DC block @1kHz:            {fund_amp(y, 1000, 16000):.0f} Hz")

    from scipy.signal import lfilter
    y2 = lfilter(d.lp_b, d.lp_a, y)
    print(f"after LP 3.4k @1kHz:             {fund_amp(y2, 1000, 16000):.0f} Hz")
    audio8 = y2[::2]
    a3 = lfilter(d.de_b, d.de_a, audio8)
    print(f"after de-emphasis @1kHz (8k):    {fund_amp(a3, 1000, 8000):.0f} Hz")
    print(f"de-emphasis |H|(1k) = {abs(1/(1+1j*2*np.pi*1000*75e-6)):.4f}")

    # stage-3 filter response around the 12.5 kHz interferer
    h3 = d.h3
    w = np.fft.rfftfreq(h3.size, 1 / (FS / 63))
    H = np.abs(np.fft.rfft(h3, 2 * h3.size)) / np.abs(np.fft.rfft(h3, 2 * h3.size)[0])
    for f in (4800, 7500, 8900, 10000, 11000, 12500, 14000):
        k = int(np.argmin(np.abs(w - f)))
        print(f"  |H3| @ {f/1000:5.1f} kHz = {20*np.log10(H[k]):8.1f} dB")

    # adjacent channel measurement
    want = fm_modulate(1800.0 * np.sin(2 * np.pi * 1000 * t), 145.0e6, FS, amp=0.5)
    intr = fm_modulate(1800.0 * np.sin(2 * np.pi * 2000 * t), 145.0125e6, FS, amp=0.5)
    out = []
    for i in range(0, n, 16384):
        a, _ = d.process((want + intr)[i:i + 16384])
        out.append(a)
    out = np.concatenate(out)
    f1 = fund_amp(out, 1000, 8000)
    f2 = fund_amp(out, 2000, 8000)
    print(f"adjacent: 1kHz={f1:.0f} Hz, 2kHz(intr)={f2:.0f} Hz, rej={20*np.log10(f1/f2):.1f} dB")


if __name__ == "__main__":
    main()
