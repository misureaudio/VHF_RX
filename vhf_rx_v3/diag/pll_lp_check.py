"""Verify the wideband PLL low-passes audio: synthesize FM tones, demod.

Convention matches test_nbfm.py: carrier at absolute f0 in the IQ,
ChannelDemod(f0, FS) with rf_center=0 (default).
"""
import sys

import numpy as np

sys.path.insert(0, "D:/Source/hermes-dir/vhf_rx_v3")

from nbfm import ChannelDemod

FS = 6.048e6


def fm_modulate(audio_hz, f_carrier, fs, amp=0.5):
    t = np.arange(audio_hz.size) / fs
    phase = 2.0 * np.pi * f_carrier * t + 2.0 * np.pi * np.cumsum(audio_hz) / fs
    return amp * np.exp(1j * phase).astype(np.complex64)


def tone_out(d, audio_in, fs_in=48000.0):
    iq = fm_modulate(audio_in, 98.0e6, FS)
    out = np.concatenate(
        [d.process(iq[i:i + (1 << 16)])[0] for i in range(0, iq.size, 1 << 16)]
    )
    # skip 0.5 s of PLL startup, measure fundamental amplitude
    out = out[int(0.5 * fs_in):]
    n = out.size
    X = np.abs(np.fft.rfft(out * np.hanning(n)))
    f = np.fft.rfftfreq(n, 1 / fs_in)
    amp = 2 * X / np.sum(np.hanning(n))
    return amp, f


for f_tone, dev in [(1000.0, 5000.0), (5000.0, 20000.0), (10000.0, 20000.0)]:
    t = np.arange(int(2 * FS)) / FS
    audio_in = dev * np.sin(2 * np.pi * f_tone * t)
    d = ChannelDemod(98.0e6, FS, fs_audio=48000, wideband=True,
                     bp_pass_hz=62000.0, bp_stop_hz=124000.0, bp_taps=128,
                     tau=50e-6, pll_f_n=25000.0)
    amp, f = tone_out(d, audio_in)
    k = int(np.argmin(np.abs(f - f_tone)))
    got = float(amp[k])
    # de-emphasis gain at f_tone (the output is de-emphasized):
    tau = 50e-6
    de = 1.0 / np.sqrt(1.0 + (2 * np.pi * f_tone * tau) ** 2)
    print(f"tone {f_tone:6.0f} Hz, dev {dev/1000:4.0f} kHz -> "
          f"out {got:8.1f} Hz (input x de-emph = {dev*de:8.1f} Hz, "
          f"ratio {got/(dev*de):.3f})")
