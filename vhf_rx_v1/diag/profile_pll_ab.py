"""A/B: discriminator vs PLL demod cost, same captured stream, back-to-back.
Captures 6 s of live RSPdx IQ once, then replays it through the full worker
DSP path with demod="disc" and demod="pll" (numba warmed up first).
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nbfm import BandSpectrum, ChannelDemod, Squelch  # noqa: E402
from sdrplay_api import AGC_50HZ, BW_5M, SDRPlay  # noqa: E402

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

    rx.start(on_iq, lambda *a: None)
    t0 = time.monotonic()
    while time.monotonic() - t0 < DUR:
        time.sleep(0.2)
    rx.stop()
    with lock:
        cap = list(cap)
    n_samples = sum(i.size for i, q in cap)
    print(f"captured {len(cap)} blocks, {n_samples} samples "
          f"({n_samples / DUR / 1e6:.3f} MSPS)")

    # pre-convert once (shared by both replays)
    iq_all = []
    for i16, q16 in cap:
        iq_all.append(((i16.astype(np.float32) + 1j * q16.astype(np.float32))
                       / prof.sample_scale).astype(np.complex64))
    big = np.concatenate(iq_all)

    # warm up numba PLL compilation (not timed)
    ChannelDemod(145.5e6, FS, demod="pll").process(big[:CHUNK])

    stream_s = n_samples / FS
    for dem in ("disc", "pll"):
        spec = BandSpectrum(FS, 145.0e6, 144.0e6, 146.0e6)
        demod = ChannelDemod(145.5e6, FS, demod=dem)
        sq = Squelch()
        t_spec = t_demod = 0.0
        t0 = time.perf_counter()
        for k in range(0, big.size - CHUNK, CHUNK):
            c = big[k:k + CHUNK]
            tc = time.perf_counter()
            spec.update(c)
            t_spec += time.perf_counter() - tc
            tc = time.perf_counter()
            audio, m = demod.process(c)
            sq.update(m["carrier_db"], audio.size * 2)
            t_demod += time.perf_counter() - tc
        t_replay = time.perf_counter() - t0
        print(f"[{dem}] replay {t_replay:.2f} s CPU for {stream_s:.2f} s "
              f"-> {t_replay / stream_s * 100:.0f}% of real time | "
              f"spectrum {t_spec / stream_s * 1000:6.1f} ms/s | "
              f"demod {t_demod / stream_s * 1000:6.1f} ms/s")


if __name__ == "__main__":
    main()
