"""Isolate the demod timing difference between profile_replay's
pending/concat loop and the prebuilt-array loop. Same capture, three
replay variants, PLL demod each:
  A) prebuilt big array, view per chunk            (profile_pll_ab style)
  B) pending/concat per chunk                      (profile_replay style)
  C) prebuilt big array, COPY per chunk            (view cost check)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from nbfm import BandSpectrum, ChannelDemod, Squelch  # noqa: E402
from sdrplay_api import (AGC_50HZ, BW_5M,  # noqa: E402
                         RSPDX_LNA_STATES_VHF_BAND3, SDRPlay)

FS = 6.048e6
I14_SCALE = 8192.0
CHUNK = 1 << 16
DUR = 6.0


def main():
    rx = SDRPlay()
    dev = rx.find_device()
    rx.get_params()
    rx.configure(fs_hz=FS, rf_hz=145.0e6, bw_khz=BW_5M,
                 lna=RSPDX_LNA_STATES_VHF_BAND3 - 1, gr_db=50,
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
    stream_s = n_samples / FS
    print(f"captured {len(cap)} blocks, {n_samples} samples")

    iq_all = []
    for i16, q16 in cap:
        iq_all.append(((i16.astype(np.float32) + 1j * q16.astype(np.float32))
                       / I14_SCALE).astype(np.complex64))
    big = np.concatenate(iq_all)

    ChannelDemod(145.5e6, FS, demod="pll").process(big[:CHUNK])  # warm

    # A) prebuilt view
    spec = BandSpectrum(FS, 145.0e6, 144.0e6, 146.0e6)
    demod = ChannelDemod(145.5e6, FS, demod="pll")
    sq = Squelch()
    t_spec = t_demod = 0.0
    for k in range(0, big.size - CHUNK, CHUNK):
        c = big[k:k + CHUNK]
        tc = time.perf_counter()
        spec.update(c)
        t_spec += time.perf_counter() - tc
        tc = time.perf_counter()
        audio, m = demod.process(c)
        sq.update(m["carrier_db"], audio.size * 2)
        t_demod += time.perf_counter() - tc
    print(f"A) prebuilt view   : spec {t_spec / stream_s * 1000:6.1f} ms/s | "
          f"demod {t_demod / stream_s * 1000:6.1f} ms/s")

    # B) pending/concat (profile_replay structure)
    spec = BandSpectrum(FS, 145.0e6, 144.0e6, 146.0e6)
    demod = ChannelDemod(145.5e6, FS, demod="pll")
    sq = Squelch()
    t_concat = t_spec = t_demod = 0.0
    pending = []
    pending_n = 0
    for iq in iq_all:
        pending.append(iq)
        pending_n += iq.size
        while pending_n >= CHUNK:
            tc = time.perf_counter()
            allc = np.concatenate(pending)
            t_concat += time.perf_counter() - tc
            tc = time.perf_counter()
            spec.update(allc[:CHUNK])
            t_spec += time.perf_counter() - tc
            tc = time.perf_counter()
            audio, m = demod.process(allc[:CHUNK])
            sq.update(m["carrier_db"], audio.size * 2)
            t_demod += time.perf_counter() - tc
            rest = allc[CHUNK:]
            pending = [rest] if rest.size else []
            pending_n = rest.size
    print(f"B) concat per chunk: concat {t_concat / stream_s * 1000:5.1f} | "
          f"spec {t_spec / stream_s * 1000:6.1f} ms/s | "
          f"demod {t_demod / stream_s * 1000:6.1f} ms/s")

    # C) prebuilt, copy per chunk
    spec = BandSpectrum(FS, 145.0e6, 144.0e6, 146.0e6)
    demod = ChannelDemod(145.5e6, FS, demod="pll")
    sq = Squelch()
    t_copy = t_spec = t_demod = 0.0
    for k in range(0, big.size - CHUNK, CHUNK):
        tc = time.perf_counter()
        c = np.array(big[k:k + CHUNK], copy=True)
        t_copy += time.perf_counter() - tc
        tc = time.perf_counter()
        spec.update(c)
        t_spec += time.perf_counter() - tc
        tc = time.perf_counter()
        audio, m = demod.process(c)
        sq.update(m["carrier_db"], audio.size * 2)
        t_demod += time.perf_counter() - tc
    print(f"C) prebuilt copy   : copy {t_copy / stream_s * 1000:5.1f} | "
          f"spec {t_spec / stream_s * 1000:6.1f} ms/s | "
          f"demod {t_demod / stream_s * 1000:6.1f} ms/s")


if __name__ == "__main__":
    main()
