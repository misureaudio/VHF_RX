"""Numeric waterfall: max per-channel power over a whole capture window.
Reports every 12.5 kHz channel that exceeded the noise floor by >= 3 dB at
any point, plus how long it stayed above the open threshold (floor + 12 dB)."""
import sys
import threading
import time

import numpy as np

sys.path.insert(0, r"C:/Users/MATTIA/source/hermes-dir/vhf_rx_v2")
from nbfm import BandSpectrum
from sdrplay_api import AGC_50HZ, BW_5M, SDRPlay

FS = 6.048e6
DUR = 90.0
LO = 145.600e6

rx = SDRPlay()
dev = rx.find_device()
prof = rx.profile
rx.get_params()
rx.configure(fs_hz=FS, rf_hz=LO, bw_khz=BW_5M, lna=prof.lna_max_2m,
             gr_db=50, agc=AGC_50HZ, agc_setpoint_dbfs=-20,
             antenna=1, bias_t=False)

spec = BandSpectrum(FS, LO, 144.0e6, 146.0e6)
chans = spec.channels
max_db = np.full(chans.size, -1e9)
open_s = np.zeros(chans.size)
lock = threading.Lock()
q = []

def on_iq(xi, xq, n, reset):
    with lock:
        q.append((np.ctypeslib.as_array(xi, shape=(n,)).copy(),
                  np.ctypeslib.as_array(xq, shape=(n,)).copy()))

rx.start(on_iq, lambda *a: None)
t0 = time.monotonic()
scale = prof.sample_scale
while time.monotonic() - t0 < DUR:
    with lock:
        blocks, q[:] = q, []
    for i16, q16 in blocks:
        iq = ((i16.astype(np.float32) + 1j * q16.astype(np.float32))
              / scale).astype(np.complex64)
        cp, _, nf = spec.update(iq)
        if cp is None:
            continue
        db = 10.0 * np.log10(np.maximum(cp, 1e-30))
        dt = iq.size / FS
        np.maximum(max_db, db, out=max_db)
        open_s += (db >= nf + 12.0) * dt
rx.stop()

floor = spec.noise_db
order = np.argsort(max_db)[::-1]
print(f"LO {LO/1e6:.3f} MHz, {DUR:.0f} s, port B, noise floor (end) {floor:.1f} dBFS")
print("channels that exceeded floor + 3 dB at ANY point in the window:")
any_hit = False
for i in order[:30]:
    if max_db[i] < floor + 3.0:
        break
    any_hit = True
    f = chans[i]
    flag = "  <== Italian repeater output (145.575-145.7875)" \
        if 145.575e6 <= f <= 145.7875e6 else ""
    print(f"  {f/1e6:9.4f} MHz  max {max_db[i]:6.1f} dBFS "
          f"({max_db[i]-floor:+5.1f} dB)  open {open_s[i]:4.1f} s{flag}")
if not any_hit:
    print("  (none — the whole window was at noise-floor level)")
