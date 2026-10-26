# diag/ — internal diagnostics

Diagnostic scripts for the RSPdx NBFM receiver. These are **not** tests:
they are exploratory/measurement tools that print what they find and exit.
The conformance test suite lives in `tests/test_nbfm.py`.

Never delete these — each one documents a specific investigation and its
measured result. Run any of them from anywhere (paths are `__file__`-
relative):

```
D:/Source/hermes-dir/.venv/Scripts/python.exe diag/<name>.py
```

Scripts marked **[HW]** need the RSPdx streaming; the rest run offline.

## Signal classification

### probe_signal.py  [HW]
Live probe: what is on a given channel? Streams the band, demodulates ONE
channel, accumulates the 16 kHz complex channel signal + 8 kHz audio, then
classifies the signal (CW/RTTY beacon = one narrow line; wideband noise-like
= flat channel spectrum; voice FM = power in the audio band).

```
python diag/probe_signal.py --band 2m --ch 144.0e6 --duration 15
python diag/probe_signal.py --band 70cm --duration 15   # strongest channel
```

Used 2026-10-25 to classify the 144.000 MHz emission: a CW beacon (single
line, ~25 kHz deviation, −24 dBFS) — see design doc §8.3.

## Demod-chain probes (offline)

### probe_chain.py
Stage-by-stage probe of the demod chain: feeds a known FM tone through
DDC → decimation → demod → LP/DC/de-emphasis and prints the amplitude at
each stage, so a discrepancy can be localized to one stage.

### probe_micro.py
Micro-probe: "where does the FM phase go?" — dumps the top spectral peaks
of the channel signal at each chain stage. Used during the decimation
design to confirm the NCO/decimation behaviour.

## Cost profiling [HW] (capture 6 s once, replay offline)

### profile_replay.py  [HW]
Captures 6 s of live RSPdx IQ, then replays it through the full worker DSP
path (convert → concat → spectrum → demod) with per-component timing in
ms/s. The per-component numbers in design doc §8.1 come from this script.

### profile_pll_ab.py  [HW]
Back-to-back A/B on ONE captured stream: full worker DSP path with
`demod="disc"` vs `demod="pll"` (numba warmed up first). Isolates the
demod-stage cost difference from machine-state run-to-run variation.
Result (2026-10-25): PLL 252.8 ms/s vs discriminator 262.9 ms/s — the PLL
is ~4 % *cheaper*; machine state (not code) explains the ±35 % variation
seen between separate replay runs.

### profile_isolate.py  [HW]
Isolates the replay-structure cost: (A) full worker path with the
interleaved int16→complex64 conversion, (B) pending/concat per chunk only,
(C) prebuilt big array with a copy per chunk (view cost check). Used to
establish that the interleaved conversion is the dominant replay overhead
(~140 ms/s of the measured demod delta) and that the array view/copy cost
is negligible.

## PLL demod validation (offline)

### validate_pll.py
Full-chain PLL validation (design doc §4.3). Five checks (2026-10-25 run,
after the PLL startup fix — see design doc §4.3 "PLL startup"):
1. **Tracking response** — 300/1000/3000 Hz tones vs
   `nbfm.pll_tracking_response` × stage gains: 300 Hz −0.1 % (the startup-kick
   slip orbit used to cost −27.4 % here), 1000 Hz −0.6 % (the gated
   conformance point, §8.4), 3000 Hz −5.8 % (full-chain effect near the loop
   crossover — the loop alone matches theory exactly, see `isolate_pll.py`).
2. **Linearity** — SND/ND at ±3/±5 kHz deviation: PLL 85.8 dB vs
   discriminator 40.4 dB at ±3 kHz (the pre-fix 51.6 dB was the startup
   transient leaking through the 8 ms DC blocker); at the out-of-spec ±5 kHz
   the loop runs in a slip state that preserves the deviation (56.1 dB) — the
   discriminator is 30.1 dB there, which is why it stays as the selectable
   fallback for comparison.
3. **Carrier offset tracking** — a +150 Hz offset leaves the tone
   unchanged (−0.00 %) and is rejected from the audio (−0.00 Hz): the
   PLL locks the offset into the NCO and the DC blocker removes the
   residual.
4. **Loop settling** — a hard 500 Hz carrier step: the loop alone
   settles to within ±10 % in 0.19 ms (f_n = 1750 Hz, ζ = 0.707); in
   the full chain the offset decays out of the audio in 19.2 ms,
   matching the DC blocker time constant τ = 1/2π/20 Hz ≈ 8 ms
   (τ·ln 10 = 18.3 ms).
5. **Noise rejection** — SND/ND on a locked carrier + noise: PLL vs
   discriminator, 13.6/23.6/33.6 dB vs 15.3/25.5/34.1 dB at 10/20/30 dB
   input SNR (the PLL sits ~1–2 dB below at low SNR, converging at
   30 dB).

### derive_tf.py
Derives the exact discrete PLL transfer function
`|T_f(e^{jw})| = |[Kp(1-u)+Ki] / [(1-u)^2 + u(Kp(1-u)+Ki)]|` (u = e^{−jw})
from the loop equations and compares it against the continuous-time
(Kp·s+Ki)/(s²+Kp·s+Ki) and the measured PLL output. This is the derivation
behind `nbfm.pll_tracking_response` — the continuous-time form is wrong
above the loop bandwidth (0.845 vs measured 1.44 at 3 kHz) and the exact
form matches to 0.1 %.

### isolate_pll.py
Feeds a pure tone directly into `nbfm._pll_block` at 16 kHz (no DDC, no
stages) and measures output amplitude vs input deviation. Localizes demod
accuracy to the PLL itself: the PLL alone matches theory exactly
(1.028/1.220/1.533 at 300/1000/3000 Hz for a 1.8 kHz dev tone). After the
startup fix (design doc §4.3 "PLL startup") the full chain matches the model
at 300/1000 Hz too (−0.1 % / −0.6 %); the remaining −5.8 % at 3000 Hz is the
full-chain effect near the loop crossover (the 6M NCO zero-order hold adds a
small feedback delay).
