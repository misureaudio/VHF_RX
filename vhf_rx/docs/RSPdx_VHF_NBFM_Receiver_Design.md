# VHF/UHF Narrowband FM Receiver on SDRPlay RSPdx — Design

**Scope:** amateur-radio receiver, bands 144–146 MHz (2 m) and 432–436 MHz (70 cm),
modulation: narrowband FM (NBFM) voice, EU standard practice.
Front end: SDRPlay RSPdx (receive-only SDR, ZIF architecture).

All regulatory/standards numbers below were verified against the primary documents
listed in §8 (downloaded and text-extracted 2026-10-24). Numbers marked *(est.)* are
engineering estimates derived from cited data, not measurements.

---

## 1. Requirements

| Item | Value |
|---|---|
| Bands | 144–146 MHz, 432–436 MHz |
| Mode | NBFM voice (analog), EU standards |
| Channel grid | 12.5 kHz (IARU R1) |
| Deviation to receive | up to ±5 kHz (worst case; standard is ±3 kHz) |
| Audio | 300–3000 Hz, 6 dB/oct de-emphasis |
| Squelch | carrier detect + optional CTCSS (and 1750 Hz tone-burst detection) |
| Throughput | whole 2 m band (2 MHz) or whole 70 cm window (4 MHz) in one tuning |
| Optional | dual-band simultaneous (needs 2× RSPdx — single tuner per unit) |

---

## 2. Standards basis (EU)

### 2.1 IARU Region 1 band plans (official, effective Dec 2020, VGC Novi Sad)

**144–146 MHz (2 m)** — FM-relevant segments:

| Segment (MHz) | Max BW | Usage |
|---|---|---|
| 144.794 – 144.9625 | 12 kHz | MGM/digital — **144.800 APRS**, DV internet-voice gateways 144.8125–144.8625 |
| 144.975 – 145.194 | 12 kHz | **FM repeater input exclusive** |
| 145.194 – 145.206 | 12 kHz | Space communication |
| 145.206 – 145.5625 | 12 kHz | FM simplex — **145.500 FM calling** |
| 145.575 – 145.7935 | 12 kHz | **FM repeater output exclusive** |
| 145.794 – 145.806 | 12 kHz | Space communication |
| 145.806 – 146.000 | 12 kHz | Satellite exclusive |

Repeater offset on 2 m: **+600 kHz** (input 144.975–145.194 ↔ output 145.575–145.7935).
Band-plan footnote (De Haan 1993, re-affirmed Tel Aviv 1996): Region 1 uses a **genuine
12.5 kHz channel spacing system** for FM repeater and simplex operation.

**430–440 MHz (70 cm)** — the 432–436 MHz window contains:

| Segment (MHz) | Max BW | Usage |
|---|---|---|
| 432.500 – 432.975 | 12 kHz | All modes — **432.500 new APRS frequency**; 432.600–432.975 = **repeater input, Region 1 standard, 25/12.5 kHz spacing, 2.0 MHz shift** |
| 433.000 – 433.375 | 12 kHz | **FM/DV repeater input, Region 1 standard, 25/12.5 kHz spacing, 1.6 MHz shift** |
| 433.400 – 433.575 | 12 kHz | FM/DV simplex — 433.400 SSTV (FM/AFSK), 433.450 DV calling, **433.500 FM calling** |
| 433.600 – 434.000 | none | All modes — 433.625–433.775 digital channels, 434.000 centre of digital experiments |
| 434.000 – 434.594 | 12 kHz | All modes / ATV — 434.450–434.575 digital channels |
| 434.594 – 434.981 | 12 kHz | All modes — **434.600–434.9875 repeater output, 12.5 kHz spacing, 1.6 or 2.0 MHz shift** |
| 435.000 – 436.000 | none | Satellite service |

So the 70 cm window carries two coexisting repeater systems (1.6 MHz and 2.0 MHz shift,
both landing on outputs 434.600–434.9875) plus simplex and digital. National
variations exist below 432 MHz (430.000–431.975 is sub-regional/national planning) —
check the national allocation for the operating country.

### 2.2 IARU R1 NBFM technical standards (VHF Handbook 10.03, 2026 ed.)

**Recommendation FM.1** (Brussels 1969, Scheveningen 1972) — *NBFM technical standard*:
- maximum modulation index **1**, audio band restricted to **3 kHz**
- maximum deviation **±3 kHz**, emission designator **12K0F3E**
- AF response **300–3000 Hz**, 12 dB/octave roll-off outside
- pre-emphasis **+6 dB/octave** (TX), de-emphasis **−6 dB/octave** (RX)

**Recommendation FM.24** — *FM repeaters, 145 MHz and 435 MHz bands*:
- vertical polarisation
- access by **CTCSS or 1750 Hz tone burst** (CTCSS preferred; since 1 Jan 2015 new FM
  repeaters shall use CTCSS on receivers as well as transmitters)
- max ERP 15 W; deviation ±3 kHz; AF 300–3000 Hz; pre/de-emphasis ±6 dB/oct
- CTCSS tone table adopted (frequencies accurate **±1 %**)

### 2.3 ETSI EN 300 086 V2.1.2 (2016-08) — EU equipment standard

"Land Mobile Service; radio equipment … 30 MHz–1 GHz; channel separations
12.5/20/25 kHz; analogue speech." The modulation limits that matter for a receiver
(Table 3, clause 7.4.3.1 — *Maximum Permissible Frequency Deviation*):

| Channel separation | MPFD |
|---|---|
| 12.5 kHz | **±2.5 kHz** |
| 20 kHz | ±4.0 kHz |
| 25 kHz | **±5.0 kHz** |

Other relevant clauses:
- normal test modulation: 1 kHz tone at **60 % of MPFD** (clause 6.1)
- audio: up to 3.0 kHz (2.55 kHz for 12.5 kHz spacing), −14 dB/oct above (Fig. 6)
- co-channel rejection (clause 8.3.3): −12…0 dB allowed at ±1.5/±3 kHz offsets for
  12.5 kHz spacing (i.e. a 12.5 kHz-channel receiver must tolerate a co-channel
  interferer up to 12 dB above the wanted signal)

### 2.4 ETSI TS 103 236 V1.1.1 (2014-01) — CTCSS / DCS

Applies to FM voice equipment per EN 300 086 / EN 300 296. 38 standard tones,
**67.0 … 250.3 Hz** (67.0, 69.3, 71.9, 74.4, 77.0, 79.7, 82.5, 85.4, 88.5, 91.5,
94.8, 97.4, 100.0, 103.5, 107.2, 110.9, 114.8, 118.8, 123.0, 127.3, 131.8, 136.5,
141.3, 146.2, 151.4, 156.7, 162.2, 167.9, 173.8, 179.9, 186.2, 192.8, 203.5, 210.7,
218.1, 225.7, 233.6, 241.8, 250.3). Recommended test set: **67 / 123 / 250.3 Hz**.
Tone accuracy ±1 % (per IARU FM.24).

### 2.5 Design implications (parameter table)

| Parameter | Standard value | Receiver design value |
|---|---|---|
| Channel spacing | 12.5 kHz | 12.5 kHz grid |
| Deviation | ±3 kHz (IARU); ±2.5 kHz (ETSI 12.5 kHz); ±5 kHz (ETSI 25 kHz, some TXs) | demod headroom for **±5 kHz** |
| Carson bandwidth (±5 kHz, 3 kHz audio) | 2(Δf+fa) = 16 kHz | channel filter passband **±8 kHz** (16 kHz) |
| Audio band | 300–3000 Hz | LPF 3.4 kHz @ 8 kHz audio rate |
| Pre/de-emphasis | 6 dB/oct | 1-pole, τ = **75 µs** (gives +6 dB/oct) |
| Squelch | CTCSS / 1750 Hz burst / carrier | carrier detect + CTCSS Goertzel bank |
| Emission | 12K0F3E | — |

---

## 3. RSPdx as the RF front end

### 3.1 Verified specifications (datasheet V1.0 + API v3 spec)

- 1 kHz–2 GHz, **ZIF** architecture (all IF bandwidths in ZIF; low-IF ≤1.536 MHz)
- IF bandwidths (3 dB): 200 / 300 / 600 kHz, 1.536 / 5 / 6 / 7 / 8 MHz
- Sample rate 2–10.66 MSPS: **14-bit ≤ 6.048 MSPS**, 12-bit ≤ 8.064, 10-bit ≤ 9.216, 8-bit above
- Antenna ports: **A** SMA 1 kHz–2 GHz 50 Ω; **B** SMA + 4.7 V bias-T (100 mA); C BNC ≤200 MHz (Hi-Z)
- Max input: **0 dBm continuous**, 10 dBm short periods
- Noise figure: 2.6 dB @100 MHz, 2.1 dB @200 MHz, 6.0 dB @340 MHz, 3.1 dB @660 MHz
  → *(est.)* ≈2.4 dB @145 MHz, ≈5 dB @435 MHz
- Internal preselector covers both bands (BP 120–250 MHz; BP 420–1000 MHz)
- **24 MHz external reference input** (GPSDO-compatible, auto-detect)
- Software AGC (50/100/500 Hz time constants) + LNA state + gain reduction (20–59 dB)
- API: **SDRplay Hardware API v3** (≥3.06 for RSPdx) — service/daemon + client library;
  IQ delivered **callback-only, int16 I/Q**. Measured delivery (RSPdx, 6.048 MSPS,
  API v3.11): average block ≈ **1 344 samples**, delivered in **bursts** (USB
  batching) — a queue of fixed depth overflows; an unbounded queue is required
  (§9.1). Measured: 6.048 MSPS sustained; noise floor −46.3 dBFS at LNA 8 +
  gain 40 (total gain 39.98 dB)

### 3.2 Antenna / front-end configuration

- 2 m: ½λ dipole or ground plane on **port A** (SMA). Port B usable if an active
  (bias-fed) indoor antenna is preferred (4.7 V / 100 mA bias-T).
- 70 cm: ¼λ ground plane or ½λ dipole on port A.
- Strong-signal environments (urban, nearby repeaters, 433 MHz ISM): put a **6 dB
  pad** in front of the RSPdx — the 0 dBm max-input limit and the internal LNA make
  this the standard linearity fix; AGC then rides on the pad.
- **433 MHz ISM coexistence (70 cm only):** EU licence-free EN 300 220-1 devices occupy
  433.05–434.79 MHz — *inside* the 432–436 window. Keyfobs/alarms/telemetry can be
  strong and will compress the front end. A bandpass filter cannot remove them;
  mitigations are the pad, RSPdx AGC, and channel-level FM capture (the digital
  channelizer + demod reject them once a channel is selected).
- **Reference:** for repeater monitoring and anything beyond casual listening, feed the
  24 MHz input from a **GPSDO**. Internal crystal ≈1 ppm → up to ~430 Hz error at
  434 MHz (marginal for a 12.5 kHz grid); GPSDO gives ~10⁻¹⁰ stability.

### 3.3 Sensitivity (estimate from datasheet NF, 12.5 kHz channel)

Thermal noise in 12.5 kHz: −174 dBm/Hz + 10·log10(12 500) = **−133.0 dBm**.

| Band | NF (datasheet, interp.) | System noise floor | 12 dB SINAD | 20 dB SNR |
|---|---|---|---|---|
| 2 m | ≈2.4 dB | ≈ −130.6 dBm | ≈ −118.6 dBm | ≈ −110.6 dBm |
| 70 cm | ≈5 dB | ≈ −128.0 dBm | ≈ −116.0 dBm | ≈ −108.0 dBm |

*(est.)* — consistent with typical handheld performance (−110…−120 dBm usable).

---

## 4. System architecture

### 4.1 Topology

The RSPdx is a **single-tuner** receiver: one band at a time.

- **Single unit:** tune 145.000 MHz (2 m) *or* 434.000 MHz (70 cm). Band switch =
  one `sdrplay_api_Update(dev, tuner, Update_Tuner_Frf, …)` call.
- **Dual-band simultaneous:** two RSPdx units (two USB ports, two API handles); each
  runs its own DSP chain. The DSP below is per-band, so the two instances are
  independent.

### 4.2 Capture configuration

| | 2 m | 70 cm |
|---|---|---|
| LO (rfFreq) | 145.000 MHz | 434.000 MHz |
| IF type | ZIF | ZIF |
| IF BW (3 dB) | 5 MHz | 8 MHz |
| fs | 6.048 MSPS (14-bit) | 6.048 MSPS (14-bit) |
| Span | 141.976–148.024 MHz | 430.976–437.024 MHz |
| Captured band | 144–146 ✓ (+2 MHz margin/side) | 432–436 ✓ (+1 MHz margin/side) |

Wider margin if desired: 12-bit @ 8.064 MSPS (±4.032 MHz) — at the cost of 2 ADC bits.
Gain: LNA at max state for the band, AGC 50 Hz enabled, set point ≈ −20 dBFS;
acknowledge `PowerOverloadChange` events (`Update_Ctrl_OverloadMsgAck`) and reduce
gain via `Update_Tuner_Gr` on overload.

### 4.3 Digital chain

```
IQ 6.048 MSPS int16
  │
  ├─(1) band-limit FIR:  2 m: pass 143.9–146.1 MHz   70 cm: pass 431.9–436.1 MHz
  │
  ├─(2) activity monitor:  8192-pt FFT every 1.35 ms
  │        → per-bin power on 12.5 kHz grid → channel activity map (hysteresis)
  │        (drives scanner UI + selects which channels to demodulate)
  │
  └─(3) per ACTIVE channel (N = 1…~20 in practice):
        NCO @ f_ch (complex DDC)
          → FIR decimate /7  → 864 kHz   (~16 taps)
          → FIR decimate /9  → 96 kHz    (~64 taps)
          → FIR channel BPF, pass ±8 kHz, stop ±12.5 kHz, /6 → 16 kHz  (~128 taps)
          → FM demod (2nd-order digital PLL, loop BW 2–3 kHz) → audio @16 kHz
          → de-emphasis 1-pole τ=75 µs
          → FIR LPF 3.4 kHz, /2 → 8 kHz audio (float32)
          → squelch (carrier detect + CTCSS Goertzel bank)
          → audio sink / logging
```

Decimation factors 7·9·6 = 378 = 6048/16 — all integer, all FIR-friendly
(6048 = 2⁵·3³·7).

**Demod reach (measured 2026-10-25):** the /7 stage's anti-alias FIR cuts at
0.40·fs/7 = **345.6 kHz**, so a channel demodulates correctly only within
**≈ ±300 kHz of the LO** (its Carson band must stay inside the passband).
Live-verified: a real +18 dB wideband emission at 144.000 MHz demodulates to
pure noise floor with the LO at 145.000 (−1.0 MHz), while the full-band
spectrum sees it fine (LO-shift test: the excess stays at 144.000 absolute
when the LO moves to 144.500). Consequence: the *full-band monitor*
(waterfall + channel grid) covers the whole IF span, but per-channel demod is
limited to the ±300 kHz window around the LO; the prototype auto-selects only
within that window and warns on explicit `--ch` outside it. Whole-band demod
in one tuning (req. §1) would need a channelized front end (polyphase filter
bank) or LO retuning on channel select — a later step (§5.2).

**Cost per channel** *(est., 2×FLOP/MAC):* ≈ 6.048M·16·2 + 864k·64·2 + 96k·128·2
≈ **0.33 GFLOP/s** → 10 active channels ≈ 3.3 GFLOP/s. Fine in C/Numba; NumPy
prototype comfortably handles 2–4 channels. The activity monitor itself (one 8192-pt
FFT / 1.35 ms) is negligible.

**Demod: 2nd-order digital PLL (implemented 2026-10-25).** With 12.5 kHz spacing the
adjacent channel is only 12.5 kHz away and EN 300 086 requires tolerance of a
co-channel interferer up to 12 dB above wanted. A 2nd-order PLL (loop BW 2–3 kHz)
integrates out wideband adjacent noise, tracks DC/IQ offset, and is *linear* in
the deviation — no sin(Δφ) folding: SND/ND 85.6 dB at ±3 kHz vs 40.4 dB for the
quadrature discriminator (kept as selectable fallback, `ChannelDemod(demod="disc")`).
Per-sample structure: wrapped phase error `e[n] = ∠x[n] − θ[n]` (arctan2 PD —
amplitude-invariant, so AGC steps don't modulate the loop gain); loop filter
`integ[n] = integ[n−1] + Ki·e[n]`, `φ[n] = Kp·e[n] + integ[n]`; NCO
`θ[n+1] = θ[n] + φ[n]`; **audio = `φ[n]·fs/2π`** (the NCO frequency — using the
integrator alone low-passes the audio and loses ~30 % of the tone). Gains:
`w_n = 2π·f_n/16000`, `Kp = 2ζw_n`, `Ki = w_n²` with `f_n = 1750 Hz`, `ζ = 0.707`
→ crossover ≈ 2.4 kHz (Jury-stable: 2Kp + Ki < 4). The loop's tracking response
`|T_f|` has a ~+1.7 dB bump near `f_n` (exact discrete form in
`nbfm.pll_tracking_response`; the continuous-time `(Kp·s+Ki)/(s²+Kp·s+Ki)` is
wrong near/above the loop BW). The PLL's PD wraps at π, so the same
deviation-aliasing limit applies (|dev| < 8 kHz at 16 kHz). At the out-of-spec
±5 kHz deviation the startup frequency error exceeds the pull-in range and the
loop runs in a slip state that preserves the deviation (audio = deviation + a DC
offset removed by the DC blocker): measured fundamental −0.8 % vs linear,
SND/ND 56.1 dB (discriminator: 30.1 dB).
Settling (measured, `diag/validate_pll.py`): a hard 500 Hz carrier step
settles within ±10 % in 0.19 ms on the loop alone; in the audio path the
offset decays out in 19.2 ms, set by the DC blocker time constant
(τ = 1/2π/20 Hz ≈ 8 ms; τ·ln 10 = 18.3 ms). A constant offset produces no
steady-state audio: the PLL locks it into the NCO and the DC blocker
removes the residual (+150 Hz offset → tone amplitude change −0.00 %,
offset level in audio −0.00 Hz).
The lock indicator (PD-error energy) is available as a free carrier-detect for
the squelch; the implemented carrier detect is the power-ratio metric of the
16 kHz baseband.

**PLL startup (measured 2026-10-25):** the decimation FIRs run from zero state,
so the first ~10 16 kHz samples of a stream are near-zero partial sums whose
phase is numerically garbage. Starting the loop there injects a one-shot phase
kick (up to ~π in the PD error); for a 300 Hz tone at 1.8 kHz deviation
(β = 6, the hardest low-frequency voice case) the kick pushed the loop out of
the capture basin into a *stable slip orbit* — the tone then came out at 72.6 %
of the tracking-response model with its harmonics at −6 dB, and the loop stayed
slipped for the entire stream (verified sample-by-sample: the slip orbit is
attractive, not a transient). The loop therefore (a) amplitude-gates the PD
(`|x| < 0.05`: hold the loop state, output 0 — also holds the loop steady
through deep fades in the live stream) and (b) initializes the NCO to the first
valid sample's phase (e₀ = 0, identical to a pre-locked start). After the fix:
the 300 Hz β = 6 tone matches the tracking-response model to −0.1 % (1820 vs
1821 Hz) with its harmonics at −119 dB; the PD error stays bounded (±0.18 rad)
and the NCO frequency tracks the input frequency to ≤ 209 Hz. The startup kick
was also corrupting the SND/ND measurements through the 8 ms DC blocker:
±3 kHz 48.5 → 85.6 dB, ±5 kHz 7.7 → 56.1 dB after the fix. Regression test:
`tests/test_nbfm.py::test_pll_low_freq_large_dev`.

**De-emphasis:** standard 6 dB/oct FM pre-emphasis ↔ τ = 75 µs. At 8 kHz audio rate:
`y[n] = α·(x[n] − x[n−1]) + (1−α)·y[n−1]`, `α = 1 − e^(−8000·75e−6) ≈ 0.811`.

**Squelch:**
1. *Carrier detect:* PLL lock indicator (or power in a 200 Hz band at DC vs. a rolling
   12.5 kHz noise-floor estimate). Open: +10 dB over floor; close: 1 s below
   threshold, 3 dB hysteresis.
2. *CTCSS (optional per channel):* Goertzel bank over the 38 TS 103 236 tones on the
   demodulated 8 kHz audio (or on baseband before demod — either works; the tones are
   sub-audible 67–250.3 Hz, ±1 % tolerance). Open after 100 ms continuous presence,
   close after 2 s absence.
3. *1750 Hz tone burst:* detect 1750 Hz ±25 Hz for ≥100 ms in the audio path (some
   repeaters use burst instead of CTCSS).

### 4.4 Channel tables (12.5 kHz grid)

**2 m (≈71 channels):**
- Repeater input: 144.975 → 145.1875 (18 ch)
- Simplex: 145.200 → 145.550 (29 ch) — 145.500 calling
- Repeater output: 145.575 → 145.7875 (18 ch)
- APRS 144.800; DV gateways 144.8125–144.8625 (5 ch)

**70 cm window 432–436 (≈87 channels):**
- Repeater input (2.0 MHz shift): 432.600 → 432.975 (31 ch)
- Repeater input (1.6 MHz shift): 433.000 → 433.375 (16 ch)
- Simplex: 433.400 → 433.575 (8 ch) — 433.500 calling
- Repeater output: 434.600 → 434.975 (32 ch)
- (digital/ATV segments 433.625–433.775, 434.000, 434.450–434.575 pass through the
  band-limit filter but need no FM demod unless you add data decoders)

### 4.5 Audio output (local sound device)

The prototype plays demodulated audio on the computer's default output device
(`--play`) via sounddevice/PortAudio (WASAPI shared mode on Windows); `--vol`
scales the level, `--device` selects the output. Pipeline (`playback.py`,
implemented 2026-10-25):

- **Rate:** 8 kHz channel audio → **48 kHz mono float32** (WASAPI-native; 6×
  integer ratio, `resample_poly` with symmetric edge padding — zero padding
  would dip the tone ~10 % at every 11 ms block edge). A bounded 2 s ring
  decouples the DSP thread from the PortAudio callback; the callback only
  copies (silence on starvation). Starved frames *while a squelch-open channel
  is active* are counted as underruns and reported at end of run; silence with
  no open channel is by design and not counted.
- **Gain:** fixed **1/5000 × vol**: 3 kHz deviation → 0.6 FS. No
  auto-normalization (it pumps with speech level changes).
- **Gating:** per-channel **5 ms linear ramp** on squelch open/close
  (click-free); only squelch-open channels are mixed; the sum is clipped at
  ±0.95.
- **Wiring:** the worker taps each channel's audio block non-destructively
  (the WAV save path is untouched) → `push(f, audio, squelch_open)` per
  channel per chunk, `flush()` at chunk end, `drop(f)` when a channel leaves
  the demod set. No output device → warn, continue recording only.

Unit-tested without hardware (`tests/test_nbfm.py::test_audio_player`): the
mix/resample/ring core is pure numpy and the PortAudio callback is called
directly — gain (1.8 kHz dev → 0.36 FS, measured 0.360), 6× resample, squelch
ramp, closed-channel silence, two-channel linear sum (0.72 FS), clip plateau
(0.95; the resample FIR rings ~2 % at the hard clip edge), callback fill +
underrun count.

---

## 5. Implementation path

1. **Prototype — Python + NumPy + ctypes (API v3).**
   No official PyPI package exists (verified); bind `sdrplay_api.dll` (Windows) with
   ctypes. Key ABI notes (measured, cost real debugging time): `sdrplay_api_Update`
   takes the value inside the param struct + bitflag (no 32-bit freq arg in v3);
   `RfFreqT` holds a `double` → `rfFreq` at offset 40 in `TunerParamsT`; `DecimationT`
   is byte-aligned → at offset 2 in `ControlParamsT`. IQ arrives callback-only as
   int16 I/Q blocks (measured: ≈1 344 samples average, bursty — §9.1); assert struct
   layouts with `ctypes.offset_of`.
   Call sequence: `Open → LockDeviceApi → GetDevices → SelectDevice (RSPdx ID = 4) →
   GetDeviceParams → SetDeviceParams (fsFreq) → SetTunerParams (ZIF, rfFreq, bwType,
   gain) → SetRspDxParams (antenna A) → StreamACbFn/EventCbFn → Init → (live: Update)
   → Uninit → ReleaseDevice → Close`.
   Deliverable: waterfall + per-channel activity + 1–4 demodulated audio channels.
2. **Full system — GNU Radio 3.9+ with `gr-sdrplay3`** (pybind11, API v3). Blocks:
   SDRPlay source → bandpass → FFT activity monitor → per-channel DDC/decimate
   (polyphase) → PLL FM demod → de-emphasis → squelch → audio sink; Qt/SDL GUI for
   waterfall, channel grid, activity and audio. Same math as §4.3; Numba or C++
   blocks for the per-channel FIRs.
3. **Alternative data path — SDRconnect** (WebSocket on TCP 5454: JSON control +
   binary IQ frames) if you prefer a daemon + thin client or want to reuse
   qt-dab-style tooling.

Performance budget: 6.048 MSPS int16 = 24.2 MB/s IQ; full-band monitor ≈ trivial;
per active channel ≈ 0.33 GFLOP/s (§4.3). A modern x86 CPU handles the whole
single-band design with headroom.

---

## 6. Verification / test plan

1. **Tone & SINAD:** sig gen at 145.000 / 434.000 MHz, 1 kHz tone at 1.8 kHz
   deviation (60 % of ±3 kHz, per EN 300 086 cl. 6.1). Sweep level; plot SINAD.
   Target: 12 dB SINAD at ≈ −118 dBm (2 m) / −116 dBm (70 cm) *(est.)*.
2. **Deviation headroom:** ±3 kHz and ±5 kHz tones → no clipping, audio 300–3000 Hz
   flat ±2 dB (2 kHz sweep), de-emphasis verified (−6 dB/oct vs. flat).
3. **Selectivity:** wanted + 12.5 kHz-adjacent equal-level signal; co-channel
   interferer at ±1.5/±3 kHz offsets up to +12 dB (EN 300 086 cl. 8.3.3 limits).
4. **Squelch:** noise-floor open/close, hysteresis, timing (1 s close-out).
5. **CTCSS:** 67 / 123 / 250.3 Hz (TS 103 236 test set) plus ±1 % edge cases
   (66.33 / 67.67 Hz) → correct open/close; 1750 Hz burst → burst flag.
6. **Field:** monitor 145.500 (FM calling) and local repeaters (both 2 m offset
   +600 kHz and 70 cm 1.6/2.0 MHz shift pairs); log activity CSV.
7. **S-meter (optional display):** IARU R1 recommendation (Brighton 1981 /
   Torremolinos 1990): 1 S-unit = 6 dB; above 30 MHz S9 = **−93 dBm**; quasi-peak
   detector, 10 ms attack, ≥500 ms decay.

---

## 7. Extensions (same DSP core)

- **DMR voice** (ETSI TS 102 361-1): 4FSK, 4.8 ksym/s, 12.5 kHz channel, 2-slot TDMA —
  dominant digital voice on EU amateur VHF/UHF; reuse the channel DDC, add a 4FSK
  slicer (e.g. libdv).
- **APRS:** 1200 baud AFSK (1300/1700 Hz) on 144.800 (2 m) / 432.500 (70 cm) —
  decode from the 8 kHz audio with a Bell-202 demod.
- **DTMF** (repeater control, per IARU FM.24): Goertzel 16-tone bank on the audio.
- **Waterfall/scan UI, CSV logging, multi-unit (dual-band) mode.**

---

## 8. Prototype status and verified performance (2026-10-25)

Implemented: `nbfm.py` (DSP core), `prototype.py` (SDRPlay → DSP → report pipeline),
`playback.py` (audio output stage, §4.5), `sdrplay_api.py` (ctypes binding, API v3.11),
`tests/test_nbfm.py` (all passing).
Demodulator is the **2nd-order digital PLL** of §4.3 (implemented 2026-10-25; the
quadrature discriminator is kept as a selectable fallback). All numbers below are
measured on the RSPdx (SerNo 201004A243) with the workspace venv
(`D:/Source/hermes-dir/.venv`, Python 3.11, numpy 2.4.6, scipy, numba).

### 8.1 Throughput — measured, before/after optimization

| Run | Processed rate | Drops (10 s) | Notes |
|---|---|---|---|
| Baseline (per-service-block Python DSP, ~1 k-sample blocks) | 0.84 MSPS (3.3× slower than real time) | 47 876 | per-block overhead dominated |
| After numba NCO + float32/complex64 + chunked worker (65 536) + shared-FFT `BandSpectrum` | 4.0–5.6 MSPS | 9 494–21 591 | worker capacity 8.5 MSPS, but fixed-depth queue overflowed on delivery bursts |
| **Final: unbounded queue** | **6.000 MSPS (2 m), 6.003 MSPS (70 cm)** | **0 / 0** | real-time, sustained 10 s, both bands |

Worker CPU budget at real time (measured breakdown, 10 s run): int16→complex64
conversion ≈ 21 % of a core, DSP (1 active channel + full-band spectrum) ≈ 47 %,
SDRPlay callback (as_array + copy + enqueue) ≈ 25 % — total ≈ 93 % of one core,
~7 % headroom. (Optimization-phase measurement; with the PLL demod the
worker-timing breakdown instead shows ≈ 46 % busy + ≈ 54 % queue-wait — §8.3.)

Per-component cost at 6.048 MSPS, 1 active channel (replay of a live capture,
`diag/profile_replay.py`):

| Stage | ms/s |
|---|---|
| int16 → complex64 conversion | 67.5 |
| chunk concatenation | 6.0 |
| full-band spectrum (8192-pt FFT + 160-channel grid + waterfall row) | 51.7 |
| channel demod (NCO, 7·9·6 decimation, PLL, DC-block, LP, de-emphasis, carrier detect) | 423.9 (433.7 in the first run) |

Per-component numbers vary with machine state between replay runs (±35 %
observed for identical code — CPU frequency/thermal state), so the
like-for-like comparison is a back-to-back A/B on one captured stream
(`diag/profile_pll_ab.py`): **PLL 252.8 ms/s vs discriminator 262.9 ms/s** —
the PLL is ~4 % *cheaper* (the numba PLL loop costs ≈1 ms/s at 16 kHz; the
discriminator's per-sample `|x|·|x[n−1]|` normalization is the expensive part).

### 8.2 What the profiling found (cost real time; keep these in mind)

- **`scipy.signal.lfilter(b, 1.0, x_complex64)` runs in complex128**: a Python
  float `a=1.0` promotes the whole recursion. Passing `np.float32(1.0)` keeps it in
  complex64 — 1.7× faster (3.34 → 1.94 ms per 65 536-sample h1 pass).
- **NCO via numpy `np.exp` is ~15× the cost of a numba geometric-recurrence fill**
  (≈250 ms/s vs ≈16 ms/s at 6.048 MSPS). Residual phase drift over one 65 536-sample
  chunk < 0.007 rad — removed by the DC blocker.
- **Per-service-block Python overhead was the original 3.3× real-time cause**: the
  service delivers ~1 344-sample blocks ≈ 4 500 times/s; chunking to 65 536 amortizes
  the decimation-chain call overhead ~50×.
- **Delivery is bursty**: bursts of a few hundred blocks arrive in ~10–20 ms, then
  silence. Queue depths 16/64/512 all overflowed (42 %/18.5 % drops); unbounded
  queue → 0 drops (worker capacity 8.5 MSPS > 6.87 MSPS average feed, so the queue
  stays bounded in practice, ≈5 MB peak).
- **API v3.11 has no `sdrplay_api_SetSampleSize`** (checked the DLL export table and
  the installed headers) — block size cannot be requested; the service decides.
- `sdrplay_api_GetDevices` leaves the device handle NULL until `SelectDevice`;
  device validity must be tested on the `valid` flag, not the handle.
- Event callback: only **event 0** carries a valid gain (dB); the power-overload
  event (1) carries gain 0.0 and must not overwrite the last real gain reading.
- Decimation must be anchored to a **global sample grid** (keep sample k iff the
  global index is divisible by 7/63/378): per-block `::7` slicing on odd-length
  blocks drifts the audio rate (measured 1 000.0 Hz only after the fix).

### 8.3 Hardware smoke test results (10 s each, 2026-10-25)

| | 2 m | 70 cm |
|---|---|---|
| LO / IF BW | 145.000 MHz / 5 MHz | 434.000 MHz / 8 MHz |
| Processed rate | 6.000 MSPS | 6.003 MSPS |
| Dropped blocks | 0 | 0 |
| Total gain (AGC 50 Hz, setpoint −20 dBFS) | 12.1 dB | 20.6 dB |
| Noise floor (median channel power) | −28.1 dBFS | −28.1 dBFS |
| Auto-selected channel | 144.000 MHz (SNR 18.8 dB over floor) | 432.000 MHz (SNR 21.8 dB over floor) |
| Squelch | closed (carrier detect 1.46 dB < 12 dB open threshold) | closed |
| Waterfall | plausible: noise floor + strong saturated interference at the upper band edge (145.9x MHz) | same pattern at 435.9x MHz |

**Follow-up probe (2026-10-25, `diag/probe_signal.py`):** the band-edge
excesses are *real external* wideband emissions, not receiver artifacts —
with the LO moved to 144.500 MHz the +18.2 dB excess stays at 144.000 MHz
absolute (432.000 likewise on 70 cm). Their demodulated audio is flat
noise (no tones, spectral flatness ≈ 0.56, carrier-detect ≈ 0 dB): not FM
voice, so the closed squelch is correct. Two root causes were found:
1. **Demod reach:** channels >±300 kHz from the LO are filtered by the /7
   anti-alias stage before demod (§4.3) — 144.000 with LO 145.000 demodulates
   to the noise floor even though the full-band spectrum sees the signal.
2. **Auto-select** could pick such unreachable channels (it picked 144.000).
   Fixed: auto-select is now restricted to the ±300 kHz window around the LO,
   and explicit `--ch` outside the window prints a warning.

Post-fix smoke runs (10 s each, 2026-10-25): 2 m — 6.000 MSPS, 0 drops,
1 overload event (the strong 144.000 emission), gain 12.1 dB, floor −27.8
dBFS, no in-window channel above the +12 dB select threshold (nothing
demodulated — correct). 70 cm — 6.000 MSPS, 0 drops, 1 overload, gain 20.6
dB, floor −28.1 dBFS, same pattern. No station was in range to verify voice
intelligibility — that is test-plan item 6 (field monitoring of 145.500 /
local repeaters); a signal inside the ±300 kHz LO window is needed to see
the squelch open live.

**PLL demod smoke runs (10 s each, 2026-10-25):** 2 m — 6.000 MSPS, 0 drops,
1 overload, gain 12.1 dB, floor −27.5 dBFS; worker-timing: queue-wait 53.7 %,
convert 21.7 %, process 24.4 %. 70 cm — 6.003 MSPS, 0 drops, 1 overload,
gain 20.6 dB, floor −28.1 dBFS; queue-wait 44.7 %, convert 31.7 %, process
23.4 %. Real-time on both bands with the worker thread ~45–55 % idle
(queue-wait) — the PLL adds no real-time pressure.

### 8.4 Unit-test suite (`tests/test_nbfm.py`, all pass, 2026-10-25)

Core DSP: ctypes struct layouts (API v3.11 ABI) · squelch open/close
hysteresis · frequency accuracy (global decimation grid) · adjacent-channel
rejection (54.3 dB at 12.5 kHz, equal level) · tone recovery (the PLL is linear
in the deviation: fundamental `D·|T_f(1 kHz)|·|H_lp|·|H_dc|·|H_de|`, `T_f` the
exact discrete loop tracking response — measured 1941 Hz vs 1954 expected) ·
PLL carrier-offset tracking (+150 Hz offset on the carrier removed, audio DC
0.00 Hz) · PLL startup at low frequency / large deviation (300 Hz @ 1.8 kHz,
β = 6: matches the tracking model to −0.1 %, no slip harmonics — regression
for the startup-kick slip orbit, §4.3) · multitone linearity (400/1000/2000 Hz
@ 800 Hz dev each, 2.4 kHz peak: each tone unchanged by the other two, IM
products ≤ −40 dB below the tones — worst −82.6 dB) · de-emphasis (−6 dB/oct,
digital bilinear at 8 kHz) · CTCSS on/off (Goertzel, ±1 % edge cases) ·
activity-monitor channel detection · audio output stage (playback.py: gain
1/5000, 6× resample, squelch ramp, mix/clip, ring + callback — §4.5).

EN 300 086 conformance (12.5 kHz spacing, cl. 6.1 normal modulation =
1 kHz @ 60 % MPFD, O.41 psophometric SND/ND):

| Test | Result |
|---|---|
| Co-channel rejection, cl. 8.3.3 (crossing must lie in −12…0 dB; unwanted raised until SND/ND = 14 dB) | PASS — crossings −8.3 dB (0 Hz), −9.1 dB (±1.5 kHz), −8.0/−8.2 dB (±3 kHz); worst −9.1 dB (discriminator: −8.3 dB) |
| Adjacent-channel selectivity, cl. 8.4 (≥ 60 dB) | PASS — ≥ +80 dB (scan ceiling), both ±12.5 kHz |
| Deviation headroom ±3 kHz / ±5 kHz: no deviation aliasing (16 kHz phase step < π), fundamental within 5 % of `D·|T_f(1 kHz)|`·stages, SND/ND ≥ 20 dB @ ±3 kHz | PASS — 67°/109° max step; 3235 Hz vs 3257 expected @ ±3 kHz (−0.7 %); SND/ND 85.6 dB @ ±3 kHz. ±5 kHz (the 25 kHz MPFD — out of spec for our 12.5 kHz grid): the loop runs in a slip state that preserves the deviation (§4.3) — fundamental 5387 Hz vs linear 5428 (−0.8 %), SND/ND 56.1 dB (info, not gated; the discriminator fallback: 30.1 dB) |

Measurement gotchas fixed along the way (kept for the record): the psoph FIR
must be trimmed of its zero-state fill-up transient (a 32 ms burst smears a
~44 Hz pedestal around the notched fundamental and cut SND/ND by 22 dB); the
co-channel scan must start quiet and raise the unwanted (a scan down from
0 dB mis-reports the crossing); the /7–/9–/6 decimation grids must be anchored
to a global sample counter (per-block slicing drifts the audio rate by 1.5 %
and smears every tone with a phase-noise skirt).

---

## 9. References (primary sources, verified 2026-10-24)

1. IARU R1 VHF band plan (effective Dec 2020, VGC Novi Sad) —
   `iaru-r1.org/wp-content/uploads/2020/12/VHF-Bandplan.pdf`
2. IARU R1 UHF band plan (effective Dec 2020) —
   `iaru-r1.org/wp-content/uploads/2021/03/UHF-Bandplan.pdf`
3. IARU R1 VHF Handbook 10.03 (2026) — FM.1 / FM.24 NBFM & repeater standards,
   S-meter recommendation — `iaru-r1.org/wp-content/uploads/2026/02/VHF_Handbook_V10_03_final.pdf`
4. ETSI EN 300 086 V2.1.2 (2016-08) — LMS 12.5/20/25 kHz analogue speech, Table 3 MPFD —
   `etsi.org/deliver/etsi_en/300001_300099/300086/02.01.02_60/en_300086v020102p.pdf`
5. ETSI TS 103 236 V1.1.1 (2014-01) — CTCSS/DCS, tone table 67.0–250.3 Hz —
   `etsi.org/deliver/etsi_ts/103200_103299/103236/01.01.01_60/ts_103236v010101p.pdf`
6. SDRplay RSPdx datasheet V1.0 — `sdrplay.nl/wp-content/uploads/RSPdxDatasheetV1.0.pdf`
7. SDRplay Hardware API v3 (≥3.06 for RSPdx) — `sdrplay.com` (API spec + `SDRplay/examples`,
   `SDRplay/ExtIO_SDRplay` reference implementation)
8. IARU R1 band-plan index page — `iaru-r1.org` → VHF & up Bandplanning (C5 committee)
