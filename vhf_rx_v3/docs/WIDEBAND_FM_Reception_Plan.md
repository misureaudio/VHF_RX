# Wideband FM Reception in 88–108 MHz — Modification Plan (vhf_rx_v2)

**Status:** plan only (no code changes made yet).
**Target:** add wideband FM (broadcast) reception to the vhf_rx v2 RSPdx NBFM
receiver prototype (`vhf_rx_v2/`), covering the 88–108 MHz FM broadcast band.

---

## 1. Goal / scope

Add an FM broadcast reception mode to the existing v2 receiver chain: tune, acquire,
demodulate and play wideband FM signals across ~88–108 MHz (EU broadcast band
87.5–108 MHz), reusing the v2 DSP core where possible.

Non-goals for the first pass: stereo / RDS decoding, digital subcarriers, multi-band
simultaneous reception.

## 2. Target signal model (broadcast FM)

| Item | Value | Source / note |
|---|---|---|
| Band | 87.5–108 MHz (user range: 88–108 MHz) | EU broadcast allocation |
| Modulation | Wideband FM (analog voice) | vs. v2 NBFM 12.5 kHz / ±3 kHz |
| Deviation | ±50 kHz (EU) / ±75 kHz (US) | **DECISION (2026-10-29): adopt **±50 kHz** (EU adherence) — supersedes the ±75 kHz design target; a ±75 kHz capability would be a reception superset but is not the EU-aligned design value |
| Stereo / RDS | Stereo: **OK** (approved); selectable **`--mono`** option; **no RDS** | |
| Channel spacing | **100 kHz** (EU) | some regions 200 kHz |
| Audio band | 30 Hz–12 kHz (mono) | stereo main +19 kHz pilot → later |
| Emission | 75K0F3E-ish; Carson BW ≈ 2(Δf + f_a) ≈ 2(75 + 12) ≈ **174 kHz** @ ±75 kHz | vs. NBFM Carson ≈ 2(5+3) ≈ 16 kHz |
| De-emphasis | 50 µs (EU broadcast) / 75 µs (US) | amateur v2 default is 75 µs |

## 3. Baseline facts in v2 that the modification must build on / overcome

From `vhf_rx_v2/` (design doc `docs/RSPdx_VHF_NBFM_Receiver_Design.md`, verified
2026-10-25):

- Capture: fs **6.048 MSPS**, IF BW **5 MHz** (2 m) / **8 MHz** (70 cm), spans ≈
  ±2 MHz around LO (design tables); RSPdx max IF BW available = **8 MHz**.
- Demod chain per channel: NCO DDC → **/7** (h1 anti-alias 0.40·fs/7 = **345.6 kHz**) →
  **/9** (h2) → **/6** (h3: channel BPF pass ±8 kHz, stop ±12.5 kHz) → **16 kHz complex**
  → FM demod (2nd-order PLL, f_n = 1750 Hz, ζ = 0.707, loop BW ≈ 2.4 kHz) → audio =
  deviation @16 kHz → DC block → LP Butterworth **3.4 kHz** → **/2** → **8 kHz** →
  de-emphasis (τ = 75 µs) → CTCSS → squelch (carrier detect) → playback (8 kHz → 48 kHz,
  fixed gain **1/5000**, 3 kHz dev → 0.6 FS, 5 ms squelch ramp).

- **Critical limit:** FM phase detector wraps at π → unambiguous deviation
  |dev| < fs_audio/2 = **±8 kHz** at 16 kHz complex. ±5 kHz verified; ±8 kHz is the hard
  ceiling. Broadcast ±50/±75 kHz is **6–9× beyond** → the final decimation to 16 kHz
  (and the ±8 kHz BPF) must be re-scaled for wideband.

- Demod reach: channels within ≈ **±300 kHz of the LO** demodulate (h1 anti-alias);
  full-band monitor (BandSpectrum) covers the captured span.

- Squelch: carrier-detect ratio (audio band vs upper band) — modulation-agnostic,
  reusable for broadcast.

## 4. Required modifications per component

### 4.1 `nbfm.py`
- `BandSpectrum`: `spacing` / `half_bw` already parameterized → use for broadcast
  (spacing **100 kHz**, half-band ≈ dev + f_a + margin ≈ **100 kHz**).
- `ChannelDemod`: parameterise the audio path for wideband (defaults keep current
  NBFM behaviour):
  - channel BPF: pass ±(dev + f_a + margin) ≈ **±100 kHz**, stopband ≥ 2× pass
    (currently firwin(128, 7500, fs=96 kHz) → ±8 kHz);
  - decimation chain: for broadcast, final complex rate must be **≥ ~2×dev ≈ ≥150 kHz**
    (design ≈ **432 kHz**): propose `/7 → BPF(±100 kHz) → /2` (i.e. /14 overall) then
    demod; the old `/9 → /6 → 16 kHz` schedule is not viable once the BPF is widened
    (96 kHz complex rate would alias ±75 kHz deviation);
  - FM demod (PLL): loop gains re-scaled to the new complex rate; loop BW should be
    re-validated/re-tuned for the wider audio band (30 Hz–12 kHz → e.g. f_n up from 1750 Hz);
  - deviation aliasing: with complex rate ≥ 2·dev, ±75 kHz demods without phase wrap;
  - audio: LP pass ≈ **12–15 kHz** (was 3.4 kHz); final audio rate ≥ 12 kHz (was /2 → 8 kHz,
    must change for broadcast, e.g. final 48 kHz or 16 kHz);
  - de-emphasis: τ = **50 µs** (broadcast) vs 75 µs (NBFM) — per-instance parameter;
  - squelch: carrier-detect only (CTCSS not used for broadcast).

### 4.2 `prototype.py`
- Add band entry (e.g. key `"wbfm"`): rf (band centre, e.g. 98 MHz), bw (BW_8M), fmin /
  fmax (88–108 MHz), spacing (100 kHz), dev, de-τ, audio params.
- `band_lna_max(profile, band)`: map new band to the RSPdx LNA-state bucket that covers
  88–108 MHz (verify bucket boundary vs the "2 m" VHF bucket).
- CLI: `--band` choices += new key; `--freq` (LO) override; `--ch` fixed channels;
  auto-select uses new spacing grid.
- Worker: BandSpectrum + ChannelDemod take band-specific parameters; **retune hook**
  (LO retune via `SDRPlay.update(rf_hz=...)`) needed for whole-band coverage:
  RSPdx max IF BW 8 MHz cannot span 88–108 MHz in one tuning → scanner/retune-on-detect.

### 4.3 `sdrplay_api.py`
- Expected: **no ABI change**. Reuse `update(rf_hz=...)`, tuner `bwType` (BW_8M), LNA,
  AGC; only new band → profile LNA mapping. Verify IF BW 8 MHz is selectable for the
  broadcast band.

### 4.4 `playback.py`
- Gain: fixed 1/5000 (3 kHz dev → 0.6 FS) → must scale with deviation / become per-band
  (75 kHz dev → ~15× → would clip at fixed gain).
- Sample rate: final audio ≥ 12 kHz (e.g. 16/48 kHz) → playback fs_in/fs_out adapt.

### 4.5 `tests/test_nbfm.py`
- New offline wideband tests: ±50/±75 kHz deviation recovery (PLL linear, no aliasing),
  SND/ND; audio LP passes 12 kHz; de-emphasis τ=50 µs; BandSpectrum 100 kHz grid detects
  broadcast channel; carrier squelch on broadcast-like FM.
- Existing NBFM tests must keep passing (defaults unchanged).

### 4.6 `diag/`
- `--band` choices flow from `BANDS`; add usage examples; probes reuse band params.

## 5. Hardware front-end gate (pre-condition, may block)

- RSPdx **internal preselector** passband (design doc §3.1: BP **120–250 MHz** / BP
  **420–1000 MHz**) vs 88–108 MHz → likely **not covered** → external bandpass
  (e.g. 88/87.5–108 MHz) mandatory before DSP work pays off.
- RSPdx **LNA state** coverage for the <108 MHz region → check LNA grouping table; may
  need external LNA or different state bucket.
- Max input 0 dBm / 10 dBm short periods; broadcast ERP → AGC + optional 6 dB pad.
- 24 MHz GPSDO reference recommended (as for repeater work).

## 6. Implementation phases (proposed)

1. **Phase 0 — decisions:** deviation target (±50/±75 kHz), final audio rate, de-τ,
   channel spacing, LO/retune strategy, playback gain handling.
2. **Phase 1 — offline DSP:** re-scaled per-channel chain (BPF, decimation, PLL, audio,
   de-emph) + BandSpectrum broadcast grid + new offline tests (no HW).
3. **Phase 2 — band plumbing:** `BANDS` entry, LNA mapping, CLI, worker/report/playback
   integration.
4. **Phase 3 — retune/scan + live validation:** LO retune-on-detect for whole band,
   auto-select, waterfall, smoke tests, probe_signal classification of a real station.
5. **Phase 4 — field validation:** SINAD, deviation, audio band, selectivity, squelch,
   field monitoring.

## 7. Verification plan (broadcast-adapted)

Tone/SINAD at broadcast deviation; deviation ±50/±75 kHz (no aliasing, SND/ND ≥ 20 dB);
audio band flat to 12 kHz; selectivity/co-channel (100 kHz spacing → easier than NBFM);
squelch; field test with real broadcast stations.

## 8. Open design decisions

1. Deviation: ±50 kHz (EU) only, or ±75 kHz (US) design target?
2. Stereo / RDS: out of scope for first pass?
3. Final audio rate for broadcast (16 kHz vs 48 kHz)?
4. Playback gain: per-band constant vs deviation-derived?
5. LO retune strategy: retune-on-detect scanner vs fixed-LO-per-station (`--freq`/`--ch`).

---

## 9. Agreed / proposed decisions (2026-10-29, vhf_rx_v3)

1. **Deviation target: ±50 kHz** (user decision — EU adherence; supersedes the ±75 kHz
   design target; ±75 kHz would be a reception capability superset but is not the
   EU-aligned design value).
2. **Stereo: approved.** Implement stereo decoding capability; add selectable
   **`--mono`** option; **no RDS**.
3. **Final audio rate (proposal):** wideband path → **48 kHz** mono by default
   (`--mono`); the NBFM path stays unchanged (final 8 kHz). Rationale: WASAPI-native,
   headroom for stereo audio (~15 kHz), playback ratio 1.0, simple gain mapping.
   Stereo = 2-channel extension later.
4. **Playback gain (proposal):** deviation-normalised — `gain = 0.6 / dev_max` per band
   (NBFM @ 3 kHz dev → 1/5000; WBFM @ 75 kHz → ≈ 1/125000). Stereo: same per-channel gain.
5. **Retune strategy (proposal):** retune-on-detect scanner for `--band wbfm`
   (LO = band centre, IF BW 8 MHz; auto-select → retune LO → demod → repeat); fixed
   `--freq`/`--ch` retained.
