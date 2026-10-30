# vhf_rx v3 — wideband FM broadcast band (88–108 MHz) band plumbing + retune

Phase 2 of the wideband FM reception plan (`docs/WIDEBAND_FM_Reception_Plan.md`);
Phase 1 (offline re-scaled wideband DSP chain) landed before this. Phase 3 adds the
retune-on-detect hook.

## What changed

### `prototype.py`
- `BANDS` += `"wbfm"`: rf = band centre (98 MHz), bw = BW_8M, fmin/fmax = 88–108 MHz,
  spacing = 100 kHz (EU), dev = 50 kHz (EU), tau = 50 µs, audio_fs = 48000, wideband = True.
- Worker: `BandSpectrum` takes the band `spacing` (100 kHz for wbfm, 12.5 kHz default).
- Worker: `_ensure_demod` builds the Phase-1 wideband `ChannelDemod` for `band == "wbfm"`
  (wideband=True, bp_pass_hz ≈ 62 kHz, bp_stop_hz ≈ 124 kHz, bp_taps 128, tau = 50 µs,
  fs_audio = 48 kHz); other bands keep the NBFM chain.
- Worker: retune-on-detect hook + retune-back-to-band-centre when the demod set empties
  (wbfm only).

### `sdrplay_api.py`
- `SDRPlay.retune(rf_hz)` — LO retune via the existing `update(rf_hz=...)` (post-Init).

### `playback.py` / CLI
- Playback gain now deviation-normalised: `gain = 0.6 / dev_max * vol`
  (NBFM @ 3 kHz → 1/5000; WBFM @ 50 kHz → ~1/83333).
- For `band == "wbfm"`, player `fs_in` = final audio rate (48 kHz).
- Report now records `spacing`, `dev_hz`, `audio_fs_hz` per band.

## Verification
- `python prototype.py --help` → `--band {2m,70cm,wbfm}`.
- Module import + `BANDS['wbfm']` / `band_lna_max('wbfm')` / `ChannelDemod` / `BandSpectrum`
  signatures verified via the project venv.
- Full suite: all NBFM + wideband + sdrplay_api tests PASS (only the pre-existing wideband
  activity-grid test remains red: SNR -0.6 dB — test-side grid sensitivity, Phase 3 follow-up).

## Not done (hardware)
- Live retune / whole-band reception on the RSPdx (external bandpass likely mandatory).

# vhf_rx v2 — antenna port selection (A/B) + explicit bias-T control

v2 = v1 + input-port selection. The connected unit is still probed
(`--radio auto/rspdx/rsp2pro`); v2 adds **which port** it listens on.
The RSPdx-only original (`vhf_rx/`, v0) and v1 (`vhf_rx_v1/`) remain
untouched.

## What changed (v1 → v2)

### `sdrplay_api.py`
- **`configure(..., antenna=0, bias_t=False)`** — the `antenna` argument
  (0 = port A, 1 = port B) already existed in v1 but nothing upstream
  passed it, so every run used port A. v2 wires it through, and adds
  **`bias_t`**: the 4.7 V / 100 mA bias feed on port B, **disabled by
  default and set explicitly** (`biasTEnable = 0` in the parameter tree
  before Init) rather than relying on the service default.
- **`update(..., bias_t=None)`** — live bias-T changes with the correct
  per-family reason bit (RSPdx ext1 0x02 / RSP2 main 0x80), matching the
  existing antenna routing (RSPdx ext1 0x04 / RSP2 main 0x200).
- `_set_bias_t()` helper mirrors `_set_antenna()`: per-family struct
  routing (`rspDxParams.biasTEnable` vs `rsp2TunerParams.biasTEnable`) —
  writing the wrong family's struct is silently ignored by the service.

### `prototype.py`
- **`--antenna {a,b}`** (default `a`) and **`--bias-t`** (flag; OFF by
  default — a passive antenna on port B must not be bias-fed).
- The startup line now states the actual front-end configuration:
  `... port B (bias-T off)`, and the JSON report records `antenna` and
  `bias_t` alongside the v1 device fields.

### `diag/probe_signal.py`
- Same `--antenna {a,b}` / `--bias-t` options, so diagnostics probe the
  port the receiver actually uses (v2 default is still `a`).

### Unchanged
- Everything else in v1: probe-based device selection, per-device sample
  scale (8192 RSPdx / 32768 RSP2pro) and LNA limits, the DSP chain
  (`nbfm.py`), the audio stage (`playback.py`), the test suite.

## Usage (typical: passive antenna on port B)

```
python prototype.py --band 2m --antenna b --play --duration 30
python prototype.py --band 70cm --antenna b --ch 434.600e6 --play
python diag/probe_signal.py --band 2m --antenna b --duration 15
```

Add `--bias-t` only for a bias-fed (active) antenna on port B.

## v1 (for reference)

v1 added the RSP2/RSP2pro as a co-equal, probe-selected front end:
`DeviceProfile`/`PROFILES` (sample scale 8192 vs 32768, LNA limits
26/20 vs 8/5, antenna path), `--radio auto/rspdx/rsp2pro`, per-device
LNA defaults, device fields in the report, `test_device_profiles()`,
and the diag scripts ported off hard-coded RSPdx constants. Full diff:
`vhf_rx_v1/CHANGES.md`.
