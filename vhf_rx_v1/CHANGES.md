# vhf_rx v1 — RSPdx + RSP2/RSP2pro dual front end

v1 is the first version that supports **both** supported SDRplay front ends in a
single codebase. The connected unit is **probed at startup** (device
enumeration by `hwVer`, never by serial prefix) and selected with `--radio`
(default `auto`). The RSPdx-only original survives untouched as `vhf_rx/`
(v0).

## What changed (v0 → v1)

### `sdrplay_api.py`
- **`DeviceProfile` + `PROFILES` + `family_of()`** — per-device parameters that
  differ between families:
  | | RSPdx | RSP2 / RSP2pro |
  |---|---|---|
  | API hwVer | 4 | 2 (no RSP2pro ID of its own; RSP2 and RSP2pro are electrically identical, datasheet v2.2) |
  | int16 full scale | 8192 (14-bit values in int16 at ≤ 6.048 MSPS) | **32768 (native 16-bit — measured** on the connected unit: max sample 18 845 > 16 383, low 4 bits all 16 patterns uniform) |
  | LNA states 2 m / 70 cm | 0–26 / 0–20 | 0–8 / 0–5 |
  | Antenna select | `rspDxParams.antennaSel` 0=A/1=B/2=C + ext1 bit 0x04 | `rsp2TunerParams.antennaSel` 5=A/6=B, `amPortSel` 1=Hi-Z/0, main update bit 0x200 |
  | Min API version | 3.06 | 3.08 |
- **`find_device(radio="auto")`** replaces the hard-coded `RSPDX_ID` default:
  - `auto`: the sole valid device is used whatever its family; if several are
    connected, the RSPdx (project baseline) is preferred — override explicitly.
  - `rspdx` / `rsp2pro`: only that family is accepted; a clear error lists
    what is actually present.
  - on return, `rx.profile` describes the selected device (scale, LNA limits,
    label).
- **Per-family antenna routing** (`_set_antenna`): writing antenna selection to
  the wrong family's struct is silently ignored by the service, so v1 routes
  through `rspDxParams` (RSPdx) or `rsp2TunerParams` (RSP2/RSP2pro) with the
  matching update-reason bit. Port B selection and bias-T now actually apply
  on the RSP2pro (v0 wrote `rspdxParams`, which the service ignores for
  hwVer 2).
- New constants: RSP2 antenna enum (5/6), `amPortSel` (1=Hi-Z/0), RSP2
  update-reason bits (0x80 bias-T / 0x100 am-port / 0x200 antenna / 0x400
  rf-notch / 0x800 ext-ref), RSP2 LNA state counts (9 / 6).
- Struct tree unchanged: the v3.15 `DevParamsT` / `RxChannelParamsT` /
  `TunerParamsT` layouts are byte-identical to what the v0 binding mirrors
  (the RSP1B / RSPdxR2 additions moved no field), so `tests/test_nbfm.py::
  test_structs` passes unmodified.

### `prototype.py`
- **`--radio {auto,rspdx,rsp2pro}`** (new; `--device` is still the *audio
  output* device).
- The hard-coded `I14_SCALE = 8192.0` is gone: the worker normalizes with the
  selected profile's `sample_scale` (8192 RSPdx / 32768 RSP2pro). A wrong
  scale silently shifts every dBFS quantity (noise floor, AGC headroom, the
  PLL's 0.05 amplitude gate, the 1/5000 playback gain) by 12 dB.
- **Band-max LNA default is now per device** (`band_lna_max(profile, band)`):
  26/20 on the RSPdx, 8/5 on the RSP2pro — v0's hard-coded RSPdx states
  would have failed `Init` with OutOfRange on the RSP2pro. `--lna` is range-
  checked against the selected profile.
- Startup prints the probed device: `API v3.15, device 1710060810 (hwVer 2)
  -> RSP2/RSP2pro (hwVer 2), int16 full scale 32768`.
- The JSON report now records `radio`, `radio_label`, `serial`, `hwver`,
  `sample_scale`.

### `tests/test_nbfm.py`
- New `test_device_profiles()` (offline): family classification by hwVer,
  per-device scales, LNA limits, min API versions, RSP2 antenna values.

### `nbfm.py`, `playback.py`
- **Unchanged** — the DSP chain and audio stage are device-independent.

## What does NOT change
- 6.048 MSPS capture plan: measured sustained **6.021 MSPS** on the RSP2pro
  (15 s average; per-second 5.81–6.19), same as the RSPdx's 6.000 MSPS.
  434 MHz / 8 MHz ZIF accepted.
- Bursty delivery (~1344-sample blocks) and the unbounded-queue design —
  identical on both devices.
- The ±300 kHz per-channel demod reach, the global decimation grid, the AGC
  setpoint range (−72…−20 dBFS below 8.064 MSPS), callback/event semantics
  (event 0 = gain change with valid gain; event 1 = power overload, gain 0.0).

## Measured on the connected RSP2pro (SerNo 1710060810, API v3.15, 2026-10-26)
- 145 MHz / 5 MHz ZIF / 6.048 MSPS, AGC off, LNA 8 + GR 0: total gain
  **53.67 dB**; 434 MHz / 8 MHz ZIF, LNA 5 + GR 0: **68.67 dB**.
- Power-overload events occur at max gain (both bands) — expected with a
  strong local signal or an open port; the prototype counts them (see the
  design doc §4.2 note on ack + gain reduction, still not wired in).

## Known follow-ups (not v1 scope)
- `ack_overload()` + automatic GR reduction on overload (design doc §4.2) —
  pre-existing gap, now more visible with the RSP2pro connected.
- If a capture ever runs at 8.064 MSPS (the wider-margin option), re-verify
  the delivered bit depth: the RSPdx's rate-dependent 14/12/10/8-bit pattern
  may apply to the RSP2pro too (at 6.048 MSPS it is 16-bit, measured).
