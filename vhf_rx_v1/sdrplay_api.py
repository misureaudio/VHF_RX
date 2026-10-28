"""ctypes bindings for the SDRplay Hardware API v3 (RSPdx + RSP2/RSP2pro).

Verified against the installed API 3.15 headers (C:\\Program Files\\SDRplay\\API\\inc)
and the v3.15 API specification PDF shipped in the install:
  * device params are obtained with sdrplay_api_GetDeviceParams() and may be
    modified BEFORE sdrplay_api_Init();
  * after Init, changes must be signalled with sdrplay_api_Update(dev, tuner,
    reason, reason_ext1) to take effect;
  * IQ arrives callback-only as int16 I/Q. The sample scale is PER DEVICE
    and must be taken from the selected DeviceProfile: RSPdx delivers
    14-bit values in 16-bit containers at <= 6.048 MSPS (full scale 8192),
    RSP2/RSP2pro deliver native 16-bit values (full scale 32768, measured).

Device selection (v1): the RSP2pro has no device ID of its own in API v3 —
it enumerates as hwVer 2 (RSP2 family; RSP2 and RSP2pro are electrically
identical, datasheet v2.2). Live devices are classified by hwVer, never by
serial prefix (serials do not reliably carry the model name).
"""
from __future__ import annotations

import ctypes as ct
import winreg
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

# --------------------------------------------------------------------------
# enums / constants (values from sdrplay_api.h, sdrplay_api_tuner.h,
# sdrplay_api_control.h, sdrplay_api_rspDx.h, sdrplay_api_callback.h)
# --------------------------------------------------------------------------

class Err:
    Success = 0
    Fail = 1
    InvalidParam = 2
    OutOfRange = 3
    GainUpdateError = 4
    RfUpdateError = 5
    FsUpdateError = 6
    HwError = 7
    AliasingError = 8
    AlreadyInitialised = 9
    NotInitialised = 10
    NotEnabled = 11
    HwVerError = 12
    OutOfMemError = 13
    ServiceNotResponding = 14
    StartPending = 15
    StopPending = 16
    InvalidMode = 17
    InvalidServiceVersion = 24

RSP1_ID = 1
RSP2_ID = 2
RSPDUO_ID = 3
RSPDX_ID = 4

TUNER_A = 1

IF_ZERO = 0
IF_450K = 450
IF_1_62M = 1620
IF_2_048M = 2048

# Bw_MHzT enum values are kHz
BW_200K = 200
BW_300K = 300
BW_600K = 600
BW_1_536M = 1536
BW_5M = 5000
BW_6M = 6000
BW_7M = 7000
BW_8M = 8000

AGC_DISABLE = 0
AGC_100HZ = 1
AGC_50HZ = 2
AGC_5HZ = 3

# sdrplay_api_ReasonForUpdateT (3.15)
UPD_DEV_FS = 0x00000001
UPD_TUNER_GR = 0x00008000
UPD_TUNER_FRF = 0x00020000
UPD_TUNER_BWTYPE = 0x00040000
UPD_TUNER_IFTYPE = 0x00080000
UPD_TUNER_LOMODE = 0x00200000
UPD_CTRL_AGC = 0x01000000
UPD_CTRL_OVERLOAD_ACK = 0x04000000

# sdrplay_api_ReasonForUpdateExtension1T
UPD1_RSPDX_HDR_ENABLE = 0x00000001
UPD1_RSPDX_BIAS_T = 0x00000002
UPD1_RSPDX_ANTENNA = 0x00000004

# RSP2/RSP2pro main update-reason bits (sdrplay_api_rsp2.h, 3.15)
UPD_RSP2_BIAS_T = 0x00000080
UPD_RSP2_AMPORTSEL = 0x00000100
UPD_RSP2_ANTENNA = 0x00000200
UPD_RSP2_RFNOTCH = 0x00000400
UPD_RSP2_EXTREF = 0x00000800

# RSP2/RSP2pro antenna / port selection (Rsp2TunerParamsT)
RSP2_ANT_A = 5        # port A (SMA 50 Ohm)
RSP2_ANT_B = 6        # port B (SMA 50 Ohm, 4.7 V bias-T)
RSP2_AMPORT_AB = 0    # RF to ports A/B
RSP2_AMPORT_HI_Z = 1  # RF to balanced Hi-Z port

# events
EVT_GAIN_CHANGE = 0
EVT_POWER_OVERLOAD = 1
EVT_DEVICE_REMOVED = 2

# LNA state counts (sdrplay_api_rspDx.h / sdrplay_api_rsp2.h, 3.15)
RSPDX_LNA_STATES_VHF_BAND3 = 27   # -> max LNAstate 26 (2 m band)
RSPDX_LNA_STATES_420MHZ = 21      # -> max LNAstate 20 (70 cm band)
RSP2_LNA_STATES_0_420MHZ = 9      # -> max LNAstate 8 (2 m band)
RSP2_LNA_STATES_420_1000MHZ = 6   # -> max LNAstate 5 (70 cm band)


# --------------------------------------------------------------------------
# device profiles (v1: RSPdx + RSP2/RSP2pro co-exist, selected by probe)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class DeviceProfile:
    """Per-device parameters that change between RSP families.

    sample_scale: int16 full scale for the IQ samples the service delivers
      (RSPdx: 14-bit values in int16; RSP2/RSP2pro: native 16-bit —
      measured, see CHANGES.md). A wrong scale silently shifts every
      dBFS quantity by 12 dB.
    lna_max_2m / lna_max_70cm: highest valid LNAstate per band (out of range
      fails Init with OutOfRange).
    min_api: minimum API version that knows this family.
    """
    key: str
    label: str
    hwvers: tuple
    sample_scale: float
    lna_max_2m: int
    lna_max_70cm: int
    min_api: float


PROFILES: dict[str, DeviceProfile] = {
    "rspdx": DeviceProfile(
        key="rspdx", label="RSPdx", hwvers=(RSPDX_ID,),
        sample_scale=8192.0,
        lna_max_2m=RSPDX_LNA_STATES_VHF_BAND3 - 1,
        lna_max_70cm=RSPDX_LNA_STATES_420MHZ - 1,
        min_api=3.06),
    # API v3 reports the RSP2pro as hwVer 2 (RSP2 family) — there is no
    # RSP2Pro_ID or RSP2pro-specific struct; it is configured through the
    # RSP2-family fields (datasheet v2.2: electrically identical to RSP2).
    "rsp2pro": DeviceProfile(
        key="rsp2pro", label="RSP2/RSP2pro (hwVer 2)", hwvers=(RSP2_ID,),
        sample_scale=32768.0,
        lna_max_2m=RSP2_LNA_STATES_0_420MHZ - 1,
        lna_max_70cm=RSP2_LNA_STATES_420_1000MHZ - 1,
        min_api=3.08),
}


def family_of(hwver: int) -> Optional[str]:
    """Classify a live device by hwVer. Match by hwVer, NOT serial prefix —
    serials do not reliably carry the model name (the connected RSP2pro
    serials as '1710060810', no prefix)."""
    for prof in PROFILES.values():
        if hwver in prof.hwvers:
            return prof.key
    return None


class SDRPlayError(RuntimeError):
    pass


# --------------------------------------------------------------------------
# struct tree (mirrors the 3.11 headers exactly; ctypes computes the same
# layout MSVC does)
# --------------------------------------------------------------------------

class FsFreqT(ct.Structure):
    _fields_ = [("fsHz", ct.c_double),
                ("syncUpdate", ct.c_ubyte),
                ("reCal", ct.c_ubyte)]

class SyncUpdateT(ct.Structure):
    _fields_ = [("sampleNum", ct.c_uint), ("period", ct.c_uint)]

class ResetFlagsT(ct.Structure):
    _fields_ = [("resetGainUpdate", ct.c_ubyte),
                ("resetRfUpdate", ct.c_ubyte),
                ("resetFsUpdate", ct.c_ubyte)]

class Rsp1aParamsT(ct.Structure):
    _fields_ = [("rfNotchEnable", ct.c_ubyte), ("rfDabNotchEnable", ct.c_ubyte)]

class Rsp2ParamsT(ct.Structure):
    _fields_ = [("extRefOutputEn", ct.c_ubyte)]

class RspDuoParamsT(ct.Structure):
    _fields_ = [("extRefOutputEn", ct.c_int)]

class RspDxParamsT(ct.Structure):
    _fields_ = [("hdrEnable", ct.c_ubyte),
                ("biasTEnable", ct.c_ubyte),
                ("antennaSel", ct.c_int),
                ("rfNotchEnable", ct.c_ubyte),
                ("rfDabNotchEnable", ct.c_ubyte)]

class DevParamsT(ct.Structure):
    _fields_ = [("ppm", ct.c_double),
                ("fsFreq", FsFreqT),
                ("syncUpdate", SyncUpdateT),
                ("resetFlags", ResetFlagsT),
                ("mode", ct.c_int),
                ("samplesPerPkt", ct.c_uint),
                ("rsp1aParams", Rsp1aParamsT),
                ("rsp2Params", Rsp2ParamsT),
                ("rspDuoParams", RspDuoParamsT),
                ("rspDxParams", RspDxParamsT)]

class GainValuesT(ct.Structure):
    _fields_ = [("curr", ct.c_float), ("max", ct.c_float), ("min", ct.c_float)]

class GainT(ct.Structure):
    _fields_ = [("gRdB", ct.c_int),
                ("LNAstate", ct.c_ubyte),
                ("syncUpdate", ct.c_ubyte),
                ("minGr", ct.c_int),
                ("gainVals", GainValuesT)]

class RfFreqT(ct.Structure):
    _fields_ = [("rfHz", ct.c_double), ("syncUpdate", ct.c_ubyte)]

class DcOffsetTunerT(ct.Structure):
    _fields_ = [("dcCal", ct.c_ubyte), ("speedUp", ct.c_ubyte),
                ("trackTime", ct.c_int), ("refreshRateTime", ct.c_int)]

class TunerParamsT(ct.Structure):
    _fields_ = [("bwType", ct.c_int),
                ("ifType", ct.c_int),
                ("loMode", ct.c_int),
                ("gain", GainT),
                ("rfFreq", RfFreqT),
                ("dcOffsetTuner", DcOffsetTunerT)]

class DcOffsetT(ct.Structure):
    _fields_ = [("DCenable", ct.c_ubyte), ("IQenable", ct.c_ubyte)]

class DecimationT(ct.Structure):
    _fields_ = [("enable", ct.c_ubyte),
                ("decimationFactor", ct.c_ubyte),
                ("wideBandSignal", ct.c_ubyte)]

class AgcT(ct.Structure):
    _fields_ = [("enable", ct.c_int),
                ("setPoint_dBfs", ct.c_int),
                ("attack_ms", ct.c_ushort),
                ("decay_ms", ct.c_ushort),
                ("decay_delay_ms", ct.c_ushort),
                ("decay_threshold_dB", ct.c_ushort),
                ("syncUpdate", ct.c_int)]

class ControlParamsT(ct.Structure):
    _fields_ = [("dcOffset", DcOffsetT),
                ("decimation", DecimationT),
                ("agc", AgcT),
                ("adsbMode", ct.c_int)]

class Rsp1aTunerParamsT(ct.Structure):
    _fields_ = [("biasTEnable", ct.c_ubyte)]

class Rsp2TunerParamsT(ct.Structure):
    _fields_ = [("biasTEnable", ct.c_ubyte),
                ("amPortSel", ct.c_int),
                ("antennaSel", ct.c_int),
                ("rfNotchEnable", ct.c_ubyte)]

class RspDuoTunerParamsT(ct.Structure):
    _fields_ = [("biasTEnable", ct.c_ubyte),
                ("tuner1AmPortSel", ct.c_int),
                ("tuner1AmNotchEnable", ct.c_ubyte),
                ("rfNotchEnable", ct.c_ubyte),
                ("rfDabNotchEnable", ct.c_ubyte)]

class RspDxTunerParamsT(ct.Structure):
    _fields_ = [("hdrBw", ct.c_int)]

class RxChannelParamsT(ct.Structure):
    _fields_ = [("tunerParams", TunerParamsT),
                ("ctrlParams", ControlParamsT),
                ("rsp1aTunerParams", Rsp1aTunerParamsT),
                ("rsp2TunerParams", Rsp2TunerParamsT),
                ("rspDuoTunerParams", RspDuoTunerParamsT),
                ("rspDxTunerParams", RspDxTunerParamsT)]

class DeviceT(ct.Structure):
    _fields_ = [("SerNo", ct.c_char * 64),
                ("hwVer", ct.c_ubyte),
                ("tuner", ct.c_int),
                ("rspDuoMode", ct.c_int),
                ("valid", ct.c_ubyte),
                ("rspDuoSampleFreq", ct.c_double),
                ("dev", ct.c_void_p)]

class DeviceParamsT(ct.Structure):
    _fields_ = [("devParams", ct.POINTER(DevParamsT)),
                ("rxChannelA", ct.POINTER(RxChannelParamsT)),
                ("rxChannelB", ct.POINTER(RxChannelParamsT))]

class ErrorInfoT(ct.Structure):
    _fields_ = [("file", ct.c_char * 256),
                ("function", ct.c_char * 256),
                ("line", ct.c_int),
                ("message", ct.c_char * 1024)]

class StreamCbParamsT(ct.Structure):
    _fields_ = [("firstSampleNum", ct.c_uint),
                ("grChanged", ct.c_int),
                ("rfChanged", ct.c_int),
                ("fsChanged", ct.c_int),
                ("numSamples", ct.c_uint)]

class GainCbParamT(ct.Structure):
    _fields_ = [("gRdB", ct.c_uint), ("lnaGRdB", ct.c_uint), ("currGain", ct.c_double)]

class EventParamsT(ct.Union):
    _fields_ = [("gainParams", GainCbParamT),
                ("powerOverloadType", ct.c_int),
                ("modeChangeType", ct.c_int)]


# callback prototypes (3.11 signatures)
STREAM_CB = ct.CFUNCTYPE(None,
                         ct.POINTER(ct.c_int16),   # xi
                         ct.POINTER(ct.c_int16),   # xq
                         ct.POINTER(StreamCbParamsT),
                         ct.c_uint,                # numSamples
                         ct.c_uint,                # reset
                         ct.c_void_p)              # cbContext

EVENT_CB = ct.CFUNCTYPE(None,
                        ct.c_int,                  # eventId
                        ct.c_int,                  # tuner
                        ct.POINTER(EventParamsT),  # params
                        ct.c_void_p)               # cbContext

class CallbackFnsT(ct.Structure):
    _fields_ = [("StreamACbFn", STREAM_CB),
                ("StreamBCbFn", STREAM_CB),
                ("EventCbFn", EVENT_CB)]


# --------------------------------------------------------------------------
# DLL location + class
# --------------------------------------------------------------------------

def find_dll() -> Path:
    with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                        r"SOFTWARE\SDRplay\Service\API") as k:
        install_dir = winreg.QueryValueEx(k, "Install_Dir")[0]
    for cand in (Path(install_dir) / "x64" / "sdrplay_api.dll",
                 Path(install_dir) / "sdrplay_api.dll"):
        if cand.exists():
            return cand
    raise SDRPlayError(f"sdrplay_api.dll not found under {install_dir}")


class SDRPlay:
    """Minimal RSPdx / RSP2-RSP2pro driver: probe -> configure -> stream
    -> live-update -> close. The selected device's DeviceProfile (sample
    scale, LNA limits, antenna path) is available as self.profile."""

    def __init__(self, dll_path: Optional[str] = None):
        dll = Path(dll_path) if dll_path else find_dll()
        self._lib = ct.CDLL(str(dll))
        L = self._lib
        L.sdrplay_api_Open.restype = ct.c_int
        L.sdrplay_api_Close.restype = ct.c_int
        L.sdrplay_api_ApiVersion.argtypes = [ct.POINTER(ct.c_float)]
        L.sdrplay_api_ApiVersion.restype = ct.c_int
        L.sdrplay_api_LockDeviceApi.restype = ct.c_int
        L.sdrplay_api_UnlockDeviceApi.restype = ct.c_int
        L.sdrplay_api_GetDevices.argtypes = [ct.POINTER(DeviceT), ct.POINTER(ct.c_uint), ct.c_uint]
        L.sdrplay_api_GetDevices.restype = ct.c_int
        L.sdrplay_api_SelectDevice.argtypes = [ct.POINTER(DeviceT)]
        L.sdrplay_api_SelectDevice.restype = ct.c_int
        L.sdrplay_api_ReleaseDevice.argtypes = [ct.POINTER(DeviceT)]
        L.sdrplay_api_ReleaseDevice.restype = ct.c_int
        L.sdrplay_api_GetErrorString.argtypes = [ct.c_int]
        L.sdrplay_api_GetErrorString.restype = ct.c_char_p
        L.sdrplay_api_GetLastError.argtypes = [ct.POINTER(DeviceT)]
        L.sdrplay_api_GetLastError.restype = ct.POINTER(ErrorInfoT)
        L.sdrplay_api_GetDeviceParams.argtypes = [ct.c_void_p, ct.POINTER(ct.POINTER(DeviceParamsT))]
        L.sdrplay_api_GetDeviceParams.restype = ct.c_int
        L.sdrplay_api_Init.argtypes = [ct.c_void_p, ct.POINTER(CallbackFnsT), ct.c_void_p]
        L.sdrplay_api_Init.restype = ct.c_int
        L.sdrplay_api_Uninit.argtypes = [ct.c_void_p]
        L.sdrplay_api_Uninit.restype = ct.c_int
        L.sdrplay_api_Update.argtypes = [ct.c_void_p, ct.c_int, ct.c_uint, ct.c_uint]
        L.sdrplay_api_Update.restype = ct.c_int

        self._check(L.sdrplay_api_Open())
        ver = ct.c_float()
        self._check(L.sdrplay_api_ApiVersion(ct.byref(ver)))
        self.api_version = ver.value
        if self.api_version < min(p.min_api for p in PROFILES.values()):
            raise SDRPlayError(f"API {self.api_version} too old "
                               f"(need >= 3.06; RSP2/RSP2pro need >= 3.08)")

        self._device = DeviceT()
        self._dev_handle: Optional[ct.c_void_p] = None
        self._params = ct.POINTER(DeviceParamsT)()
        self.profile: Optional[DeviceProfile] = None   # set by find_device()
        self._cbs: list = []          # keep ctypes callback objects alive
        self._cb_ctx: int = 0
        self.on_iq: Optional[Callable] = None
        self.on_event: Optional[Callable] = None
        self.total_gain_db: float = float("nan")

    # ------------------------------------------------------------------
    def _check(self, err: int) -> None:
        if err != Err.Success:
            msg = self._lib.sdrplay_api_GetErrorString(ct.c_int(err))
            raise SDRPlayError(f"SDRplay API error {err}: {msg}")

    def error_string(self, err: int) -> str:
        return self._lib.sdrplay_api_GetErrorString(ct.c_int(err)).decode()

    # ------------------------------------------------------------------
    def find_device(self, radio: str = "auto") -> DeviceT:
        """Probe the bus and select a device.

        radio: "auto" (default), "rspdx" or "rsp2pro" (case-insensitive).
          * auto: the sole valid device is used whatever its family; if
            several are connected, the RSPdx (this project's baseline)
            is preferred over the RSP2/RSP2pro — override with an
            explicit choice.
          * rspdx / rsp2pro: only a device of that family is accepted.
        Selection matches by hwVer (family_of), never by serial prefix.
        On return, self.profile describes the selected device.
        """
        if isinstance(radio, int):          # backward compat: raw id
            radio = {RSPDX_ID: "rspdx", RSP2_ID: "rsp2pro"}.get(radio, str(radio))
        radio = str(radio).lower()
        if radio != "auto" and radio not in PROFILES:
            raise SDRPlayError(f"unknown radio '{radio}' (choose auto/rspdx/rsp2pro)")

        self._check(self._lib.sdrplay_api_LockDeviceApi())
        try:
            devices = (DeviceT * 16)()
            n = ct.c_uint()
            self._check(self._lib.sdrplay_api_GetDevices(devices, ct.byref(n), 16))
        finally:
            self._check(self._lib.sdrplay_api_UnlockDeviceApi())
        # GetDevices leaves dev=NULL; the handle is filled by SelectDevice
        valid = [d for d in devices if d.valid]
        if not valid:
            raise SDRPlayError("no SDRplay devices visible (is the unit connected & service running?)")
        classified: list[tuple[str, DeviceT]] = []
        for d in valid:
            fam = family_of(d.hwVer)
            classified.append((fam if fam is not None else f"unknown-hwVer{d.hwVer}", d))

        if radio == "auto":
            known = [c for c in classified if c[0] in PROFILES]
            if len(valid) == 1:
                pick = known[0] if known else None
                if pick is None:
                    raise SDRPlayError(
                        f"single device is hwVer {valid[0].hwVer} (serial {self._ser(valid[0])}) — "
                        f"not a supported family (RSPdx=4, RSP2/RSP2pro=2)")
            else:
                # several devices: prefer the project baseline (RSPdx)
                pick = next((c for c in known if c[0] == "rspdx"),
                            next((c for c in known if c[0] == "rsp2pro"), None))
                if pick is None:
                    raise SDRPlayError(
                        f"several devices connected, none supported: "
                        f"{[(self._ser(d), h) for h, d in classified]}")
                print(f"[probe] {len(valid)} devices connected, auto-selecting "
                      f"{PROFILES[pick[0]].label} — use --radio to override")
        else:
            pick = next((c for c in classified if c[0] == radio), None)
            if pick is None:
                present = [(self._ser(d), h) for h, d in classified]
                raise SDRPlayError(f"no {PROFILES[radio].label} connected; present: {present}")

        fam, dev = pick
        prof = PROFILES[fam]
        if self.api_version < prof.min_api:
            raise SDRPlayError(f"API {self.api_version} too old for {prof.label} (need >= {prof.min_api})")
        self.profile = prof
        self._device = dev
        self._check(self._lib.sdrplay_api_SelectDevice(ct.byref(self._device)))
        self._dev_handle = ct.c_void_p(self._device.dev)
        return self._device

    @staticmethod
    def _ser(d: DeviceT) -> str:
        return d.SerNo.split(b"\x00")[0].decode(errors="replace")

    @staticmethod
    def _device_id(d: DeviceT) -> int:
        # hwVer encodes the family for RSP devices; the service also exposes
        # the id via the device list ordering — use the serial prefix heuristic
        # plus hwVer. RSPdx hwVer is 0x04 (4) in API v3.
        return int.from_bytes(d.hwVer.to_bytes(1, "little"), "little")

    # ------------------------------------------------------------------
    def get_params(self) -> None:
        self._check(self._lib.sdrplay_api_GetDeviceParams(
            self._dev_handle, ct.byref(self._params)))
        self.dev_params: DevParamsT = self._params.contents.devParams.contents
        self.rx_a: RxChannelParamsT = self._params.contents.rxChannelA.contents

    @property
    def tuner(self) -> TunerParamsT:
        return self.rx_a.tunerParams

    @property
    def ctrl(self) -> ControlParamsT:
        return self.rx_a.ctrlParams

    @property
    def rspdx(self) -> RspDxParamsT:
        return self.dev_params.rspDxParams

    @property
    def rsp2_tuner(self) -> Rsp2TunerParamsT:
        return self.rx_a.rsp2TunerParams

    # ------------------------------------------------------------------
    def _set_antenna(self, antenna: int) -> None:
        """antenna: 0 = port A, 1 = port B (per-family struct routing —
        writing to the wrong family's struct is silently ignored by the
        service)."""
        if self.profile is None:
            raise SDRPlayError("no device selected (call find_device first)")
        if self.profile.key == "rspdx":
            self.rspdx.antennaSel = int(antenna)      # 0=A, 1=B, 2=C
        else:
            if antenna == 0:
                self.rsp2_tuner.antennaSel = RSP2_ANT_A
            elif antenna == 1:
                self.rsp2_tuner.antennaSel = RSP2_ANT_B
            else:
                raise SDRPlayError(
                    f"antenna {antenna} not available on {self.profile.label} (0=port A, 1=port B)")

    # ------------------------------------------------------------------
    def configure(self, fs_hz: float, rf_hz: float, bw_khz: int,
                  lna: int, gr_db: int, agc: int = AGC_50HZ,
                  agc_setpoint_dbfs: int = -20,
                  antenna: int = 0) -> None:
        """Modify the parameter tree (call before start())."""
        self.dev_params.fsFreq.fsHz = float(fs_hz)
        t = self.tuner
        t.rfFreq.rfHz = float(rf_hz)
        t.ifType = IF_ZERO
        t.bwType = int(bw_khz)
        t.gain.gRdB = int(gr_db)
        t.gain.LNAstate = int(lna)
        c = self.ctrl
        c.agc.enable = int(agc)
        c.agc.setPoint_dBfs = int(agc_setpoint_dbfs)
        c.dcOffset.DCenable = 1
        c.dcOffset.IQenable = 1
        self._set_antenna(antenna)

    def start(self, on_iq: Callable, on_event: Optional[Callable] = None) -> None:
        """Register callbacks and start streaming (Init)."""
        self.on_iq = on_iq
        self.on_event = on_event
        self._cb_ctx = 1

        def _stream(xi, xq, params, num_samples, reset, ctx):
            try:
                n = int(num_samples)
                iq = ct.POINTER(ct.c_int16)
                # fast path: view the raw buffers, no copy
                i = ct.cast(xi, ct.POINTER(ct.c_int16))
                q = ct.cast(xq, ct.POINTER(ct.c_int16))
                self.on_iq(i, q, n, int(reset))
            except Exception:
                import traceback; traceback.print_exc()

        def _event(event_id, tuner, params, ctx):
            try:
                if event_id == EVT_GAIN_CHANGE and params:
                    self.total_gain_db = float(params.contents.gainParams.currGain)
                if self.on_event:
                    self.on_event(int(event_id), int(tuner),
                                  float(params.contents.gainParams.currGain) if params else float("nan"))
            except Exception:
                import traceback; traceback.print_exc()

        self._stream_cb = STREAM_CB(_stream)
        self._event_cb = EVENT_CB(_event)
        self._cbs = [self._stream_cb, self._event_cb]
        fns = CallbackFnsT(self._stream_cb, self._stream_cb, self._event_cb)
        self._check(self._lib.sdrplay_api_Init(self._dev_handle, ct.byref(fns), ct.c_void_p(self._cb_ctx)))

    # ------------------------------------------------------------------
    def update(self, rf_hz: Optional[float] = None, gr_db: Optional[int] = None,
               lna: Optional[int] = None, bw_khz: Optional[int] = None,
               antenna: Optional[int] = None,
               agc_setpoint_dbfs: Optional[int] = None,
               fs_hz: Optional[float] = None) -> None:
        """Apply live parameter changes (post-Init)."""
        reason = 0
        ext1 = 0
        if fs_hz is not None:
            self.dev_params.fsFreq.fsHz = float(fs_hz)
            reason |= UPD_DEV_FS
        t = self.tuner
        if rf_hz is not None:
            t.rfFreq.rfHz = float(rf_hz)
            reason |= UPD_TUNER_FRF
        if gr_db is not None or lna is not None:
            if gr_db is not None:
                t.gain.gRdB = int(gr_db)
            if lna is not None:
                t.gain.LNAstate = int(lna)
            reason |= UPD_TUNER_GR
        if bw_khz is not None:
            t.bwType = int(bw_khz)
            reason |= UPD_TUNER_BWTYPE
        if agc_setpoint_dbfs is not None:
            self.ctrl.agc.setPoint_dBfs = int(agc_setpoint_dbfs)
            reason |= UPD_CTRL_AGC
        if antenna is not None:
            self._set_antenna(antenna)
            if self.profile is not None and self.profile.key == "rspdx":
                ext1 |= UPD1_RSPDX_ANTENNA
            else:
                reason |= UPD_RSP2_ANTENNA
        if reason or ext1:
            self._check(self._lib.sdrplay_api_Update(
                self._dev_handle, TUNER_A, reason, ext1))

    def ack_overload(self) -> None:
        self._check(self._lib.sdrplay_api_Update(
            self._dev_handle, TUNER_A, UPD_CTRL_OVERLOAD_ACK, 0))

    # ------------------------------------------------------------------
    def stop(self) -> None:
        if self._dev_handle:
            try:
                self._lib.sdrplay_api_Uninit(self._dev_handle)
            except Exception:
                pass
            self._lib.sdrplay_api_ReleaseDevice(ct.byref(self._device))
            self._dev_handle = None
        self._lib.sdrplay_api_Close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.stop()
