from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .external import discover_executable, discover_sextractor_config_dir, stilts_jar


FILTERS_DEFAULT = [
    "U", "F378", "F395", "F410", "F430", "G",
    "F515", "R", "F660", "I", "F861", "Z",
]

DETECTION_FILTERS_DEFAULT = ["G", "R", "I", "Z"]

# Exact PHOT_APERTURES from SPLUSCALIB SExsplus.sex.
APERTURES_SPLUSCALIB = [
    1.81818, 3.63636, 5.45455, 7.27273, 10.90909, 14.54545,
    18.18182, 21.09091, 24.00000, 26.90909, 29.81818, 32.72727,
    35.63636, 38.54545, 41.45455, 44.36364, 47.27273, 50.18182,
    53.09091, 56.00000, 58.90909, 61.81818, 64.72727, 67.63636,
    70.54545, 73.45455, 76.36364, 79.27273, 82.18182, 85.09091,
    88.00000, 90.90909,
]

SYNTHETIC_MAG_COLUMNS = {
    "U": "USER_CTIO_S_PLUS_mag_u",
    "F378": "USER_CTIO_S_PLUS_mag_F378",
    "F395": "USER_CTIO_S_PLUS_mag_F395",
    "F410": "USER_CTIO_S_PLUS_mag_F410",
    "F430": "USER_CTIO_S_PLUS_mag_F430",
    "G": "USER_CTIO_S_PLUS_mag_g",
    "F515": "USER_CTIO_S_PLUS_mag_F515",
    "R": "USER_CTIO_S_PLUS_mag_r",
    "F660": "USER_CTIO_S_PLUS_mag_F660",
    "I": "USER_CTIO_S_PLUS_mag_i",
    "F861": "USER_CTIO_S_PLUS_mag_F861",
    "Z": "USER_CTIO_S_PLUS_mag_z",
}


def read_conf(path: str | Path) -> dict[str, str]:
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    out: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" in line:
            left, right = line.split("=", 1)
            if left.strip() and " " not in left.strip():
                out[left.strip()] = right.strip()
                continue
        parts = line.split(None, 1)
        out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


def _bool(v: str | None, default: bool) -> bool:
    if v is None:
        return default
    x = v.strip().lower()
    if x in {"1", "true", "yes", "y", "on"}:
        return True
    if x in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"Invalid boolean: {v}")


def _list(v: str | None, default: list[str]) -> list[str]:
    if not v:
        return list(default)
    return [x.strip() for x in v.split(",") if x.strip()]


def _float_list(v: str | None, default: list[float]) -> list[float]:
    if not v:
        return list(default)
    return [float(x.strip()) for x in v.split(",") if x.strip()]


def _expand(v: list[float], n: int, name: str) -> list[float]:
    if len(v) == 1:
        return v * n
    if len(v) != n:
        raise ValueError(f"{name}: expected one value or {n}, got {len(v)}")
    return v


@dataclass(frozen=True)
class Settings:
    config_file: Path
    input_path: Path
    output_path: Path
    gaiaxp_path: Path
    image_pattern: str
    weight_pattern: str
    filters: list[str]
    detection_filters: list[str]
    use_weight: bool
    inst_zp: float
    apercorr_diameter: float
    apercorr_s2ncut: tuple[float, float]
    apercorr_starcut: float
    apercorr_max_aperture: float
    gaia_mag_low: list[float]
    gaia_mag_up: list[float]
    gaia_filename_glob: str
    gaia_healpix_level: int
    gaia_footprint_margin_arcmin: float
    swarp_path: str
    sextractor_path: str
    java_path: str
    funpack_path: str | None
    sextractor_config_dir: Path
    stilts_path: Path
    keep_work_images: bool
    make_plots: bool

    @classmethod
    def from_file(cls, path: str | Path) -> "Settings":
        p = Path(path).expanduser().resolve()
        raw = read_conf(p)

        def required(k: str) -> str:
            v = raw.get(k, "").strip()
            if not v:
                raise ValueError(f"Missing required configuration parameter: {k}")
            return v

        filters = _list(raw.get("filters"), FILTERS_DEFAULT)
        unknown = [f for f in filters if f not in SYNTHETIC_MAG_COLUMNS]
        if unknown:
            raise ValueError(f"Unknown filters: {unknown}")

        det = _list(raw.get("detection_filters") or raw.get("detection_image"), DETECTION_FILTERS_DEFAULT)
        for f in det:
            if f not in filters:
                raise ValueError(f"Detection filter {f} is not in filters")

        low = _expand(
            _float_list(raw.get("gaiaxpsp_zp_fitting_mag_low") or raw.get("zp_mag_low"), [14.0]),
            len(filters),
            "gaiaxpsp_zp_fitting_mag_low",
        )
        up = _expand(
            _float_list(raw.get("gaiaxpsp_zp_fitting_mag_up") or raw.get("zp_mag_up"), [19.0]),
            len(filters),
            "gaiaxpsp_zp_fitting_mag_up",
        )
        s2n = _float_list(raw.get("apercorr_s2ncut"), [30.0, 1000.0])
        if len(s2n) != 2:
            raise ValueError("apercorr_s2ncut must have exactly two values")

        funpack: str | None
        try:
            funpack = discover_executable(raw.get("funpack_path", "auto"), ["funpack"], "funpack")
        except FileNotFoundError:
            funpack = None

        return cls(
            config_file=p,
            input_path=Path(required("input_path")).expanduser().resolve(),
            output_path=Path(required("output_path")).expanduser().resolve(),
            gaiaxp_path=Path(required("gaiaxp_path")).expanduser(),
            image_pattern=raw.get("image_pattern", "{field}_{filter}.fits.fz"),
            weight_pattern=raw.get("weight_pattern", "{field}_{filter}.weight.fits.fz"),
            filters=filters,
            detection_filters=det,
            use_weight=_bool(raw.get("use_weight"), True),
            inst_zp=float(raw.get("inst_zp", "20")),
            apercorr_diameter=float(raw.get("apercorr_diameter", "5.45454")),
            apercorr_s2ncut=(float(s2n[0]), float(s2n[1])),
            apercorr_starcut=float(raw.get("apercorr_starcut", "0.9")),
            apercorr_max_aperture=float(raw.get("apercorr_max_aperture", "36.363636")),
            gaia_mag_low=low,
            gaia_mag_up=up,
            gaia_filename_glob=raw.get("gaia_filename_glob", "XpContinuousMeanSpectrum_*_synthetic.fits"),
            gaia_healpix_level=int(raw.get("gaia_healpix_level", "8")),
            gaia_footprint_margin_arcmin=float(raw.get("gaia_footprint_margin_arcmin", "3.0")),
            swarp_path=discover_executable(raw.get("swarp_path", "auto"), ["SWarp", "swarp"], "SWarp"),
            sextractor_path=discover_executable(raw.get("sextractor_path", "auto"), ["source-extractor", "sextractor", "sex"], "SExtractor"),
            java_path=discover_executable(raw.get("java_path", "auto"), ["java"], "Java"),
            funpack_path=funpack,
            sextractor_config_dir=discover_sextractor_config_dir(raw.get("sextractor_config_dir", "auto")),
            stilts_path=stilts_jar(),
            keep_work_images=_bool(raw.get("keep_work_images"), True),
            make_plots=_bool(raw.get("make_plots"), True),
        )

    def mag_limits(self, filt: str) -> tuple[float, float]:
        i = self.filters.index(filt)
        return self.gaia_mag_low[i], self.gaia_mag_up[i]
