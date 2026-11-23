from __future__ import annotations

import re
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import ICRS, SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy_healpix import HEALPix

from .config import SYNTHETIC_MAG_COLUMNS, Settings
from .log import Logger
from .paths import FieldPaths

_RANGE_RE = re.compile(r"_(\d+)-(\d+)_synthetic\.fits$", re.IGNORECASE)


def parse_catalog_range(path: Path) -> tuple[int, int] | None:
    m = _RANGE_RE.search(path.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def _field_cone(detection: Path, margin_arcmin: float) -> tuple[SkyCoord, u.Quantity]:
    with fits.open(detection, memmap=True) as hdul:
        hdu = next(h for h in hdul if h.data is not None and getattr(h.data, "ndim", 0) == 2)
        ny, nx = hdu.data.shape
        w = WCS(hdu.header)
        center = w.pixel_to_world((nx - 1) / 2.0, (ny - 1) / 2.0).icrs
        corners = w.calc_footprint(axes=(nx, ny))
    cc = SkyCoord(corners[:, 0] * u.deg, corners[:, 1] * u.deg, frame="icrs")
    radius = np.max(center.separation(cc)) + margin_arcmin * u.arcmin
    return center, radius


def _needed_level8_pixels(center: SkyCoord, radius: u.Quantity) -> np.ndarray:
    hp = HEALPix(nside=2**8, order="nested", frame=ICRS())
    return np.asarray(hp.cone_search_skycoord(center, radius), dtype=np.int64)


#====================
#Etapa GaiaXP substituida integralmente: o SITH 2.0 NAO baixa espectros e NAO
#executa GaiaXPy. Ele le a fotometria sintetica S-PLUS ja calculada nos FITS
#XpContinuousMeanSpectrum_*_synthetic.fits indicados por gaiaxp_path no .conf.
#Todo o procedimento posterior de crossmatch e zero point volta a seguir o
#SPLUSCALIB.
#M. Agenor
#====================
def select_gaia_files(settings: Settings, detection: Path, log: Logger) -> tuple[list[Path], SkyCoord, u.Quantity]:
    if settings.gaia_healpix_level != 8:
        raise ValueError(
            "The filename ranges XpContinuousMeanSpectrum_XXXXXX-YYYYYY_synthetic.fits "
            "are interpreted as Gaia source_id HEALPix level 8 ranges; set gaia_healpix_level 8."
        )
    center, radius = _field_cone(detection, settings.gaia_footprint_margin_arcmin)
    pixels = _needed_level8_pixels(center, radius)
    selected: list[Path] = []
    for p in sorted(settings.gaiaxp_path.glob(settings.gaia_filename_glob)):
        r = parse_catalog_range(p)
        if r is None:
            continue
        lo, hi = r
        if np.any((pixels >= lo) & (pixels <= hi)):
            selected.append(p)
    if not selected:
        raise RuntimeError(
            f"No Gaia synthetic filename range overlaps field pixels {pixels.min()}..{pixels.max()}"
        )
    log.write(
        f"Gaia synthetic selection: center=({center.ra.deg:.6f},{center.dec.deg:.6f}), "
        f"radius={radius.to_value(u.deg):.3f} deg, level8 pixels={pixels.min()}..{pixels.max()}"
    )
    for p in selected:
        log.write(f"  Gaia synthetic file: {p.name}")
    return selected, center, radius


def _table_hdu(path: Path, needed: set[str]):
    hdul = fits.open(path, memmap=True)
    for hdu in hdul:
        names = set(getattr(getattr(hdu, "data", None), "names", []) or [])
        if needed.issubset(names):
            return hdul, hdu
    hdul.close()
    raise KeyError(f"No table HDU in {path.name} contains required columns: {sorted(needed)}")


def _read_subset(path: Path, center: SkyCoord, radius: u.Quantity, filters: list[str]) -> Table | None:
    needed = {"source_id", "ra", "dec", *[SYNTHETIC_MAG_COLUMNS[f] for f in filters]}
    hdul, hdu = _table_hdu(path, needed)
    try:
        data = hdu.data
        ra = np.asarray(data["ra"], dtype=float)
        dec = np.asarray(data["dec"], dtype=float)
        sky = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
        mask = sky.separation(center) <= radius
        idx = np.flatnonzero(mask)
        if len(idx) == 0:
            return None
        out = Table()
        out["source_id"] = np.asarray(data["source_id"])[idx]
        # Names deliberately mimic the GaiaXPy product consumed by SPLUSCALIB.
        out["RA_ICRS"] = ra[idx]
        out["DE_ICRS"] = dec[idx]
        for filt in filters:
            out[f"SPLUS_{filt}_mod"] = np.asarray(data[SYNTHETIC_MAG_COLUMNS[filt]], dtype=float)[idx]
        return out
    finally:
        hdul.close()


def build_gaia_synthetic_subset(settings: Settings, paths: FieldPaths, detection: Path, log: Logger, force: bool = False) -> Path:
    out = paths.gaiaxp / "gaiaXPSP_phot_precomputed.fits"
    if out.is_file() and not force:
        return out
    files, center, radius = select_gaia_files(settings, detection, log)
    chunks: list[Table] = []
    for p in files:
        log.write(f"Reading precomputed synthetic photometry: {p.name}")
        t = _read_subset(p, center, radius, settings.filters)
        if t is not None and len(t):
            chunks.append(t)
    if not chunks:
        raise RuntimeError("Selected Gaia synthetic files contain no sources in the field region")
    tab = chunks[0] if len(chunks) == 1 else vstack(chunks, metadata_conflicts="silent")
    _, idx = np.unique(np.asarray(tab["source_id"]), return_index=True)
    tab = tab[np.sort(idx)]
    tab.write(out, overwrite=True)
    log.write(f"Gaia synthetic subset: {len(tab)} sources -> {out}")
    return out
