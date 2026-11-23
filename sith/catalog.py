from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .config import Settings
from .paths import FieldPaths
from .zeropoint import ZeroPoint


def write_calibrated_master(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    master_file: Path,
    zps: dict[str, ZeroPoint],
) -> Path:
    """
    Apply the same net SPLUSCALIB calibration to PStotal magnitudes.

    SPLUSCALIB final.zp = inst_zp + GaiaXP_offset and then applies
    (final.zp - sex_mag_zp), hence the net correction is GaiaXP_offset.
    """
    tab = Table.read(master_file)
    for filt in settings.filters:
        col = np.asarray(tab[f"SPLUS_{filt}"], dtype=float).copy()
        good = col != 99
        col[good] += zps[filt].final_zp - settings.inst_zp
        tab[f"SPLUS_{filt}"] = col
    tab.meta["SITHVER"] = "2.0.0"
    tab.meta["FIELD"] = field
    tab.meta["INSTZP"] = settings.inst_zp
    out = paths.calibration / f"{field}_SITH2_calibrated_master.fits"
    tab.write(out, overwrite=True)
    return out


def write_calibrated_filter_catalogs(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    corrected_catalogs: dict[str, Path],
    zps: dict[str, ZeroPoint],
) -> dict[str, Path]:
    """Keep the full SExtractor/aperture-corrected catalog and append calibrated PStotal."""
    outdir = paths.calibration / "dual"
    outdir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}
    for filt in settings.filters:
        hdul = fits.open(corrected_catalogs[filt])
        data = hdul[1].data
        offset = zps[filt].final_zp - settings.inst_zp
        cols = [c for c in data.columns]

        mag = np.array(data["MAG_PStotal"], dtype=float, copy=True)
        good = mag != 99
        mag[good] += offset
        cols.append(fits.Column(name="MAG_PStotal_CAL", format="1E", array=mag))

        # Same net zero-point correction used by SPLUSCALIB for the standard
        # SExtractor magnitude estimators; appended as calibrated copies so the
        # original instrumental columns remain available for diagnostics.
        for name in ["MAG_AUTO", "MAG_PETRO", "MAG_ISO"]:
            arr = np.array(data[name], dtype=float, copy=True)
            ok = arr != 99
            arr[ok] += offset
            cols.append(fits.Column(name=f"{name}_CAL", format="1E", array=arr))

        aper = np.array(data["MAG_APER"], dtype=float, copy=True)
        ok = aper != 99
        aper[ok] += offset
        cols.append(fits.Column(name="MAG_APER_CAL", format=f"{aper.shape[1]}E", array=aper))
        out = outdir / f"{field}_{filt}_dual.fits"
        fits.BinTableHDU.from_columns(cols).writeto(out, overwrite=True)
        hdul.close()
        result[filt] = out
    return result
