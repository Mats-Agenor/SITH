from __future__ import annotations

from pathlib import Path

from .config import Settings
from .log import Logger
from .paths import FieldPaths


def _stilts(settings: Settings, *args: str) -> list[str]:
    return [settings.java_path, "-jar", str(settings.stilts_path), *args]


def crossmatch_splus_gaiaxp(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    master_file: Path,
    gaia_file: Path,
    log: Logger,
    force: bool = False,
) -> Path:
    out = paths.gaiaxp / f"{field}_SPLUS_dual_gaiaXPSP_phot.fits"
    if out.exists() and force:
        out.unlink()
    if out.exists():
        return out

    # Exact SPLUSCALIB Gaia XPSP matcher:
    # tmatch2, matcher=sky, params=3 arcsec, join=1and2.
    args = [
        "tmatch2",
        f"ifmt1=fits",
        f"in1={master_file}",
        "values1=RAJ2000 DEJ2000",
        f"ifmt2=fits",
        f"in2={gaia_file}",
        "values2=RA_ICRS DE_ICRS",
        f"out={out}",
        "ofmt=fits",
        "matcher=sky",
        "params=3",
        "join=1and2",
    ]
    log.write("Crossmatching S-PLUS PStotal x precomputed GaiaXP with original 3 arcsec STILTS rule")
    log.command(_stilts(settings, *args))
    if not out.is_file():
        raise RuntimeError(f"STILTS did not create {out}")
    return out
