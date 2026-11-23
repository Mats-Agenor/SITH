from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.table import Table

from .config import Settings
from .log import Logger
from .paths import FieldPaths


#====================
#Estrutura substituida para restaurar o workflow do SITH 1. Como a fotometria
#foi executada em dual mode, NUMBER e a identidade da fonte e e usado
#diretamente para combinar os 12 filtros. O tmatchn espacial entre catalogos
#dual foi removido; nao se redescobre por RA/DEC uma identidade que o dual mode
#ja estabeleceu.
#M. Agenor
#====================
def build_master(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    corrected: dict[str, Path],
    log: Logger,
    force: bool = False,
) -> Path:
    out = paths.master / f"{field}_master_photometry_only_dual.fits"
    if out.is_file() and not force:
        log.write(f"Instrumental PStotal master already exists: {out}")
        return out

    ref_filter = settings.filters[0]
    ref = Table.read(corrected[ref_filter])
    ref_number = np.asarray(ref["NUMBER"])
    n = len(ref)

    master = Table()
    master["calib_ID"] = np.asarray(
        [f"calib_SITH_{field}_dual_{i:07d}" for i in range(1, n + 1)],
        dtype="U50",
    )
    master["RAJ2000"] = np.asarray(ref["ALPHA_J2000"], dtype=np.float64)
    master["DEJ2000"] = np.asarray(ref["DELTA_J2000"], dtype=np.float64)
    master["X"] = np.asarray(ref["X_IMAGE"], dtype=np.float32)
    master["Y"] = np.asarray(ref["Y_IMAGE"], dtype=np.float32)
    master["CLASS_STAR"] = np.asarray(ref["CLASS_STAR"], dtype=np.float32)

    ndet = np.zeros(n, dtype=np.int16)

    for filt in settings.filters:
        tab = Table.read(corrected[filt])
        number = np.asarray(tab["NUMBER"])

        if len(number) != n or not np.array_equal(number, ref_number):
            raise RuntimeError(
                f"Dual NUMBER mismatch while building master: {ref_filter} versus {filt}. "
                "Spatial re-matching is intentionally disabled."
            )

        master[f"{filt}_SExNumber"] = number.astype(np.int64, copy=False)

        mag = np.asarray(tab["MAG_PStotal"], dtype=float).copy()
        err = np.asarray(tab["MAGERR_PStotal"], dtype=float).copy()

        good = np.isfinite(mag) & (mag != 99)
        ndet[good] += 1
        mag[~np.isfinite(mag)] = 99.0
        err[~np.isfinite(err)] = 99.0

        master[f"SPLUS_{filt}"] = mag.astype(np.float32)
        master[f"SPLUS_{filt}_err"] = err.astype(np.float32)

    master["NDET"] = ndet
    master.meta["FIELD"] = field
    master.meta["JOINMODE"] = "DUAL_NUMBER"
    master.write(out, overwrite=True)

    log.write(
        f"Instrumental PStotal master built directly by dual NUMBER: {out} "
        f"({len(master)} sources)"
    )
    return out
