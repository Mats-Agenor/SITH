from __future__ import annotations

import shutil
from pathlib import Path

from .aperture import run_aperture_correction
from .catalog import write_calibrated_filter_catalogs, write_calibrated_master
from .config import Settings
from .crossmatch import crossmatch_splus_gaiaxp
from .detection import build_detection
from .gaiaxp import build_gaia_synthetic_subset
from .inputs import prepare_images, validate_inputs
from .log import Logger
from .master import build_master
from .paths import FieldPaths
from .photometry import run_dual_photometry
from .zeropoint import estimate_zero_points


def run(field: str, config_file: str | Path, force: bool = False, check_only: bool = False) -> Path | None:
    settings = Settings.from_file(config_file)
    paths = FieldPaths.create(settings.output_path, field)
    log = Logger(paths.logs / "sith2.log")

    log.write("SITH 2.0 — S-PLUS Integrated Tool for Homogeneous Photometry")
    log.write("=" * 76)
    log.write(f"Field      : {field}")
    log.write(f"Input      : {settings.input_path / field}")
    log.write(f"Output     : {paths.root}")
    log.write(f"GaiaXP FITS: {settings.gaiaxp_path}")
    log.write(f"SExtractor : {settings.sextractor_path}")
    log.write(f"SWarp      : {settings.swarp_path}")
    log.write(f"Java/STILTS: {settings.java_path} / {settings.stilts_path}")
    log.write("=" * 76)

    validate_inputs(settings, field, log)
    if check_only:
        log.write("CHECK OK: inputs and dependencies are available.")
        return None

    shutil.copy2(settings.config_file, paths.logs / settings.config_file.name)

    log.write("[1/9] Preparing reduced FITS inputs")
    images = prepare_images(settings, field, paths, log, force=force)

    log.write("[2/9] Building SPLUSCALIB-equivalent GRIZ SWarp detection image")
    detection, _ = build_detection(settings, field, paths, images, log, force=force)

    log.write("[3/9] Running SPLUSCALIB-equivalent SExtractor dual photometry")
    catalogs = run_dual_photometry(settings, field, paths, images, detection, log, force=force)

    log.write("[4/9] Running original SPLUSCALIB aperture-correction algorithm")
    _, corrected = run_aperture_correction(settings, field, paths, catalogs, log)

    log.write("[5/9] Building original SPLUSCALIB dual PStotal master (STILTS 1 arcsec)")
    master = build_master(settings, field, paths, corrected, log, force=force)

    log.write("[6/9] Reading precomputed GaiaXP S-PLUS synthetic photometry")
    gaia = build_gaia_synthetic_subset(settings, paths, detection, log, force=force)

    log.write("[7/9] Original SPLUSCALIB GaiaXP crossmatch rule (STILTS 3 arcsec)")
    cross = crossmatch_splus_gaiaxp(settings, field, paths, master, gaia, log, force=force)

    log.write("[8/9] Original SPLUSCALIB GaiaXP KDE zero-point estimator")
    zps = estimate_zero_points(settings, field, paths, cross, log)

    log.write("[9/9] Applying final PStotal calibration")
    final = write_calibrated_master(settings, field, paths, master, zps)
    write_calibrated_filter_catalogs(settings, field, paths, corrected, zps)

    log.write(f"DONE: {final}")
    log.write(f"Zero points: {paths.calibration / (field + '_final.zp')}")
    log.write(f"Diagnostics: {paths.diagnostics}")

    if not settings.keep_work_images:
        shutil.rmtree(paths.root / "work", ignore_errors=True)
        log.write("Removed decompressed work images (keep_work_images=False)")
    return final
