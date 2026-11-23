from __future__ import annotations

import os
import shutil
from pathlib import Path

from astropy.io import fits

from .config import Settings
from .fitsutil import image_hdu_index, image_header, seeing_arcsec
from .log import Logger
from .paths import FieldPaths


def source_file(settings: Settings, field: str, filt: str, weight: bool = False) -> Path:
    pattern = settings.weight_pattern if weight else settings.image_pattern
    return settings.input_path / field / pattern.format(field=field, filter=filt)


def validate_inputs(settings: Settings, field: str, log: Logger) -> None:
    field_dir = settings.input_path / field
    if not field_dir.is_dir():
        raise FileNotFoundError(f"Input field directory not found: {field_dir}")

    missing_sci: list[str] = []
    missing_wei: list[str] = []
    log.write("Checking configured input images:")
    for filt in settings.filters:
        sci = source_file(settings, field, filt)
        wei = source_file(settings, field, filt, True)
        sci_ok = sci.is_file()
        wei_ok = wei.is_file()
        log.write(f"  {filt:<5} science={'OK' if sci_ok else 'MISSING':<7} weight={'OK' if wei_ok else 'MISSING':<7}")
        if not sci_ok:
            missing_sci.append(filt)
        if settings.use_weight and not wei_ok:
            missing_wei.append(filt)
        if sci_ok:
            h = image_header(sci)
            for key in ("CRVAL1", "CRVAL2", "GAIN", "SATURATE"):
                if key not in h:
                    raise KeyError(f"Required SPLUSCALIB header keyword {key} missing in {sci}")
            # Requires the real FWHMSEXT quantity; no invented seeing fallback.
            seeing_arcsec(sci)
    if missing_sci:
        raise FileNotFoundError(f"Missing science filters: {', '.join(missing_sci)}")
    if missing_wei:
        raise FileNotFoundError(f"Missing weight filters: {', '.join(missing_wei)}")

    if not settings.gaiaxp_path.is_dir():
        raise FileNotFoundError(
            f"GaiaXP directory is not accessible: {settings.gaiaxp_path}\n"
            "If this is GVFS/SMB, mount the share in the current graphical login first."
        )
    gaia_files = list(settings.gaiaxp_path.glob(settings.gaia_filename_glob))
    if not gaia_files:
        raise FileNotFoundError(
            f"No Gaia synthetic FITS matching {settings.gaia_filename_glob!r} in {settings.gaiaxp_path}"
        )
    log.write(f"Gaia synthetic catalogs visible: {len(gaia_files)}")
    log.write(f"Exact SExtractor convolution: {settings.sextractor_config_dir / 'tophat_3.0_3x3.conv'}")
    log.write(f"Bundled STILTS: {settings.stilts_path}")


def _astropy_decompress(src: Path, dst: Path) -> None:
    idx = image_hdu_index(src)
    with fits.open(src, memmap=False) as hdul:
        hdu = hdul[idx]
        fits.PrimaryHDU(data=hdu.data, header=hdu.header.copy()).writeto(dst, overwrite=True)


def _materialize(src: Path, dst: Path, settings: Settings, log: Logger, force: bool) -> None:
    if dst.exists() and not force:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()

    lower = src.name.lower()
    if lower.endswith(".fits.fz") or lower.endswith(".fz"):
        if settings.funpack_path:
            log.command([settings.funpack_path, "-O", str(dst), str(src)])
        else:
            log.write(f"funpack unavailable; Astropy decompressing {src.name}")
            _astropy_decompress(src, dst)
    elif lower.endswith(".fits") or lower.endswith(".fit"):
        try:
            os.symlink(src.resolve(), dst)
        except OSError:
            shutil.copy2(src, dst)
    else:
        raise ValueError(f"Unsupported FITS extension: {src}")


def prepare_images(settings: Settings, field: str, paths: FieldPaths, log: Logger, force: bool = False) -> dict[str, tuple[Path, Path | None]]:
    """Materialize user inputs using the internal filenames expected by SPLUSCALIB."""
    images: dict[str, tuple[Path, Path | None]] = {}
    for filt in settings.filters:
        sci = paths.work_images / f"{field}_{filt}_swp.fits"
        _materialize(source_file(settings, field, filt), sci, settings, log, force)
        wei: Path | None = None
        if settings.use_weight:
            wei = paths.work_images / f"{field}_{filt}_swpweight.fits"
            _materialize(source_file(settings, field, filt, True), wei, settings, log, force)
        images[filt] = (sci, wei)
    return images
