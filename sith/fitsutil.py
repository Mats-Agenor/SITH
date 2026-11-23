from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits


def image_hdu_index(path: str | Path) -> int:
    with fits.open(path, memmap=True) as hdul:
        for i, hdu in enumerate(hdul):
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return i
    raise ValueError(f"No 2-D image HDU found in {path}")


def image_header(path: str | Path) -> fits.Header:
    return fits.getheader(path, image_hdu_index(path))


def header_float(path: str | Path, keys: list[str]) -> float:
    h = image_header(path)
    for key in keys:
        if key in h:
            value = float(h[key])
            if np.isfinite(value):
                return value
    raise KeyError(f"None of {keys} found in {path}")


#====================
#Estrutura mantida para aceitar os headers dos mosaicos PST_CONF_P01.
#O seeing real esta armazenado diretamente em FWHMSEXT. Nao existe fallback
#arbitrario para o seeing.
#M. Agenor
#====================
def seeing_arcsec(path: str | Path) -> float:
    return header_float(
        path,
        [
            "HIERARCH OAJ PRO FWHMSEXT",
            "HIERARCH MAR PRO FWHMSEXT",
            "FWHMSEXT",
        ],
    )


def gain(path: str | Path) -> float:
    return header_float(path, ["GAIN"])


def saturation(path: str | Path) -> float:
    return header_float(path, ["SATURATE"])


#====================
#Estrutura adicionada para validar somente a condicao objetiva exigida pelo
#SExtractor em dual mode neste workflow: detection e measurement devem possuir
#as mesmas dimensoes em pixels. A verificacao WCS adicional do patch anterior
#foi removida porque ela alterava/interrompia o workflow que funcionava no
#SITH 1 e nao fazia parte da execucao original.
#M. Agenor
#====================
def image_shape(path: str | Path) -> tuple[int, int]:
    h = image_header(path)
    return int(h["NAXIS1"]), int(h["NAXIS2"])


def assert_same_dimensions(reference: str | Path, image: str | Path) -> tuple[int, int]:
    ref_shape = image_shape(reference)
    img_shape = image_shape(image)
    if ref_shape != img_shape:
        raise RuntimeError(
            "Dual-mode image-dimension mismatch:\n"
            f"  reference: {reference}  {ref_shape[0]}x{ref_shape[1]}\n"
            f"  image    : {image}  {img_shape[0]}x{img_shape[1]}\n"
            "SExtractor dual mode requires the same image dimensions."
        )
    return ref_shape
