from __future__ import annotations

from pathlib import Path

from .config import Settings
from .fitsutil import assert_same_dimensions, image_header
from .log import Logger
from .paths import FieldPaths


#====================
#Estrutura restaurada a partir do workflow do SITH 1, que produziu a
#fotometria correta para os mosaicos PST_CONF_P01. A imagem de deteccao GRIZ
#usa CENTER_TYPE MANUAL e o centro da imagem de referencia.
#
#Correcao adicional: o SITH 2.0 nao possui os atributos swarp_image_size,
#swarp_nthreads e pixel_scale na classe Settings. Estes valores agora sao
#obtidos diretamente do header FITS da imagem de referencia quando aplicavel:
#NAXIS1/NAXIS2 -> IMAGE_SIZE e PIXSCALE -> PIXEL_SCALE. NTHREADS permanece 2,
#que era o valor utilizado no workflow do SITH 1.
#M. Agenor
#====================
def write_swarp_config(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    images: dict[str, tuple[Path, Path | None]],
) -> Path:
    ref = images[settings.detection_filters[0]][0]
    h = image_header(ref)

    center = f"{float(h['CRVAL1'])},{float(h['CRVAL2'])}"

    # O SITH 2.0 nao armazena mais essas opcoes em Settings. Para os mosaicos
    # reduzidos, a propria imagem de referencia e a fonte correta desses dados.
    naxis1 = int(h['NAXIS1'])
    naxis2 = int(h['NAXIS2'])
    image_size = f"{naxis1},{naxis2}"

    if 'PIXSCALE' in h:
        pixel_scale = float(h['PIXSCALE'])
    else:
        # Compatibilidade com os mosaicos S-PLUS originais, cuja escala usada
        # pelo SPLUSCALIB/SITH 1 e 0.55 arcsec/pixel. Nao altera os pixels.
        pixel_scale = 0.55

    swarp_nthreads = 2

    image_out = paths.detection / f"{field}_detection.fits"
    weight_out = paths.detection / f"{field}_detectionweight.fits"
    xml_out = paths.detection / "swarp_output.xml"
    config_file = paths.detection / f"{field}.swarp"

    combine_type = "WEIGHTED" if settings.use_weight else "AVERAGE"
    weight_type = "MAP_WEIGHT" if settings.use_weight else "NONE"

    text = f"""IMAGEOUT_NAME        {image_out}
WEIGHTOUT_NAME       {weight_out}
RESAMPLE_DIR         {paths.detection}
XML_NAME             {xml_out}
COMBINE_TYPE         {combine_type}
WEIGHT_TYPE          {weight_type}
BACK_DEFAULT         0.0
BACK_FILTERSIZE      5
BACK_FILTTHRESH      0.0
BACK_SIZE            576
BACK_TYPE            AUTO
BLANK_BADPIXELS      N
CELESTIAL_TYPE       EQUATORIAL
CENTER_TYPE          MANUAL
CENTER               {center}
COMBINE               Y
COMBINE_BUFSIZE       1024
COPY_KEYWORDS         OBJECT,TELESCOP,INSTRUME
DELETE_TMPFILES       Y
FSCALASTRO_TYPE       VARIABLE
FSCALE_DEFAULT        1.0
FSCALE_KEYWORD        FLXSCALE
GAIN_DEFAULT          0.0
GAIN_KEYWORD          GAIN
HEADER_ONLY           N
HEADER_SUFFIX         .head
IMAGE_SIZE            {image_size}
INTERPOLATE           N
MEM_MAX               2048
NNODES                1
NODE_INDEX            0
NTHREADS              {swarp_nthreads}
OVERSAMPLING          0
PIXELSCALE_TYPE       MANUAL
PIXEL_SCALE           {pixel_scale}
PROJECTION_ERR        0.001
PROJECTION_TYPE       TAN
RESAMPLE              Y
RESAMPLE_SUFFIX       resamp.fits
RESAMPLING_TYPE       LANCZOS3
RESCALE_WEIGHTS       Y
SATLEV_DEFAULT        60000.0
SATLEV_KEYWORD        SATURATE
SUBTRACT_BACK         Y
VERBOSE_TYPE          NORMAL
WEIGHT_SUFFIX         weight.fits
WRITE_FILEINFO        N
WRITE_XML             Y
"""
    config_file.write_text(text, encoding="utf-8")
    return config_file


def build_detection(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    images: dict[str, tuple[Path, Path | None]],
    log: Logger,
    force: bool = False,
) -> tuple[Path, Path | None]:
    detection = paths.detection / f"{field}_detection.fits"
    dweight = paths.detection / f"{field}_detectionweight.fits"

    if detection.exists() and (not settings.use_weight or dweight.exists()) and not force:
        # Mantem apenas a verificacao objetiva de dimensoes usada pelo workflow
        # corrigido. Nao introduz uma nova exigencia numerica de igualdade WCS.
        for filt in settings.filters:
            assert_same_dimensions(detection, images[filt][0])
        log.write(f"Detection image already exists: {detection}")
        return detection, dweight if settings.use_weight else None

    image_list = paths.detection / "combine_images.txt"
    image_list.write_text(
        "\n".join(str(images[f][0]) for f in settings.detection_filters) + "\n",
        encoding="utf-8",
    )
    cfg = write_swarp_config(settings, field, paths, images)

    log.write(
        f"Building detection image from {', '.join(settings.detection_filters)} "
        f"with SITH 1 workflow..."
    )
    log.command([settings.swarp_path, f"@{image_list}", "-c", str(cfg)])

    if not detection.is_file():
        raise RuntimeError(f"SWarp finished without producing {detection}")
    if settings.use_weight and not dweight.is_file():
        raise RuntimeError(f"SWarp finished without producing {dweight}")

    for filt in settings.filters:
        assert_same_dimensions(detection, images[filt][0])

    return detection, dweight if settings.use_weight else None
