from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.table import Table

from .config import Settings
from .fitsutil import assert_same_dimensions, gain, saturation, seeing_arcsec
from .log import Logger
from .paths import FieldPaths


def _sex_template() -> str:
    return (Path(__file__).resolve().parent / "resources" / "SExsplus.sex").read_text(encoding="utf-8")


def _param_file() -> Path:
    return Path(__file__).resolve().parent / "resources" / "SExsplus.param"


#====================
#Estrutura mantida para usar exatamente os parametros cientificos do
#spluscalib no Source Extractor. O filtro tophat_3.0_3x3.conv e o default.nnw
#continuam obrigatorios; somente os caminhos locais sao substituidos.
#M. Agenor
#====================
def write_dual_config(
    settings: Settings,
    catalog: Path,
    measurement_image: Path,
    detection: Path,
    out: Path,
) -> None:
    template = _sex_template()
    conv = settings.sextractor_config_dir / "tophat_3.0_3x3.conv"
    nnw = settings.sextractor_config_dir / "default.nnw"
    if not conv.is_file() or not nnw.is_file():
        raise FileNotFoundError("Exact SPLUSCALIB SExtractor resources are missing")

    template = template.replace("{path_to_sex}/config/tophat_3.0_3x3.conv", str(conv))
    template = template.replace("{path_to_sex}/config/default.nnw", str(nnw))

    text = template.format(
        catalog_file=catalog,
        param_file=_param_file(),
        path_to_sex="",
        inst_zp=settings.inst_zp,
        satur=saturation(measurement_image),
        seeing=seeing_arcsec(measurement_image),
        gain=gain(measurement_image),
    )

    if settings.use_weight:
        measurement_weight = Path(str(measurement_image).replace(".fits", "weight.fits"))
        detection_weight = Path(str(detection).replace(".fits", "weight.fits"))
        if not measurement_weight.is_file() or not detection_weight.is_file():
            raise FileNotFoundError(
                f"Dual weight image missing: detection={detection_weight}, measurement={measurement_weight}"
            )
        text += (
            "WEIGHT_TYPE MAP_WEIGHT, MAP_WEIGHT\n"
            f"WEIGHT_IMAGE {detection_weight}, {measurement_weight}\n"
        )

    out.write_text(text, encoding="utf-8")


#====================
#Estrutura restaurada para o workflow do SITH 1: o dual mode e executado
#diretamente para os 12 filtros e NUMBER e tratado como identidade da fonte.
#Aqui mantem-se somente a igualdade #das dimensoes e a verificacao rigorosa
#de NUMBER apos a fotometria.
#M. Agenor
#====================
def run_dual_photometry(
    settings: Settings,
    field: str,
    paths: FieldPaths,
    images: dict[str, tuple[Path, Path | None]],
    detection: Path,
    log: Logger,
    force: bool = False,
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}

    for filt in settings.filters:
        image = images[filt][0]
        cat = paths.catalogs / f"sex_{field}_{filt}_dual.fits"
        cfg = paths.sexconf / f"{field}_{filt}_dual.sex"
        outputs[filt] = cat

        assert_same_dimensions(detection, image)

        if cat.is_file() and not force:
            log.write(f"Dual SExtractor catalog already exists: {cat.name}")
            continue

        write_dual_config(settings, cat, image, detection, cfg)
        log.write(f"SExtractor dual mode: {filt}")
        log.command([settings.sextractor_path, f"{detection},{image}", "-c", str(cfg)])
        if not cat.is_file():
            raise RuntimeError(f"SExtractor did not create {cat}")

    verify_dual_catalogs(outputs)
    return outputs


def verify_dual_catalogs(catalogs: dict[str, Path]) -> None:
    ref_filter = next(iter(catalogs))
    ref_tab = Table.read(catalogs[ref_filter])
    ref = np.asarray(ref_tab["NUMBER"])

    for filt, path in catalogs.items():
        tab = Table.read(path)
        number = np.asarray(tab["NUMBER"])
        if len(ref) != len(number) or not np.array_equal(ref, number):
            raise RuntimeError(
                "Dual-mode identity mismatch: NUMBER is not identical between "
                f"{ref_filter} and {filt}. The master must not be spatially re-matched."
            )
