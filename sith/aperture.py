from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import linregress

from .config import Settings
from .log import Logger
from .paths import FieldPaths


def get_apertures_from_sexconf(sexconf: Path) -> np.ndarray:
    with sexconf.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if parts[0] == "PHOT_APERTURES":
                value = "".join(parts[1:])
                return np.array(value.split(","), dtype=float)
    raise ValueError(f"Could not find PHOT_APERTURES in {sexconf}")


# The next scientific functions are direct ports of the corresponding
# SPLUSCALIB algorithms. Their cuts and estimators are intentionally unchanged.
def star_selector(
    catalog: Table,
    s2ncut: tuple[float, float] = (30, 1000),
    starcut: float = 0.9,
    mag_partition: float = 2,
) -> tuple[Table, float]:
    s2n = catalog["FLUX_AUTO"] / catalog["FLUXERR_AUTO"]
    conditions = (
        (catalog["CLASS_STAR"] > starcut)
        & (s2n > s2ncut[0])
        & (s2n < s2ncut[1])
        & (catalog["FLAGS"] == 0)
    )
    select = catalog[conditions]
    inferior, median_fwhm, superior = np.percentile(select["FWHM_WORLD"], [16, 50, 84])
    conditions2 = (select["FWHM_WORLD"] > inferior) & (select["FWHM_WORLD"] < superior)
    select = select[conditions2]
    select.sort("MAG_AUTO")
    select = select[0:int(len(select) / mag_partition)]
    return select, float(median_fwhm)


def obtain_aperture_correction(
    catalog: Path,
    filt: str,
    sexconf: Path,
    aperture: float,
    s2ncut: tuple[float, float],
    starcut: float,
    max_aperture: float = 72.72727,
    mag_partition: float = 2,
    convergence_slope: float = 1e-2,
    check_convergence: bool = True,
) -> tuple[float, float, float, str | None, int]:
    tab = Table.read(catalog)
    apertures = get_apertures_from_sexconf(sexconf)
    round_list = np.around(apertures, 3)
    round_aper = np.around(aperture, 3)
    matches = np.where(round_list == round_aper)[0]
    if len(matches) == 0:
        raise ValueError(f"Aperture {aperture} not present in {sexconf}")
    aperture_id = int(matches[0])

    select, _ = star_selector(tab, s2ncut=s2ncut, starcut=starcut, mag_partition=mag_partition)
    magnitude_corr_list = []
    for star in select:
        magnitude_corr_list.append(star["MAG_APER"] - star["MAG_APER"][aperture_id])

    mincorr, mediancorr, maxcorr = np.percentile(magnitude_corr_list, [16, 50, 84], axis=0)

    correction = np.nan
    correction_low = np.nan
    correction_up = np.nan
    slope = np.nan
    for i in range(len(apertures)):
        if apertures[i] <= max_aperture:
            correction = mediancorr[i]
            correction_low = mincorr[i]
            correction_up = maxcorr[i]
            j = len(apertures) - i
            if check_convergence and i > 2:
                slope_radii = apertures[-(j + 3):-j]
                slope_corrs = mediancorr[-(j + 3):-j]
                slope = linregress(slope_radii, slope_corrs)[0]

    convergence: str | None = None
    if check_convergence:
        if abs(slope) > convergence_slope:
            convergence = "Not_converged"
        else:
            convergence = "Converged"

    return float(correction), float(correction_low), float(correction_up), convergence, len(select)


def growth_curve(catalog: Path, s2ncut=(30, 1000), starcut=0.9, mag_partition=2):
    tab = Table.read(catalog)
    select, median_fwhm = star_selector(tab, s2ncut=s2ncut, starcut=starcut, mag_partition=mag_partition)
    mlist = []
    for star in select:
        mlist.append(np.diff(star["MAG_APER"]))
    percentile = np.percentile(mlist, [2.5, 16, 50, 84, 97.5], axis=0)
    return percentile, median_fwhm, len(select)


def plot_growth_curve(
    catalog: Path,
    filt: str,
    sexconf: Path,
    save_file: Path,
    aperture: float,
    max_aperture: float,
    s2ncut: tuple[float, float],
    starcut: float,
) -> None:
    import matplotlib.pyplot as plt

    apertures = get_apertures_from_sexconf(sexconf)
    result, median_fwhm, starcount = growth_curve(catalog, s2ncut=s2ncut, starcut=starcut, mag_partition=2)
    minmag2, minmag, medianmag, maxmag, maxmag2 = result
    diameter = [(a + b) / 2 for a, b in zip(apertures[:], apertures[1:])]

    plt.figure(figsize=(5, 4))
    maxy, miny = 0.1, -0.5
    plt.plot(diameter, medianmag, color="red")
    plt.fill_between(diameter, minmag, maxmag, color="orange", alpha=0.7)
    plt.fill_between(diameter, minmag2, maxmag2, color="orange", alpha=0.3)
    plt.plot([median_fwhm, median_fwhm], [miny, maxy], color="darkslategray", label="Median FWHM")
    plt.plot([max_aperture, max_aperture], [miny, maxy], "-", color="purple", label=f"{max_aperture} pix")
    plt.plot([aperture, aperture], [miny, maxy], "-", color="blue", label=f"{aperture:.3f}-diameter")
    plt.plot([0, max(diameter)], [0, 0], color="blue")
    plt.fill_between([0, max(diameter)], [-1e-2, -1e-2], [1e-2, 1e-2], color="blue", alpha=0.3)
    plt.ylim(miny, maxy)
    plt.xlim(0, max(diameter))
    plt.legend()
    plt.xlabel("Aperture Diameter (pix)")
    plt.ylabel(r"$m_{k+1} - m_{k}$")
    plt.title(f"Magnitude growth curve in {filt}, {starcount} stars")
    plt.savefig(save_file, bbox_inches="tight")
    plt.close()


def write_apercorr_file(settings: Settings, field: str, paths: FieldPaths, catalogs: dict[str, Path], log: Logger) -> Path:
    out = paths.aperture / f"{field}_apercorr_dual.zp"
    lines: list[str] = []
    for filt in settings.filters:
        sexconf = paths.sexconf / f"{field}_{filt}_dual.sex"
        corr, low, up, conv, nstars = obtain_aperture_correction(
            catalog=catalogs[filt],
            filt=filt,
            sexconf=sexconf,
            aperture=settings.apercorr_diameter,
            s2ncut=settings.apercorr_s2ncut,
            starcut=settings.apercorr_starcut,
            max_aperture=settings.apercorr_max_aperture,
        )
        suffix = f" {conv}" if conv is not None else ""
        lines.append(f"SPLUS_{filt} {corr} {low} {up}{suffix}\n")
        log.write(f"Aperture correction SPLUS_{filt}: {corr:+.6f} (N={nstars}, {conv})")
        if settings.make_plots:
            plot_growth_curve(
                catalogs[filt], filt, sexconf,
                paths.diagnostics / f"growth_curve_{field}_{filt}_dual.png",
                settings.apercorr_diameter,
                settings.apercorr_max_aperture,
                settings.apercorr_s2ncut,
                settings.apercorr_starcut,
            )
    out.write_text("".join(lines), encoding="utf-8")
    return out


def read_apercorr(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        p = line.split()
        out[p[0]] = float(p[1])
    return out


def apply_aperture_correction(
    catalog: Path,
    filt: str,
    sexconf: Path,
    aperture: float,
    apercorr_file: Path,
    save_file: Path,
) -> None:
    cat = fits.open(catalog)
    cat_data = cat[1].data
    corr_dict = read_apercorr(apercorr_file)
    new_data = [col for col in cat_data.columns]

    apertures = get_apertures_from_sexconf(sexconf)
    round_list = np.around(apertures, 3)
    round_aper = np.around(aperture, 3)
    aper_id = int(np.where(round_list == round_aper)[0][0])

    aper_corr = corr_dict[f"SPLUS_{filt}"]
    flux_corr = 10 ** (-aper_corr / 2.5)

    col = copy.deepcopy(cat_data.columns["MAG_APER"].array[:, aper_id])
    f = col != 99
    col[f] = col[f] + aper_corr
    new_data.append(fits.Column(name="MAG_PStotal", format="1E", array=col))

    col = copy.deepcopy(cat_data.columns["MAGERR_APER"].array[:, aper_id])
    new_data.append(fits.Column(name="MAGERR_PStotal", format="1E", array=col))

    col = copy.deepcopy(cat_data.columns["FLUX_APER"].array[:, aper_id] * flux_corr)
    new_data.append(fits.Column(name="FLUX_PStotal", format="1E", array=col))

    col = copy.deepcopy(cat_data.columns["FLUXERR_APER"].array[:, aper_id])
    new_data.append(fits.Column(name="FLUXERR_PStotal", format="1E", array=col))

    fits.BinTableHDU.from_columns(new_data).writeto(save_file, overwrite=True)
    cat.close()


def run_aperture_correction(settings: Settings, field: str, paths: FieldPaths, catalogs: dict[str, Path], log: Logger) -> tuple[Path, dict[str, Path]]:
    zp = write_apercorr_file(settings, field, paths, catalogs, log)
    corrected: dict[str, Path] = {}
    for filt in settings.filters:
        out = paths.aperture / f"sex_{field}_{filt}_dual_apercorr.fits"
        apply_aperture_correction(
            catalogs[filt], filt, paths.sexconf / f"{field}_{filt}_dual.sex",
            settings.apercorr_diameter, zp, out,
        )
        corrected[filt] = out
    return zp, corrected
