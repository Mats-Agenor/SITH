from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table
from sklearn.neighbors import KernelDensity

from .config import Settings
from .log import Logger
from .paths import FieldPaths


@dataclass(frozen=True)
class ZeroPoint:
    filt: str
    offset: float
    final_zp: float
    mag_low: float
    mag_up: float
    nstars: int


def mean_robust(x, low=3, high=3):
    """Direct port of spluscalib mean_robust."""
    x = np.array(x)
    mean_x = np.nanmean(x)
    std_x = np.nanstd(x)
    x = x[(x > (mean_x - low * std_x)) & (x < (mean_x + high * std_x))]
    return np.mean(x)


def get_filter_zeropoint(obs_mag_array: np.ndarray, model_mag_array: np.ndarray, cut=(14, 19)) -> float:
    """Direct mathematical port of spluscalib get_filter_zeropoint()."""
    f = (model_mag_array >= cut[0]) & (model_mag_array <= cut[1])
    f = f & (~np.isnan(obs_mag_array))
    delta_array = model_mag_array[f] - obs_mag_array[f]
    delta_array = np.asarray(delta_array).reshape(-1, 1)
    if len(delta_array) < 1:
        raise RuntimeError(f"No stars available for zero point in magnitude range {cut}")
    kde_dens = KernelDensity(kernel="gaussian", bandwidth=0.05).fit(delta_array)
    x = np.arange(-10, 10, 0.0001)
    y = np.exp(kde_dens.score_samples(x.reshape(-1, 1)))
    mode = x[y == np.max(y)][0]
    return float(mode)


def _write_zp(path: Path, values: dict[str, float], filters: list[str]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for filt in filters:
            key = f"SPLUS_{filt}"
            fh.write(f"{key:s} {values[key]:.5f}\n")


def _plot_diagnostic(
    settings: Settings,
    paths: FieldPaths,
    filt: str,
    obs: np.ndarray,
    model: np.ndarray,
    cut: tuple[float, float],
    zp: float,
) -> None:
    if not settings.make_plots:
        return
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    # Same plot preparation as spluscalib: discard 99-like instrumental values
    # for visualization, then select the model magnitude interval.
    valid = np.asarray(obs) < 50
    x = np.asarray(model, dtype=float)[valid]
    y = x - np.asarray(obs, dtype=float)[valid]
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    selection = (x >= cut[0]) & (x <= cut[1])
    delta = y[selection].reshape(-1, 1)
    if len(delta) == 0:
        return

    kde = KernelDensity(kernel="gaussian", bandwidth=0.05).fit(delta)
    y_dens = np.arange(-10, 10, 0.001)
    x_dens = np.exp(kde.score_samples(y_dens.reshape(-1, 1)))
    mode = float(y_dens[x_dens == np.max(x_dens)][0])
    mu = float(np.mean(y[selection]))
    mu_robust = float(mean_robust(y[selection]))

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.04)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharey=ax1)

    #====================
    #Estrutura visual adicionada para diagnostico:
    #azul claro = matches fora do corte; laranja = estrelas usadas no ZP.
    #Isto altera apenas a figura, nao o estimador do spluscalib.
    #M. Agenor
    #====================
    ax1.scatter(x[~selection], y[~selection], c="#6baed6", s=28, alpha=0.16, edgecolors="none", zorder=1)
    ax1.scatter(x[selection], y[selection], c="#ff7f0e", s=28, alpha=0.60, edgecolors="none", zorder=2)
    ax1.plot([cut[0], cut[1]], [mode, mode], color="#FF3219", linestyle="-", zorder=6, label=f"mode: {mode:.4f}")
    ax1.axhline(zp, color="#000000", linestyle="--", zorder=5, label=f"zp: {zp:.4f}")
    ax1.axhline(mu, color="#1932FF", linestyle="--", zorder=4, label=f"mean: {mu:.4f}")
    ax1.axhline(mu_robust, color="#19DD32", linestyle=":", zorder=4, label=f"mean_robust: {mu_robust:.4f}")

    xlim = [cut[0] - 2, cut[1] + 4]
    ylim = [mode - 1.5, mode + 1.5]
    ax1.text(cut[0] - 1.5, mode + 1.3, f"SPLUS_{filt}", fontsize=20)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xlabel("Model")
    ax1.set_ylabel("Model - Instrumental")
    ax1.legend(fontsize=9, loc="upper right")

    ax2.plot(x_dens / np.max(x_dens), y_dens, color="#1f77b4", zorder=3)
    ax2.plot([0, 2], [zp, zp], color="#000000", linestyle="--", zorder=3)
    ax2.plot([0, 2], [mode, mode], color="#FF3219", linestyle="-", zorder=2)
    ax2.plot([0, 2], [mu, mu], color="#1932FF", linestyle="--", zorder=1)
    ax2.plot([0, 2], [mu_robust, mu_robust], color="#19DD32", linestyle=":", zorder=1)
    ax2.set_xlim([0, 1.1])
    ax2.set_ylim(ylim)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)

    for i in np.arange(-10, 10, 0.5):
        ax1.plot([0, 30], [i, i], color="#666666", alpha=0.3, linewidth=0.5, zorder=-5)
        ax2.plot([0, 30], [i, i], color="#666666", alpha=0.3, linewidth=0.5, zorder=-5)
    for i in np.arange(-10, 10, 0.1):
        ax1.plot([0, 30], [i, i], color="#aaaaaa", alpha=0.12, linewidth=0.3, zorder=-6)
        ax2.plot([0, 30], [i, i], color="#aaaaaa", alpha=0.12, linewidth=0.3, zorder=-6)
    for i in np.arange(0, 30, 1):
        ax1.plot([i, i], [-10, 10], color="#666666", alpha=0.3, linewidth=0.5, zorder=-5)

    fig.savefig(paths.diagnostics / f"zeropoint_{filt}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def estimate_zero_points(settings: Settings, field: str, paths: FieldPaths, crossmatch_file: Path, log: Logger) -> dict[str, ZeroPoint]:
    tab = Table.read(crossmatch_file)
    result: dict[str, ZeroPoint] = {}
    offsets: dict[str, float] = {}
    finals: dict[str, float] = {}

    for filt in settings.filters:
        obs = np.asarray(tab[f"SPLUS_{filt}"], dtype=float)
        model = np.asarray(tab[f"SPLUS_{filt}_mod"], dtype=float)
        low, up = settings.mag_limits(filt)
        zp = get_filter_zeropoint(obs, model, cut=(low, up))
        selection = (model >= low) & (model <= up) & (~np.isnan(obs))
        key = f"SPLUS_{filt}"
        offsets[key] = zp
        finals[key] = settings.inst_zp + zp
        result[filt] = ZeroPoint(filt, zp, settings.inst_zp + zp, low, up, int(np.sum(selection)))
        log.write(f"{key}: GaiaXP offset={zp:+.5f}; final ZP={settings.inst_zp + zp:.5f}; Ncut={np.sum(selection)}")
        _plot_diagnostic(settings, paths, filt, obs, model, (low, up), zp)

    # Same role as SPLUSCALIB gaiaXPSP temporary/final offset and final.zp.
    _write_zp(paths.calibration / f"{field}_gaiaxpsp.zp", offsets, settings.filters)
    _write_zp(paths.calibration / f"{field}_final.zp", finals, settings.filters)
    return result
