# SITH 

**SITH — S-PLUS Integrated Tool for Homogeneous Photometry**

## Scientific chain

1. Materialize the reduced science/weight images under the internal `_swp` names.
2. Build the GRIZ detection image with the original `swarpsplus.swarp` template.
3. Run Source Extractor in dual mode using the original `SExsplus.sex`,
   `SExsplus.param`, 32 apertures, top-hat convolution, thresholds, background,
   Kron/Petrosian settings and instrumental ZP=20.
4. Calculate aperture correction with direct ports of SPLUSCALIB
   `star_selector()` and `obtain_aperture_correction()`.
5. Apply the correction exactly as SPLUSCALIB constructs `MAG_PStotal`.
6. Build the dual PStotal master with the original STILTS `tmatchn` sky matcher
   at **1 arcsec** and `multimode=group`.
7. Read precomputed `XpContinuousMeanSpectrum_*_synthetic.fits` files from
   `gaiaxp_path`; no Gaia download and no GaiaXPy are executed.
8. Crossmatch S-PLUS and the precomputed synthetic catalogue using the original
   STILTS `tmatch2`, `matcher=sky`, **params=3**, `join=1and2` rule.
9. Estimate each GaiaXP offset with the original Gaussian KDE bandwidth 0.05 and
   the original `np.arange(-10, 10, 0.0001)` grid.
10. Write `FIELD_gaiaxpsp.zp`, `FIELD_final.zp`, diagnostic figures and the
    calibrated PStotal master.

## Install

```bash
python3 -m pip install -e . --no-build-isolation
```

System programs required:

- Source Extractor (`source-extractor` is supported)
- SWarp
- Java (STILTS is bundled with SITH)
- `funpack` is optional; Astropy can decompress `.fits.fz` when unavailable

The exact Source Extractor files `default.nnw` and `tophat_3.0_3x3.conv` are
required.  On Ubuntu/Debian they are commonly in `/usr/share/source-extractor`.

## Check before running

```bash
sith SGR107 ./JEDIS.conf --check
```

## Run

```bash
sith SGR107 ./JEDIS.conf
```

To overwrite existing intermediate products:

```bash
sith SGR107 ./JEDIS.conf --force
```

## Important outputs

For an output root `/.../PST_CONF_P01_sith2` and field `SGR107`:

```text
SGR107/
├── detection/
├── photometry/
│   ├── catalogs/
│   ├── sexconf/
│   └── master/SGR107_master_photometry_only_dual.fits
├── aperture_correction/
│   ├── SGR107_apercorr_dual.zp
│   └── sex_SGR107_*_dual_apercorr.fits
├── gaiaxp/
│   ├── gaiaXPSP_phot_precomputed.fits
│   └── SGR107_SPLUS_dual_gaiaXPSP_phot.fits
├── diagnostics/
│   ├── growth_curve_*.png
│   └── zeropoint_*.png
└── calibration/
    ├── SGR107_gaiaxpsp.zp
    ├── SGR107_final.zp
    ├── SGR107_SITH2_calibrated_master.fits
    └── dual/
```

In zero-point diagnostics, pale blue sources are matches outside the fitting
magnitude interval and orange sources are the stars actually used in the KDE
zero-point fit.
