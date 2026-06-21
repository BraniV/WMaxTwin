# WMaxTwin Doppler / Nason-type wavelet cross-validation in MATLAB

This package contains a self-contained MATLAB implementation of Example II from the WMaxTwin paper: Nason-type wavelet cross-validation for the Doppler signal.

The main file is:

```text
WMaxTwin_Doppler_Nason_CV_MATLAB_Notebook.m
```

It is written as a MATLAB notebook-style script using `%%` sections. Open it in MATLAB and run the sections interactively, or run the full script from the command line.

## Main statistical point

Nason's even/odd split is grid-compatible. MaxTwin and WMaxTwin masks are generally not dyadic subgrids. Therefore the wavelet fit must be written as a masked wavelet penalized least-squares problem on the original grid, rather than as an ordinary DWT of the retained observations.

The script compares:

1. Noisy data.
2. Levelwise SURE.
3. Nason even/odd CV with a fixed global threshold.
4. Location-only MaxTwin with a fixed global threshold.
5. WMaxTwin with high-frequency multiscale geometry and a ramp threshold family.
6. Oracle ramp, included only as an unattainable reference.

## MATLAB dependencies

Base MATLAB only. The script implements periodic `sym4` and `db2` DWT/IDWT directly from filter coefficients. No Wavelet Toolbox is required.

## Default run

The default settings are:

```matlab
N          = 1024;
SNR        = 7;
SIGMA      = 1;
WAVELET    = 'sym4';
J0         = 3;
NREP       = 18;
SEED       = 2026;
N_ITER     = 30;
WEIGHT_MODE = 'paper';
```

The explicit WMaxTwin feature weights for levels `j = 3,...,9` are

```text
0.015, 0.036, 0.074, 0.124, 0.182, 0.248, 0.321
```

These weights emphasize high-frequency levels while retaining a small contribution from lower levels.

## Outputs

When run, the script writes outputs to `outputs_matlab/`:

```text
doppler_replicates_MATLAB.csv
doppler_amse_summary_MATLAB.csv
fig01_doppler_signal.png
fig02_feature_weights.png
fig03_one_replication.png
fig04_split_mask.png
fig05_thresholds.png
fig06_cv_curve.png
fig07_amse_bar.png
fig08_slope_diagnostics.png
```

## Reference outputs

The directory `reference_outputs_from_python/` contains reference outputs from the executed Python notebook. MATLAB and Python random number generators differ, so exact replicate-level values will not be identical. The expected qualitative ordering is the same: WMaxTwin ramp should be close to the oracle ramp and should improve on the fixed-threshold MaxTwin construction.

## Live Script note

MATLAB's native notebook format is the Live Script format `.mlx`. This package provides the plain-text `.m` version because it is better for GitHub review, diffs, and reproducibility. In MATLAB, open the `.m` file and use **Save As > Live Script** if a native `.mlx` version is desired.
