# WMaxTwin claw-density MATLAB notebook

This directory contains a self-contained MATLAB implementation of Example I from the WMaxTwin paper:

**WMaxTwin for the Marron-Wand claw density**.

The main file is:

- `WMaxTwin_Claw_Density_MATLAB_Notebook.m`

It is both a standalone MATLAB script and a notebook-style script. The `%%` headings define MATLAB code sections. Open it in MATLAB and run the sections sequentially, or run the entire file.

## Native MATLAB notebooks

MATLAB's native notebook format is the **Live Script**, with extension `.mlx`. A `.mlx` file is convenient for interactive display, but it is a binary/zipped format and is less GitHub-friendly than a plain `.m` script. For repository use, the sectioned `.m` file is preferable because it is readable and version-control friendly.

To convert this script to a Live Script:

1. Open `WMaxTwin_Claw_Density_MATLAB_Notebook.m` in MATLAB.
2. Choose **Save As > Live Script**.
3. MATLAB will save a native `.mlx` notebook.

## Requirements

The implementation uses only base MATLAB functions. It does not require the Wavelet Toolbox or Statistics and Machine Learning Toolbox. Haar atoms and Gaussian densities are implemented directly in the script.

## What the script does

The script compares three balanced half-sample split geometries:

1. Random splitting.
2. MaxTwin+ using strengthened nonwavelet local-rank features.
3. WMaxTwin+ using the same local-rank features plus Haar scale-location atoms evaluated at the sample locations.

The same sample, bandwidth, half-sample KDE estimator, and AMSE loss are used for all three split geometries. Thus AMSE differences are attributable to split geometry.

## Outputs

Running the script creates a directory named `wmaxtwin_claw_matlab_output` containing:

- `wmaxtwin_plus_amse_table_matlab.csv`
- `wmaxtwin_plus_replicates_matlab.csv`
- `fig_wmaxtwin_plus_split_geometry_matlab.png`
- `fig_wmaxtwin_plus_density_estimates_matlab.png`
- `fig_wmaxtwin_plus_amse_bars_matlab.png`
- `fig_wmaxtwin_plus_amse_differences_matlab.png`

Because MATLAB and Python use different random-number generators, the exact Monte Carlo numbers will not be bitwise identical to the Python version. The expected pattern is the same: Random has the largest AMSE, MaxTwin+ improves, and WMaxTwin+ improves further.
