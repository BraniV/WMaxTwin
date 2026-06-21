# WMaxTwin Tecator Functional Regression, MATLAB version

This folder contains a standalone MATLAB notebook-style script for Example III in the WMaxTwin paper:

**Multiscale MaxTwin: Data Splitting as Multiscale Feature Geometry**

The example studies scalar-on-function regression for the StatLib Tecator near-infrared spectroscopy data. The response is fat percentage. The script compares internal validation split geometries inside the 172-sample calibration set and evaluates the selected wavelet-ridge model on the 43-sample external test set.

## Main file

- `WMaxTwin_Tecator_Functional_Regression_MATLAB_Notebook.m`

The script is organized using MATLAB `%%` sections, so it can be run cell by cell like a notebook. In MATLAB, you may also save it as a Live Script (`.mlx`) from the editor.

## Dependencies

The script is designed to avoid specialized toolboxes. It implements directly:

- Tecator parsing
- interpolation to 128 wavelengths
- calibration-only standardization
- orthonormal Haar transform matrix
- response-informed WMaxTwin scale weights
- Random, SPlit, Twinning, DUPLEX, MaxTwin, and WMaxTwin split construction
- ridge regression with internal validation over `jmax` and ridge penalty
- custom simple boxplots, so MATLAB `boxplot` is not required

## Outputs

Running the script creates:

- `tecator_wmaxtwin_outputs_MATLAB/tecator_split_results_MATLAB.csv`
- `tecator_wmaxtwin_outputs_MATLAB/tecator_split_summary_MATLAB.csv`
- `tecator_wmaxtwin_outputs_MATLAB/tecator_scale_weights_MATLAB.csv`
- `tecator_wmaxtwin_outputs_MATLAB/tecator_wmaxtwin_table_MATLAB.tex`
- eight PNG figures matching the Python/R notebook structure

## Repetitions

The default value is

```matlab
MATLAB_REPS = 8;
```

Increase this to 50 for a paper-scale Monte Carlo run. Because MATLAB and Python have different random-number streams and because tie-breaking in the split constructors is random, numerical values will not be bitwise identical to the Python notebook. The expected qualitative pattern should match: raw MaxTwin can select a deceptively low validation error but worse external test error, while moderate WMaxTwin perturbations stabilize the selected wavelet-ridge model.

## Reference outputs

The folder `reference_outputs_from_python/` contains the checked Python outputs used in the Jupyter notebook. They are included only for comparison.

## Data

The file `tecator_statlib.txt` is included. If it is removed, the MATLAB script will attempt to download it from StatLib.
