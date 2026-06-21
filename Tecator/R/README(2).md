# WMaxTwin Tecator Functional Regression, R version

This GitHub-ready package provides an annotated R implementation of the Tecator functional-regression example from the WMaxTwin paper.

## Files

- `WMaxTwin_Tecator_Functional_Regression_R_Notebook.Rmd` -- RStudio notebook.
- `WMaxTwin_Tecator_Functional_Regression_R_Notebook.ipynb` -- Jupyter notebook using an R kernel.
- `wmaxtwin_tecator_functional_regression_R.R` -- standalone annotated R script.
- `tecator_statlib.txt` -- local copy of the Tecator data file used by the notebook.
- `reference_outputs_from_python/` -- figures and tables from the executed Python notebook, for comparison.

## Requirements

The standalone R script uses only base R, `stats`, `graphics`, and `utils`.  These ship with standard R.  To knit the Rmd, install `knitr` and `rmarkdown`.  To run the `.ipynb`, install an R kernel such as `IRkernel`.

```r
install.packages(c("knitr", "rmarkdown", "IRkernel"))
IRkernel::installspec()
```

## Run from command line

```bash
Rscript wmaxtwin_tecator_functional_regression_R.R 8
```

The optional number is the number of repeated internal splits.  Use 8 for a quick reproduction and 50 or more for a fuller Monte Carlo run.

## Output

Running the script creates `tecator_wmaxtwin_outputs_R/` with:

- `tecator_split_results_R.csv`
- `tecator_split_summary_R.csv`
- `tecator_scale_weights_R.csv`
- `tecator_wmaxtwin_table_R.tex`
- `fig1_tecator_spectra.png`
- `fig2_quartile_mean_spectra.png`
- `fig3_response_scale_weights.png`
- `fig4_validation_curves_rep0.png`
- `fig5_test_rmse_boxplot.png`
- `fig6_selected_jmax_boxplot.png`
- `fig7_selected_alpha_boxplot.png`
- `fig8_gamma_path_rmse.png`

## Note on numerical reproducibility

R and NumPy use different random-number generators and different random tie-breaking in the split construction.  The R results should therefore be compared to the Python outputs by checking the qualitative pattern, the calibration-only scale weights, and the nested WMaxTwin path rather than by expecting bit-for-bit identical split orientations.
