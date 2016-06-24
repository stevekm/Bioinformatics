Instructions for use:

- Update the file `peaks-per-sample_report/report/peaks-per-sample_report.Rnw` with the following items:
- - title and project ID information in the `report_setup` code chunk (update all relevant file paths)
- - pipeline results and script information in the `run_pipeline` code chunk
- Compile RNW report
- - use the `knitr` package in R to compile the RNW; `knit("peaks-per-sample_report/report/peaks-per-sample_report.Rnw")`. This will produce TEX output file.
- - compile the TEX with `pdflatex` to produce a PDF report; `pdflatex peaks-per-sample_report/report/peaks-per-sample_report.tex`. Run this 2 - 3 times to complete compilation.
