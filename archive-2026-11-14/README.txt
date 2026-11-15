Eldor et al. (2023) rebuttal: code and outputs from 14 Nov 2026.

code/      Python and R scripts. The R scripts need the KRLS package (source: github.com/cran/KRLS, version 1.7-1).
           To re-run in a fresh workspace: apt-get update; apt-get install -y --no-install-recommends r-base-core;
           git clone https://github.com/cran/KRLS; R CMD INSTALL KRLS.
           The digitise scripts need the paper's PDF: figures are extracted with `pdfimages -png` (pages 5, 7, 8, 9, 10, 11, 13, 14).
logs/      Raw output of every simulation run reported in the project notes.
digitised/ The digitised lines and dots.
What each script produced is listed in the project documents eldor-krls-note (section 7) and eldor-figure-digitisation (section 9).
