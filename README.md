# Mixed Distribution Flood Frequency Analysis code for Roop-Eckart et al. (2018)

# Citation

Roop-Eckart, K. J., Lee, B. S., Spence, C., Russo, T., Klaus, K. (2018), "Effects of Mixed Distribution Statistical Flood Frequency Models on Dam Safety Assessments: A Case Study of the Pueblo Dam". Github Repository. https://github.com/Joelroopeckart/mixed-distributions

The thesis, titled "Analysis of Mixed Distribution Statistical Flood Frequency Models and Implications for Dam Safety Assessments" is published in the Pennsylvania State University library system, and will be publicly released two years after July 17th, 2018

Paper coming soon

# Overview

This study quantifies the "Effects of Mixed Distribution Statistical Flood Frequency Models on Dam Safety Assessments" through "A Case Study of the Pueblo Dam". Neglecting mixed distribution of peak flows may lead to large underestimations of flood hazard and dam failure.

This code is documentation of the methods used in Roop-Eckart et al. (Coming soon)

# Requirements

The R scripts were written and tested in R v3.4.0 (2017-04-21) and require the following packages:

DEoptim, fExtremes, triangle, RColorBrewer, and zipfR

The R code automatically installs and loads the necessary packages. 
To install and load them manually see the following example:

install.packages("DEoptim")
library(DEoptim)

# Instructions

To run this script:
1) Download the "mixed_distribution_FFA" folder
2) Open the "Code" folder
3) Open "Master_script_thesis.R" in R or Rstudio and change the working directories to your folder structure
IMPORTANT NOTE: the working directory must be changed in three places, line 20 and 70 of the "Master_script_thesis.R", and line 22 of the "hydrograph_overtopping_script.R" code.
4) If you are using windows, run the "Master_script_thesis.R". The run should take roughly 24 core hours.
5) If you are using Mac/Linux, make the necessary edits to the code before running
    In the "LoadingData_Uncertainty.R" code, change "method" from "wininet" to "curl on line 5
    In the "hydrograph_overtopping_script.R" code, uncomment out lines 37-39 and change flrout-kjr-edited.exe to flrout-kjr-edited.out on line 42
    
IMPORTANT NOTE: The results of this study were calculated on Windows, non-Windows execution should, but may not, perfectly reproduce the results.

# Authorship

K. Joel Roop-Eckart wrote the R code, conducted the analysis, and wrote the paper.
Ben Seiyon Lee contributed significantly to the likelihood functions and provided invaluable statistical coding guidance.
Caitlin Spence contributed to the problem framing and direction of the project.
Tess Russo contributed to the conception of the project.
Klaus Keller contributed to the conception of the project, the conceptual guidance on the analysis, and co-drafted the paper.

# Acknowledgements

We would like to thank John F. England for sharing the TREX modeled hydrograph data from (England et al., 2006) and for sharing the source code for the United States Bureau of Reclamation Flood Routing for Dams (FLROUT) program, citations and FLROUT link below.

England, J.F. Jr., Klawon, J.E., Klinger, R.E. and Bauer, T.R. (2006) Flood Hazard Study, Pueblo Dam, Colorado, Final Report, Bureau of Reclamation, Denver, CO, June, 160 p. and seven appendices.

United States Bureau of Reclamation, Technical Services Center, Waterways and Concrete Dams Group. (1993). Flood Routing for Dams (FLROUT). https://sites.google.com/a/alumni.colostate.edu/jengland/resources

# Compiling FLROUT

The source code, examples, the user manual, and compiling instructions for FLROUT can be found at (https://sites.google.com/a/alumni.colostate.edu/jengland/resources). This repository contains the uncompiled and compiled forms of the edited version used for the analysis.

To compile FLROUT in Windows for this analysis, using Gfortran, execute the following command in command line

gfortran flrout-kjr-edited.for -o flrout-kjr-edited.exe

Note: include or previously specify the proper file paths to gfortran and each FLROUT file location.

# Contact

K. Joel Roop-Eckart, E-mail: kjr30@psu.edu

Corresponding author:
Klaus Keller, E-mail: kzk10@psu.edu
