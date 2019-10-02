CODE and DATA for paper "Bayesian Scalar on Image Regression with Non-ignorable Non-response"

This is the ReadMe document for running the data analysis presented in the paper.

*************** 
Platform: The implementation was done in R (R version 3.6.1).

*************** 
"BSOINN" package (package version 1.0):

Main functions are included in the R package "BSOINN", please install the package first to implement the demo codes.

Please note that the brain images are too large to be included in the package or in the folder of demo code for the real data analysis. We therefore require the eigenscores from the functional principal component analysis (FPCA) on the images as inputs for the functions in the “BSOINN” package. The FPCA in the paper can be easily implemented using matlab or the “fast.svd” function in the R `corpcor' package (package version 1.6.9), which is demonstrated in the demo code for the simulation study.  We also include a demo code that uses the “fast.svd” function to perform FPCA on images. 

****************
Installing "BSOINN":

on Windows:

1. Download and install R software (R version 3.6.1) from http://cran.r-project.org.
2. Download and install Rtools (Rtools version 34) from http://cran.r-project.org/bin/windows/Rtools/. During the installation process, please check the box provided to edit the system PATH, so that the C++ compiler included in Rtools can be used by R.
3. Download and install Rstudio software (Rstudio version 1.2.5001) from https://www.rstudio.com/.
4. Install packages “Rcpp” (package version 1.0.2) and “RcppArmadillo” (package version 0.9.700.2.0) from CRAN inside the R software.
5. Install package “BSOINN” from the local package archive “BSOINN_1.0.tar.gz”.

on Mac OS (versions of software and packages are the same as above):

1. Download and install R software.
2. Download and install C++ Toolchain, Xcode, from Mac “App Store”. After the Xcode is installed, you need to open Xcode once to accept the license agreement.
3. Download and install Rstudio software.
4. Run the following code in R to put the path of “libgfortran”" into “FLIBS”:

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("FLIBS = -L`gfortran -print-file-name=libgfortran.dylib | xargs dirname`", file = M, sep = "\n", append = TRUE)

5. Install packages “Rcpp” and “RcppArmadillo” from CRAN inside the R software.
6. Install package “BSOINN” from the local package archive “BSOINN_1.0.tar.gz”.

****************
Democode-simulation:

Demo code of parameter estimation for the simulation study is in the folder "Democode-simulation-estimation".
*Please also install the package "corpcor" from CRAN to use the function "fast.svd"
*The program "simulation_example_estimation.R" generates simulated data sets and conducts the analysis.
*When implementing the code, the convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.
*The approximate runtime for the demo code is 85 seconds for one replication. (We evaluate the approximate runtimes for the demo codes on a computer with a 64-bit six-core 2.2 GHz Intel Core i7-8750H processor with 16 GB of main memory)

Demo code of out-of-sample prediction for the simulation study is in the folder "Democode-simulation-prediction".
*Please also install the package "corpcor" from CRAN to use the function "fast.svd"
*The program "simulation_example_prediction.R" generates simulated data sets and conducts the out-of-sample prediction.
*The approximate runtime for the demo code is 25 seconds for one replication.

****************
Democode-realdata:

Demo code of parameter estimation for the real data analysis is in the folder "Democode-realdata-estimation".
*The program "realdata_example_estimation.R" reads data (based on the ADNI study) from files :
--- Y.txt: learning scores of the subjects at the 36th month, with missingness
--- COV.txt:  gender, age, education, race, whether ever married, apoe4,  whether MCI or AD   
--- XI.txt: five eigenscores corresponding to the first five eigenimages of the RAVENS data.
--- R.txt: 0, observed, 1, missing
*When implementing the code, the convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.
*The approximate runtime for the demo code is 45 seconds.

Demo code of out-of-sample prediction for the real data analysis is in the folder "Democode-realdata-prediction".
*The program "realdata_example_prediction.R" reads data (based on the ADNI study) from files :
--- Y.txt: learning scores of the subjects at the 36th month, with missingness
--- Y_true.txt: true learning scores of the subjects at the 36th month, with non-responses being filled using the most recently observed learning scores of the subjects.
--- COV.txt:  gender, age, education, race, whether ever married, apoe4,  whether MCI or AD   
--- XI.txt: five eigenscores corresponding to the first five eigenimages of the RAVENS data.
--- R.txt: 0, observed, 1, missing
*The approximate runtime for the demo code is 10 seconds for one random partition analysis.

****************
Democode-FPCA-on-images:

Demo code of reading RAVENS images and performing FPCA on the images is in the folder "Democode-FPCA-on-images".
*The program "Democode_FPCA_fastsvd.R" demonstrates how to read images and use "fast.svd" function to perform FPCA on images.
*Note that the program should be adjusted before use. 

