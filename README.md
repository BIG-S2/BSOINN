# BSOLNN: Bayesian Scalar on Image Regression with Non-ignorable Non-response

# Table of contents
 [BSOINN-package]
 
    i. [Overview]
    
    ii. [Installing BSOINN]
    
    iii. [Democode-FPCA-on-images]
    
    iv. [Democode-simulation]
    
    v. [Democode-realdata]
    
    vi. [Implementing Democode]
    

# BSOINN-package <a name="The packageod"></a>

"BSOINN" package (package version 1.0) conducts Bayesian inference for a scalar on image regression model with non-ignorable non-response, which implements the method proposed in the paper "Bayesian scalar on image regression with non-ignorable non-response".

# Overview  <a name="Overview"></a>

This package contains three functions for the Bayesian scalar on image (BSOI) regression model. 

The three functions of the package are as follows:

1. The function "BSOIFull" conducts Bayesian inference for a BSOI with fully-observed response.

2. The function "BSOIIN" conducts Bayesian inference for a BSOI with ignorable non-response.

3. The function "BSOINN" conducts Bayesian inference for a BSOI with non-ignorable non-response.

The three functions are based on the functional principal component analysis (FPCA) for imaging data, which can be implemented using the "fast.svd" function in the package of "corpcor" (package version 1.6.9). Details of the FPCA can be found in the paper and the democode for the simulation study. A demo code that uses the ¡°fast.svd¡± function to perform FPCA on images is also introduced below. 

# Installing BSOINN  <a name="Installing BSOINN"></a>

* Installing BSOINN on Windows:
1. Download and install R software (R version 3.6.1) from http://cran.r-project.org.
2. Download and install Rtools (Rtools version 34) from  http://cran.r-project.org/bin/windows/Rtools/. During the installation process, please check the box provided to edit the system PATH, so that the C++ compiler included in Rtools can be used by R.
3. Download and install Rstudio software (Rstudio version 1.2.5001) from https://www.rstudio.com/.
4. Install packages "Rcpp" (package version 1.0.2) and "RcppArmadillo" (package version 0.9.700.2.0) from CRAN inside the R software.
5. Install package "BSOINN" from the local package archive "BSOINN_1.0.tar.gz".

* Installing BSOINN on Mac OS (versions of software and packages are the same as above):
1. Download and install R software.
2. Download and install C++ Toolchain, Xcode, from Mac "App Store". After the Xcode is installed, you need to open Xcode once to accept the license agreement.
3. Download and install Rstudio software.
4. Run the following code in R to put the path of "libgfortran"" into "FLIBS": 
```{r, eval=FALSE, warning=FALSE,message=FALSE}
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("FLIBS = -L`gfortran -print-file-name=libgfortran.dylib | xargs dirname`", file = M, sep = "\n", append = TRUE)
```
5. Install packages "Rcpp" and "RcppArmadillo" from CRAN inside the R software.
6. Install package "BSOINN" from the local package archive "BSOINN_1.0.tar.gz".

# Democode-FPCA-on-images <a name="Democode-FPCA-on-images"></a>

This demo code is for reading imaging data and performing FPCA on the images. Please refer to Section 2.1 of the manuscript for the description of FPCA.

# Democode-simulation <a name="Democode-simulation"></a>

The democode of parameter estimation for the simulation study from the paper can be found in the file "simulation_example_estimation.R" from the folder "Democode-simulation-estimation". 

The democode of out-of-sample prediction for the simulation study from the paper can be found in the file "simulation_example_prediction.R" from the folder "Democode-simulation-prediction". 

After installing the package "BSOINN" (this package) and "corpcor" (can be installed from CRAN), these R programs can be implemented. 

Please refer to the R files for detailed codes and explanations.

# Democode-realdata <a name="Democode-realdata"></a>

The democode of parameter estimation for the real data analysis using ADNI data set from the paper can be found in the file "realdata_example_estimation.R" from the folder "Democode-realdata-estimation". 

The democode of out-of-sample prediction for the real data analysis using ADNI data set from the paper can be found in the file "realdata_example_prediction.R" from the folder "Democode-realdata-prediction". 

Data from the ADNI study are included:

* "Y.txt": learning scores of the subjects at the 36th month, with non-response.

* "Y_true.txt": true learning scores of the subjects at the 36th month, with non-responses being filled using the most recently observed learning scores of the subjects.

* "COV.txt": gender, age, education (the instrument),  race, whether ever married, apoe gene, whether MCI or AD.

* "XI.txt": five eigenscores corresponding to the first five eigenimages of the imaging data.

* "R.txt": 0 denotes response and 1 denotes non-response.

After installing the package "BSOINN" (this package), these R programs can be implemented. 

Please refer to the R files for detailed codes and explanations.

# Implementing Democode <a name="Implementing Democode"></a>

1. Democode-FPCA-on-images:

* the code should be adjusted before use. The guideline is provided in the code. 

2. Democode-simulation-estimation:

* set the folder "Democode-simulation-estimation" as the working directory in R.

*  implementing the file "simulation_example_estimation.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("simulation_example_estimation.R")
```

*  The convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.

*  The approximate runtime for the demo code is 85 seconds for one replication. (We evaluate the approximate runtimes for the demo codes on a computer with a 64-bit six-core 2.2 GHz Intel Core i7-8750H processor with 16 GB of main memory)

3. Democode-simulation-prediction:

* set the folder "Democode-simulation-prediction" as the working directory in R.

* implementing the file "simulation_example_prediction.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("simulation_example_prediction.R")
```

* The out-of-sample prediction results are automatically shown on the R console.

* The approximate runtime for the demo code is 25 seconds for one replication.


4. Democode-realdata-estimation:

* set the folder "Democode-realdata-estimation" as the working directory in R.

* implementing the file "realdata_example_estimation.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("realdata_example_estimation.R")
```

* The convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.

* The approximate runtime for the demo code is 45 seconds.


5. Democode-realdata-prediction:

* set the folder "Democode-realdata-prediction" as the working directory in R.

* implementing the file "realdata_example_prediction.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("realdata_example_prediction.R")
```

* The out-of-sample prediction results are automatically shown on the R console.

* The approximate runtime for the demo code is 10 seconds for one random partition analysis.



