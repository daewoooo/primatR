<img src="https://github.com/daewoooo/primatR/raw/master/primatR_logo.png" />
=========================================================================

# primatR
An R package of useful functions to process and analyze Strand-seq data.

Installation
------------

### Development version from Github
To install the development version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>=3.5.0) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   #### To install required packages  
   if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
   install("GenomicRanges")
	 install.packages("devtools")
	 library(devtools)  

4. To install breakpointR package from github	 
	 #### Option1
	 install_github("daewoooo/primatR", force=TRUE)  
	 #### Option2 
	 install_git("git://github.com/daewoooo/primatR.git", branch = "master")  

Report Errors
-------------

If you encounter errors of any kind, please report an [issue here](https://github.com/daewoooo/primatR/issues/new).