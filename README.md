gislkit
====


## Installation ##

We can install Gislkit using the following commands:
```r
#install.packages("devtools")
library(devtools)
install_github("th86/gislkit")
```

If you encounter the error messages about RCurl, installing the following packages may fix the problem.
```r
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

## Usage ##

To find the attractors in a gene expression matrix ge using CENPA, PTPRC, COL3A1 as seeds:
```r
attractorList<-attractorSearch(ge, c("CENPA","PTPRC","COL3A1"))
```
The function returns a list of converged attractors.


To find all the genomically localized attractors in a gene expression matrix:
```r
data(grch37.geneSymbol)
GLattractorList<-GLattractorSearch(ge,genome=grch37_genesymbol)
```
The function returns a list of converged genomically localized attractors.


