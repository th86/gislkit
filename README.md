gislkit
====


## Installation ##

We can install Gislkit using the following commands:
```r
#install.packages("devtools")
library(devtools)
install_github("th86/gislkit")
```

A Windows user may need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) and use the following command to install the package locally:
```r
install.packages("gislkit-master", repo=NULL, type="source")
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
