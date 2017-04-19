gislkit
====

Download this package from https://github.com/th86/gislkit

Uncompress this package gislkit-master.zip as a directory

Install the package using the following command

```r
install.packages("gislkit-master", repo=NULL, type="source")
```

Windows users need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) beforehand.

## Usage ##

Find the attractors in a gene expression matrix ge using CENPA, PTPRC, COL3A1 as seeds.
```r
attractorList<-attractorSearch(ge, c("CENPA","PTPRC","COL3A1"))
```
The function returns a list of converged attractors.


Find all the genomically localized attractors in a gene expression matrix
```r
data(grch37.geneSymbol)
GLattractorList<-GLattractorSearch(ge,genome=grch37_genesymbol)
```
The function returns a list of converged genomically localized attractors.
