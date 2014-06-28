#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

double 	exactci_c(const double* , const double*, R_len_t);
SEXP 	exactciR2C(SEXP, SEXP);

double 	equalci_c(const double* , const double*, R_len_t);
SEXP 	equalciR2C(SEXP, SEXP);