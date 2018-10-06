#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <R.h>

float log2f(float);
double log2d(double);
void xToZ(const double*, double*, int, int, int, double, double);
double mean(double*, int);
double std(double*,int);
double meani(int*, int);
double stdi(int*,int);
/* void SplineKnots(int*,int,int); */
void knotVector(double*, int, int);
void findWeights(const double *, const double *, double *, int, int, int, double, double);
double entropy1(const double*, int, int);
double entropy2(const double*, const double*, int, int);
double mi2(const double*, const double*, int, int, int, int, int);


// export R function
void mi2R(const double *, const double *, int *, int *, int *, double *, int *, int *);
void entropy1R(const double *, int *, int *, int *, double *);
void entropy2R(const double *, const double *, int *, int *, int *, double *);
void centropy2R(const double *, const double *, int *, int *, int *, double *);
void mi2DiffBins(const double *, const double *, int *, int *, int *, int *, int *, double *, int *, int *);
void mi2vs1R(const double *, const double *, const double *, int *, int *, int *, double *, int *);
void getAllMIWz(const double *, const double* , double *, int , int , int , int , int , int );
void getAllMIWz_R(const double *, const double* , double *, int *, int *, int *, int *, int *, int *);
void getAllMI3Wz_R(const double *, const double* , const double* , double *, int *, int *, int *, int *, int *, int *);
void mi3(const double *x, const double *, const double *, int *, int *, int *, double *);
