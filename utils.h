#ifndef __UTILS_H__
#define __UTILS_H__

#include "sislin.h"

void imprimeVetor(double *v, int n);

void multiplicaDiagVetor(double *matrizDiag, double *v, double *dest, int n);

void calculaCInv(SL *sl, double *cInv);

void multiplicaMatrizVetor(SL *sl, double *v, double *dest);


double multiplicaVetores(double *v1, double *v2, int n);


void calculaResiduo(SL *sl, double *x, double *r);


double calculaNormaEuclidiana(double *v, int n);
double calculaNormaMax(double *v1, double *v2, int n);
double calculaNormaMaxRelativa(double *v1, double *v2, int n);

#endif