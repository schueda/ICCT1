#ifndef __UTILS_H__
#define __UTILS_H__

#include "sislin.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void imprimeVetor(double *v, int n, FILE *fp);
void imprimeVetorTeste(double *v, int n);

void multiplicaDiagVetor(double *matrizDiag, double *v, double *dest, int n);

void calculaMInv(SL *sl, double *cInv);

void multiplicaMatrizVetor(SL *sl, double *v, double *dest);

SL *simetrizaSL(SL *slA, SL *slB);

double multiplicaVetores(double *v1, double *v2, int n);

void obtemMatrizTransposta(SL *sl);

void calculaResiduo(SL *sl, double *x, double *r);


double calculaNormaEuclidiana(double *v, int n);
double calculaNormaMax(double *v1, double *v2, int n);
double calculaNormaMaxRelativa(double *v1, double *v2, int n);

double timestamp(void);

#endif