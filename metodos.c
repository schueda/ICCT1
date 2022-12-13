#include <stdlib.h>
#include <stdio.h>
#include "metodos.h"
#include "utils.h"

void preCondicionado(SL *sl, double *x, double erro, int maxIt) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta, norma;

    double *r = (double *) malloc(sl->n * sizeof(double));
    double *rProx = (double *) malloc(sl->n * sizeof(double));
    double *z = (double *) malloc(sl->n * sizeof(double));
    double *zProx = (double *) malloc(sl->n * sizeof(double));
    double *d = (double *) malloc(sl->n * sizeof(double));
    double *Ad = (double *) malloc(sl->n * sizeof(double));
    double *cInv = (double *) malloc(sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    calculaCInv(sl, cInv);

    multiplicaDiagVetor(cInv, r, z, sl->n);
    copiaVetor(d, z, sl->n);

    while (k < maxIt) {
        multiplicaMatrizVetor(sl, d, Ad);

        alpha = multiplicaVetores(z, r, sl->n) / multiplicaVetores(d, Ad, sl->n);

        for (i = 0; i < sl->n; i++) {
            x[i] += alpha * d[i];
            rProx[i] = r[i] - alpha * Ad[i];
        }

        norma = calculaNormaMax(rProx, sl->n);
        if (norma < erro) {
            break;
        }

        multiplicaDiagVetor(cInv, r, zProx, sl->n);

        beta = multiplicaVetores(zProx, rProx, sl->n) / multiplicaVetores(z, r, sl->n);

        for (i = 0; i < sl->n; i++) {
            d[i] = zProx[i] + beta * d[i];
            r[i] = rProx[i];
            z[i] = zProx[i];
        }

        k++;
    }

    printf("norma: %f\n", norma);

    free(r);
    free(rProx);
    free(z);
    free(zProx);
    free(d);
    free(Ad);
    free(cInv);
}


void gradienteConjugado(SL *sl, double *x, double erro, int maxIt) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double *r = (double *) malloc(sl->n * sizeof(double));
    double *p = (double *) malloc(sl->n * sizeof(double));
    double *Ap = (double *) malloc(sl->n * sizeof(double));
    double *rProx = (double *) malloc(sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    double norma = calculaNormaEuclidiana(r, sl->n);
    if (norma < erro) {
        return;
    }

    copiaVetor(p, r, sl->n);

    while (k < maxIt && norma > erro) {
        multiplicaMatrizVetor(sl, p, Ap);

        alpha = multiplicaVetores(r, r, sl->n) / multiplicaVetores(p, Ap, sl->n);

        for (i = 0; i < sl->n; i++) {
            x[i] += alpha * p[i];
            rProx[i] = r[i] - alpha * Ap[i];
        }

        norma = calculaNormaEuclidiana(rProx, sl->n);
        if (norma < erro) {
            break;
        }

        beta = multiplicaVetores(rProx, rProx, sl->n) / multiplicaVetores(r, r, sl->n);

        for (i = 0; i < sl->n; i++) {
            r[i] = rProx[i];
            p[i] = r[i] + beta * p[i];
        }

        k++;
    }

    printf("norma: %f\n", norma);

    free(r);
    free(p);
    free(Ap);
    free(rProx);
}