#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "metodos.h"
#include "utils.h"

void gradienteConjugado(SL *sl, double *x, double erro, int maxIt) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double *xProx = (double *) malloc(sl->n * sizeof(double));

    double *r = (double *) malloc(sl->n * sizeof(double));
    double *p = (double *) malloc(sl->n * sizeof(double));
    double *Ap = (double *) malloc(sl->n * sizeof(double));
    double *rProx = (double *) malloc(sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    memcpy(p, r, sl->n * sizeof(double));

    while (k < maxIt) {
        multiplicaMatrizVetor(sl, p, Ap);

        alpha = multiplicaVetores(r, r, sl->n) / multiplicaVetores(p, Ap, sl->n);

        for (i = 0; i < sl->n; i++) {
            xProx[i] += alpha * p[i];
            rProx[i] = r[i] - alpha * Ap[i];
        }

        
        if (calculaNormaMaxRelativa(xProx, x, sl->n) < erro) {
            break;
        }
        memcpy(x, xProx, sl->n * sizeof(double));

        beta = multiplicaVetores(rProx, rProx, sl->n) / multiplicaVetores(r, r, sl->n);

        memcpy(r, rProx, sl->n * sizeof(double));
        for (i = 0; i < sl->n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        k++;
    }

    free(r);
    free(p);
    free(Ap);
    free(rProx);
}

void preCondicionado(SL *sl, double *x, double erro, int maxIt) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double *xProx = (double *) malloc(sl->n * sizeof(double));

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
    memcpy(d, z, sl->n * sizeof(double));

    while (k < maxIt) {
        multiplicaMatrizVetor(sl, d, Ad);

        alpha = multiplicaVetores(z, r, sl->n) / multiplicaVetores(d, Ad, sl->n);

        for (i = 0; i < sl->n; i++) {
            xProx[i] += alpha * d[i];
            rProx[i] = r[i] - alpha * Ad[i];
        }

        if (calculaNormaMaxRelativa(xProx, x, sl->n) < erro) {
            break;
        }
        memcpy(x, xProx, sl->n * sizeof(double));

        multiplicaDiagVetor(cInv, r, zProx, sl->n);

        beta = multiplicaVetores(zProx, rProx, sl->n) / multiplicaVetores(z, r, sl->n);

        for (i = 0; i < sl->n; i++) {
            d[i] = zProx[i] + beta * d[i];
        }
        
        memcpy(r, rProx, sl->n * sizeof(double));
        memcpy(z, zProx, sl->n * sizeof(double));


        k++;
    }

    free(r);
    free(rProx);
    free(z);
    free(zProx);
    free(d);
    free(Ad);
    free(cInv);
}