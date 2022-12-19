#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "metodos.h"
#include "utils.h"

void gradienteConjugado(SL *sl, double *x, double erro, int maxIt, FILE *fp) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double norma;

    double timeStampI, timeStampF, timeSum = 0.0;

    double *xProx = (double *) calloc(sl->n, sl->n * sizeof(double));

    double *r = (double *) calloc(sl->n, sl->n * sizeof(double));
    double *p = (double *) calloc(sl->n, sl->n * sizeof(double));
    double *Ap = (double *) calloc(sl->n, sl->n * sizeof(double));
    double *rProx = (double *) calloc(sl->n, sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    memcpy(p, r, sl->n * sizeof(double));

    while (k < maxIt) {
        timeStampI = timestamp();
        multiplicaMatrizVetor(sl, p, Ap);

        alpha = multiplicaVetores(r, r, sl->n) / multiplicaVetores(p, Ap, sl->n);

        for (i = 0; i < sl->n; i++) {
            xProx[i] = x[i] + alpha * p[i];
            rProx[i] = r[i] - alpha * Ap[i];
        }

        
        if (calculaNormaMaxRelativa(xProx, x, sl->n) < erro) {
            break;
        }

        beta = multiplicaVetores(rProx, rProx, sl->n) / multiplicaVetores(r, r, sl->n);

        memcpy(r, rProx, sl->n * sizeof(double));
        for (i = 0; i < sl->n; i++) {
            p[i] = r[i] + beta * p[i];
        }
        
        timeStampF = timestamp();

        k++;

        fprintf(fp, "# iter %d: %.15g\n", k, calculaNormaMax(xProx, x, sl->n));
        memcpy(x, xProx, sl->n * sizeof(double));

        timeSum += timeStampF - timeStampI;
    }

    timeStampI = timestamp();
    norma = calculaNormaEuclidiana(r, sl->n);
    timeStampF = timestamp();
    fprintf(fp, "# residuo: %.15g\n", norma);

    fprintf(fp, "# Tempo PC: %.15g\n", 0.0);
    fprintf(fp, "# Tempo iter: %.15g\n", timeSum/k);
    fprintf(fp, "# Tempo residuo: %.15g\n", timeStampF - timeStampI);

    free(xProx);

    free(r);
    free(p);
    free(Ap);
    free(rProx);
}

void preCondicionado(SL *sl, double *x, double erro, int maxIt, FILE *fp) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double norma;

    double timeStampI, timeStampF, timeSum = 0.0, tempoPC;

    double *xProx = (double *) malloc(sl->n * sizeof(double));

    double *r = (double *) malloc(sl->n * sizeof(double));
    double *rProx = (double *) malloc(sl->n * sizeof(double));
    double *z = (double *) malloc(sl->n * sizeof(double));
    double *zProx = (double *) malloc(sl->n * sizeof(double));
    double *p = (double *) malloc(sl->n * sizeof(double));
    double *Ap = (double *) malloc(sl->n * sizeof(double));
    double *mInv = (double *) malloc(sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    timeStampI = timestamp();
    calculaMInv(sl, mInv);
    timeStampF = timestamp();
    tempoPC = timeStampF - timeStampI;

    multiplicaDiagVetor(mInv, r, z, sl->n);
    memcpy(p, z, sl->n * sizeof(double));

    while (k < maxIt) {
        timeStampI = timestamp();
        multiplicaMatrizVetor(sl, p, Ap);

        alpha = multiplicaVetores(r, z, sl->n) / multiplicaVetores(p, Ap, sl->n);

        for (i = 0; i < sl->n; i++) {
            xProx[i] = x[i] + alpha * p[i];
            rProx[i] = r[i] - alpha * Ap[i];
        }

        if (calculaNormaMaxRelativa(xProx, x, sl->n) < erro) {
            break;
        }

        multiplicaDiagVetor(mInv, rProx, zProx, sl->n);

        beta = multiplicaVetores(rProx, zProx, sl->n) / multiplicaVetores(r, z, sl->n);

        for (i = 0; i < sl->n; i++) {
            p[i] = zProx[i] + beta * p[i];
        }

        memcpy(r, rProx, sl->n * sizeof(double));
        memcpy(z, zProx, sl->n * sizeof(double));
        
        timeStampF = timestamp();

        k++;

        fprintf(fp, "# iter %d: %.15g\n", k, calculaNormaMaxRelativa(xProx, x, sl->n));
        memcpy(x, xProx, sl->n * sizeof(double));

        timeSum += timeStampF - timeStampI;
    }

    timeStampI = timestamp();
    norma = calculaNormaEuclidiana(r, sl->n);
    timeStampF = timestamp();
    fprintf(fp, "# residuo: %.15g\n", norma);

    fprintf(fp, "# Tempo PC: %.15g\n", tempoPC);
    fprintf(fp, "# Tempo iter: %.15g\n", timeSum/k);
    fprintf(fp, "# Tempo residuo: %.15g\n", timeStampF - timeStampI);

    free(r);
    free(rProx);
    free(z);
    free(zProx);
    free(p);
    free(Ap);
    free(mInv);
}