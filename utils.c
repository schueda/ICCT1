#include <stdio.h>
#include <math.h>
#include "utils.h"

void imprimeVetor(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%f ", v[i]);
    }
    printf("\n");
}

void multiplicaDiagVetor(double *matrizDiag, double *v, double *dest, int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = matrizDiag[i] * v[i];
    }
}

void calculaCInv(SL *sl, double *cInv) {
    int posDiag = 0;

    int i;
    for (i = 0; i < floor(sl->k/2); i++) {
        cInv[i] = 1/sl->linhas[i][posDiag];
        posDiag++;
    }

    for (i = floor(sl->k/2); i < sl->n; i++) {
        cInv[i] = 1/sl->linhas[i][posDiag];
    }
}

void multiplicaMatrizVetor(SL *sl, double *v, double *dest) {
    int tam_linha = ceil((double) sl->k/2);
    int i, j, offset;
    for (i = 0; i < floor(sl->k/2); i++) {
        dest[i] = 0;
        for (j = 0; j < tam_linha; j++) {
            dest[i] += sl->linhas[i][j] * v[j];
        }
        tam_linha++;
    }

    offset = 0;
    for (i = floor(sl->k/2); i < sl->n - floor(sl->k/2); i++) {
        dest[i] = 0;
        for (j = 0; j < sl->k; j++) {
            dest[i] += sl->linhas[i][j] * v[j+offset];
        }
        offset++;
    }

    tam_linha--;
    for (i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        dest[i] = 0;
        for (j = 0; j < tam_linha; j++) {
            dest[i] += sl->linhas[i][j] * v[j+offset];
        }
        tam_linha--;
        offset++;
    }
}

void copiaVetor(double *dest, double *orig, int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = orig[i];
    }
}

double multiplicaVetores(double *v1, double *v2, int n) {
    double soma = 0;
    for (int i = 0; i < n; i++) {
        soma += v1[i] * v2[i];
    }
    return soma;
}

/*!
 \brief Função que calcula o residuo de um sistema linear k-diagonal
 \param sl Sistema linear
 \param x Vetor solução do sistema linear
 \param r Vetor que receberá os valores do resíduo do sistema linear
*/
void calculaResiduo(SL *sl, double *x, double *r) {
    int tam_linha = ceil((double) sl->k/2), i, j, offset;

    for (i = 0; i < floor(sl->k/2); i++) {
        r[i] = sl->b[i];
        for (j = 0; j < tam_linha; j++) {
            r[i] -= sl->linhas[i][j] * x[j];
        }
        tam_linha++;
    }

    offset = 0;
    for (i = floor(sl->k/2); i < sl->n - floor(sl->k/2); i++) {
        r[i] = sl->b[i];
        for (j = 0; j < sl->k; j++) {
            r[i] -= sl->linhas[i][j] * x[j + offset];
        }
        offset++;
    }

    tam_linha--;
    for (i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        r[i] = sl->b[i];
        for (j = 0; j < tam_linha; j++) {
            r[i] -= sl->linhas[i][j] * x[j + offset];
        }
        tam_linha--;
        offset++;
    }
}

double calculaNormaEuclidiana(double *v, int n) {
    double soma = 0;
    for (int i = 0; i < n; i++) {
        soma += v[i] * v[i];
    }
    return sqrt(soma);
}

double calculaNormaMax(double *v, int n) {
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(v[i]) > max) {
            max = fabs(v[i]);
        }
    }
    return max;
}