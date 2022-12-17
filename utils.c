#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"

void imprimeVetor(double *v, int n, FILE *fp) {
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f ", v[i]);
    }
    fprintf(fp, "\n");
}

void imprimeVetorTeste(double *v, int n) {
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

void calculaMInv(SL *sl, double *cInv) {
    int posDiag = 0;

    int i;
    for (i = 0; i < floor(sl->k/2); i++) {
        cInv[i] = 1/sl->A[i][posDiag];
        posDiag++;
    }

    for (i = floor(sl->k/2); i < sl->n; i++) {
        cInv[i] = 1/sl->A[i][posDiag];
    }
}

double valorNaMatriz(SL *sl, int i, int j) {
    int offset;

    int tam_linha = ceil((double) sl->k/2);
    tam_linha = MIN(tam_linha + i, sl->k);
    tam_linha = MIN(tam_linha, sl->n-1 - i + ceil((double) sl->k/2));


    offset = MAX(i - floor(sl->k/2), 0);

    if (j - offset >= tam_linha) return 0;

    if (j - offset < 0) return 0;

    return sl->A[i][j-offset];
}

void obtemColuna(SL *sl, int j, int ini, int fim, double *dest) {
    int i;
    for (i = ini; i <= fim; i++) {
        dest[i-ini] = valorNaMatriz(sl, i, j);;
    }
}

SL *simetrizaSL(SL *slA, SL *slB) {
    int tam_linha, tam_coluna, i, j;
    int offset, resultOffset = 0;

    SL *slDest = alocaSL(slA->n, slA->k + (int) 2 * floor((double) slA->k/2));
    double *col = (double *) calloc(slB->k, slB->k * sizeof(double));

    tam_linha = ceil((double) slDest->k/2);
    tam_coluna = ceil((double) slB->k/2);

    for (i = 0; i < floor(slDest->k/2); i++) {
        for (j = 0; j < tam_linha; j++) {
            obtemColuna(slB, j, 0, tam_coluna-1, col);

            slDest->A[i][j] = multiplicaVetores(slA->A[i], col, tam_linha);
        }
        tam_linha++;
        tam_coluna++;
    }

    for (i = floor(slDest->k/2); i < slDest->n - floor(slDest->k/2); i++) {
        offset = i - (int) floor((double) slA->k/2);
        for (j = 0; j < slDest->k; j++) {
            obtemColuna(slB, j + resultOffset, offset, offset + slA->k-1, col);

            slDest->A[i][j] = multiplicaVetores(slA->A[i], col, slA->k);
        };
        resultOffset++;
    }

    tam_linha--;
    tam_coluna--;
    for (i = slDest->n - floor(slDest->k/2); i < slDest->n; i++) {
        offset = i - (int) floor((double) slA->k/2);
        for (j = 0; j < tam_linha; j++) {
            obtemColuna(slB, j + resultOffset, offset, offset + tam_coluna-1, col);

            slDest->A[i][j] = multiplicaVetores(slA->A[i], col, tam_coluna);
        }
        tam_linha--;
        tam_coluna--;
        resultOffset++;
    }

    multiplicaMatrizVetor(slA, slB->b, slDest->b);

    free(col);

    return slDest;
}

void multiplicaMatrizVetor(SL *sl, double *v, double *dest) {
    int tam_linha = ceil((double) sl->k/2);
    int i, j, offset;
    for (i = 0; i < floor(sl->k/2); i++) {
        dest[i] = 0;
        for (j = 0; j < tam_linha; j++) {
            dest[i] += sl->A[i][j] * v[j];
        }
        tam_linha++;
    }

    offset = 0;
    for (i = floor(sl->k/2); i < sl->n - floor(sl->k/2); i++) {
        dest[i] = 0;
        for (j = 0; j < sl->k; j++) {
            dest[i] += sl->A[i][j] * v[j+offset];
        }
        offset++;
    }

    tam_linha--;
    for (i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        dest[i] = 0;
        for (j = 0; j < tam_linha; j++) {
            dest[i] += sl->A[i][j] * v[j+offset];
        }
        tam_linha--;
        offset++;
    }
}

double multiplicaVetores(double *v1, double *v2, int n) {
    double soma = 0;
    for (int i = 0; i < n; i++) {
        soma += v1[i] * v2[i];
    }
    return soma;
}

double *itemNaMatriz(SL *sl, int i, int j) {
    int offset;

    if (i < floor(sl->k/2)) {
        return &sl->A[i][j];
    }

    offset = i - floor(sl->k/2);
    return &sl->A[i][j-offset];
}

void obtemMatrizTransposta(SL *sl) {
    double *sup, *inf;

    double aux;
    int i, j;
    int dec = 1;


    for (i = 0; i < sl->n - floor(sl->k/2); i++) {
        for (j = i + 1; j < i+1 + floor(sl->k/2); j++) {
            sup = itemNaMatriz(sl, i, j);
            inf = itemNaMatriz(sl, j, i);

            aux = *sup;
            *sup = *inf;
            *inf = aux;
        }
    }

    for (i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        for (j = i+1; j < i+1 + floor(sl->k/2) - dec; j++) {
            sup = itemNaMatriz(sl, i, j);
            inf = itemNaMatriz(sl, j, i);

            aux = *sup;
            *sup = *inf;
            *inf = aux;
        }
        dec++;
    }

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
            r[i] -= sl->A[i][j] * x[j];
        }
        tam_linha++;
    }

    offset = 0;
    for (i = floor(sl->k/2); i < sl->n - floor(sl->k/2); i++) {
        r[i] = sl->b[i];
        for (j = 0; j < sl->k; j++) {
            r[i] -= sl->A[i][j] * x[j + offset];
        }
        offset++;
    }

    tam_linha--;
    for (i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        r[i] = sl->b[i];
        for (j = 0; j < tam_linha; j++) {
            r[i] -= sl->A[i][j] * x[j + offset];
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

double calculaNormaMax(double *v1, double *v2, int n) {
    double max = 0.0;
    double aux;
    for (int i = 0; i < n; i++) {
        aux = fabs(v1[i] - v2[i]);
        if (aux > max) {
            max = aux;
        }
    }
    return max;
}

double calculaNormaMaxRelativa(double *v1, double *v2, int n) {
    double max = 0.0;
    double aux;
    for (int i = 0; i < n; i++) {
        aux = fabs(v1[i] - v2[i]) / fabs(v1[i]);
        if (aux > max) {
            max = aux;
        }
    }
    return max;
}


double timestamp(void) {
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
}