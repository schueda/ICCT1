// Autores: André Schueda Menezes e Marcus Augusto Ferreira Dudeque

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"


/*!
 \brief Função que simplifica o acesso ao valor de um elemento da matriz A,
        permitindo que a matriz seja acessada como se fosse uma matriz normal.
        Para valores fora do escopo da matriz é retornado o valor 0.

 \param sl Ponteiro para o sistema linear que contém a matriz A.
 \param i Coordenada i do elemento a ser acessado.
 \param j Coordenada j do elemento a ser acessado.

 \return Valor do elemento Aij.
*/
double valorNaMatriz(SL *sl, int i, int j) {
    int offset;

    int tam_linha = ceil((double) sl->k/2);
    tam_linha = MIN(tam_linha + i, sl->k);
    tam_linha = MIN(tam_linha, sl->n-1 - i + ceil((double) sl->k/2));


    offset = MAX(i - (int) floor((double) sl->k/2), 0);

    if (j - offset >= tam_linha) return 0;

    if (j - offset < 0) return 0;

    return sl->A[i*sl->k + j-offset];
}


/*!
 \brief Função que simplifica o acesso ao endereço de um elemento da matriz A,
        permitindo que a matriz seja acessada como se fosse uma matriz normal.
        Para valores fora do escopo da matriz é retornado NULL.

 \param sl Ponteiro para o sistema linear que contém a matriz A.
 \param i Coordenada i do elemento a ser acessado.
 \param j Coordenada j do elemento a ser acessado.

 \return Endereço do elemento Aij.
*/
double *itemNaMatriz(SL *sl, int i, int j) {
    int offset;

    int tam_linha = ceil((double) sl->k/2);
    tam_linha = MIN(tam_linha + i, sl->k);
    tam_linha = MIN(tam_linha, sl->n-1 - i + ceil((double) sl->k/2));


    offset = MAX(i - (int) floor((double) sl->k/2), 0);

    if (j - offset >= tam_linha) return NULL;

    if (j - offset < 0) return NULL;

    return &sl->A[i*sl->k + j-offset];
}


/*!
 \brief Função que copia uma linha da matriz A do sl ao vetor dest.

 \param sl Ponteiro para o sistema linear que contém a matriz A.
 \param i Coordenada i da linha a ser copiada.
 \param dest Vetor que receberá a linha.
*/
void obtemLinha(SL *sl, int i, double *dest) {
    int j;
    for (j = 0; j < sl->k; j++) {
        dest[j] = sl->A[i*sl->k + j];
    }
}


/*!
 \brief Função que copia uma coluna da matriz A do sl ao vetor dest.

 \param sl Ponteiro para o sistema linear que contém a matriz A.
 \param j Coordenada j da coluna a ser copiada.
 \param ini Coordenada i inicial da coluna a ser copiada.
 \param fim Coordenada i final da coluna a ser copiada.
 \param dest Vetor que receberá a coluna.
*/
void obtemColuna(SL *sl, int j, int ini, int fim, double *dest) {
    int i;
    for (i = ini; i <= fim; i++) {
        dest[i-ini] = valorNaMatriz(sl, i, j);
    }
}



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


void calculaMInv(SL *sl, double *mInv) {
    int j = 0;

    int i;
    for (i = 0; i < (int) floor((double) sl->k/2); i++) {
        mInv[i] = 1/sl->A[i*sl->k + j];
        j++;
    }

    for (i = (int) floor((double) sl->k/2); i < sl->n; i++) {
        mInv[i] = 1/sl->A[i*sl->k + j];
    }
}


void multiplicaMatrizVetor(SL *sl, double *v, double *dest) {
    int tam_linha = ceil((double) sl->k/2);
    int i, j, offset;

    for (i = 0; i < sl->n; i++) {
        offset = MAX(i - (int) floor((double) sl->k/2), 0);
        dest[i] = 0;
        for (j = 0; j < sl->k; j++) {
            dest[i] += sl->A[i*sl->k + j] * v[j+offset];
        }
    }
}


void obtemMatrizTransposta(SL *sl) {
    double *sup, *inf;

    double aux;
    int i, j;
    int tam_linha;

    for (i = 0; i < sl->n; i++) {
        tam_linha = MIN(floor((double) sl->k/2), sl->n-1 - i);

        for (j = i+1; j < i+1 + tam_linha; j++) {
            sup = itemNaMatriz(sl, i, j);
            inf = itemNaMatriz(sl, j, i);

            if (sup != NULL && inf != NULL) {
                aux = *sup;
                *sup = *inf;
                *inf = aux;
            }
        }
    }
}


SL *multiplicaMatrizSL(SL *slA, SL *slB) {
    int tam_linha, i, j;
    int offset, resultOffset;

    int k = slA->k;
    int n = slA->n;
    int newK = 2 * k - 1;

    SL *slDest = alocaSL(n, newK);

    double *lin = (double *) calloc(k, sizeof(double));
    double *col = (double *) calloc(k, sizeof(double));


    tam_linha = ceil((double) newK/2);

    for (i = 0; i < n; i++) {
        offset = MAX(0, i - (int) floor((double) k/2));
        resultOffset = MAX(0, i - (int) floor((double) newK/2));

        for (j = 0; j < newK; j++) {
            obtemLinha(slA, i, lin);
            obtemColuna(slB, j + resultOffset, offset, offset + k-1, col);

            slDest->A[i*newK + j] = multiplicaVetores(lin, col, k);
        }
    }

    multiplicaMatrizVetor(slA, slB->b, slDest->b);

    free(lin);
    free(col);

    return slDest;
}


double multiplicaVetores(double *v1, double *v2, int n) {
    double soma = 0;
    for (int i = 0; i < n; i++) {
        soma += v1[i] * v2[i];
    }
    return soma;
}


void calculaResiduo(SL *sl, double *x, double *r) {
    int i, j;
    int offset;

    for (i = 0; i<sl->n; i++) {
        offset = MAX(i - (int) floor((double) sl->k/2), 0);
        r[i] = sl->b[i];
        for (j = 0; j < sl->k; j++) {
            r[i] -= sl->A[i*sl->k + j] * x[j+offset];
        }
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