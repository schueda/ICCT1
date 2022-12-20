#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "sislin.h"
#include "utils.h"

/*!
 \brief Função que gera os coeficientes de um sistema linear k-diagonal
 \param i Coordenada i do elemento a ser calculado 0<=i
 \param j Coordenada j do elemento a ser calculado j<n
 \param k Numero de diagonais da matriz A

 \return Valor do elemento Aij
*/
double generateRandomA(unsigned int i, unsigned int j, unsigned int k) {
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

/*!
 \brief Função que gera os termos independentes de um sistema linear k-diagonal
 \param k Numero de diagonais da matriz A

 \return Valor do elemento Bi
 */
double generateRandomB(unsigned int k) {
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k<<2) * (double)rand() * invRandMax;
}

/*!
 \brief Função que aloca o sistema linear k-diagonal
 \param n Tamanho do sistema linear
 \param k Numero de diagonais da matriz A

 \return Ponteiro para o SL criado
 */
SL *alocaSL(int n, int k) {
    int i;
    int tam_linha = ceil((double) k/2);

    SL *sl = (SL *) malloc(sizeof(SL));
    sl->n = n;
    sl->k = k;
    sl->A = (double *) calloc(n * k, n * k * sizeof(double *));
    sl->b = (double *) malloc(n * sizeof(double));
    
    return sl;
}

void populaSL(SL *sl) {
    int tam_linha = ceil((double) sl->k/2);
    int offset = 0;
    int i, j;

    for (i = 0; i < sl->n; i++) {
        tam_linha = MIN(tam_linha + i, sl->k);
        tam_linha = MIN(tam_linha, sl->n-1 - i + ceil((double) sl->k/2));

        offset = MAX(0, i - (int) floor((double) sl->k/2));

        for (j = 0; j < tam_linha; j++) {
            sl->A[i*sl->k + j] = generateRandomA(i, j + offset, sl->k);
        }
    }
    
    for (int i = 0; i < sl->n; i++) {
        sl->b[i] = generateRandomB(sl->k);
    }
}

void copiaSL(SL *slDest, SL *slOrigin) {
    slDest->n = slOrigin->n;
    slDest->k = slOrigin->k;

    memcpy(slDest->A, slOrigin->A, slDest->n * slDest->k * sizeof(double));
    memcpy(slDest->b, slOrigin->b, slDest->n * sizeof(double));
}


/*!
 \brief Função que imprime o sistema linear k-diagonal
 \param sl Ponteiro para o SL a ser impresso
 */
void imprimeSL(SL *sl) {
    for (int i = 0; i < sl->n; i++) {
        for (int j = 0; j < sl->k; j++) {
            printf("%f ", sl->A[i*sl->k + j]);
        }
        printf("\n");
    }

    for (int i = 0; i < sl->n; i++) {
        printf("%f ", sl->b[i]);
    }

    printf("\n\n");
}

/*!
 \brief Função que desaloca o sistema linear k-diagonal
 \param sl Ponteiro para o SL a ser desalocado
 */
void destroiSL(SL *sl) {
    free(sl->A);
    free(sl->b);
    free(sl);
}