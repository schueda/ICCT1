#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "sislin.h"

/*!
 \brief Função que gera os coeficientes de um sistema linear k-diagonal
 \param i Coordenada i do elemento a ser calculado 0<=i
 \param j Coordenada j do elemento a ser calculado j<n
 \param k Numero de diagonais da matriz A

 \return Valor do elemento Aij
*/
double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

/*!
 \brief Função que gera os termos independentes de um sistema linear k-diagonal
 \param k Numero de diagonais da matriz A

 \return Valor do elemento Bi
 */
double generateRandomB( unsigned int k )
{
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k<<2) * (double)rand() * invRandMax;
}

/*!
 \brief Função que gera e popula o sistema linear k-diagonal
 \param n Tamanho do sistema linear
 \param k Numero de diagonais da matriz A

 \return Ponteiro para o SL criado
 */
SL *criaSL(int n, int k) {
    int tam_linha = ceil((double) k/2);
    int offset = 0;
    int i;

    SL *sl = (SL *) malloc(sizeof(SL));
    sl->n = n;
    sl->k = k;
    sl->linhas = (double **) malloc(n * sizeof(double *));


    for (i = 0; i < floor(k/2); i++) {
        sl->linhas[i] = (double *) malloc(tam_linha * sizeof(double));

        for (int j = 0; j < tam_linha; j++) {
            sl->linhas[i][j] = generateRandomA(i, j, k);
        }
        tam_linha++;
    }

    for (i = floor(k/2); i < n-floor(k/2); i++) {
        sl->linhas[i] = (double *) malloc(k * sizeof(double));

        for (int j = 0; j < k; j++) {
            sl->linhas[i][j] = generateRandomA(i, j+offset, k);
        }
        offset++;
    }

    tam_linha--;
    for (i=n-floor(k/2); i < n; i++) {
        sl->linhas[i] = (double *) malloc(tam_linha * sizeof(double));

        for (int j = 0; j < tam_linha; j++) {
            sl->linhas[i][j] = generateRandomA(i, j+offset, k);
        }
        tam_linha--;
        offset++;
    }

    sl->b = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        sl->b[i] = generateRandomB(k);
    }

    return sl;
}

/*!
 \brief Função que imprime o sistema linear k-diagonal
 \param sl Ponteiro para o SL a ser impresso
 */
void imprimeSL(SL *sl) {
    int tam_linha = ceil((double) sl->k/2);
    for (int i = 0; i < floor(sl->k/2); i++) {
        for (int j = 0; j < tam_linha; j++) {
            printf("%f ", sl->linhas[i][j]);
        }
        printf("\n");
        tam_linha++;
    }

    for (int i = floor(sl->k/2); i < sl->n - floor(sl->k/2); i++) {
        for (int j = 0; j < sl->k; j++) {
            printf("%f ", sl->linhas[i][j]);
        }
        printf("\n");
    }

    tam_linha--;
    for (int i = sl->n - floor(sl->k/2); i < sl->n; i++) {
        for (int j = 0; j < tam_linha; j++) {
            printf("%f ", sl->linhas[i][j]);
        }
        printf("\n");
        tam_linha--;
    }

    printf("\n");
    for (int i = 0; i < sl->n; i++) {
        printf("%f ", sl->b[i]);
    }

    printf("\n");
}

/*!
 \brief Função que desaloca o sistema linear k-diagonal
 \param sl Ponteiro para o SL a ser desalocado
 */
void destroiSL(SL *sl) {
    for (int i = 0; i < sl->n; i++) {
        free(sl->linhas[i]);
    }

    free(sl->linhas);
    free(sl->b);
    free(sl);
}