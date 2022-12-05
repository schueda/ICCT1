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
 ***********************/
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

// Cria sistema linear
SL *criaSL(int n, int k) {
    SL *sl = (SL *) malloc(sizeof(SL));
    sl->n = n;
    sl->k = k;
    sl->diagonais = (double **) malloc(k * sizeof(double *));

    int sub_tam = (int) k/2;
    for (int i = 0; i < k; i++) {
        sl->diagonais[i] = (double *) malloc((n - abs(sub_tam)) * sizeof(double));
        for (int j = 0; j < n - sub_tam; j++) {
            sl->diagonais[i][j] = generateRandomA(i, j, k);
        }
        sub_tam--;
    }

    sl->b = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        sl->b[i] = generateRandomB(k);
    }

    return sl;
}

// Imprime sistema linear
void imprimeSL(SL *sl) {
    int sub_tam = (int) sl->k/2;
    for (int i=0; i<sl->k; i++) {
        for (int j=0; j<sl->n - abs(sub_tam); j++) {
            printf("%f ", sl->diagonais[i][j]);
        }
        printf("\n");
        sub_tam--;
    }
    printf("\n");
    for (int i=0; i<sl->n; i++) {
        printf("%f ", sl->b[i]);
    }
}


// Destroi sistema linear
void destroiSL(SL *sl) {
    for (int i = 0; i < sl->k; i++) {
        free(sl->diagonais[i]);
    }
    free(sl->diagonais);
    free(sl->b);
    free(sl);
}