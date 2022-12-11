#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "sislin.h"

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

double calculaNormaL2(double *r, int n) {
    double soma = 0;
    for (int i = 0; i < n; i++) {
        soma += r[i] * r[i];
    }
    return sqrt(soma);
}

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

void gradienteConjugado(SL *sl, double *x, double erro, int maxIt) {
    int k = 0, tam_linha, offset, i, j;

    double alpha, beta;

    double *r = (double *) malloc(sl->n * sizeof(double));
    double *p = (double *) malloc(sl->n * sizeof(double));
    double *Ap = (double *) malloc(sl->n * sizeof(double));
    double *rProx = (double *) malloc(sl->n * sizeof(double));

    calculaResiduo(sl, x, r);

    double norma = calculaNormaL2(r, sl->n);
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

        norma = calculaNormaL2(rProx, sl->n);
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

    free(r);
    free(p);
    free(Ap);
    free(rProx);
}



void imprimeVetor(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%f ", v[i]);
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