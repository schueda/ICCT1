// Autores: André Schueda Menezes e Marcus Augusto Ferreira Dudeque

#ifndef __UTILS_H__
#define __UTILS_H__

#include "sislin.h"

/*!
 \brief Macro que retorna o maior valor entre dois valores.

 \param x Primeiro valor.
 \param y Segundo valor.

 \return Maior valor entre x e y.
*/
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


/*!
 \brief Macro que retorna o menor valor entre dois valores.

 \param x Primeiro valor.
 \param y Segundo valor.

 \return Menor valor entre x e y.
*/
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


/*!
 \brief Função que imprime os valores de um vetor em um arquivo.

 \param v Ponteiro para o vetor a ser impresso.
 \param n Tamanho do vetor.
 \param fp Ponteiro para o arquivo no qual o vetor será impresso.
*/
void imprimeVetor(double *v, int n, FILE *fp);


/*!
 \brief Função que imprime os valores de um vetor na saída padrão.

 \param v Ponteiro para o vetor a ser impresso.
 \param n Tamanho do vetor.
*/
void imprimeVetorTeste(double *v, int n);


/*!
 \brief Função que multiplica uma matriz diagonal com um vetor. o resultado é
        armazenado no vetor dest.

 \param matrizDiag Ponteiro para a matriz a ser multiplicada.
 \param v Ponteiro para o vetor a ser multiplicado.
 \param dest Ponteiro para o vetor que receberá o resultado da multiplicação.
 \param n Tamanho da matriz e do vetor.
*/
void multiplicaDiagVetor(double *matrizDiag, double *v, double *dest, int n);


/*!
 \brief Função que calcula o pré-condicionador de Jacobi, que consiste em
        calcular a inversa da matriz diagonal. 

 \param sl Ponteiro para o sistema linear de onde será extraída a matriz diagonal.
 \param mInv Ponteiro para o vetor que receberá o pré-condicionador.
*/
void calculaMInv(SL *sl, double *mInv);


/*!
 \brief Função que multiplica uma matriz k-diagonal com um vetor. o resultado é
        armazenado no vetor dest.

 \param sl Ponteiro para o sistema linear da matriz a ser multiplicada.
 \param v Ponteiro para o vetor a ser multiplicado.
 \param dest Ponteiro para o vetor que receberá o resultado da multiplicação.
*/
void multiplicaMatrizVetor(SL *sl, double *v, double *dest);



/*!
 \brief Função que transforma a matriz k-diagonal de um sistema linear em sua
        transposta. A operação é realizada trocando os valores das diagonais 
        superiores da matriz pelos valores das diagonais inferiores.

 \param sl Ponteiro para o sistema linear da matriz a ser multiplicada.
 \param v Ponteiro para o vetor a ser multiplicado.
 \param dest Ponteiro para o vetor que receberá o resultado da multiplicação.
*/
void obtemMatrizTransposta(SL *sl);


/*!
 \brief Função que multiplica um sistema linear por uma matriz k-diagonal.

 \param slA Ponteiro para o sistema linear de onde será obtida a matriz a ser 
            multiplicada.
 \param slB Ponteiro para o sistema linear que será multiplicado pela matriz.

 \return Ponteiro para o sistema linear resultante da multiplicação.
*/
SL *multiplicaMatrizSL(SL *slA, SL *slB);


/*!
 \brief Função que realiza o produto escalar de dois vetores.

 \param v1 Ponteiro para o primeiro vetor do produto.
 \param v2 Ponteiro para o segundo vetor do produto.
 \param n Tamanho dos vetores.

 \return Resultado do produto escalar.
*/
double multiplicaVetores(double *v1, double *v2, int n);


/*!
 \brief Função que calcula o residuo de um sistema linear. O residuo é calculado
        pela diferença entre o vetor b e o produto da matriz A pelo vetor x.

 \param sl Ponteiro para o sistema linear.
 \param x Ponteiro para o vetor solução.
 \param r Ponteiro para o vetor que receberá o residuo.
*/
void calculaResiduo(SL *sl, double *x, double *r);


/*!
 \brief Função que calcula a norma euclidiana, também conhecida como norma L2,
        de um vetor. A norma euclidiana é calculada pela raiz quadrada da soma
        dos quadrados dos elementos do vetor.

 \param v Ponteiro para o vetor.
 \param n Tamanho do vetor.

 \return Norma euclidiana do vetor.
*/
double calculaNormaEuclidiana(double *v, int n);


/*!
 \brief Função que calcula a norma máxima absoluta da diferença de dois vetores.
        A norma máxima absoluta é calculada pela maior diferença absoluta entre
        dois valores de mesma posição nos vetores.

 \param v1 Ponteiro para o primeiro vetor.
 \param v2 Ponteiro para o segundo vetor.
 \param n Tamanho dos vetores.

 \return Norma máxima da diferença dos vetores.
*/
double calculaNormaMax(double *v1, double *v2, int n);


/*!
 \brief Função que calcula a norma máxima relativa da diferença de dois vetores.
    A norma máxima relativa é calculada pela maior diferença relativa entre
    dois valores de mesma posição nos vetores.

 \param v1 Ponteiro para o primeiro vetor.
 \param v2 Ponteiro para o segundo vetor.
 \param n Tamanho dos vetores.

 \return Norma máxima relativa da diferença dos vetores.
*/
double calculaNormaMaxRelativa(double *v1, double *v2, int n);


/*!
 \brief Função que calcula um timestamp, para que a diferença entre dois
        timestamps seja o tempo de execução de um trecho de código.

 \return O timestamp.
*/
double timestamp(void);

#endif