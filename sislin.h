// Autores: André Schueda Menezes e Marcus Augusto Ferreira Dudeque

#ifndef __SISLIN_H__
#define __SISLIN_H__


/*!
    \brief Estrutura que representa um sistema linear k-diagonal.

    Para representar o sistema linear, foi utilizada essa estrutura de dados,
    de forma que a matriz k-diagonal é armazenada em um vetor unidimensional
    mas que representa uma matriz de banda, bidimensional. Apesar de o acesso
    a memória ser menos intuitivo, pois é utilizada a aritmética de ponteiros,
    essa estrutura de dados é mais eficiente, pois todos os elementos da matriz
    estão contiguos na memória, o que permite que o acesso a memória seja mais
    rápido, devido a localidade espacial da memória cache.

    \param A Ponteiro para o vetor A, que é a matriz k-diagonal.
    \param b Ponteiro para o vetor b, que é o vetor de termos independentes.
    \param n Tamanho do sistema linear.
    \param k Numero de diagonais da matriz A.
*/
typedef struct {
    double *A;
    double *b;
    int n;
    int k;
} SL;


/*!
 \brief Função que aloca o sistema linear k-diagonal.

 \param n Tamanho do sistema linear.
 \param k Numero de diagonais da matriz A.

 \return Ponteiro para o SL criado.
*/
SL *alocaSL(int n, int k);


/*!
 \brief Função que popula o sistema linear k-diagonal com valores aleatórios que
        satisfazem a condição de diagonal dominante.

 \param sl Ponteiro para o sistema linear a ser populado.
*/
void populaSL(SL *sl);


/*!
 \brief Função que copia os dados de um sistema linear para outro.

 \param slDest Ponteiro para o sistema linear destino.
 \param slOrigin Ponteiro para o sistema linear origem.
*/
void copiaSL(SL *slDest, SL *slOrigin);


/*!
 \brief Função que imprime o sistema linear k-diagonal.

 \param sl Ponteiro para o SL a ser impresso.
*/
void imprimeSL(SL *sl, FILE *fp);


/*!
 \brief Função que desaloca o sistema linear k-diagonal.

 \param sl Ponteiro para o SL a ser desalocado.
*/
void destroiSL(SL *sl);

#endif