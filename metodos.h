#ifndef __METODOS_H__
#define __METODOS_H__

#include "sislin.h"

/*!
 \brief Função que resolve o sistema linear pelo método do gradiente conjugado
        sem pré-condicionador.

 \param sl Ponteiro para o sistema linear a ser resolvido.
 \param x Ponteiro para o vetor que receberá a solução do sistema linear.
 \param erro Erro máximo permitido. Essa valor será considerado como critério de 
             parada quando comparado com o valor da norma máxima do vetor x e sua
             próxima iteração, de forma que, quando a diferença entre eles for
             suficientemente pequena, o resultado convergiu.
 \param maxIt Número máximo de iterações. Caso o número de iterações seja
              excedido, o programa será encerrado.
 \param fp Ponteiro para o arquivo no qual o resultado será impresso.
*/
void gradienteConjugado(SL *sl, double *x, double erro, int maxIt, FILE *fp);

/*!
 \brief Função que resolve o sistema linear pelo método do gradiente conjugado
        com o pré-condicionador de Jacobi, que consiste na inversa da diagonal
        principal da matriz A do sistema a ser solucionado.

 \param sl Ponteiro para o sistema linear a ser resolvido.
 \param x Ponteiro para o vetor que receberá a solução do sistema linear.
 \param erro Erro máximo permitido. Essa valor será considerado como critério de 
             parada quando comparado com o valor da norma máxima do vetor x e sua
             próxima iteração, de forma que, quando a diferença entre eles for
             suficientemente pequena, o resultado convergiu.
 \param maxIt Número máximo de iterações. Caso o número de iterações seja
              excedido, o programa será encerrado.
 \param fp Ponteiro para o arquivo no qual o resultado será impresso.
*/
void preCondicionado(SL *sl, double *x, double erro, int maxIt, FILE *fp);

#endif