#ifndef __METODOS_H__
#define __METODOS_H__

#include "sislin.h"

// Resolve sistema linear pelo método do gradiente conjugado com pré condicionador
void preCondicionado(SL *sl, double *x, double erro, int maxIt);

// Resolve sistema linear pelo método do gradiente conjugado
void gradienteConjugado(SL *sl, double *x, double erro, int maxIt);

#endif