#ifndef __SISLIN_H__
#define __SISLIN_H__

typedef struct {
    double **linhas;
    double *b;
    int n;
    int k;
} SL;

// Cria sistema linear
SL *criaSL(int n, int k);

// Imprime sistema linear
void imprimeSL(SL *sl);

// Destroi sistema linear
void destroiSL(SL *sl);

#endif