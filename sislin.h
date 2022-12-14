#ifndef __SISLIN_H__
#define __SISLIN_H__

typedef struct {
    double **A;
    double *b;
    int n;
    int k;
} SL;

// Aloca sistema linear
SL *alocaSL(int n, int k);

void populaSL(SL *sl);

void copiaSL(SL *slDest, SL *slOrigin);

// Imprime sistema linear
void imprimeSL(SL *sl);

// Destroi sistema linear
void destroiSL(SL *sl);

#endif