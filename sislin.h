typedef struct {
    double **diagonais;
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