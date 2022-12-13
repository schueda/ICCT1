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
void imprimeVetor(double *v, int n);

// Resolve sistema linear pelo método do gradiente conjugado com pré condicionador
void preCondicionado(SL *sl, double *x, double erro, int maxIt);

// Resolve sistema linear pelo método do gradiente conjugado
void gradienteConjugado(SL *sl, double *x, double erro, int maxIt);

// Destroi sistema linear
void destroiSL(SL *sl);