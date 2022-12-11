#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sislin.h"

int main(int argc, char *argv[]) {
    int c;
    int n = -1;
    int k = -1;
    int p = -1;
    int i = -1;
    double e = 0.0;
    FILE *fp = NULL;

    while((c = getopt(argc, argv, "n:k:p:i:e:o:")) != -1) {
        switch(c) {
            case 'n':
                n = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'p':
                p = atoi(optarg);
                break;
            case 'i':
                i = atoi(optarg);
                break;
            case 'e':
                e = atof(optarg);
                break;
            case 'o':
                fp = fopen(optarg, "w");
                if (fp == NULL) {
                    fprintf(stderr, "Não foi possível abrir o arquivo.");
                    exit(1);
                }
                break;
            default:
                fprintf (stderr, "Uso: %s -n <tamanho do sistema linear> -k <numero de diagonais> -p <pre-condicionador> -i <máximo de iterações> -e <erro aproximado máximo absoluto> -o <arquivo_saida>\n", argv[0]);
                exit(1);
        }
    }

    if (n == -1) {
        fprintf(stderr, "A flag obrigatória -n não foi passada.\n");
        exit(1);
    }
    if (n <= 10) {
        fprintf(stderr, "O tamanho do sistema linear deve ser maior que 10\n");
        exit(1);
    }

    if (k == -1) {
        fprintf(stderr, "A flag obrigatória -k não foi passada.\n");
        exit(1);
    }
    if (k <= 1 || k % 2 == 0) {
        fprintf(stderr, "O número de diagonais deve ser ímpar e maior que 1\n");
        exit(1);
    }

    if (p == -1) {
        fprintf(stderr, "A flag obrigatória -p não foi passada.\n");
        exit(1);
    }

    if (i <= 0) {
        fprintf(stderr, "A flag obrigatória -i é inválida.\n");
        exit(1);
    }

    if (fp == NULL) {
        fprintf(stderr, "A flag obrigatória -o não foi passada.\n");
        exit(1);
    }

    srand(20222);

    SL *sl = criaSL(n, k);
    imprimeSL(sl);
    double *x = (double *) malloc(n * sizeof(double));
    gradienteConjugado(sl, x, e, i);
    free(x);
    destroiSL(sl);

    fclose(fp);
    return 0;
}
