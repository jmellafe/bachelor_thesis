//Ejecucion del algoritmo de Gillespie para procesos estocasticos a tiempo continuo sobre
//2D Lattice cuadrada. El proceso estocástico utilizado será e voter model, la distribución de
// probabilidad temporal del cambio será de poisson, todos los procesos tendran el mismo
// coeficiente (aunque podriamos ponerlos diferentes sin ningun problema, el programa lo
//contemplaria)
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>


int L = 16;

static double poisConst = 1.;

double poisMean2 = 0;


static size_t mcPas = 10000;

void fillPoisLat(double **poisLat) {

    int i, j;

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            poisLat[i][j] = poisConst;
        }
    }

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            poisMean2 += 1 / (poisLat[i][j] * poisLat[i][j]);
        }
    }

    poisMean2 = 1 / poisMean2;


}

double meanMtrx(double **mtrx, int lenRow, int lenCol) {

    int i, j;
    double mean = 0.;

    for (i = 0; i < lenRow; i++) {
        for (j = 0; j < lenCol; j++) {
            mean += mtrx[i][j];
        }
    }

    mean /= (double) (lenRow * lenCol);

    return mean;
}

double meanMtrxInt(int **mtrx, int lenRow, int lenCol) {

    int i, j;
    double mean = 0.;

    for (i = 0; i < lenRow; i++) {
        for (j = 0; j < lenCol; j++) {
            mean += mtrx[i][j];
        }
    }

    mean /= (double) (lenRow * lenCol);

    return mean;
}


double rand_expo(double param) {

    double u = (double) rand() / RAND_MAX;


    return -log(1 - u) / param;

}


double rand_time(double param) {

    double u = (double) rand() / RAND_MAX;


    return sqrt(-2. * poisMean2 * log(1 - u));

}

size_t rand_process(double **poisLat, int numRow, int numCol) {

    if (1) {

        return (size_t) numRow * numCol * rand() / RAND_MAX;
    } else {
        return 0;
    }
}

int transPos(int pos) {

    if (pos == L) {
        return 0;
    }
    if (pos == -1) {
        return L - 1;
    }
    return pos;


}


int main() {

//Definimos la lattice de estado que se irá actualizando y la lattice con los coeficientes
//poisson


    size_t i, j, k, m;

    double **poisLat;
    int **stateLat;

    stateLat = (int **) malloc(L * sizeof(int *));
    poisLat = (double **) malloc(L * sizeof(double *));

    srand(time(NULL));

    printf("%f\n", log(2.7));


    for (i = 0; i < L; i++) {
        poisLat[i] = malloc(L * sizeof(double));
        stateLat[i] = malloc(L * sizeof(int));

        for (j = 0; j < L; j++) {
            int r = rand();

            if ((double) r / RAND_MAX < 0.5) {
                stateLat[i][j] = 1;
                printf("1 ");
            } else {
                stateLat[i][j] = 0;

                printf("0 ");
            }
        }
        printf("\n");
    }


    fillPoisLat(poisLat);
//
//    double poisMed = meanMtrx(poisLat, L, L);
//
//
//    poisMed = L*L*poisMed;

    double tiempo = 0., **allData;
    size_t proc;


    allData = (double **) malloc(mcPas * sizeof(double *));


    int vecino, row, col;
    size_t maxMc = mcPas;


    for (i = 0; i < mcPas; i++) {

        for (j = 0; j < L * L; j++) {

//            Calculamos el tiempo hasta el proximo cambio y sumamos al tiempo global,
//            y vemos cual es el proceso cambiado

            tiempo += rand_time(poisMean2);

            proc = rand_process(poisLat, L, L);

            row = (int) proc / L;

            col = (int) proc % L;

//            El nodo seleccionado copia a un vecino al azar

            vecino = (int) 4 * (double) rand() / RAND_MAX;

            if (vecino == 0) {
                stateLat[col][row] = stateLat[col][transPos(row - 1)];
            } else if (vecino == 1) {
                stateLat[col][row] = stateLat[col][transPos(row + 1)];
            } else if (vecino == 2) {
                stateLat[col][row] = stateLat[transPos(col + 1)][row];
            } else {
                stateLat[col][row] = stateLat[transPos(col - 1)][row];
            }

        }

        allData[i] = malloc(2 * sizeof(double));

        allData[i][0] = tiempo;
        allData[i][1] = meanMtrxInt(stateLat, L, L);

        if (fabs(allData[i][1] - 1.) < 0.1 / (double) (L * L) || fabs(allData[i][1]) < 0.1 / (double) (L * L)) {
            maxMc = i;
            printf("Se ha alcanzado el equilibrio (todos opinan lo mismo) a tiempo \nT=%f \n\n", tiempo);
            break;
        }
    }

    printf("Estos numeros deberían ser similares\n");
    printf("%f\n", maxMc * L * L * sqrt(M_PI * poisMean2 / 2.));
    printf("%f\n", tiempo);
    printf("%d\n", (int) maxMc);


    FILE *fout = fopen("/home/alex/CLionProjects/tfg/results/voter_model_gillespie.dat", "w");
    fprintf(fout, "# t, Media Opinión\n");
    for (i = 0; i < maxMc; i++)
        fprintf(fout, "%f %f\n", allData[i][0], allData[i][1]);

    fclose(fout);


    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            if (stateLat[i][j] == 1) {
                printf("1 ");
            } else {
                printf("0 ");
            }

        }
        printf("\n");


    }


    for (i = 0; i < L; i++) {
        free(poisLat[i]);
        free(stateLat[i]);
    }
    for (i = 0; i < maxMc; i++) {
        free(allData[i]);
    }

    free(poisLat);
    free(stateLat);
    free(allData);


    poisLat = NULL;
    stateLat = NULL;
    allData = NULL;


    return 0;
}