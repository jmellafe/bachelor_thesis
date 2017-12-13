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


int L = 10;

static double poisConst = 1.;

double poisMean2 = 0;

static size_t mcPas = 10000000, numProm = 100, maxL = 22;

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

void fillStateLat(int **stateLat) {

    int i, j;

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            if ((double) rand() / RAND_MAX < 0.5) {
                stateLat[i][j] = 1;
            } else {
                stateLat[i][j] = 0;
            }
        }
    }


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


long double rand_expo(double param) {

    long double u = (long double) rand() / RAND_MAX;


    return -logl(1. - u) / param;

}


double rand_time(double param) {

    double u = (double) rand() / RAND_MAX;


    return sqrt(-2. * poisMean2 * log(1 - u));

}

size_t rand_process(double **poisLat, int numRow, int numCol) {

    if (1) {

        return (size_t) rand()%(numRow*numCol);
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
    unsigned int seed = (unsigned int) time(NULL);

    srand(seed++);

    size_t i, j, k, m;
    double mediaEstado;
//
//    double poisMed = meanMtrx(poisLat, L, L);
//
//
//    poisMed = L*L*poisMed;

    long double tiempo = 0.;
    size_t proc;


    int vecino, row, col;

    long double **allData;
    size_t numIters = (size_t) (maxL-L)/10+1;

    allData = (long double **) malloc(numIters * sizeof(long double *));

    double **poisLat;
    int **stateLat;

    for (m = 0; L < maxL; L += 10, m++) {
        printf("L=%d \n", L);

        double itersMax = 0.;

        allData[m] = malloc(3 * sizeof(long double));


//        Inicializamos las matrices


        stateLat = (int **) malloc(L * sizeof(int *));
        poisLat = (double **) malloc(L * sizeof(double *));


        for (i = 0; i < L; i++) {
            poisLat[i] = malloc(L * sizeof(double));
            stateLat[i] = malloc(L * sizeof(int));

        }


        fillPoisLat(poisLat);

        size_t maxMc = mcPas;
        long double consTime = 0., consTimeDesv = 0.;
        for (k = 0; k < numProm; k++) {
            printf("%d \n", (int) k);


            srand(seed++);
            fillStateLat(stateLat);

            size_t sumStat = 0;

//            fill weights

            for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {

                    sumStat += stateLat[i][j];
                }
            }

            maxMc = mcPas;
            tiempo = 0.;


            for (i = 0; i < mcPas; i++) {

                for (j = 0; j < L * L; j++) {

//            Calculamos el tiempo hasta el proximo cambio y sumamos al tiempo global,
//            y vemos cual es el proceso cambiado

                    tiempo += rand_time(poisMean2);

                    proc = rand_process(poisLat, L, L);

                    row = (int) proc / L;

                    col = (int) proc % L;

//            El nodo seleccionado copia a un vecino al azar

                    vecino = rand() % 4;

                    sumStat -= stateLat[row][col];

                    if (vecino == 0) {
                        stateLat[row][col] = stateLat[row][transPos(col - 1)];
                    } else if (vecino == 1) {
                        stateLat[row][col] = stateLat[row][transPos(col + 1)];
                    } else if (vecino == 2) {
                        stateLat[row][col] = stateLat[transPos(row + 1)][col];
                    } else {
                        stateLat[row][col] = stateLat[transPos(row - 1)][col];
                    }
                    sumStat += stateLat[row][col];


                    if (sumStat == L*L || sumStat == 0) {
                        maxMc = i;
                        consTime += tiempo / (double) numProm;
                        printf("tiempo prom %d \n", (int) maxMc);
                        itersMax += (double) maxMc / (double) numProm;
                        consTimeDesv += tiempo * tiempo / (double) numProm;
                        goto endMC;
                    }
                }
            }

            printf("No se ha llegado al consenso a L=%d con %d pasos monte carlo\n", L, (int) mcPas);

            endMC:;

        }

        printf("iters medias consenso %f \n", itersMax);

        allData[m][0] = L * L;
        allData[m][1] = consTime;
        allData[m][2] = consTimeDesv;

//        Liberamos
        for (i = 0; i < L; i++) {
            free(poisLat[i]);
            free(stateLat[i]);
        }
        poisLat = NULL;
        stateLat = NULL;
    }


    FILE *fout = fopen("/home/alex/CLionProjects/tfg/results/tiempo_consenso.dat", "w");
    fprintf(fout, "#L, tiempo consenso, desv\n");
    for (i = 0; i < numIters; i++)
        fprintf(fout, "%Le %Le %Le\n", allData[i][0], allData[i][1],
                sqrtl(allData[i][2] - allData[i][1] * allData[i][1]));

    fclose(fout);


    for (i = 0; i < numIters; i++) {
        free(allData[i]);
    }
    free(allData);
    allData = NULL;

    return 0;
}