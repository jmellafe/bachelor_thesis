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


int L = 20;

static double poisConst = 1.;

double poisMean2 = 0;

static size_t mcPas = 10000000, numProm = 10000, maxL = 22;

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


double **histogram(double *data, size_t numBars, size_t len){
    double **histo, max, min, h;

    max = data[0];
    min = data[0];
    for(int i=1;i<len;i++){
        if(min > data[i]){
            min = data[i];
        }
        else if(max<data[i]){
            max = data[i];
        }
    }

    h = (max-min)/(double) numBars;

    histo = (double **)malloc(3*sizeof(double*));

    for(int i=0;i<2;i++){
        histo[i] = malloc(numBars*sizeof(double));
    }
    for(int i=0;i<numBars;i++){
        histo[0][i] = min+(i+1./2.)*h;
        histo[1][i] = 0.;
    }

    for(int i=0;i<len;i++){
        histo[1][(int)((data[i]-min)/h)]++;
    }

    for(int i=0;i<numBars;i++){
        histo[1][i] /= (double)(h*len);
    }

    double norma = 0.;
    for(int i=0;i<numBars;i++){
        norma += h*histo[1][i];
    }

    return histo;




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
    double *allTimes, coefCons = 0., coefConsLog = 0.;
    allTimes = malloc(numProm*sizeof(double));

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
                        allTimes[k] = (double)tiempo;
                        coefConsLog += log((double)tiempo);
                        coefCons += tiempo;
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


    size_t numBars = 100;
    double **histo;

    histo = histogram(allTimes, numBars, numProm);

    double s = log(coefCons/(double)numProm) - coefConsLog/(double)numProm;

    double sigma = (3-s+sqrt(pow((s-3),2.) + 24.*s))/(12.*s);

    printf("Coeficiente sigma %f \n", sigma );
    printf("Coeficiente otro %f \n", coefCons/(double)(numProm*sigma));


    FILE *fout = fopen("/home/alex/CLionProjects/tfg/results/histogram.dat", "w");
    fprintf(fout, "#L, tiempo consenso, desv\n");
    for (i = 0; i < numBars; i++)
        fprintf(fout, "%f %f\n", histo[0][i], histo[1][i]);

    fclose(fout);


    for (i = 0; i < numIters; i++) {
        free(allData[i]);
    }
    for(i=0;i<2;i++){
        free(histo[i]);
    }
    free(allData);
    free(histo);
    allData = NULL;
    histo = NULL;

    return 0;
}