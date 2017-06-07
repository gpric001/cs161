#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void mulijk(double*, double*, double*, int);
void mulikj(double*, double*, double*, int);
void mulkji(double*, double*, double*, int);
double getError(double*, double*, int);

const int N = 4000;

int main() {
    srand(time(0));
    double *a, *b, *c1, *c2, *c3;
    time_t start, end;
    double ijk_time = 0;
    double ikj_time = 0;
    double kji_time = 0;

    a = (double *) malloc(N*N*sizeof(double));
    b = (double *) malloc(N*N*sizeof(double));
    c1 = (double *) calloc(N*N, sizeof(double));
    c2 = (double *) calloc(N*N, sizeof(double));
    c3 = (double *) calloc(N*N, sizeof(double));

    int i;
    for (i = 0; i < N*N; ++i) {
        a[i] = rand() % 1000000 / 1000.0;
        b[i] = rand() % 1000000 / 1000.0;
    }
    
    start = clock();
    mulijk(a, b, c1, N);
    end = clock();
    ijk_time = (double)(end - start) / CLOCKS_PER_SEC;
    start = clock();
    mulikj(a, b, c2, N);
    end = clock();
    ikj_time = (double)(end - start) / CLOCKS_PER_SEC;
    start = clock();
    mulkji(a, b, c3, N);
    end = clock();
    kji_time = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("Time for ijk version: %.6lf\n", ijk_time);
    printf("Time for ikj version: %.6lf\n", ikj_time);
    printf("Time for kji version: %.6lf\n", kji_time);

    double error1 = getError(c1, c2, N);
    double error2 = getError(c1, c3, N);

    printf("Error for c2: %.6lf\n", error1);
    printf("Error for c3: %.6lf\n", error2);

    return 0;
}

double getError(double* model, double* source, int n) {
    double error = 0;
    int i;
    for (i = 0; i < N * N; ++i) {
        double derror = model[i] - source[i];
        error += derror * derror;
    }
    return error;
}

void mulijk(double* a, double* b, double* c, int n) {
    int i, j, k;
    for (i = 0; i < n; ++i) {
       for (j = 0; j < n; ++j) {
          for (k = 0; k < n; ++k) {
             c[i * n + j] = c[i * n + j] + a[i * n + k] * b[k * n + j];
          }
       }
    }
} 

void mulikj(double* a, double* b, double* c, int n) {
    int i, j, k;
    for (i = 0; i < n; ++i) {
       for (k = 0; k < n; ++k) {
          for (j = 0; j < n; ++j) {
             c[i * n + j] = c[i * n + j] + a[i * n + k] * b[k * n + j];
          }
       }
    }
} 

void mulkji(double* a, double* b, double* c, int n) {
    int i, j, k;
    for (k = 0; k < n; ++k) {
       for (j = 0; j < n; ++j) {
          for (i = 0; i < n; ++i) {
             c[i * n + j] = c[i * n + j] + a[i * n + k] * b[k * n + j];
          }
       }
    }
} 
