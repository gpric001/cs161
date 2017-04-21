#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double* mm1(double* A, double* B, int n);
double* mm2(double* A, double* B, int n);
double* mm3(double* A, double* B, int n);
void printMatrix(double* M, int n);
double* fillMatrix(int n);
double getTimeUsed(clock_t start, clock_t end);

const int NUM_OPS = 3;

int main(void) {
    double* (*mmOps[NUM_OPS]) (double* A, double* B, int n);
    mmOps[0] = mm1;
    mmOps[1] = mm2;
    mmOps[2] = mm3;
    clock_t start, end;
    double cpuTimes[NUM_OPS];
    int n, currentOp;
	double* A = fillMatrix(n);
	double* B = fillMatrix(n);
    for(n = 1000; n <= 2500; n += 500) {
        for(currentOp = 0; currentOp < NUM_OPS; ++currentOp) {
            start = clock();
            double* C = (*mmOps[currentOp]) (A, B, n);
            end = clock();
            cpuTimes[currentOp] += getTimeUsed(start, end);
            free(C);
        }
        printf("For size %d by %d matrix\n", n, n);
    }
    printf("Unoptimized matrix multiplication: %f [seconds]\n", cpuTimes[0]);
    printf("Optimized matrix multiplication v1: %f [seconds]\n", cpuTimes[1]);
    printf("Optimized matrix multiplication v2: %f [seconds]\n", cpuTimes[2]);
	free(A);
	free(B);
	return 0;
}

double getTimeUsed(clock_t start, clock_t end) {
    return ((double)(end - start)) / CLOCKS_PER_SEC;
}

double* fillMatrix(int n) {
	double* M = (double*)malloc(sizeof(double) * n * n);
	int row, col;
	for (row = 0; row < n; ++row) {
		for (col = 0; col < n; ++col) {
			M[row * n + col] = row + col + 0.1;
		}
	}
	return M;
}

double* mm1(double* A, double* B, int n) {
	double* C = (double*)calloc(n * n, sizeof(double));
	int i, j, k;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			for (k = 0; k < n; ++k) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}
		}
	}
	return C;
}

double* mm2(double* A, double* B, int n) {
	double* C = (double*)calloc(n * n, sizeof(double));
	int i, j, k;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			register double t = C[i * n + j];
			for (k = 0; k < n; ++k) {
				t += A[i * n + k] * B[k * n + j];
			}
			C[i * n + j] = t;
		}
	}
	return C;
}

double* mm3(double *A, double*B, int n) {
    double* C = (double*)calloc(n * n, sizeof(double));
    int i, j, k;
    for (i = 0; i < n; i += 2) {
        for (j = 0; j < n; j += 2) {
            register double c11 = C[i * n + j];
            register double c12 = C[i * n + j + 1];
            register double c21 = C[(i + 1) * n + j];
            register double c22 = C[(i + 1) * n + j + 1];
            for (k = 0; k < n; k += 2) {
                register double a11 = A[i * n + k];
                register double a12 = A[i * n + k + 1];
                register double a21 = A[(i + 1) * n + k];
                register double a22 = A[(i + 1) * n + k + 1];
                register double b11 = B[k * n + j];
                register double b12 = B[k * n + j + 1];
                register double b21 = B[(k + 1) * n + j];
                register double b22 = B[(k + 1) * n + j + 1];
                c11 += (a11 * b11) + (a12 * b21);
                c12 += (a11 * b12) + (a12 * b22);
                c21 += (a21 * b11) + (a22 * b21);
                c22 += (a21 * b12) + (a22 * b22);
            }
            C[i * n + j] = c11;
            C[i * n + j + 1] = c12;
            C[(i + 1) * n + j] = c21;
            C[(i + 1) * n + j + 1] = c22;
        }
    }
    return C;
}

void printMatrix(double* M, int n) {
	int row, col;
	for (row = 0; row < n; ++row) {
		for (col = 0; col < n; ++col) {
			printf("%f ", M[row * n + col]);
		}
		printf("\n");
	}
	printf("\n");
}
