#include <stdio.h>

double* mm1(double* A, double* B, int n);
double* mm2(double* A, double* B, int n);
void printMatrix(double* M, int n);
double* fillMatrix(int n);

int main(void) {
	int n = 1500;
	double* A = fillMatrix(n);
	double* B = fillMatrix(n);
	//double* C = mm1(A, B, n);
	double* C = mm2(A, B, n);
	//printf("Matrix A\n");
	//printMatrix(A, n);
	//printf("Matrix B\n");
	//printMatrix(B, n);
	//printf("Result of A * B\n");
	//printMatrix(C, n);
	free(A);
	free(B);
	free(C);
	return 0;
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
	double* C = (double*)malloc(sizeof(double) * n * n);
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
	double* C = (double*)malloc(sizeof(double) * n * n);
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
