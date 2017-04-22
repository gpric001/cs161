// Created by Guthrie Price
// 4/22/2017
// For CS 161 at UCR Spring 2017

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double* mm1(double* A, double* B, int n);
double* mm2(double* A, double* B, int n);
double* mm3(double* A, double* B, int n);
void randMatrix(double* M, int n);
double getTimeUsed(clock_t start, clock_t end);
void freeMatricies(double* matricies[], int rows, int cols);
void outputResults(double cpuTimes[], double* matricies[], int rows, int cols);
void outputMatrix(double* M, int n);
void reportStats(double* model, double* sample, int n);

const int NUM_OPS = 3; // Number of different types of matrix multiplication
const int NUM_TESTS = 4; // Number of different test sizes
const int MIN_SIZE = 1000; // The minimum dimension for the matrix
const int STEP_SIZE = 500; // The amount to increase the dimension per test

int main(void) {
    srand(time(NULL));

    // Array of pointers to the different matrix operations
    double* (*mmOps[NUM_OPS]) (double* A, double* B, int n);
    mmOps[0] = mm1;
    mmOps[1] = mm2;
    mmOps[2] = mm3;

    clock_t start, end;
    int n, currentOp;
    double cpuTimes[NUM_TESTS * NUM_OPS];
    double* cMatricies[NUM_TESTS * NUM_OPS];

    // For each test case compute the multiplication of random matricies
    // using different matrix muliplication procedures. Record the resulting
    // matrix as well as the time in took for each multiplication.
    for(n = 0; n < NUM_TESTS; ++n) {
        int matrixSize = n * STEP_SIZE + MIN_SIZE;
        double* A = malloc(sizeof(double) * matrixSize * matrixSize);
        double* B = malloc(sizeof(double) * matrixSize * matrixSize);
        randMatrix(A, matrixSize);
        randMatrix(B, matrixSize);
        for(currentOp = 0; currentOp < NUM_OPS; ++currentOp) {
            start = clock();
            double* C = (*mmOps[currentOp]) (A, B, matrixSize);
            end = clock();
            cpuTimes[NUM_OPS * n + currentOp] = getTimeUsed(start, end);
            cMatricies[NUM_OPS * n + currentOp] = C;
        }
        free(A);
        free(B);
    }

    // Output the results and free memory
    outputResults(cpuTimes, cMatricies, NUM_TESTS, NUM_OPS);
    freeMatricies(cMatricies, NUM_TESTS, NUM_OPS);
	return 0;
}

//=============================================================================
// Reports the maximum and mean squared error between the model and sample
// matrix.
// Both matricies are assumed to be square and of the same size.
//
// INPUTS
// model:  The matrix you are comparing against.
// sample: The matrix you are computing the error of.
// n:      The number of columns or number of rows.
//
// OUTPUTS
// None
//
// SIDE EFFECTS
// Prints the maximum error and mean squared error to standard output.
void reportStats(double* model, double* sample, int n) {
    int i, totalElements = n * n;
    double maxError = 0, meanSqError = 0;
    for (i = 0; i < totalElements; ++i) {
        double error = fabs(model[i] - sample[i]);
        double sqError = error * error;
        if (error > maxError) maxError = error;
        meanSqError += sqError;
    }
    printf("Max error: %f\n", maxError);
    printf("Mean Squared error: %f\n", (meanSqError / totalElements));
}

//=============================================================================
// Reports the time taken for each of the matrix multipliations.
// Also reports the errors between the last and second to last matricies
// genrated by the matrix operations (see reportStats).
// cpuTImes and matricies are assumed to be of the same size.
//
// INPUTS
// cpuTimes:  An array of CPU times.
// matricies: An array of square matricies.
// rows:      The number of rows in cpuTimes and matricies.
// cols:      The number of cols in cpuTimes and matricies.
//
// OUTPUTS
// None
//
// SIDE EFFECTS
// Outputs CPU Times and errors to standard output.
void outputResults(double cpuTimes[], double* matricies[], int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; ++i) {
        int dim = i * STEP_SIZE + MIN_SIZE;
        printf("Test results for matrix of size %d by %d\n", dim, dim);
        for (j = 0; j < cols; ++j) {
            printf("Time of mm%d: %f [seconds]\n", j + 1, cpuTimes[cols * i + j]);
            if (j == cols - 1) {
                reportStats(matricies[cols * i + j - 1], 
                            matricies[cols * i + j], 
                            dim);
            }
        }
    }
}

//=============================================================================
// Calculates the time used by the CPU in a time interval.
//
// INPUTS
// start: The beginning of the time interval.
// end:   The end of the time interval.
//
// OUTPUTS
// The CPU time used between start and end.
//
// SIDE EFFECTS
// None
double getTimeUsed(clock_t start, clock_t end) {
    return ((double) (end - start)) / CLOCKS_PER_SEC;
}

//=============================================================================
// Fills a square matrix with random values between 0 and 1.
//
// INPUTS
// M: The square matrix to be filled.
// n: The number of rows or number of columns of M.
//
// OUTPUTS
// None
//
// SIDE EFFECTS
// M is filled with random values between 0 and 1. Whatever data was in M is
// lost.
void randMatrix(double* M, int n) {
	int row, col;
	for (row = 0; row < n; ++row) {
		for (col = 0; col < n; ++col) {
			M[row * n + col] = ((double) rand()) / RAND_MAX;
		}
	}
}

//=============================================================================
// Free all memory pointed to be the values in an array of pointers.
//
// INPUTS
// matricies: An array of pointers that point to memeory allocated on the heap.
// rows:      The number of rows of matricies.
// cols:      The number of cols of matricies.
//
// OUTPUTS
// None
//
// SIDE EFFECTS
// All the data pointed to be the pointers in matricies is lost
void freeMatricies(double* matricies[], int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            free(matricies[cols * i + j]);
        }
    }   
}

//=============================================================================
// A naive implementation of matrix multiplication of two square matricies.
//
// INPUTS
// A: The first matrix to be multiplied.
// B: The second matrix to be multiplied.
// n: The number or rows or number of columns of A and B.
//
// OUTPUTS
// C: The result of the matrix multiplication A * B.
//
// SIDE EFFECTS
// Memory on the heap has been allocated for C. Remember to free it.
double* mm1(double* A, double* B, int n) {
	double* C = (double*) calloc(n * n, sizeof(double));
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

//============================================================================
// A slightly optimized version of matrix multiplication of two square
// matricies. Only loads and stores the memory location of an element of C
// once. Calculates the value of an element of C in a register.
//
// INPUTS
// A: The first matrix to be multiplied.
// B: The second matrix to be multiplied.
// n: The number or rows or number of columns of A and B.
//
// OUTPUTS
// C: The result of the matrix multiplication A * B.
//
// SIDE EFFECTS
// Memory on the heap has been allocated for C. Remember to free it.=
double* mm2(double* A, double* B, int n) {
	double* C = (double*) calloc(n * n, sizeof(double));
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

//=============================================================================
// A more optimized version of matrix multiplication of two square matricies.
// Calculates the matrix C by breaking the multiplication into many 2 x 2 
// matrix multiplications. where each 2 x 2 matrix multiplication is done in
// registers.
//
// INPUTS
// A: The first matrix to be multiplied.
// B: The second matrix to be multiplied.
// n: The number or rows or number of columns of A and B.
//
// OUTPUTS
// C: The result of the matrix multiplication A * B.
//
// SIDE EFFECTS
// Memory on the heap has been allocated for C. Remember to free it.double* mm3(double* A, double* B, int n) {
double* mm3(double* A, double* B, int n) {
    double* C = (double*) calloc(n * n, sizeof(double));
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

//=============================================================================
// Outputs the elements of a matrix to the standard output.
//
// INPUTS
// M: The square matrix to be output.
// n: The number of rows or number of columns of M.
//
// OUTPUTS
// None
//
// SIDE EFFECTS
// The elements of M will be output to the standard output.
void outputMatrix(double* M, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%f ", M[n * i + j]);
        }
        printf("\n");
    }
    printf("\n");
}
