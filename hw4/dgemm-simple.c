#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int N = 2048;

int main()
{
    int n = N;
    int i, j, k;
    double *a, *b, *c;
    time_t start, end;

    /* malloc */
    a = (double *) malloc(N*N*sizeof(double));
    b = (double *) malloc(N*N*sizeof(double));
    c = (double *) malloc(N*N*sizeof(double));
    if(!a || !b || !c)
    {
	printf("allocate failed!\n");
	exit(-1);
    }

    /* initialization */
    for(i = 0; i < N*N; i++)
    {
	a[i] = rand()%1000000/1000.0;
	b[i] = rand()%1000000/1000.0;
	c[i] = rand()%1000000/1000.0;
    }

    /* ijk version */
    start = clock();
    for(i = 0; i < n; i++)
	 for(j = 0; j < n; j++)
	    for(k = 0; k < n; k++)
		 c[i*n+j] += a[i*n+k] * b[k*n+j];
    end = clock();
    printf("Time for ijk version: %.6lf\n", (double)(end-start)/CLOCKS_PER_SEC);

}
