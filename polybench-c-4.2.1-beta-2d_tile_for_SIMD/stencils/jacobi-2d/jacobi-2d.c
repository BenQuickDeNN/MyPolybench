/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* jacobi-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "jacobi-2d.h"

#define SIZE_L1_DCACHE 32768 // L1 data cache size (byte)
#define SIZE_L2_CACHE 262144 // L2 cache size (byte)
#define NUM_ARRAY 2 // the number of arrays referenced

int TI, TI1; // tile size for i dimension
int TJ, TJ1; // tile size for j dimension

/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
	B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_jacobi_2d(int tsteps,
			    int n,
			    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
			    DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int t, ii, i1, i, jj, j1, j;
  int I1; // index for rest part
	int J1;
	const int end_i = _PB_N - 2;
	const int end_ii = end_i - ((_PB_N - 1) % TI) + 1;	
	const int end_j = end_i;
	const int end_jj = end_jj;

#pragma scop
	if (TI > 1) // tiling is valuable only if tile size is larger than 1
	{
		for (t = 0; t < _PB_TSTEPS; t++)
		{
			for (ii = 1; ii < end_ii; ii+=TI)
			{
				for (jj = 1; jj < end_jj; jj+=TJ)
				{
					if(ii==1)
					{
						TI1 = TI;
						i1 = 1;
					}
					else
					{
						TI1 = TI - 1;
						i1 = ii + 1;
					}
					if(jj==1)
					{
						TJ1 = TJ;
						j1 = 1;
					}
					else
					{
						TJ1 = TJ - 1;
						j1 = jj + 1;
					}
					/* Compute tiles of array B */
					for (i = i1; i < i1 + TI1; i++)
					{
						/* Compute tiles of array B */
						for (j = j1; j < j1 + TJ1; j++)
							B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
						/* Precompute down boundary points of array B */
						j = j1 + TJ1;
						if (j < end_j)
							B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
					}
					/* Precompute right boundary points of array B */
					i = i1 + TI1;
					if (i < end_i)
					{
						for (j = 1; j < _PB_N - 1; j++)
							B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
					}
					for (i = i1; i < i1 + TI1; i++)
					{
						/* Compute tiles of array A */
						for (j = j1; j < j1 + TJ1; j++)
							A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
						/* Precompute down boundary points of array A */
						j = j1 + TJ1;
						if (j < end_j)
							A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
					}
					/* Precompute right boundary points of array A */
					i = i1 + TI1;
					if (i < end_i)
					{
						for (j = 1; j < _PB_N - 1; j++)
							B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
					}
				}		
			}
		}
	}
  for (t = 0; t < _PB_TSTEPS; t++)
    {
      for (i = 1; i < _PB_N - 1; i++)
	for (j = 1; j < _PB_N - 1; j++)
	  B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
      for (i = 1; i < _PB_N - 1; i++)
	for (j = 1; j < _PB_N - 1; j++)
	  A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
    }
#pragma endscop

}


int main(int argc, char** argv)
{
	TI = floor(sqrt((SIZE_L1_DCACHE + SIZE_L2_CACHE) / sizeof(DATA_TYPE)));
	TI = TI>((N-3)/2)?floor((N-3)/2):TI; // TI should not be larger than half compute bound
	TJ = TI; // square tiling
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;
	printf("sizeof data_type = %d\r\n", sizeof(DATA_TYPE));
	printf("TI=%d\r\n", TI); // print tile size

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;
	/* Set up timer */
  clock_t start, finish;
  long duration;
  start = clock();
  /* Run kernel. */
  kernel_jacobi_2d(tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));
	/* Stop timer */
  finish = clock();
  duration = (long)(finish - start);
  /* Stop and print timer. */
  polybench_stop_instruments;
  //polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
	printf("duration = %ld us\r\n", duration);
  return 0;
}
