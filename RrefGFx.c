
/*
% function [G,rank] = RrefGFx(A,q)
% Returns reduced row echelon form of A and rank of A;

% matrix reduction is done in GF(q) algebra.
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#define LOG_SESSION 0
#define INF         1000

void reduceMatrix();
#if LOG_SESSION
void printMatrix(int rows, int columns);
FILE *logFile;
#endif
FILE *gfxFile;

#ifdef MATLAB_MEX_FILE
double *pAin;			// double precision input matrix from Matlab workspace
double *pAout;			// double precision reduced output matrix to Matlab workspace
double *pRank;			// double precision rank
int **A;			// input matrix, represented as a 2D array of q-ary elements
int **pA;
#else
unsigned int A[10] = { 0xFF, 0x7F, 0x3F, 0x1F, 0xF, 0x7, 0x3, 0x1 };
#endif

int m, n;			// number of rows, columns of A
int *jb;
int rank;
int *pA1;			// input matrix pointer
int *pA2;			// input matrix pointer
unsigned char vectors[256];	// Galois field LUTs
unsigned char symbols[256];
int q;
int prevQ = 0;

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int i = 0, j = 0, k;

	int element;
	char s[20];

	if (nrhs != 2) {
		mexErrMsgTxt("Not enough input arguments!");
		return;
	}

	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	q = (int)mxGetScalar(prhs[1]);

	A = (int **)calloc(m, sizeof(int *));
	if (A == NULL) {
		printf("Alloc failed (row pointers to A)\n");
		return;
	}
	pA = A;
	for (i = 0; i < m; i++) {
		pA1 = (int *)calloc(n, sizeof(int));
		if (pA1 == NULL) {
			printf("Alloc failed (A)\n");
			return;
		}
		*pA++ = pA1;
	}

// Read input matrix from Matlab workspace and convert to binary arrays.
// Matrices are stored in an m*n array as follows:
// A[0]                = [row 1, columns 1-32]
// A[1]                = [row 1, columns 33-64]
// A[2:(n-1)]           =[     etc     ]
// A[n]                = [row 2, columns 1-32]
// A[n+1]              = [row 2, columns 33-64]
// A[n:(2*n-1)] = [        etc     ]
// A[2*n:m*n]           =[     etc     ]

	for (i = 0; i < m; i++) {
		pAin = mxGetPr(prhs[0]) + i;
		for (j = 0; j < n; j++) {
			// Read elements from matrix and convert from symbols to powers
			element = (int)floor(*pAin - 0.5);
			if (element < 0) {

				element = -INF;	// approximate -infinity
			}
			pAin += m;
			A[i][j] = element;
		}
	}

	if (q != prevQ) {
		sprintf(s, "gf%d_8x.vec", q);
		gfxFile = fopen(s, "rb");
		if (gfxFile == NULL) {
			printf("File open failed (%s)\n", s);
			return;
		}
		fread(vectors, 1, q, gfxFile);
		fclose(gfxFile);
		sprintf(s, "gf%d_8x.sym", q);
		gfxFile = fopen(s, "rb");
		if (gfxFile == NULL) {
			printf("File open failed (%s)\n", s);
			return;
		} else {
			printf("Read new tables for q=%d\n", q);
			prevQ = q;
		}
		fread(symbols, 1, q, gfxFile);
		fclose(gfxFile);
	}

	reduceMatrix();

	plhs[0] = mxCreateDoubleMatrix(rank, n, mxREAL);
	for (i = 0; i < rank; i++) {
		pAout = mxGetPr(plhs[0]) + i;
		for (j = 0; j < n; j++) {
			if (A[i][j] < 0) {
				*pAout = (double)0;
			} else {
				*pAout = (double)(A[i][j] + 1);
			}

			pAout += rank;
		}
	}
	pRank = mxGetPr(plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL));
	*pRank = (double)rank;

	pA = A;
	for (i = 0; i < m; i++) {
		pA1 = *pA++;
		free(pA1);
	}
	free(A);
	return;
}

#endif

void main()
{
	// DEBUG: SET INPUTS
	m = 8;			// rows; length of array
	n = 8;			// columns; width of each word in bits

	reduceMatrix();
}

void reduceMatrix()
{
	int i = 0, j = 0, k, l, jj;
	int a, b, c, aa, a1, b1;
	int kmax;
	int temp;
	int maxElement;
	int loopCount = 0;

	rank = 0;
	jb = calloc(n, sizeof(int));
	if (jb == NULL) {
		printf("Alloc failed (jb)\n");
		return;
	}
#if LOG_SESSION
	logFile = fopen("RrefGFx.log", "w");

	if (logFile == NULL) {
		printf("Error opening log file\n");
	}
	printMatrix(m, n);
	fprintf(logFile, "\n\n");
#endif

	// Loop over the entire matrix.
	i = 0;
	j = 0;

	// Matrices with >32 columns are stored in an m*nw array as follows:
	// A[0]               = [row 1, columns 1-32]
	// A[1]               = [row 1, columns 33-64]
	// A[2:(nw-1)]        =[       etc     ]
	// A[nw]              = [row 2, columns 1-32]
	// A[nw+1]                 = [row 2, columns 33-64]
	// A[nw:(2*nw-1)] = [          etc     ]
	// A[2*nw:m*nw] = [            etc     ]

	while ((i < m) & (j < n)) {
#if LOG_SESSION
		fprintf(logFile, "row %d, column %d\n", i, j);
#endif
		// Find index of largest element in the remainder of column j.
		// [p,k] = max(abs(A(i:m,j)));
		maxElement = -INF;
		for (k = i; k < m; k++) {
			if (A[k][j] > maxElement) {
				maxElement = A[k][j];
				kmax = k;
			}
			loopCount++;	// DEBUG
		}

		if (maxElement < 0) {
			// The column is negligible, zero it out
			for (k = i; k < m; k++) {
				A[k][j] = -INF;
			}
		} else {
			k = kmax;

#if LOG_SESSION
			fprintf(logFile, "Pivot row k = %d\n", k);
#endif

			// Remember column index
			// jb = [jb j];
			jb[rank++] = j;
#if LOG_SESSION
			fprintf(logFile, "Rank = %d\n", rank);

			// Swap i-th and k-th rows.
			// A([i k],j:n) = A([k i],j:n);
			fprintf(logFile, "Swap rows k,i (%d,%d)\n", i, k);
#endif
			for (l = 0; l < n; l++) {
				temp = A[i][l];
				A[i][l] = A[k][l];
				A[k][l] = temp;
			}

			// Divide the pivot row by the pivot element.
			// A(i,j:n) = A(i,j:n)/A(i,j);
			a = A[i][j];
			if (a < 0) {
				printf("Error: divide by zero\n");
			} else {
				for (jj = j; jj < n; jj++) {
					if (A[i][jj] >= a) {
						// Subtract powers to divide elements
						A[i][jj] =
						    (A[i][jj] - a) % (q - 1);
					} else if (A[i][jj] >= 0) {
						// Subtract powers to divide elements
						// Note: must wrap negative result into range [0,q-2]
						A[i][jj] =
						    (A[i][jj] - a + q -
						     1) % (q - 1);
					}
				}
			}

			// Subtract multiples of the pivot row from all the other rows.
			// for k = [1:i-1 i+1:m]
			for (k = 0; k < m; k++) {

				if (k != i) {
					aa = A[k][j];
					if (aa >= 0) {
						//A(k,j:n) = A(k,j:n) - A(k,j)*A(i,j:n);
						for (jj = j; jj < n; jj++) {
							if (A[i][jj] >= 0) {
								a1 = (aa +
								      A[i][jj])
								    % (q - 1);
							} else {
								a1 = -INF;
							}

							if (a1 >= 0) {
								if (a1 >= q - 1) {
									printf
									    ("Error: invalid symbol (a1)\n");
									return;
								}
								// Note: Matlab version adds 2 due to index starting at 1
								a = (int)
								    vectors[a1 +
									    1];
							} else {
								a = 0;
							}

							b1 = A[k][jj];
							if (b1 >= 0) {
								if (b1 >= q - 1) {
									printf
									    ("Error: invalid symbol (b1)\n");
									return;
								}
								b = (int)
								    vectors[b1 +
									    1];
							} else {
								b = 0;
							}

							c = a ^ b;
							if (c > 0) {

								if (c >= q) {
									printf
									    ("Error: Symbol index too large\n");
									return;
								}
								A[k][jj] =
								    (int)
								    symbols[c] -
								    1;
							} else {
								A[k][jj] = -INF;
							}
						}	// for jj=j:n
					}	// if aa >= 0
				}	// if(k!=i)
			}	// for k = [1:i-1 i+1:m]

			i++;
		}		// if(maxElementFound)

#if LOG_SESSION
		printMatrix(m, n);
#endif
		j++;
	}			// while((i < m) & (j < n))

	free(jb);

#if LOG_SESSION
	// Print out results
	fprintf(logFile, "\n------------\n");
	printMatrix(rank, n);
	for (i = 0; i < rank; i++) {
		fprintf(logFile, "jb[%d] = %d\n", i, jb[i]);
	}
	fprintf(logFile, "\nRank = %d\n\n", rank);
	fprintf(logFile, "Loop count = %d\n", loopCount);
	fclose(logFile);
#endif
}

#if LOG_SESSION

void printMatrix(int rows, int columns)
{
	int i, j;

	fprintf(logFile, "A[%d,%d]:\n", rows, columns);
	for (i = 0; i < rows; i++) {
		fprintf(logFile, "%2d: ", i);
		for (j = 0; j < columns; j++) {
			fprintf(logFile, "%d\t", A[i][j]);
		}
		fprintf(logFile, "\n");
	}
	fprintf(logFile, "\n");
}

#endif
APPENDIX C:Matlab Data Generation Source Code
