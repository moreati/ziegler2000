
/*
% function [G,rank] = RrefGF2(A)
% Returns reduced row echelon form of A and rank of A;
% matrix reduction is done in GF(2) algebra.
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#define LOG_SESSION            0

void reduceMatrix();
#if LOG_SESSION
void printMatrix(int rows);
FILE *logFile;
#endif

#ifdef MATLAB_MEX_FILE
double *pAin;			// double precision input matrix from Matlab workspace
double *pAout;			// double precision reduced output matrix to Matlab workspace
double *pRank;			// double precision rank
unsigned long *A;		// input matrix, represented as a vector of binary numbers
#else
unsigned int A[10] = { 0xFF, 0x7F, 0x3F, 0x1F, 0xF, 0x7, 0x3, 0x1 };
#endif

int m, n;			// number of rows, columns of A
int nw;				// number of parallel words concatenated to form all columns of A
int *jb;
int rank;
unsigned long *pA1;		// input matrix pointer
unsigned long *pA2;		// input matrix pointer

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int i = 0, j = 0, k;
	unsigned int bit, word;
	int maxColumn;

	if (nrhs != 1) {
		mexErrMsgTxt("Not enough input arguments!");
		return;
	}

	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	nw = (int)ceil((double)n / (double)32);
#if LOG_SESSION
	printf("A[%d,%d]:\n", m, n);
	printf("nw = %d\n", nw);
#endif

	A = calloc(m * nw, sizeof(unsigned long));
	if (A == NULL) {
		printf("Alloc failed\n");
		return;
	}
	pA1 = A;

	// Read input matrix from Matlab workspace and convert to binary arrays.
	if (n <= 32) {
		for (i = 0; i < m; i++) {
			pAin = mxGetPr(prhs[0]) + i;
			word = 0;
			for (j = 0; j < n; j++) {
				bit = (int)floor(*pAin + 0.5);
				pAin += m;
				if (bit != 0 && bit != 1) {
					printf("A[%d,%d] = %f\n", i, j,
					       *(pAin - 1));
					mexErrMsgTxt("Nonbinary matrix!");
					return;
				}
				word |= bit << (n - j - 1);
			}
			*pA1++ = (unsigned long)word;
		}
	}
	// Matrices with >32 columns are stored in an m*nw array as follows:
	// A[0]               = [row 1, columns 1-32]
	// A[1]               = [row 1, columns 33-64]
	// A[2:(nw-1)]        =[      etc     ]

//     A[nw]          = [row 2, columns 1-32]
//     A[nw+1]             = [row 2, columns 33-64]
//     A[nw:(2*nw-1)] = [      etc     ]
//     A[2*nw:m*nw] = [        etc     ]

	else {
		for (i = 0; i < m; i++) {
			pAin = mxGetPr(prhs[0]) + i;
			for (k = 0; k < nw; k++) {
				word = 0;
				if (n < (k + 1) * 32) {
					maxColumn = n - k * 32;
				} else {
					maxColumn = 32;
				}

				for (j = 0; j < maxColumn; j++) {
					bit = (int)floor(*pAin + 0.5);
					pAin += m;
					if (bit != 0 && bit != 1) {
						printf("A[%d,%d] = %f\n", i, j,
						       *(pAin - m));
						mexErrMsgTxt
						    ("Nonbinary matrix!");
						return;
					}
					word |= bit << (maxColumn - j - 1);
				}
				*pA1++ = (unsigned long)word;
			}
		}
	}

	reduceMatrix();

	plhs[0] = mxCreateDoubleMatrix(rank, n, mxREAL);
	if (n <= 32) {
		for (i = 0; i < rank; i++) {
			pAout = mxGetPr(plhs[0]) + i;
			for (j = n - 1; j >= 0; j--) {

				bit = (A[i] >> j & 0x1);
				*pAout = (double)bit;
				pAout += rank;
			}
		}
	} else {
		pA1 = A;
		for (i = 0; i < rank; i++) {
			pAout = mxGetPr(plhs[0]) + i;
			for (k = 0; k < nw; k++) {
				if (n < (k + 1) * 32) {
					maxColumn = n - k * 32;
				} else {
					maxColumn = 32;
				}

				for (j = maxColumn - 1; j >= 0; j--) {
					bit = (*pA1 >> j & 0x1);
					*pAout = (double)bit;
					pAout += rank;
				}
				pA1++;
			}
		}
	}
	pRank = mxGetPr(plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL));
	*pRank = (double)rank;
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
	int i = 0, j = 0, k, l;
	unsigned int columnMask, temp;
	int leadingOneFound = 0;
	int loopCount = 0;
	int j32 = 0;
	int jnw = 0;

	rank = 0;
	jb = calloc(n, sizeof(int));
	if (jb == NULL) {
		printf("Alloc failed\n");
	}
#if LOG_SESSION
	logFile = fopen("RrefGF2.log", "w");
	if (logFile == NULL) {
		printf("Error opening log file\n");
	}
	fprintf(logFile, "A[%d,%d]:\n", m, n);
	printMatrix(m);
	fprintf(logFile, "\n\n");
#endif

	// Loop over the entire matrix.
	i = 0;
	j = 0;
	if (n <= 32) {
		while ((i < m) & (j < n)) {
#if LOG_SESSION
			fprintf(logFile, "row %d, column %d\n", i, j);
#endif
			// Find index of largest element in the remainder of column j (i.e., leading one)
			// [p,k] = max(abs(A(i:m,j)));
			leadingOneFound = 0;
			columnMask = 0x1 << (n - j - 1);
			for (k = i; k < m; k++) {
				if (A[k] & columnMask) {
#if LOG_SESSION
					fprintf(logFile, "Pivot row k = %d\n",
						k);

#endif
					leadingOneFound = 1;
					break;
				}
				loopCount++;	// DEBUG
			}
			if (leadingOneFound) {
				// Remember column index
				// jb = [jb j];
				jb[rank++] = j;
#if LOG_SESSION
				fprintf(logFile, "Rank = %d\n", rank);
				// Swap i-th and k-th rows.
				// A([i k],j:n) = A([k i],j:n);
				fprintf(logFile, "Swap rows k,i (%d,%d)\n", i,
					k);
#endif
				temp = A[i];
				A[i] = A[k];
				A[k] = temp;

				// Subtract multiples of the pivot row from all the other rows.
				// for k = [1:i-1 i+1:m]
				for (k = 0; k < i; k++) {
					// A(k,j:n) = rem(A(k,j:n) + A(k,j)*A(i,j:n), q);
					if (A[k] & columnMask) {
						A[k] ^= A[i];
					}
					loopCount++;	// DEBUG
				}
				for (k = i + 1; k < m; k++) {
					// A(k,j:n) = rem(A(k,j:n) + A(k,j)*A(i,j:n), q);
					if (A[k] & columnMask) {
						A[k] ^= A[i];
					}
					loopCount++;	// DEBUG
				}

				i++;
			}
#if LOG_SESSION
			printMatrix(m);
#endif
			j++;
		}

	}

	else {
		// Matrices with >32 columns are stored in an m*nw array as follows:
		// A[0]               = [row 1, columns 1-32]
		// A[1]               = [row 1, columns 33-64]
		// A[2:(nw-1)]        =[      etc     ]
		// A[nw]              = [row 2, columns 1-32]
		// A[nw+1]            = [row 2, columns 33-64]
		// A[nw:(2*nw-1)] = [         etc     ]
		// A[2*nw:m*nw] = [           etc     ]

		j32 = 0;
		jnw = 0;
		while ((i < m) & (j < n)) {
#if LOG_SESSION
			fprintf(logFile, "row %d, column %d\n", i, j);
#endif
			// Find index of largest element in the remainder of column j.(i.e., leading one)
			// [p,k] = max(abs(A(i:m,j)));
			leadingOneFound = 0;
			if (jnw < nw - 1) {
				columnMask = 0x1 << (31 - j32);
			} else {
				columnMask = 0x1 << (n - jnw * 32 - j32 - 1);
			}

			pA1 = A + i * nw + jnw;
			for (k = i; k < m; k++) {
				if (*pA1 & columnMask) {
#if LOG_SESSION
					fprintf(logFile, "Pivot row k = %d\n",
						k);
#endif
					leadingOneFound = 1;
					break;
				}
				pA1 += nw;
				loopCount++;	// DEBUG
			}

			if (leadingOneFound) {
				// Remember column index
				// jb = [jb j];
				jb[rank++] = j32 + jnw * 32;
#if LOG_SESSION
				fprintf(logFile, "Rank = %d\n", rank);
				// Swap i-th and k-th rows.
				// A([i k],j:n) = A([k i],j:n);
				fprintf(logFile, "Swap rows k,i (%d,%d)\n", i,
					k);
#endif
				pA1 = A + i * nw;
				pA2 = A + k * nw;
				for (l = 0; l < nw; l++) {
					temp = *pA1;
					*pA1++ = *pA2;
					*pA2++ = temp;
				}

				// Subtract multiples of the pivot row from all the other rows.
				// for k = [1:i-1 i+1:m]
				for (k = 0; k < i; k++) {
					pA2 = A + k * nw;
					// A(k,j:n) = rem(A(k,j:n) + A(k,j)*A(i,j:n), q);
					if (*(pA2 + jnw) & columnMask) {
						pA1 = A + i * nw;
						for (l = 0; l < nw; l++) {
							*(pA2 + l) ^= *pA1++;
							loopCount++;	// DEBUG
						}
					}
				}

				for (k = i + 1; k < m; k++) {
					pA2 = A + k * nw;
					// A(k,j:n) = rem(A(k,j:n) + A(k,j)*A(i,j:n), q);
					if (*(pA2 + jnw) & columnMask) {
						pA1 = A + i * nw;
						for (l = 0; l < nw; l++) {
							*(pA2 + l) ^= *pA1++;
							loopCount++;	// DEBUG
						}

					}
				}

				i++;
			}
#if LOG_SESSION
			printMatrix(m);
#endif
			j++;
			j32++;
			if (j32 >= 32) {
				j32 = 0;
				jnw++;
			}
		}
	}

#if LOG_SESSION
	// Print out results
	fprintf(logFile, "\n------------\n");
	printMatrix(rank);
	for (i = 0; i < rank; i++) {
		fprintf(logFile, "jb[%d] = %d\n", i, jb[i]);
	}
	fprintf(logFile, "\nRank = %d\n\n", rank);
	fprintf(logFile, "Loop count = %d\n", loopCount);
	fclose(logFile);
#endif
}

#if LOG_SESSION

void printMatrix(int rows)
{
	int i, j, k;
	int maxColumn;

	if (rows > 10) {
		return;
	}

	fprintf(logFile, "A:\n");

	if (n <= 32) {
		for (i = 0; i < rows; i++) {
			fprintf(logFile, "%2d: ", i);
			for (j = n - 1; j >= 0; j--) {
				fprintf(logFile, "%d", A[i] >> j & 0x1);
			}
			fprintf(logFile, "\n");
		}
	} else {
		pA1 = A;
		for (i = 0; i < rows; i++) {
			fprintf(logFile, "%2d: ", i);
			for (k = 0; k < nw; k++) {
				if (n < (k + 1) * 32) {
					maxColumn = n - k * 32;
				} else {
					maxColumn = 32;
				}

				for (j = maxColumn - 1; j >= 0; j--) {
					fprintf(logFile, "%d", *pA1 >> j & 0x1);
				}
				pA1++;
			}
			fprintf(logFile, "\n");
		}
	}
	fprintf(logFile, "\n");
}

#endif
