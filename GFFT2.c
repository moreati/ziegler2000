
/*
% function [X] = GFFT2(x,WordCount,n,qm,mx)
%     Galois Field Fourier Transform:
%     Time-domain vector x is over GF(q), elements are q-ary
%     Frequency-domain vector X is over GF(q^m), elements are q^m-ary
%     For BCH: n = (q^m-1)
%     For RS: n = q-1, (m=1), hence x and X are both in GF(q)
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#define LOG_SESSION           0

void calcGFFT();
#if LOG_SESSION
void printMatrix(int *A);
FILE *logFile;
#endif


#ifdef MATLAB_MEX_FILE
double *pAin;          // double precision input matrix from Matlab workspace
double *pVector;       // GF(q^m) vectors input from Matlab workspace
double *pSymbol;       // GF(q^m) symbols input from Matlab workspace
double *pAout;         // double precision reduced output matrix to Matlab workspace
int *x;                // input matrix, represented as a vector of binary numbers
int *vectors,*pVectors;         // GF(q^m) vectors
int *symbols,*pSymbols;         // GF(q^m) symbols
#else
// DEBUG C CODE
int x[7] = {1,1,0,0,1,0,0};

                                                  121




int vectors[8] = {0,1,2,4,3,6,7,5};
int symbols[8] = {0,1,2,4,3,7,5,6};
#endif

int *X;                // input matrix, represented as a vector of binary numbers
int m,n;               // number of rows, columns of x
int q,qm;              // alphabet size q, Galois field q^m
int *px,*pX;           // input/output matrix pointer




#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i=0,j=0,k;

    if (nrhs < 5)
    {
         mexErrMsgTxt("Not enough input arguments!");
         return;
    }

     m = mxGetM(prhs[0]);
     n = mxGetN(prhs[0]);
     q = (int)mxGetScalar(prhs[3]);
     qm = (int)mxGetScalar(prhs[4]);
#if LOG_SESSION
     printf("x[%d,%d]:\n",m,n);
#endif

    x = calloc(m*n,sizeof(int));
    if(x==NULL)
    {
         printf("x alloc failed\n");
         return;
    }
    px = x;
    X = calloc(m*n,sizeof(int));
    if(X==NULL)
    {
         printf("X alloc failed\n");
         return;
    }
    pX = X;
    vectors = calloc(qm,sizeof(int));
    if(vectors==NULL)
    {
         printf("vectors alloc failed\n");




     return;
}
symbols = calloc(qm,sizeof(int));
if(symbols==NULL)
{
     printf("symbols alloc failed\n");
     return;
}

// Read input matrix from Matlab workspace; store columns
// such that code words are sequential
for(i=0; i<m; i++)
{
     pAin = mxGetPr(prhs[0]) + i;
     for(j=0; j<n; j++)
     {
          // Convert symbols {0,1,2=alpha,...} to powers {-inf=0,0=1,1=alpha,...}
          // by subtracting 1 (-1 represents 0 = alpha^-inf)
          *px++ = (int)*pAin - 1;
          pAin += m;
     }
}
// Read Galois field vectors and symbols from Matlab workspace
pVector = mxGetPr(prhs[1]);
pVectors = vectors;
for(i=0; i<qm; i++)
{
     *pVectors++ = (int)*pVector++;
}
pSymbol = mxGetPr(prhs[2]);
pSymbols = symbols;
for(i=0; i<qm; i++)
{
     *pSymbols++ = (int)*pSymbol++;
}

// Take GFFT of each row of input matrix
calcGFFT();

plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
pX = X;
for(i=0; i<m; i++)
{
     pAout = mxGetPr(plhs[0]) + i;
     for(j=n-1; j>=0; j--)
     {
          *pAout = (double)*pX++;
          pAout += m;
     }
}




    return;
}

#endif




void main()
{
    int i;
    int *x2;

    // DEBUG: SET INPUTS
    m = 1;        // rows
    n = 7;        // columns
    q = 2;
    qm = 8;
    x2 = calloc(m*n,sizeof(int));
    if(x2==NULL)
    {
         printf("x2 alloc failed\n");
         return;
    }
    for(i=0; i<n; i++)
    {
         x2[i] = x[n-i-1]-1;          // convert to powers of alpha
    }
    for(i=0; i<n; i++)
    {
         x[i] = x2[i];
    }

    X = calloc(m*n,sizeof(int));
    if(X==NULL)
    {
        printf("X alloc failed\n");
        return;
    }
    pX = X;

    calcGFFT();
}


void calcGFFT()
{
    int K1,K2,qm1;
    int i,j,k,l,ij;




    int xij,Xij;

#if LOG_SESSION
     logFile = fopen("Gfft2.log","w");
     if(logFile==NULL)
     {
          printf("Error opening log file\n");
     }
     fprintf(logFile,"q = %d, q^m = %d\n",q,qm);
     printMatrix(x);
     fprintf(logFile,"\n\n");
#endif

    qm1 = qm-1;
    K1 = (int)(qm1/n);
    K2 = (int)(qm1/(q-1));
    pX = X;

    for(l=0; l<m; l++)
    {
         for(j=0; j<n; j++)        // frequency index
         {
              Xij = 0;
              px = x + n*l;
              for(i=0; i<n; i++)   // time index
              {
                   // Convert GFFT coefficients from powers of gamma to alpha:
                   // power of gamma, where gamma is an nth root of unity in GF(q^m)
                   ij = (i*j)%n;
                   // power of alpha, where alpha is an (q^m-1) root of unity in GF(q^m)
                   ij = (ij*K1)%qm1;

                   // Convert input coefficients from powers of beta to alpha:
                   // power of beta, where beta is an (q-1) root of unity in GF(q)
                   k = *px++;
                   if(k >= 0) // only add if input is nonzero (power > -inf)
                   {
                        // power of alpha, where alpha is an (q^m-1) root of unity in GF(q^m)
                        k = (k*K2)%qm1;
                        xij = (ij+k)%qm1;
                        Xij ^= vectors[xij+1];      // XOR is GF(qm) addition
                   }
              }
              *pX++ = symbols[Xij];

#if LOG_SESSION
           fprintf(logFile,"X[%d] = %d\n",j,symbols[Xij]);
           printMatrix(X);
#endif
       }




    }

#if LOG_SESSION
     // Print out results
     fprintf(logFile,"\n------------\n");
     printMatrix(m);
     fprintf(logFile,"Loop count = %d\n",loopCount);
     fclose(logFile);
#endif
}



#if LOG_SESSION

void printMatrix(int *A)
{
    int i,j,k;
    int maxColumn;
    int *pA;

    if(m>10)
    {
         return;
    }
    fprintf(logFile,"A[%d,%d]:\n",m,n);

    pA = A;
    for(i=0; i<m; i++)
    {
         fprintf(logFile,"%2d: ",i);
         for(j=0; j<n; j++)
         {
              fprintf(logFile,"%d ",*pA++);
         }
         fprintf(logFile,"\n");
    }
    fprintf(logFile,"\n");
}

#endif




