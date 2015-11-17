#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>
#include "my_functions.h"


#ifndef PI
#define PI 3.14159265358979323846
#endif

int main(int argc, char * argv[]) {
    int M = 200; // Problem size, assume square case
    int threadCount = -1;

    int c;
    extern char* optarg;
    const char* optstring = "n:p:";
    while ((c = getopt(argc, argv, optstring)) !=-1){
        switch (c) {
            case 'n': M           = atoi(optarg); break;
            case 'p': threadCount = atoi(optarg); break;
        }
    }

    double alpha       = 2.0;
    double waveNumber  = 1.0;
    double alphaX      = alpha;     // Laplace operator-x prefactor
    double alphaY      = alpha;     // Laplace operator-y prefactor
    double waveNumberX = waveNumber;  // test problem x-coordinate wave number
    double waveNumberY = waveNumber;  // test problem y-coordinate wave number

    int   xPanel = M;  // X panel count
    int   yPanel = M;  // Y panel count

    double xMin = 0.0; double xMax = 1.0;
    double yMin = 0.0; double yMax = 1.0;
    double hx = (xMax-xMin)/(double)xPanel;
    double hy = (yMax-yMin)/(double)yPanel;

    double dt_0 = 4*hx*hx/(alpha*PI*PI);  // Relaxation timestep

    // Echo input parameters
    printf("\n==== Douglas ADI Program Start ====\n");
    printf("** alpha_x/alpha_y : %f/%f\n", alphaX, alphaY);
    printf("** X/Y Panel Count : %d/%d\n", xPanel, yPanel);

    int N = M+1; // vector size
    double* restrict f __attribute__((aligned(64))) = (double*) _mm_malloc(N*N*sizeof(double),64);
    double* restrict uk __attribute__((aligned(64))) = (double*) _mm_malloc(N*N*sizeof(double),64);
    double* restrict uLast __attribute__((aligned(64))) = (double*) _mm_malloc(N*N*sizeof(double),64);
    memset(uk, 0.0, N*N);
    double tol = 1.0e-6;  // Stopping tolerance

    __assume_aligned(f, 64);
    __assume_aligned(uk, 64);
    __assume_aligned(uLast, 64);
    /*
        -> y    (j)
      ---------------------
     |
     || x
     |v
     |
     | (i)
     |

     */
    double x;
    double y;

    for(int j = 0; j < N; j++){
        y =  yMin + j*hy;
        for(int i = 0; i < N; i++){
            x = xMin + i*hx;
            f[j*N + i] = cos(2.0*PI*x*waveNumberX)*cos(2.0*PI*y*waveNumberY);
        }
    }
    int index = 0;
    for(int i = 0; i < N; i++){ // Update x
        x   =   xMin + i*hx;
        y   =   yMin + index*hy;
        uk[i]  = test_init(x, y, waveNumberX, waveNumberY, alphaX, alphaY);
        f [i]  = 0.0;
    }
    index   = M;
    for(int i = 0; i < N; i++){
        x   = xMin + i*hx;
        y   = yMin + index*hy;
        uk[index*N + i]  = test_init(x, y, waveNumberX, waveNumberY, alphaX, alphaY);
        f [index*N + i]  = 0.0;
    }
    index = 0;
    for(int i = 0; i < N; i++){
        x = xMin + index*hx;
        y = yMin + i*hy;
        uk[index + i*N]  = test_init(x, y, waveNumberX, waveNumberY, alphaX, alphaY);
        f [index + i*N]  = 0.0;
    }
    index = M;
    for(int i = 0; i < N; i++){
        x = xMin + index*hx;
        y = yMin + i*hy;
        uk[index + i*N]  = test_init(x, y, waveNumberX, waveNumberY, alphaX, alphaY);
        f [index + i*N]  = 0.0;
    }

   // Set up OpenMP parameters
#ifdef _OPENMP
    if(threadCount > omp_get_max_threads() || threadCount <= 0) omp_set_num_threads(omp_get_max_threads());
    else omp_set_num_threads(threadCount);
    printf("\n==>Using OpenMP With %d Threads\n\n",omp_get_max_threads());
#endif

    // Initialize solver variables
    const int maxSweepSize = (int)ceil(log(M)/log(2)) + 1;
    double dt = dt_0;

    double ax = dt*alphaX/(hx*hx);
    double ay = dt*alphaY/(hy*hy);
    double bx = -2.0 * ax;
    double by = -2.0 * ay;

    double diffNorm = 2*tol;
    int   iterMax  = 200;
    int   iter     = 0;
    double ukNorm = 0;

    double* restrict fstar __attribute__((aligned(64))) = (double*) _mm_malloc(maxSweepSize*N*N*sizeof(double),64);
    for (int i = 0; i < maxSweepSize; i++) {
        int offset = i*N*N;
        for (int j = 0; j < N*N; j++) {
            fstar[offset + j] = f[j]*dt;
        }
        dt *= 2*2;
    }

    double* restrict scratch __attribute__((aligned(64))) = (double*) _mm_malloc(N*sizeof(double),64);

    printf("\n");

    double t0 = omp_get_wtime();
    while((diffNorm > tol)&&(iter < iterMax)){
        memcpy(uLast, uk, N*N*sizeof(double)); // This can be avoided by swapping the roles for each iteration
        for (int j = 0; j < maxSweepSize; j++) {
            relaxOperation(uk, fstar+j*N*N, scratch, ax, bx, N);
            ax *= 2*2;
            bx *= 2*2;
        }
        for (int j = maxSweepSize-1; j > -1 ; j--) {
            ax /= 2*2;
            bx /= 2*2;
            relaxOperation(uk, fstar+j*N*N, scratch, ax, bx, N);
        }
        ukNorm = normInf(uk, N);    // Check difference between iterates
        for ( int i = 0; i < N*N; i++) {
            uLast[i] -= uk[i];
        }
        diffNorm = normInf(uLast, N) / ukNorm; // uLast.normInf();
        iter++; // Update iterate
    }
    double t1 = omp_get_wtime();

    _mm_free(scratch);
    _mm_free(fstar);
    _mm_free(f);
    _mm_free(uk);
    _mm_free(uLast);

    printf("==== Multi-scale Timstep Relaxation Output XXXX ====\n");
    printf("** Difference between iterates : \t\t\t\t%3.10f\n",diffNorm);
    printf("** Iteration Count for multi-scale time scheme :\t %d \n\n",iter*2*maxSweepSize);
    // Timing results
    printf("** Time it takes to solve the multi-scale timestep solution is %e ms\n", t1-t0);

    return 0;
}
