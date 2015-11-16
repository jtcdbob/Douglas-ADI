//
//  main.c
//  DouglasADI
//
//  Created by Bob Chen on 11/12/15.
//  Copyright © 2015 Bob. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#define PI 3.14159265358979323846

void transpose(float* restrict u_t, const float* restrict u, const int n){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u_t[j*n+i] = u[i*n+j];
        }
    }
}

//void apply_tri(float* restrict vOut, const float* restrict vIn,
//               const float* restrict diag,  const float* restrict udiag,
//               const float* restrict ldiag, const int n){
//    vOut[0] = diag[0] * vIn[0] + udiag[0] * vIn[1];
//    vOut[n-1] = diag[n-1] * vIn[n-1] + ldiag[n-2] * vIn[n-2];
//    for(int i = 1; i < n-1; i++){
//        vOut[i] = ldiag[i-1]*vIn[i-1] + diag[i]*vIn[i] + udiag[i]*vIn[i+1];
//    }
//}

//void solve_tri(float* restrict x, const float* restrict ldiag,
//               const float* restrict diag, const float* restrict udiag,
//               const int n) {
//    // Allocate scratch space.
//    float* cprime = (float*) malloc(sizeof(float) * n);
//    cprime[0] = udiag[0] / diag[0];
//    x[0] = x[0] / diag[0];
//    // loop from 1 to N - 1 inclusive
//    for (int i = 1; i < n; i++) {
//        float m = 1.0 / (diag[i] - ldiag[i-1] * cprime[i - 1]);
//        cprime[i] = udiag[i] * m;
//        x[i] = (x[i] - ldiag[i-1] * x[i - 1]) * m;
//    }
//    // loop from N - 2 to 0 inclusive, safely testing loop end condition
//    for (int i = n - 1; i-- > 0; )
//        x[i] = x[i] - cprime[i] * x[i + 1];
//    // free scratch space
//    free(cprime);
//}

void apply_tri_special(float* restrict u, float* restrict u_temp, const float a,  const float b, const int n){
    /*
     This operator is designed to apply a tridiagonal matrix multiplication with the following
     structure:

     0  0   0   0   0   0
     a  b   a   0   0   0
     0  a   b   a   0   0
     ......
     0  0   0   a   b   a
     0  0   0   0   0   0

     There is little room for generalizing its application but it should yield great performance boost
     */
//    float * u_copy = (float*) malloc(n * sizeof(float));
    memcpy(u_temp, u, n*sizeof(float));
    for(int i = 1; i < n-1; i++){
        u[i] = a*(u_temp[i-1]+u_temp[i+1]) + b * u_temp[i];
    }
    u[0] = 0.0;
    u[n-1] = 0.0;
//    free(u_copy);
}
void apply_tri_special_plus(float* restrict u, float* restrict u_temp, const float a,  const float b, const int n){
    /*
     This operator is designed to apply a tridiagonal matrix multiplication with the following
     structure:

     1  0   0   0   0   0
     a  b+1 a   0   0   0
     0  a   b+1 a   0   0
     ......
     0  0   0   a   b+1 a
     0  0   0   0   0   1

     There is little room for generalizing its application but it should yield great performance boost
     */
//    float* scratch = (float*) malloc(sizeof(float) * n);
    memcpy(u_temp, u, n*sizeof(float));
    for(int i = 1; i < n-1; i++){
        u[i] += a*(u_temp[i-1]+u_temp[i+1]) + b * u_temp[i];
    }
//    free(scratch);
}

void solve_tri_special(float* restrict x, const float a, const float b, const int n){
//void solve_tri_special(float* restrict x, float* restrict cprime, const float a, const float b, const int n){
     /*
     This solve is designed to solve a tridiagonal linear system using thomas' algorithm with the following structure:

     1  0   0   0   0   0
     -a 1-b -a  0   0   0
     0  -a  1-b -a  0   0
     ......
     0  0   0   -a  1-b -a
     0  0   0   0   0   1

     There is little room for generalizing its application but it should yield great performance boost
     */

    // Allocate scratch space.
    float* restrict cprime __attribute__((aligned(64))) = (float*) _mm_malloc(sizeof(float) * n, 64);
    cprime[0] = 0;
    float astar = -a;
    float bstar = 1-b;
    float m;
    // loop from 1 to N - 2 inclusive
    for (int i = 1; i < n-1; i++) {
        m = 1.0 / (bstar - astar * cprime[i - 1]);
        cprime[i] = astar * m;
        x[i] = (x[i] - astar * x[i - 1]) * m; // I'm not sure about this step, might step into the last one
    }
    cprime[n-1] = 0;
    // loop from N - 2 to 0 inclusive, safely testing loop end condition
    for (int i = n - 1; i-- > 0; )
        x[i] = x[i] - cprime[i] * x[i + 1];

    _mm_free(cprime);
}

void relaxOperation(float * restrict u, const float * restrict fstar, float* scratch, float a, float b, int n){

    float* ustar = (float*) malloc(n*n*sizeof(float));
    float* u_t = (float*) malloc(n*n*sizeof(float));
    memcpy(ustar, u, n*n*sizeof(float));
    transpose(u_t, ustar, n);
    // keep two copies and then do it.
    for (int i = 1; i < n-1; i++) {
        apply_tri_special_plus(u+i*n, scratch, a/2, b/2, n);
    }
    for (int i = 1; i < n-1; i++) {
        apply_tri_special(u_t+i*n, scratch, a, b, n);
    }
    memset(u_t, 0.0, n*sizeof(float));
    memset(u_t+(n-1)*n, 0.0, n*sizeof(float));
    for (int i = 0; i < n*n; i++) {
        u[i] -= fstar[i];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u[i*n + j] += u_t[j*n+i];
        }
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 1; i < n-1; i++) {
        solve_tri_special(u+i*n, a/2, b/2, n);
    }
    transpose(ustar, u, n);
    for (int i = 0; i < n*n; i++) {
        ustar[i] -= 0.5*u_t[i];
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 1; i < n-1; i++) {
        //solve_tri_special(ustar+i*n, scratch, a/2, b/2, n);
        solve_tri_special(ustar+i*n, a/2, b/2, n);
    }
    transpose(u, ustar, n);

    free(ustar);
    free(u_t);
}
float normInf(const float* restrict u, const int n){
    float max = u[0];
    for (int i = 1; i < n*n; i++) {
        max = (max < fabs(u[i]))? fabs(u[i]):max;
    }
    return max;
}
float test_init(const float x, const float y, const float wx, const float wy, const float ax, const float ay){
    float d2Xfact = (2.0*PI*wx)*(2.0*PI*wx);
    float d2Yfact = (2.0*PI*wy)*(2.0*PI*wy);
    return -(cos(2.0*PI*x*wx)*cos(2.0*PI*y*wy))/(ax*d2Xfact + ay*d2Yfact);
}



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

    float alpha       = 2.0;
    float waveNumber  = 1.0;
    float alphaX      = alpha;     // Laplace operator-x prefactor
    float alphaY      = alpha;     // Laplace operator-y prefactor
    float waveNumberX = waveNumber;  // test problem x-coordinate wave number
    float waveNumberY = waveNumber;  // test problem y-coordinate wave number

    int   xPanel = M;  // X panel count
    int   yPanel = M;  // Y panel count

    float xMin = 0.0; float xMax = 1.0;
    float yMin = 0.0; float yMax = 1.0;
    float hx = (xMax-xMin)/(float)xPanel;
    float hy = (yMax-yMin)/(float)yPanel;

    float dt_0 = 4*hx*hx/(alpha*PI*PI);  // Relaxation timestep

    // Echo input parameters
    printf("\n==== Douglas ADI Program Start ====\n");
    printf("** alpha_x/alpha_y : %f/%f\n", alphaX, alphaY);
    printf("** X/Y Panel Count : %d/%d\n", xPanel, yPanel);

    int N = M+1; // vector size
    //float* f      = (float*) malloc(N*N*sizeof(float));
    //float* uk     = (float*) calloc(N*N, sizeof(float)); // initialize to 0
    //float* uLast  = (float*) malloc(N*N*sizeof(float));
    float* restrict f __attribute__((aligned(64))) = (float*) _mm_malloc(N*N*sizeof(float),64);
    float* restrict uk __attribute__((aligned(64))) = (float*) _mm_malloc(N*N*sizeof(float),64);
    float* restrict uLast __attribute__((aligned(64))) = (float*) _mm_malloc(N*N*sizeof(float),64);
    memset(uk, 0.0, N*N);
    float tol = 1.0e-6;  // Stopping tolerance


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
    float x;
    float y;
//    for (int i = 0; i < N*N; i++) {
//        uk[i] = i;
//    }

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
    float dt = dt_0;

    float ax = dt*alphaX/(hx*hx);
    float ay = dt*alphaY/(hy*hy);
    float bx = -2.0 * ax;
    float by = -2.0 * ay;

    float diffNorm = 2*tol;
    int   iterMax  = 200;
    int   iter     = 0;
    float ukNorm = 0;

//    ax = 1;
//    bx = 2;
//    float a[5] = {4,18,3,1,1};
//    float b[5];
//    solve_tri_special(a, ax, bx, 5);
//    apply_tri_special_plus(a, ax, bx, 5);
//    b[0]=b[0];

    //float* fstar = (float*) malloc(maxSweepSize*N*N*sizeof(float));
    float* restrict fstar __attribute__((aligned(64))) = (float*) _mm_malloc(maxSweepSize*N*N*sizeof(float),64);
    for (int i = 0; i < maxSweepSize; i++) {
        int offset = i*N*N;
        for (int j = 0; j < N*N; j++) {
            fstar[offset + j] = f[j]*dt;
        }
        dt *= 2*2;
    }

    float* restrict scratch __attribute__((aligned(64))) = (float*) _mm_malloc(N*sizeof(float),64);

    printf("\n");

//    relaxOperation(uk, fstar, scratch, ax, bx, N);
//    for (int i = 0 ; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            printf("%f \t", uk[i*N+j]);
//        }
//        printf("\n");
//    }

    float t0 = omp_get_wtime();
    while((diffNorm > tol)&&(iter < iterMax)){
        memcpy(uLast, uk, N*N*sizeof(float)); // This can be avoided by swapping the roles for each iteration
        for (int j = 0; j < maxSweepSize; j++) {
            relaxOperation(uk, fstar+j*N*N, scratch, ax, bx, N);
//            dt *= 2*2;
            ax *= 2*2;
            bx *= 2*2;
        }
//        printf("print final values:\n");
//        for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                printf("%f\t", uk[i*N + j]);
//            }
//            printf("\n");
//        }
//
        for (int j = maxSweepSize-1; j > -1 ; j--) {
//            dt /= 2*2;
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
    float t1 = omp_get_wtime();

    _mm_free(scratch);
    _mm_free(fstar);
    _mm_free(f);
    _mm_free(uk);
    _mm_free(uLast);

    printf("==== Multi-scale Timstep Relaxation Output XXXX ====\n");
    printf("** Difference between iterates : \t\t\t\t%3.10f\n",diffNorm);
    printf("** Iteration Count for multi-scale time scheme :\t %d \n\n",iter*2*maxSweepSize);
    // Timing results
//    printf("** Time it takes to initialize the multi-scale timestep solution is %e ms\n", elapsed_1);
    printf("** Time it takes to solve the multi-scale timestep solution is %e ms\n", t1-t0);
//    printf("** The total time is is %e ms\n", elapsed_1+elapsed_2);

    return 0;
}
