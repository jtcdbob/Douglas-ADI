#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

// New transpose function that deals with cache miss problem
// Exact implementation is modified from
// http://stackoverflow.com/questions/5200338/a-cache-efficient-matrix-transpose-program
void transpose_cache(double* restrict dst,const double* restrict src, size_t n){
    __assume_aligned(dst, 64);
    __assume_aligned(src, 64);
    size_t block=64;
    for(size_t i = 0;i < n-block; i += block){
        for(size_t j=0; j < n; ++j ){
            for(size_t b = 0; b < block ; ++b){
                dst[j*n+b]=src[b*n+j];
            }
        }
    }
    size_t start = n - n%block;
    for(size_t j=0; j < n; ++j ){
        for(size_t b = start; b < n; ++b){
            dst[j*n+b]=src[b*n+j];
        }
    }
}

void vec_subtract(double* restrict u, const double* restrict v, size_t n){
    __assume_aligned(u, 64);
    __assume_aligned(v, 64);
    size_t block=128;
    for(size_t i = 0;i < n-block; i += block){
        for(size_t j=0; j < n; ++j ){
            for(size_t b = 0; b < block ; ++b){
                u[j*n+b] -= v[j*n+b];
            }
        }
    }
    size_t start = n - n%block;
    for(size_t j=0; j < n; ++j ){
        for(size_t b = start; b < n; ++b){
            u[j*n+b] -= v[j*n+b];
        }
    }
}

void vec_transpose_add(double* restrict u, const double* restrict v, size_t n){
    __assume_aligned(u, 64);
    __assume_aligned(v, 64);
    size_t block=128;
    for(size_t i = 0;i < n-block; i += block){
        for(size_t j=0; j < n; ++j ){
            for(size_t b = 0; b < block ; ++b){
                u[b*n+j] += v[j*n+b];
                //u[j*n+b] += v[b*n+j];
            }
        }
    }
    size_t start = n - n%block;
    for(size_t j=0; j < n; ++j ){
        for(size_t b = start; b < n; ++b){
            u[b*n+j] += v[j*n+b];
           //u[j*n+b] += v[b*n+j];
        }
    }
}

void vec_subtract_half(double* restrict u, const double* restrict v, size_t n){
    __assume_aligned(u, 64);
    __assume_aligned(v, 64);
    size_t block=128;
    for(size_t i = 0;i < n-block; i += block){
        for(size_t j=0; j < n; ++j ){
            for(size_t b = 0; b < block ; ++b){
                u[j*n+b] -= 0.5*v[j*n+b];
            }
        }
    }
    size_t start = n - n%block;
    for(size_t j=0; j < n; ++j ){
        for(size_t b = start; b < n; ++b){
            u[j*n+b] -= 0.5*v[j*n+b];
        }
    }
}

void apply_tri_special(double* restrict u, double* restrict u_temp, const double a,  const double b, const int n){
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
    __assume_aligned(u, 64);
    __assume_aligned(u_temp, 64);

    memcpy(u_temp, u, n*sizeof(double));
    for(int i = 1; i < n-1; i++){
        u[i] = a*(u_temp[i-1]+u_temp[i+1]) + b * u_temp[i];
    }
    u[0] = 0.0;
    u[n-1] = 0.0;
}
void apply_tri_special_plus(double* restrict u, double* restrict u_temp, const double a,  const double b, const int n){
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
    __assume_aligned(u, 64);
    __assume_aligned(u_temp, 64);

    memcpy(u_temp, u, n*sizeof(double));
    for(int i = 1; i < n-1; i++){
        u[i] += a*(u_temp[i-1]+u_temp[i+1]) + b * u_temp[i];
    }
}

void solve_tri_special(double* restrict x, const double a, const double b, const int n){
//void solve_tri_special(double* restrict x, double* restrict cprime, const double a, const double b, const int n){
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
    double* restrict cprime __attribute__((aligned(64))) = (double*) _mm_malloc(sizeof(double) * n, 64);
    cprime[0] = 0;
    double astar = -a;
    double bstar = 1-b;
    double m;
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

void relaxOperation(double * restrict u, const double * restrict fstar, double* scratch, double a, double b, int n){

    double* restrict ustar __attribute__((aligned(64))) = (double*) _mm_malloc(n*n*sizeof(double),64);
    double* restrict u_t __attribute__((aligned(64))) = (double*) _mm_malloc(n*n*sizeof(double),64);
    memcpy(ustar, u, n*n*sizeof(double));
    transpose_cache(u_t, ustar, n);
    // keep two copies and then do it.
    for (int i = 1; i < n-1; i++) {
        apply_tri_special_plus(u+i*n, scratch, a/2, b/2, n);
    }
    for (int i = 1; i < n-1; i++) {
        apply_tri_special(u_t+i*n, scratch, a, b, n);
    }
    memset(u_t, 0.0, n*sizeof(double));
    memset(u_t+(n-1)*n, 0.0, n*sizeof(double));
    vec_subtract(u, fstar, n);
    vec_transpose_add(u, u_t, n);
#ifdef _OPENMP
//#pragma omp parallel for
#endif
    for (int i = 1; i < n-1; i++) {
        solve_tri_special(u+i*n, a/2, b/2, n);
    }
    transpose_cache(ustar, u, n);
    vec_subtract_half(ustar, u_t, n);
#ifdef _OPENMP
//#pragma omp parallel for
#endif
    for (int i = 1; i < n-1; i++) {
        solve_tri_special(ustar+i*n, a/2, b/2, n);
    }
    transpose_cache(u, ustar, n);

    _mm_free(ustar);
    _mm_free(u_t);
}
double normInf(const double* restrict u, const int n){
    double max = u[0];
    for (int i = 1; i < n*n; i++) {
        max = (max < fabs(u[i]))? fabs(u[i]):max;
    }
    return max;
}
double test_init(const double x, const double y, const double wx, const double wy, const double ax, const double ay){
    double d2Xfact = (2.0*PI*wx)*(2.0*PI*wx);
    double d2Yfact = (2.0*PI*wy)*(2.0*PI*wy);
    return -(cos(2.0*PI*x*wx)*cos(2.0*PI*y*wy))/(ax*d2Xfact + ay*d2Yfact);
}

// Functions that are for more general purposes, not in use right now
//
//void apply_tri(double* restrict vOut, const double* restrict vIn,
//               const double* restrict diag,  const double* restrict udiag,
//               const double* restrict ldiag, const int n){
//    vOut[0] = diag[0] * vIn[0] + udiag[0] * vIn[1];
//    vOut[n-1] = diag[n-1] * vIn[n-1] + ldiag[n-2] * vIn[n-2];
//    for(int i = 1; i < n-1; i++){
//        vOut[i] = ldiag[i-1]*vIn[i-1] + diag[i]*vIn[i] + udiag[i]*vIn[i+1];
//    }
//}

//void solve_tri(double* restrict x, const double* restrict ldiag,
//               const double* restrict diag, const double* restrict udiag,
//               const int n) {
//    // Allocate scratch space.
//    double* cprime = (double*) malloc(sizeof(double) * n);
//    cprime[0] = udiag[0] / diag[0];
//    x[0] = x[0] / diag[0];
//    // loop from 1 to N - 1 inclusive
//    for (int i = 1; i < n; i++) {
//        double m = 1.0 / (diag[i] - ldiag[i-1] * cprime[i - 1]);
//        cprime[i] = udiag[i] * m;
//        x[i] = (x[i] - ldiag[i-1] * x[i - 1]) * m;
//    }
//    // loop from N - 2 to 0 inclusive, safely testing loop end condition
//    for (int i = n - 1; i-- > 0; )
//        x[i] = x[i] - cprime[i] * x[i + 1];
//    // free scratch space
//    free(cprime);
//}

//void transpose(double* restrict u_t, const double* restrict u, const int n){
//    __assume_aligned(u_t, 64);
//    __assume_aligned(u, 64);
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            u_t[j*n+i] = u[i*n+j];
//        }
//    }
//}
