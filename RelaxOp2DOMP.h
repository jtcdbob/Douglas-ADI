#ifndef _70A_HW8_RelaxOp2D_h
#define _70A_HW8_RelaxOp2D_h

#include "TriSolver.h"
#include "GridFun2D.h"
class RelaxOp2D{
public:
    double dt;
    TriSolver triSolX, triSolY;
    TriOperator triOpX, triOpY;
    GridFun2D uStar, Fstar, uTemp;
    std::vector<double> uX_Zero, uY_Zero;
    void initialize(double, double, double, const GridFun2D&);
    void apply(GridFun2D&, GridFun2D&);
};

void RelaxOp2D::initialize(double timestep, double alphaX, double alphaY, const GridFun2D& f){
    dt = timestep;
    uStar = f; uTemp = f;   // Instantiation of temporary variables
    Fstar = f; Fstar *= dt; // Initialize Fstar*dt

    long M = f.xPanel;      // Define system X size from the problem

    uX_Zero.resize(M+1);    // Square matrix assumed
    uY_Zero.resize(M+1);

    for (long i = 1; i < M+1; i++) {
        uX_Zero[i] = 0.0;
    }
    for (long i = 1; i < M+1; i++) {
        uY_Zero[i] = 0.0;
    }

    std::vector<double> loDiag(M);
    std::vector<double> upDiag(M);
    std::vector<double> diag(M+1);

    double hx = f.hx;       // Define grid size for x;
    double a = alphaX*dt/(2*hx*hx), b = -2.0 * a; // Interior Points
    for(long i = 0; i < M; i++){
        loDiag[i] = a;
        upDiag[i] = a;
        diag[i] = b;
    }
    diag[0] = 0.0; diag[M] = 0.0;   // B.C.s
    upDiag[0] = 0.0; loDiag[M-1] = 0.0;

    triOpX = TriOperator(M+1, loDiag, diag, upDiag); // Forward Solver in the x direction


    for(long i = 0; i < M ; i++){ // Change the elements for backsolver
        loDiag[i] = -a;
        upDiag[i] = -a;
        diag[i] = 1 - b;
    }
    diag[0] = 1.0; diag[M] = 1.0;   // B.C.s
    upDiag[0] = 0.0; loDiag[M-1] = 0.0;

    triSolX = TriSolver(M+1, loDiag, diag, upDiag); // Backsolver in the x direction


    long N = f.yPanel; //Define system Y size from the problem
    loDiag.resize(N);
    upDiag.resize(N);
    diag.resize(N+1);

    double hy = f.hy; // Define grid size for x;
    a = alphaY*dt/(hy*hy); b = -2.0 * a; // Interior Points
    for(long i = 0; i < N; i++){
        loDiag[i] = a;
        upDiag[i] = a;
        diag[i] = b;
    }
    diag[0] = 0.0; diag[N] = 0.0;   // B.C.s
    upDiag[0] = 0.0; loDiag[N-1] = 0.0;

    triOpY = TriOperator(N+1, loDiag, diag, upDiag);    // Forward Solver in the y direction

    a = alphaY*dt/(2*hy*hy); b = -2.0 * a;
    for(long i = 0; i < N; i++){
        loDiag[i] = - a;
        upDiag[i] = - a;
        diag[i] = 1 - b;
    }
    diag[0] = 1.0; diag[N] = 1.0; // B.C.s
    upDiag[0] = 0.0; loDiag[N-1] = 0.0;

    triSolY = TriSolver(N+1, loDiag, diag, upDiag);     // Backsolver in the y direction
}

void RelaxOp2D::apply(GridFun2D& uIn, GridFun2D& uOut){
    std::vector<double> uXTemp, uYTemp;
    std::vector<double> uXTempNew, uYTempNew;
    uXTemp.resize(uStar.values.getIndex2Size()); uXTempNew.resize(uStar.values.getIndex2Size());
    uYTemp.resize(uStar.values.getIndex1Size()); uYTempNew.resize(uStar.values.getIndex1Size());

    long i;
    uOut = uIn; // Apply identity operator
#pragma omp parallel for private(i) firstprivate(uXTemp,uXTempNew) schedule(static) default(shared)
    for (i = 1; i < uIn.yPanel; i++) { // Apply Forward X operator
        uIn.extractXslice(i, uXTemp);
        triOpX.apply(uXTemp, uXTempNew);
        uStar.insertXslice(i, uXTempNew);
    }
    uStar.insertXslice(0, uX_Zero); //Don't touch the boundary
    uStar.insertXslice(uIn.yPanel, uX_Zero);

#pragma omp parallel for private(i) firstprivate(uYTemp,uYTempNew) schedule(static) default(shared)
    for (i = 1; i < uIn.xPanel; i++) { // Apply Forward Y operator
        uIn.extractYslice(i, uYTemp);
        triOpY.apply(uYTemp, uYTempNew);
        uTemp.insertYslice(i, uYTempNew);
    }
    uTemp.insertXslice(0, uY_Zero); //Don't touch the boundary
    uTemp.insertXslice(uIn.xPanel, uY_Zero);

    uOut += uStar; // Complete arithmetic in the x direction
    uOut += uTemp;

    uOut -= Fstar;
#pragma omp parallel for private(i) firstprivate(uXTemp,uXTempNew) schedule(static) default(shared)
    for (i = 1; i < uIn.yPanel; i++) { // Apply Backsolver in the x direction
        uOut.extractXslice(i, uXTemp);
        triSolX.apply(uXTemp, uXTempNew);
        uStar.insertXslice(i, uXTempNew);
    }
    uOut.extractXslice(0, uXTemp); //Don't touch the boundary
    uStar.insertXslice(0, uXTemp);
    uOut.extractXslice(uIn.xPanel, uXTemp);
    uStar.insertXslice(uIn.xPanel, uXTemp);

    uTemp *= 0.5; // Divide by 2
    uStar -= uTemp; // Complete arithmetics in the y direction


#pragma omp parallel for private(i) firstprivate(uYTemp,uYTempNew) schedule(static) default(shared)
    for (i = 1; i < uIn.xPanel; i++) { // Apply Backsolver in the y direction
        uStar.extractYslice(i, uYTemp);
        triSolY.apply(uYTemp, uYTempNew);
        uOut.insertYslice(i, uYTempNew);
    }
    uStar.extractYslice(0, uYTemp); //Don't touch the boundary
    uOut.insertYslice(0, uYTemp);
    uStar.extractYslice(uIn.xPanel, uYTemp);
    uOut.insertYslice(uIn.xPanel, uYTemp);
}
////################# DEBUG #######################
//    std::cout << "######### After Y BackSolve ###########" << std::endl;
//    for (long j = 0; j < uIn.yPanel + 1; j++) {
//        for (long i = 0; i < uIn.xPanel + 1; i++) {
//            std::cout << uOut.values(j,i) << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "####################" << std::endl;

#endif
