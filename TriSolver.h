//
//  TriSolver.h
//  270A_HW2
//
//  Created by Bob Chen on 1/21/15.
//  Copyright (c) 2015 Bob. All rights reserved.
//

#ifndef _TriSolver_h
#define _TriSolver_h

class TriOperator{
    long systemSize;
    std::vector<double> loDiag;
    std::vector<double> diag;
    std::vector<double> upDiag;
public:
    TriOperator(){initialize();}; //Null Constructor
    TriOperator(const TriOperator&); //Copy constructor
    TriOperator(long, std::vector<double>&, std::vector<double>&, std::vector<double>&);
    virtual ~TriOperator(){}; //Default destructor
    void initialize(); //Null initializer, resizes the internal arrays to zero length
    void initialize(const TriOperator& ); //Copy initializer.
    void initialize(long, std::vector<double>&, std::vector<double>&, std::vector<double>&);
    void apply(std::vector<double>&, std::vector<double>&); //Apply vOut = T*vIn
};
TriOperator::TriOperator(const TriOperator& T){initialize(T);}; //Copy constructor
TriOperator::TriOperator(long systemSize, std::vector<double>& loDiag,
                         std::vector<double>& diag, std::vector<double>& upDiag){
    initialize(systemSize,loDiag,diag,upDiag);
}
void TriOperator::initialize(){
    systemSize = 0;
    loDiag.clear();
    diag.clear();
    upDiag.clear();
};
void TriOperator::initialize(const TriOperator& T){
    systemSize = T.systemSize;
    loDiag     = T.loDiag;
    diag       = T.diag;
    upDiag     = T.upDiag;
};
void TriOperator::initialize(long systemSize, std::vector<double>& loDiag,
                             std::vector<double>& diag, std::vector<double>& upDiag){
    this->systemSize = systemSize;
    this->loDiag     = loDiag;
    this->diag       = diag;
    this->upDiag     = upDiag;
}
void TriOperator::apply(std::vector<double>& vIn, std::vector<double>& vOut){
    if(systemSize == 1){ //Trivial case
        vOut[0] = diag[0]*vIn[0];
        return;
    }
    vOut[0] = diag[0]*vIn[0] + upDiag[0]*vIn[1];
    long i;
    for(i = 1; i < systemSize-1; i++){
        vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i] + upDiag[i]*vIn[i+1];
    }
    i = systemSize-1;
    vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i];
}

class TriSolver{
    long systemSize;
    std::vector<double> loDiag;
    std::vector<double> diag;
    std::vector<double> upDiag;
public:
    TriSolver(){initialize();};     // Null Constructor
    TriSolver(const TriSolver&);    // Copy Constructor
    TriSolver(long, std::vector<double>&, std::vector<double>&, std::vector<double>&);
    virtual ~TriSolver(){};         // Default destructor
    void initialize();              // Null initializer, resizes the internal arrays to zero length
    void initialize(const TriSolver&); // Copy initializer
    void initialize(long, std::vector<double>&, std::vector<double>& , std::vector<double>&);
    void apply(std::vector<double>, std::vector<double>&); // Back solve for the system

};

TriSolver::TriSolver(const TriSolver& T){initialize(T);}; // Copy constructor
TriSolver::TriSolver(long systemSize, std::vector<double>& loDiag,
                     std::vector<double>& diag, std::vector<double>& upDiag){
    initialize(systemSize,loDiag,diag,upDiag);
}

void TriSolver::initialize(){
    systemSize = 0;
    loDiag.clear();
    diag.clear();
    upDiag.clear();
};
void TriSolver::initialize(const TriSolver& T){
    systemSize = T.systemSize;
    loDiag     = T.loDiag;
    diag       = T.diag;
    upDiag     = T.upDiag;
};
void TriSolver::initialize(long systemSize, std::vector<double>& loDiag,
                           std::vector<double>& diag, std::vector<double>& upDiag){
    this->systemSize = systemSize;
    this->loDiag     = loDiag;
    this->diag       = diag;
    this->upDiag     = upDiag;
}
void TriSolver::apply(std::vector<double> f, std::vector<double>& u){
    std::vector<double> upDiag_temp = upDiag;
    std::vector<double> diag_temp = diag;
    if(systemSize == 1){ //Trivial case
        u[0] = f[0]/diag_temp[0];
        return;
    }
    //LU factorization.
    long i;
    upDiag_temp[0] = upDiag_temp[0]/diag_temp[0];
    for (i = 1; i < systemSize; i++){
        diag_temp[i] = diag_temp[i] - loDiag[i-1]*upDiag_temp[i-1];
        upDiag_temp[i] = upDiag_temp[i]/diag_temp[i];
    }
    // Forward Solve Lw = f
    f[0] = f[0]/diag_temp[0];
    for (i = 1; i < systemSize; i++) {
        f[i] = (f[i] - loDiag[i-1]*f[i-1])/diag_temp[i];
    }
    // Backward Solve Uu = w
    u[systemSize - 1] = f[systemSize-1];
    for (i = systemSize - 2; i > -1; i--) {
        u[i] = f[i] - upDiag_temp[i]*u[i+1];
    }
}
#endif
