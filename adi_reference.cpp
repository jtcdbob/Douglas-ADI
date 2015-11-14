
class TestProblem2D{
public:
    TestProblem2D(double alphaX, double waveNumberX,double xMin, double xMax,
                  double alphaY, double waveNumberY,double yMin, double yMax){
        this->alphaX      = alphaX;
        this->alphaY      = alphaY;
        this->waveNumberX = waveNumberX;
        this->waveNumberY = waveNumberY;
        this->pi          = 3.1415926535897932;
    }
    double operator()(double x, double y){
        double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
        double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
        return -(cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY))/(alphaX*d2Xfact + alphaY*d2Yfact);
    }
    double normInf(){
        double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
        double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
        return 1.0/(alphaX*d2Xfact + alphaY*d2Yfact);
    }
    double pi;
    double alphaX;
    double alphaY;
    double waveNumberX;
    double waveNumberY;
};

// This function returns the maximal value of the residual of
// (alpha_x*Delta_x + alpha_y*Delta_y) u = f at interior points
double evaluateResidualNormInf(double alpha_x, double alpha_y, GridFun2D& u, GridFun2D& f){
    double residualNorm = 0.0;
    double residualValue;

    double hy = u.hy;
    double hx = u.hx;
    for(long i = 1; i <  u.xPanel; i++){
        for(long j = 1; j <  u.yPanel; j++){
            residualValue  = alpha_x*((u.values(i+1,j) - 2.0*u.values(i,j) + u.values(i-1,j))/(hx*hx))
            + alpha_y*((u.values(i,j-1) - 2.0*u.values(i,j) + u.values(i,j+1))/(hy*hy)) - f.values(i,j);
            if(fabs(residualValue) > residualNorm){residualNorm = fabs(residualValue);}
        }
    }
    return residualNorm;
}

int main(int argc, char** argv){

    TestProblem2D testSoln(alphaX, waveNumberX, xMin, xMax, alphaY, waveNumberY, yMin, yMax); // Instantiate the test problem solution
    GridFun2D f    (xPanel,xMin,xMax,yPanel,yMin,yMax); // Instantiate 2D grid functions
    GridFun2D uk   (xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D ukp1 (xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uLast(xPanel,xMin,xMax,yPanel,yMin,yMax);

    // Instantiate and initialize the relaxation operator for multi-scale scheme

    ClockIt clock1;
    clock1.start();
    const int maxSweepSize = (int)ceil(log(systemSizeM)/log(2)) + 1;
    double dt = dt_0;
    std::vector<RelaxOp2D>  relaxOpVec(maxSweepSize);
    for (long i = 0; i < maxSweepSize; i++) {
        relaxOpVec[i].initialize(dt, alphaX, alphaY, f);
        dt *= 2*2;
    }
    clock1.stop();

    // Relaxation loop
    double diffNorm = 2*tol;
    long iterMax  = 4000;
    long iter     = 0;
    double ukNorm;

    ClockIt clock2;
    clock2.start();

    long j;
    while((diffNorm > tol)&&(iter < iterMax)){
        uLast = uk;
        for (j = 0; j < maxSweepSize; j++) {
            relaxOpVec[j].apply(uk, ukp1);
            uk = ukp1;
        }
        for (j = 1; j < maxSweepSize + 1; j++) {
            relaxOpVec[maxSweepSize - j].apply(uk, ukp1);
            uk = ukp1;
        }
        ukNorm = ukp1.normInf();    // Check difference between iterates
        uLast -= ukp1;
        diffNorm = uLast.normInf()/ukNorm; // uLast.normInf();
        iter++; // Update iterate
    }
    clock2.stop();

    // Results for multiscale
    std::cout << "==== Multi-scale Timstep Relaxation Output XXXX ====" << std::endl;
    std::cout << "** Difference between iterates : " << diffNorm << std::endl;
    std::cout << "** Iteration Count for multi-scale time scheme        : " << iter*2*maxSweepSize     << std::endl << std::endl;

    // Timing results
    std::cout << "** Time it takes to initialize the multi-scale timestep solution is "<< clock1.getMilliSecElapsedTime() <<" ms\n";
    std::cout << "** Time it takes to solve the multi-scale timestep solution is "<< clock2.getMilliSecElapsedTime() <<" ms\n";
    std::cout << "** The total time is is "<< clock1.getMilliSecElapsedTime() + clock2.getMilliSecElapsedTime() <<" ms\n";

    return 0;
}
