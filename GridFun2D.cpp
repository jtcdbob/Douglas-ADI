#include "GridFun2D.h"

GridFun2D::GridFun2D(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax){
    initialize(xPanel, xMin, xMax, yPanel, yMin, yMax);
}
GridFun2D::GridFun2D(const GridFun2D& G){
    initialize(G);
}
void GridFun2D::initialize(){
    hx = 0.0; xMin = 0.0; xMax = 0.0; xPanel = 0;
    hy = 0.0; yMin = 0.0; yMax = 0.0; yPanel = 0;
    values.initialize();
}
void GridFun2D::initialize(const GridFun2D& G){
    hx = G.hx; xMin = G.xMin; xMax = G.xMax; xPanel = G.xPanel;
    hy = G.hy; yMin = G.yMin; yMax = G.yMax; yPanel = G.yPanel;
    values = G.values;
}
void GridFun2D::initialize(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax){
    this->xPanel = xPanel; this->xMin = xMin; this->xMax = xMax;
    this->yPanel = yPanel; this->yMin = yMin; this->yMax = yMax;
    hx = (xMax-xMin)/xPanel; hy = (yMax-yMin)/yPanel;
    values.initialize(xPanel+1,yPanel+1);
}
void GridFun2D::operator=(const GridFun2D& G){
    this->xPanel = G.xPanel; this->xMin = G.xMin; this->xMax = G.xMax; this->hx = G.hx;
    this->yPanel = G.yPanel; this->yMin = G.yMin; this->yMax = G.yMax; this->hy = G.hy;
    this->values = G.values;
}
void GridFun2D::operator+=(const GridFun2D& G){
    long i,j;
    for(i = 0; i < values.getIndex1Size(); i++){
        for(j = 0; j < values.getIndex2Size(); j++){
            values(i,j) += G.values(i,j);
        }
    }
}
void GridFun2D::operator-=(const GridFun2D& G){
    for(long i = 0; i < values.getIndex1Size(); i++){
        for(long j = 0; j < values.getIndex2Size(); j++){
            values(i,j) -= G.values(i,j);
        }
    }
}
void GridFun2D::operator*=(const double alpha){
    long i,j;
    for( i = 0; i < values.getIndex1Size(); i++){
        for( j = 0; j < values.getIndex2Size(); j++){
            values(i,j) *= alpha;
        }
    }
}
void GridFun2D::operator/=(const double alpha){
    for(long i = 0; i < values.getIndex1Size(); i++){
        for(long j = 0; j < values.getIndex2Size(); j++){
            values(i,j) /= alpha;
        }
    }
}
void GridFun2D::setToValue(const double d){
    for(long i = 0; i < values.getIndex1Size(); i++){
        for(long j = 0; j < values.getIndex2Size(); j++){
            values(i,j) = d;
        }
    }
}
double GridFun2D::normInf(){
    double maximum = fabs(values(0,0));
    for(long i = 0; i < values.getIndex1Size(); i++){
        for(long j = 0; j < values.getIndex2Size(); j++){
            if(fabs(values(i,j)) > maximum){
                maximum = fabs(values(i,j));
            }
        }
    }
    return maximum;
}
std::ostream& operator<<(std::ostream& outStream, const GridFun2D V){
    for (long k = 0; k < V.xPanel + 1; k++) {
        for (long l = 0; l < V.yPanel + 1; l++) {
            outStream << V.values(l,k) << "\t" << std::endl;
        }
    }
    return outStream;
}

void GridFun2D::extractXslice(const long index, std::vector<double>& uOut){
    for (long j = 0; j < values.getIndex2Size(); j++) uOut[j] = values(j,index);
}
void GridFun2D::insertXslice(const long index, std::vector<double>& uIn){
    for (long j = 0; j < values.getIndex2Size(); j++) values(j,index) = uIn[j];
}
void GridFun2D::extractYslice(const long index, std::vector<double>& uOut){
    for (long j = 0; j < values.getIndex1Size(); j++) uOut[j] = values(index,j);
}
void GridFun2D::insertYslice(const long index, std::vector<double>& uIn){
    for (long j = 0; j < values.getIndex1Size(); j++) values(index,j) = uIn[j];
}
