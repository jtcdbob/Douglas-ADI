#ifndef _DoubleArray2D_
#define _DoubleArray2D_

#include <iostream>
#include <iomanip>

#include <cassert> // Toggle off asserts if not in debug mode
#ifndef  _DEBUG
#define _NDEBUG
#endif

namespace Math270A{
    //####################################################################
    //                    Math270A_DoubleArray2D.h
    //####################################################################
    /**
     Provides a very "light weight" two dimensional array structure
     with initialization capabilities and optional bounds checking.

     <pre>
     The beginning index is 0               : (C convention)
     Data for the array is stored by ROWS   : (C convention)
     Access using (*,*), e.g. A(i,j) for (i,j)th element.: (NOT C convention)
     </pre>

     The copy constructor creates a duplicate instance.<p>

     Bounds checking is toggled on by specifying the compiler pre-processor
     define _DEBUG
     <p>
     Created for use in Math 270A<p>

     @author Chris Anderson (C) UCLA
     @version  Jan. 26, 2015
     */
    //#####################################################################
    //
    //#####################################################################
    class DoubleArray2D{
        public:
        //###################################################################
        //                 Constructors/Initialization
        //###################################################################
        DoubleArray2D(){dataPtr = 0; index1Size = 0; index2Size = 0;};

        DoubleArray2D(long m, long n){
            dataPtr = 0; index1Size = 0; index2Size = 0;
            initialize(m,n);
        };

        DoubleArray2D(const DoubleArray2D& D){
            dataPtr = 0; index1Size = 0; index2Size = 0;
            initialize(D);
        };
        virtual ~DoubleArray2D(){
            if(dataPtr !=  0) delete [] dataPtr;
        }

        void initialize(){
            if(dataPtr != 0) delete [] dataPtr;
            dataPtr = 0;
            index1Size = 0;
            index2Size = 0;
        }
        void initialize(const DoubleArray2D& D){
            if(D.dataPtr == 0){ // Initialization of null object
                initialize();
                return;
            }
            if((index1Size !=  D.index1Size)||(index2Size != D.index2Size)){
                if(dataPtr != 0) delete [] dataPtr;
                index1Size    = D.index1Size;
                index2Size    = D.index2Size;
                dataPtr       = new double[index1Size*index2Size];
            }
            for(long i = 0; i < index1Size*index2Size; i++){
                dataPtr[i] = D.dataPtr[i];
            }
        }
        void initialize(long m, long n){
            if((index1Size != m)||(index2Size != n)){
                if(dataPtr != 0) delete [] dataPtr;
                index1Size = m;
                index2Size = n;
                dataPtr = new double[index1Size*index2Size];
            }
            for(long i = 0; i  < index1Size*index2Size; i++){
                dataPtr[i] = 0.0;
            }
        }
        void operator=(const DoubleArray2D& D){
            initialize(D);
        }
        //
        //###################################################################
        //      Element Access with bounds checking toggled
        //      using _DEBUG
        //###################################################################
        //
#ifdef _DEBUG
        double&  operator()(long i1, long i2){
            assert(boundsCheck(i1, 0, index1Size-1,1));
            assert(boundsCheck(i2, 0, index2Size-1,2));
            return *(dataPtr +  i2 + i1*index2Size);
        };
        const double&  operator()(long i1, long i2) const{
            assert(boundsCheck(i1, 0, index1Size-1,1));
            assert(boundsCheck(i2, 0, index2Size-1,2));
            return *(dataPtr +   i2  + i1*index2Size);
        };
#else
        inline double&  operator()(long i1, long i2){
            return *(dataPtr +  i2 + i1*index2Size);
        };
        inline const double&  operator()(long i1, long i2) const{
            return *(dataPtr +  i2  + i1*index2Size);
        };
#endif
        //
        //###################################################################
        //                Array Structure Access Functions
        //###################################################################
        //
        /// Returns the number of indices in the 1st coordinate direction
        long getIndex1Size()  const {return index1Size;}
        /// Returns the number of indices in the 2nd coordinate direction
        long getIndex2Size()  const {return index2Size;}
        /// Returns the total number of data values
        long getDataSize()    const {return index1Size*index2Size;}
        /// Returns a pointer to the double array containing the data values
        double* getDataPointer() const {return dataPtr;}
        //  Input/Output
        //
        //  Prints out values as as if they were in the first Cartesian
        //  quadrant --- not in matrix indexing.
        //
        //
        friend std::ostream&  operator <<(std::ostream& outStream, const DoubleArray2D& A){
            for(long j = A.index2Size-1; j >=  0; j--){
                for(long i = 0; i <=  A.index1Size-1; i++){
                    outStream <<  "\t" << A(i,j) << " ";
}
                outStream << std::endl;
            }
            return outStream;
        }
        //
        //###################################################################
        //                      Class Data Members
        //###################################################################
        //
        protected :
        double*      dataPtr;     // data pointer
        long      index1Size;     // coordinate 1 size
        long      index2Size;     // coordinate 2 size
        //
        //###################################################################
        //                      Bounds Checking
        //###################################################################
        //
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const{
            if((i < begin)||(i  > end)){
                cerr << "Array index " << coordinate << " out of bounds " << endl;
                cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << endl;
                return false;
            }
            return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif
    };
}       /* Namespace Math270A */
#endif  /* _DoubleArray2D_    */
