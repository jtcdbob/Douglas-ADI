#ifndef _ClockIt_
#define _ClockIt_
// This version of Class ClockIt can be used on either a Linux/Unix or MS system
// to time portions of code.
// The class depends upon calls to the gettimeofday(..) (Linux: in <sys/time.h> ) or the
// QueryPerformanceCounter (MS : in <Windows.h>) functions.
// Author: Chris Anderson. Version date: Mon 16 Feb 2015 09:03:58 AM PST
// Usage Snippet:
/*
 int main()
 {
 ClockIt clock1;              // declare ClockIt instance
 clock1.start();              // start the instance (reads the clock)
 // === Code to be timed === //  // carry out computation
 clock1.stop();               // stop the clock (reads the clock again)
 cout << "Time elapsed in Milli-seconds" << clock1.getMilliSecElapsedTime() << endl;
 }
 */

#ifndef _MSC_VER
//
// This version of Class ClockIt can be used on a Linux/Unix system to time portions of code. It
// depends upon calls to the gettimeofday(..) function in sys/time.h to read the system
// clock.
//
#include <sys/time.h>

class ClockIt{
public:
    ClockIt(){elapsedTime = 0;}
    ClockIt(const ClockIt& clock){
        elapsedTime = clock.elapsedTime;
        t1          = clock.t1;
        t2          = clock.t2;
    }

    // Reads the clock
    inline void start() {gettimeofday(&t1, 0);}
    // Reads the clock again, and computes the
    // the elapsed time.
    // Returns the time in milliseconds since start() called
    inline double stop(){
        gettimeofday(&t2, 0);// Compute elapsed time in milliseconds
        elapsedTime  = (t2.tv_sec  - t1.tv_sec)  * 1000.0;    // seconds to milliseconds
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;    // microseconds to milliseconds
        return elapsedTime;
    }
    double getMilliSecElapsedTime(){return elapsedTime;}
    double getMicroSecElapsedTime(){return elapsedTime*1000.0;}
    double getSecElapsedTime(){return elapsedTime/1000.0;}

    timeval t1, t2;
    double elapsedTime;
};
#else
// This version of Class ClockIt can be used on a MS system to time portions of code. It
// depends upon calls to the QPC timer function to read the system
// clock.
#include <Windows.h>
class ClockIt{
public:
    ClockIt(){elapsedTime = 0;}
    ClockIt(const ClockIt& clock){
        elapsedTime = clock.elapsedTime;
        t1 = clock.t1;
        t2 = clock.t2;
    }
    // Reads the clock
    inline void start(){
        QueryPerformanceFrequency(&Frequency);
        QueryPerformanceCounter(&StartingTime);
    }
    // Reads the clock again, and computes the
    // the elapsed time.
    // Returns the time in milliseconds since start() called
    inline double stop(){
        QueryPerformanceCounter(&EndingTime);
        ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
        // Compute elapsed time in milliseconds
        ElapsedMicroseconds.QuadPart *= 1000000;
        ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
        elapsedTime = (double)(ElapsedMicroseconds.QuadPart)/1000.0;
        return elapsedTime;
    }
    double getMilliSecElapsedTime(){return elapsedTime;}
    double getMicroSecElapsedTime(){return elapsedTime*1000.0;}
    double getSecElapsedTime(){return elapsedTime / 1000.0;}

    timeval t1, t2;
    double elapsedTime;

    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
};

#endif

#endif
