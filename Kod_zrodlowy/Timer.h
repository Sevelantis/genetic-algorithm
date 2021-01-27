//
// Created by oskro on 12/19/2020.
//

#ifndef CODE_TIMER_H
#define CODE_TIMER_H


#include<cstdint>
#include<chrono>
#include<iostream>

using namespace std;

struct Timer{
    chrono::high_resolution_clock::time_point start_time;
    chrono::high_resolution_clock::time_point stop_time;
    uint64_t duration = -1;

    Timer(){}

    void start()
    {
        start_time = chrono::high_resolution_clock::now();
    }

    void stop()
    {
        stop_time = chrono::high_resolution_clock::now();
        duration =  chrono::duration_cast<chrono::nanoseconds>(stop_time-start_time).count();
    }

    void info()
    {
        cout << "time: " << float(duration)/1000000000.0 << "s" <<endl;
    }

    uint64_t getTime()
    {
        stop_time = chrono::high_resolution_clock::now();
        duration =  chrono::duration_cast<chrono::nanoseconds>(stop_time-start_time).count();
        return duration;
    }
};


#endif //CODE_TIMER_H
