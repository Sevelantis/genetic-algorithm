//
// Created by oskro on 12/19/2020.
//

#ifndef CODE_SW_H
#define CODE_SW_H

#define BOLTZMANN_CONSTANT 1.3806503e-23

#include<iostream>
#include<vector>
#include <climits>
#include "Timer.h"

using namespace std;

class SW
{

private:
    std::vector<std::vector<int>> cost_matrix;//input data
    int stop;
    double temp, temp_start;
    double cool_rate;
    double temp_min;
    Timer timer;

private:
    void swap_greedy(vector<int> &path);
    void _swap(vector<int> path, vector<int>& path_dest, int i, int j);
    void _inverse(vector<int> path, vector<int>& path_dest, int i, int j);
    void _insert(vector<int> path, vector<int>& path_dest, int i, int j);
    double initial_temperature();
    void randomize(vector<int>& path);

public:
    SW(vector<vector<int>> matrix,
       int stop,
       double cool_rate);
    ~SW();
    void find_best_path();
    void print_best_path();
    void print_path(vector<int> path);
    int cost(vector<int> path);

    vector<int> path_best;
    uint64_t time_found;
};


#endif //CODE_SW_H
