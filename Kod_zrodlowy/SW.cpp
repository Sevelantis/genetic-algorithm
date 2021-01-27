//
// Created by oskro on 12/19/2020.
//
#define TEMP_MIN .01
#define M_E		2.7182818284590452354

#include <numeric>
#include <cmath>
#include <cfloat>
#include <random>
#include <algorithm>
#include "SW.h"

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> dist(0, std::nextafter(1, DBL_MAX));
int count_of_vectors = 1;// powinno byc 1, proba ulepszenia funkcji wybierajacej najlepszego sasiada
int dim_g;
vector<vector<int>> paths_swap, paths_insert, paths_inverse;
vector<vector<int>> paths[3];

void SW::find_best_path()
{
    //deklaracje
    int dim = cost_matrix.size() - 1;
    vector<int> path_neigh = vector<int>(dim);//without 0
    iota(path_neigh.begin(), path_neigh.end(), 1);// path_neigh<=1,2,3,4...
    randomize(path_neigh);
    vector<int> path_curr = path_best = path_neigh;
    //prepare for greedy
    dim_g = dim;
    paths_swap.clear();
    paths_insert.clear();
    paths_inverse.clear();
    paths_swap.resize(count_of_vectors, vector<int>(dim_g));
    paths_insert.resize(count_of_vectors, vector<int>(dim_g));
    paths_inverse.resize(count_of_vectors, vector<int>(dim_g));
    paths[0] = paths_swap;
    paths[1] = paths_insert;
    paths[2] = paths_inverse;

    int delta;
    int L =  dim * (dim-1); // epoka
    int it = 0;
    bool heat = false;
//    if(dim < 100) L *= 1.25;
    if(dim>=100 and dim <=200) L /= 2;
    if(dim > 200)  L /= 6;
    //dziala przez okreslony czas[s] lub po osiagnieciu temp.min.
    timer.start();
    while(timer.getTime()/1000000000 < stop)
    {
        //for L times temperature is stable, doesnt change
        for (int i = 0; i < L; i++)
        {
            //search for solutions near current solution

            path_neigh = path_curr;

            //wyznaczanie losowo sasiedniego rozwiazania
            swap_greedy(path_neigh);

            //czy aktualne rozwiazanie jest najlepsze?
            if(cost(path_neigh)<cost(path_best))
            {
                time_found = timer.getTime();
                path_best=path_neigh;
                cout << cost(path_best) << endl;
                timer.info();
            }

            //czy rozwiazanie sasiednie jest lepsze od aktualnego?
            delta = cost(path_neigh) - cost(path_curr);
            if(delta < 0 || dist(mt) < (double)pow(M_E, (-delta/temp)))
                path_curr = path_neigh;
        }
        //after L times change temperature
        it++;//cauchy iterator
        if(!heat)
            //Geometric
//            temp *= cool_rate;
            //Cauchy
            temp = temp_start / (cool_rate * it);
        if(heat)
        {
            //geometric
//            temp = 16*cost(path_best) / sqrt(path_best.size());
            it = sqrt(cost(path_best))*2;
            heat = false;
        }
        if(temp<temp_min) heat = true;
    }
}


void SW::swap_greedy(vector<int> &path)
{
    int i,j,min=INT32_MAX;
    for (int k = 0; k < count_of_vectors; ++k)
    {
        i = dist(mt) * (dim_g-1) + 1;
        j = dist(mt) * (dim_g-1) + 1;

        _swap(path, paths[0][k], i, j);
        _insert(path, paths[1][k], i, j);
        _inverse(path, paths[2][k], i, j);

        for (int l = 0; l < 3; ++l)
        {
            if(min > cost(paths[l][k]))
            {
                min = cost(paths[l][k]);
                path = paths[l][k];
            }
        }
    }
}

int SW::cost(vector<int> path) {
    int length = cost_matrix[0][path[0]];
    for (int i = 0; i < path.size()-1; ++i)
    {
        length += cost_matrix[path[i]][path[i + 1]];
    }
    return length + cost_matrix[path[path.size()-1]][0];
}

void SW::print_path(vector<int> path)
{
    cout<<endl;
    for (std::vector<int>::const_iterator i = path.begin(); i != path.end(); ++i)
        std::cout << *i << ' ';
    cout<<endl;
}

void SW::print_best_path()
{
    cout << path_best.size()<<" = size \nbest path:\n0 - ";
    for (int i = 0; i < path_best.size(); ++i)
    {
        cout << path_best[i] << " - ";
    }
    cout << "0 ,\tcost: " << cost(path_best) << endl;
}


void SW::randomize(vector<int> &path)
{
    for (int i = path.size(); i > 0; --i)
        swap(path[i - 1], path[rand() % i]);
}

double SW::initial_temperature()
{
    double initial_temp = 0.;
    double tmp_diff;

    int dim = cost_matrix.size() - 1;
    vector<int> permutation_first(dim);
    vector<int> permutation_second(dim);
    iota(permutation_first.begin(),permutation_first.end(),1);
    permutation_second = permutation_first;

    for (int i = 0; i < 1000; i++)
    {
        randomize(permutation_first);
        randomize(permutation_second);
        tmp_diff = abs(cost(permutation_first) - cost(permutation_second));
        if(tmp_diff > initial_temp)
            initial_temp = tmp_diff;
    }
    return initial_temp;
}

void SW::_swap(vector<int> path, vector<int>& path_dest, int i, int j)
{
    int dim = path.size();
    copy(path.begin(),path.end(),path_dest.begin());
    if(i > dim || j > dim || i == j) return;
    path_dest[i] = path[j];
    path_dest[j] = path[i];
}

void SW::_inverse(vector<int> path, vector<int> &path_dest, int i, int j)
{
    int dim = path.size();
    if (i > dim || j > dim || i == j)
    {
        std::copy(path.begin(), path.begin() + dim, path_dest.begin());
        return;
    }
    if (i > j) { // Ustawienie w rosnacej kolejnosci
        int copy_v = i;
        i = j;
        j = copy_v;
    }
    std::copy(path.begin(), path.begin() + dim, path_dest.begin());
    std::reverse(path_dest.begin() + i, path_dest.begin() + j + 1);
}

void SW::_insert(vector<int> path, vector<int> &path_dest, int i, int j)
{
    int dim = path.size();
    if (i < j && i<dim && j < dim)
    {
        std::copy(path.begin(), path.begin() + i, path_dest.begin());
        //przesuniecie w lewo
        std::copy(path.begin()+i+1, path.begin() + j + 1, path_dest.begin() + i);
        path_dest[j] = path.begin()[i];
        std::copy(path.begin() + j + 1, path.begin() + dim, path_dest.begin() + j + 1 );
    }
    else if (i > j && i < dim && j < dim)
    {
        std::copy(path.begin(), path.begin() + j, path_dest.begin());
        //przesuniecie w prawo
        std::copy(path.begin() + j, path.begin() + i, path_dest.begin() + j + 1);
        path_dest[j] = path.begin()[i];
        std::copy(path.begin() + i + 1, path.begin() + dim, path_dest.begin() + i + 1);
    }
    else
        std::copy(path.begin(), path.begin() + dim, path_dest.begin());
}

SW::SW(vector<vector<int>> matrix,
       int stop,
       double cool_rate
)
{
    this->cost_matrix = std::move(matrix);
    this->stop = stop;
    this->cool_rate = cool_rate;
    this->temp_start = this->temp = initial_temperature();
    this->temp_min = TEMP_MIN;
}

SW::~SW()
{
    cost_matrix.clear();
    path_best.clear();
    cout <<"SW deleted\n";
}
