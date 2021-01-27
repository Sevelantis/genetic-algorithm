//
// Created by oskro on 12/19/2020.
//

#ifndef CODE_ALGORITHM_CPP
#define CODE_ALGORITHM_CPP

#define SWAP 1
#define INSERT 2
#define INVERSE 3

#include<iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <fstream>
#include <cfloat>
#include <random>
#include <iomanip>
#include <string>

#include "Timer.h"
#include <numeric>

#define POP_SIZE_MIN 800
#define POP_SIZE_MID 8000
#define POP_SIZE_MAX 20000

using namespace std;

/*
 * this class contains genetical algorithm- GA
 */

class Algorithm
{
public:
    enum CrossoverType{
        PMX, OX
    };

private:

    std::vector<std::vector<int>> cost_matrix;//input data
    int stop;
    Timer timer;
    int size;

    double fMutation;
    double fCrossover;
    int fPopulation;

    CrossoverType crossoverType;

    inline double rnd()
    {
        static random_device rd;
        static std::mt19937 mt(rd());
        static std::uniform_real_distribution<double> dist(0, std::nextafter(1, DBL_MAX));
        return dist(mt);
    }

private:

    void _swap(vector<int>& path_dest ,int i, int j)
    {
        if(i >size || j > size || i == j) return; // w razie w
        int tmp = path_dest[i];
        path_dest[i] = path_dest[j];
        path_dest[j] = tmp;
    }

    void _inverse(vector<int>& path_dest ,int i, int j)
    {
        if (i > size || j > size || i == j) // w razie w
            return;

        if (i > j)
        { // Ustawienie w rosnacej kolejnosci
            int copy_v = i;
            i = j;
            j = copy_v;
        }
        std::reverse(path_dest.begin() + i, path_dest.begin() + j + 1);
    }

    void _insert(vector<int>& pathDest, int i, int j)
    {
        if (i < j && i<size && j < size)
        {
            //przesuniecie w lewo
            int tmp = pathDest[i];
            std::copy(pathDest.begin() + i + 1, pathDest.begin() + j + 1, pathDest.begin() + i);
            pathDest[j] = tmp;
            std::copy(pathDest.begin() + j + 1, pathDest.end(), pathDest.begin() + j + 1 );
        }
        else if (i > j && i < size && j < size) // prawdopodobnie nigdy sie nie wywola
        {
            //przesuniecie w prawo
            int tmp = pathDest[i];
            std::copy(pathDest.begin() + j, pathDest.begin() + i, pathDest.begin() + j + 1);
            pathDest[j] = tmp;
            std::copy(pathDest.begin() + i + 1, pathDest.end(), pathDest.begin() + i + 1);
        }
    }

    void randomize(vector<int>& path)
    {
        for (int i = path.size(); i > 0; --i)
            swap(path[i - 1], path[rand() % i]);
    }

    struct Single // osobnik
    {
        Single() = default;
        int cost{};
        vector<int> path;
    };

    struct SingleComparator : public std::binary_function<Single, Single, bool>
    {
        bool operator()(const Single& s1, const Single& s2) const
        {
            return s1.cost < s2.cost; // Porównywarka osobnikow - mniejszy koszt na gorze
        }
    };

    struct Pair
    {
        Single s1, s2;
    };

    struct Item
    {
        Item() = default;
        int r{};
        int s{};
        int index{};

        bool has(int &candidate, int &childNr)
        {
            return (childNr == 1) ? (candidate == r) : (candidate == s);
        }

        int map(int &candidate)
        {
            if(candidate == r)
                return s;
            if(candidate == s)
                return r;
        }
    };

    void initPopulation(vector<Single> &newPopulation)
    {
        vector<int> path = bestSolution.path;
        Single single;
        for (int i = 0; i < fPopulation; ++i)
        {
            randomize(path);

            single.path = path;
            single.cost = cost(path);
            newPopulation.push_back(single);
        }
    }

    void mutate(vector<Single> &population, bool mutate = false)
    {
        // generator pseudolosowy

        for(auto &single : population)
        {
            double mutationProbability = rnd(); //double (0,1)
            if(mutationProbability <= fMutation || mutate)
            {
                // sekcja dopasowania z przedzialu 1 - 'size' - losowanie
                int left = rnd() * (size - 2) + 1;
                int right = rnd() * (size - 1) + 1;
                while (left >= right) right = rnd() * (size - 1) + 1;

                const int MODE = rnd() * 3 + 1;

                if (MODE == INSERT)
                    _insert(single.path, left, right);
                else if (MODE == SWAP)
                    _swap(single.path, left, right);
                else if (MODE == INVERSE)
                    _inverse(single.path, left, right);
            }
        }

    }

    int avgPopulationCost(vector<Single> &population)
    {
        int sum = 0;
        for(auto &single : population)
            sum += single.cost;
        return sum / fPopulation;
    }

    void crossover(vector<Single> &parentPopulation)
    {
        if(crossoverType == OX)
            crossOX(parentPopulation);

        else if (crossoverType == PMX)
            crossPMX(parentPopulation);
    }

    void crossOX(vector<Single> &parentPopulation)
    {
        // stworzenie pustych miast juz niedostepnych dla krzyzowanych osobnikow
        vector<bool> tabuS1(size, false);
        vector<bool> tabuS2(size, false);

        vector<Single>childPopulation;
        Pair pairParent;
        Pair pairChild;
        // rand osobnikow krzyzowac dopoki nie powstanie cala populacja
        while(childPopulation.size() != fPopulation)//
        {
            int ii = rnd() * parentPopulation.size();
            int jj = rnd() * parentPopulation.size();

            pairParent.s1.path = pairChild.s1.path = parentPopulation[ii].path;
            pairParent.s2.path = pairChild.s2.path = parentPopulation[jj].path;


            double crossoverProbability = rnd(); // zwraca double z przedzialu (0,1)
            if(crossoverProbability <= fCrossover)
            {
                // sekcja dopasowania z przedzialu 1 - 'size' - losowanie
                int left = rnd() * (size - 1) + 1;  //CHECK size - 2
                int right = rnd() * (size - 1) + 1;

                while (left >= right)
                {
                    left = rnd() * (size - 2) + 1;  //CHECK
                    right = rnd() * (size - 1) + 1;
                }

                Item item;
                //krzyzowanie wstepne - sekcja dopasowania
                int j;
                for (j = left; j <= right; ++j)
                {
                    // zdefiniowanie odwzorowania
                    item.s = pairParent.s1.path[j];  //upper R
                    item.r = pairParent.s2.path[j];  //lower S

                    //krzyzowanie w zakresie <left,right> ---|left,right|---
                    pairChild.s1.path[j] = item.r;
                    pairChild.s2.path[j] = item.s;

                    // dodanie do listy juz niedostepnych miast dla krzyzowanych osobnikow
                    tabuS1[item.r - 1] = true;
                    tabuS2[item.s - 1] = true;
                }

                // utworzenie kandydata
                int candidate;

                int tabuOffset1 = 0;
                int tabuOffset2 = 0;
                // prawa strona zakresu
                for (j = right + 1; j < size; ++j)
                {
                    putCandidateOX(tabuS1,pairParent.s1, pairChild.s1, j,tabuOffset1);
                    putCandidateOX(tabuS2,pairParent.s2, pairChild.s2, j,tabuOffset2);
                }

                // krzyzowanie po lewej stronie zakresu - analogicznie
                for (j = 0; j < left; ++j)
                {
                    putCandidateOX(tabuS1,pairParent.s1, pairChild.s1, j,tabuOffset1);
                    putCandidateOX(tabuS2,pairParent.s2, pairChild.s2, j,tabuOffset2);
                }

                // czyszczenie przed nastepna petla7
                std::fill(tabuS1.begin(), tabuS1.end(), false);
                std::fill(tabuS2.begin(), tabuS2.end(), false);
            }// end if

            // dodanie potomstwa do nowej populacji
            pairChild.s1.cost = cost(pairChild.s1.path);
            pairChild.s2.cost = cost(pairChild.s2.path);
            childPopulation.push_back(pairChild.s1);
            childPopulation.push_back(pairChild.s2);
        }//end while

        parentPopulation = childPopulation;
    }// end crossOX

    inline void putCandidateOX(vector<bool> &tabu, Single &parent, Single &child, int j, int &tabuOffset)
    {
        int candidate = parent.path[j];
        int candidatePos = j + tabuOffset;
        bool first  = true;
        while(tabu[candidate - 1])
        {
            first = false;
            candidatePos++;
            if(candidatePos >= size)
                candidatePos %= size;
            candidate = parent.path[candidatePos];
        }
        if(!first)
            tabuOffset++;
        //wpisanie kandydata
        child.path[j] = candidate;
        tabu[candidate - 1] = true;
    }

    void crossPMX(vector<Single> &parentPopulation)
    {
        // stworzenie pustych miast juz niedostepnych dla krzyzowanych osobnikow
        vector<bool> tabuSingle1(size,false);
        vector<bool> tabuSingle2(size,false);

        vector<Single>childPopulation;
        Pair p;
        // randomowych osobnikow krzyzowac dopoki nie powstanie cala populacja
        while(childPopulation.size() != fPopulation)//
        {
            // wylosowac indeksy osobnikow do skrzyzowania

            int ii = rnd() * parentPopulation.size();
            int jj = rnd() * parentPopulation.size();

            p.s1.path = parentPopulation[ii].path;
            p.s2.path = parentPopulation[jj].path;

            double crossoverProbability = rnd(); // zwraca double z przedzialu (0,1)
            if(crossoverProbability <= fCrossover)
            {
                // sekcja dopasowania z przedzialu 1 - 'size' - losowanie
                int left = rnd() * (size - 1) + 1;  //CHECK size - 2
                int right = rnd() * (size - 1) + 1;

                while (left >= right)
                {
                    left = rnd() * (size - 2) + 1;  //CHECK
                    right = rnd() * (size - 1) + 1;
                }

                // stworzenie pustej mapy odwzorowan
                vector<Item> map;
                Item item;

                //krzyzowanie wstepne - sekcja dopasowania
                for (int j = left; j <= right; ++j)
                {
                    // zdefiniowanie odwzorowania
                    item.s = p.s1.path[j];
                    item.r = p.s2.path[j];

                    //krzyzowanie w zakresie <left,right> ---|left,right|---
                    p.s1.path[j] = item.r;
                    p.s2.path[j] = item.s;

                    // dodanie elementu do mapy odwzorowan
                    map.push_back(item);

                    // dodanie do listy juz niedostepnych miast dla krzyzowanych osobnikow
                    tabuSingle1[item.r - 1] = true;
                    tabuSingle2[item.s - 1] = true;
                }

                // krzyzowanie po prawej stronie zakresu
                for (int j = right + 1; j < size; ++j)
                {
                    // dziecko 1
                    putCandidatePMX(map, tabuSingle1, p.s1.path[j], 1);
                    // dziecko 2 - analogicznie
                    putCandidatePMX(map, tabuSingle2, p.s2.path[j], 2);
                }

                // krzyzowanie po lewej stronie zakresu - analogicznie
                for (int j = 0; j < left; ++j)
                {
                    // dziecko 1
                    putCandidatePMX(map, tabuSingle1, p.s1.path[j], 1);
                    // dziecko 2
                    putCandidatePMX(map, tabuSingle2, p.s2.path[j], 2);
                }

                // czyszczenie przed nastepna petla7
                std::fill(tabuSingle1.begin(), tabuSingle1.end(), false);
                std::fill(tabuSingle2.begin(), tabuSingle2.end(), false);
            }// end if

            // dodanie potomstwa do nowej populacji
            p.s1.cost = cost(p.s1.path);
            p.s2.cost = cost(p.s2.path);
            childPopulation.push_back(p.s1);
            childPopulation.push_back(p.s2);

        }//end while
        parentPopulation = childPopulation;
    }//end crossPMX

    inline void putCandidatePMX(vector<Item> &map, vector<bool> &tabu, int &candidate, int childNr)
    {
        // dopoki kandydat powtarza sie
        while(tabu[candidate - 1])
        {
            // znalezc nowego kandydata na podstawie mapy odwzorowan
            for (auto & el : map)    //mozna zoptymalizowac
                if(el.has(candidate, childNr))
                {
                    candidate = el.map(candidate);
                    break;
                }
        }
        tabu[candidate-1] = true;
    }

    void selectPopulation(vector<Single> &population)
    {
        // SELEKCJA osobnikow - najlepsi sa na gorze
        sort(population.begin(), population.end(), SingleComparator());
        vector<Single> selectedPopulation;

        //selekcja elitarna
        int avgCost = avgPopulationCost(population);
        Single single;
        for(int i = 0; i < fPopulation * fCut; i++)
        {
            single = population[i];
            if(single.cost < avgCost) // pozbywam sie gorszych osobnikow (wiekszy koszt)
                selectedPopulation.push_back(single);
            else if(rnd() <= fWeakProb) // akceptuje gorszych z pewnym prawdopodobienstwem
                selectedPopulation.push_back(single);
        }

        for(int i = fPopulation * fCut; i < fPopulation ; i++)
            if(rnd() <= fWeakProb*.2) // akceptuje gorszych z pewnym mniejszym prawdopodobienstwem
                selectedPopulation.push_back(population[i]);

        // kiedy zaden nie zostanie zaakceptowany
        if(selectedPopulation.empty())
        {
            mutate(lastPopulation, true);
            selectedPopulation = lastPopulation;
        }

        population = selectedPopulation;
    }

private:
    double fWeakProb = .175, fWeakProbCpy;
    double fCut = .175, fCutCpy;

    vector<Single> lastPopulation;

    void setupParams()              // wspolczynniki wyznaczano dla krzyzowania PMX
    {
        if(size < 80)                       //ftv47
        {
//            fPopulation = POP_SIZE_MIN;

            if(fPopulation == POP_SIZE_MIN) // 30s,stuck
            {
                fWeakProb = .335;   //.025
                fCut = .15;         //.98
            }else if(fPopulation == POP_SIZE_MID) // 15s, stuck
            {
                fWeakProb = .0195;   //.125
                fCut = .75;         //.75
            }else if(fPopulation == POP_SIZE_MAX) //15s,stuck
            {
                fWeakProb = .185;
                fCut = .75;
            }
        }
        else if(size >= 80 && size <= 200) // ftv170
        {
//            fPopulation = POP_SIZE_MIN;

            if(fPopulation == POP_SIZE_MIN) //60s, stuck
            {
                fWeakProb = .024;
                fCut = .98;
            }else if(fPopulation == POP_SIZE_MID) // 160s, stuck
            {
                fWeakProb = .125;   //.125
                fCut = .75;         //.75
            }else if(fPopulation == POP_SIZE_MAX)  // 160s, stuck
            {
                fWeakProb = .175;
                fCut = .7;
            }
        }
        else if(size > 200)                 //rbg403
        {
//            fPopulation = POP_SIZE_MIN;

            if(fPopulation == POP_SIZE_MIN) // po 60s stuck
            {
                fWeakProb = .024;   //.025
                fCut = .98;         //.98
            }else if(fPopulation == POP_SIZE_MID)// 120s moglby dalej
            {
                fWeakProb = .125;
                fCut = .25;
            }else if(fPopulation == POP_SIZE_MAX) // 120s moglby dalej
            {
                fWeakProb = .175;
                fCut = .175;
            }
        }
    }

public:
    int it = 0;
    int itWithNoUpdate = 0;
    Single bestSolution;
    uint64_t time_found;

    void GA()
    {
        //setup parameters (probabilities boundaries) to maximize profit
        setupParams();

        bestSolution.path = vector<int>(size);
        iota(bestSolution.path.begin(), bestSolution.path.end(), 1);
        randomize(bestSolution.path);
        bestSolution.cost = cost(bestSolution.path);

        vector<Single> population;
        initPopulation(population);
        lastPopulation = population;
        fCutCpy = fCut;
        fWeakProbCpy = fWeakProb;

        timer.start();
        while(timer.getTime()/1000000000 < stop)
        {
            it++;
            selectPopulation(population);
            crossover(population);
            mutate(population);
            // sprawdzenie najlepszego rozwiazania
            Single singleCandidate = bestSingle(population);
            if(singleCandidate.cost < bestSolution.cost)
            {
                this->bestSolution = singleCandidate;
                time_found = timer.getTime();
                cout << bestSolution.cost <<" at " << it << " at ";
                timer.info();

                // when stuck
                itWithNoUpdate=0;
                fCut = fCutCpy;
                fWeakProb = fWeakProbCpy;
                lastPopulation = population;
            }else
                itWithNoUpdate++; // count loops when stuck

            if(itWithNoUpdate > 125) // execute on stuck
            {
                fCut = rnd();
                fWeakProb = rnd();
                cout << "zdarzenie krytyczne proba poprawy\n";
                timer.info();
                itWithNoUpdate = 0;
                mutate(population, true);
            }
        }
    }

    Single bestSingle(vector<Single> &population)
    {
        Single s = *min_element(population.begin(), population.end(), SingleComparator());
        return s;
    }
    int fPopulationCpy;
    Algorithm(vector<vector<int>> matrix, int stop,
              double fCrossover, double fMutation, int fPopulation, CrossoverType crossoverType)
    {
        this->cost_matrix = std::move(matrix);
        this->stop = stop;
        this->fCrossover = fCrossover;
        this->fMutation = fMutation;
        this->fPopulation = fPopulationCpy = fPopulation;
        this->crossoverType = crossoverType;

        this->size = cost_matrix.size() - 1;

        this->timer = Timer();
    }

    void setParams(double fCut, double fWeakProb)
    {
        this->fCut = fCut;
        this->fWeakProb = fWeakProb;
    }

    void print_best_path()
    {
        cout << bestSolution.path.size() << " = size \nbest path:\n0 - ";
        for (int i = 0; i < bestSolution.path.size(); ++i)
        {
            cout << bestSolution.path[i] << " - ";
        }
        cout << "0 ,\tcost: " << cost(bestSolution.path) << endl;
    }

    void print_path(vector<int> path)
    {
        cout<<endl;
        for (std::vector<int>::const_iterator i = path.begin(); i != path.end(); ++i)
            std::cout << *i << ' ';
        cout<<endl;
    }

    int cost(const vector<int>& path)
    {
        int length = cost_matrix[0][path[0]];
        for (int i = 0; i < path.size() - 1; ++i)
            length += cost_matrix[path[i]][path[i + 1]];

        return length + cost_matrix[path[path.size() - 1]][0];
    }

    int val(const vector<int>& pA, const vector<int>& pB)
    {
        int len = cost_matrix[0][pB[0]] - cost_matrix[0][pA[0]];
        for (int i = 0; i < pA.size() - 1; ++i)
            len += cost_matrix[pB[i]][pB[i + 1]] - cost_matrix[pA[i]][pA[i + 1]];

        len += cost_matrix[pB[pB.size() - 1]][0] - cost_matrix[pA[pA.size() - 1]][0];
        return len;
    }

    static std::vector<std::vector<int>> read_from_file(const char* datafile)
    {
        ifstream f;
        vector<vector<int>> matrix;

        f.open(datafile);

        if(f.is_open())
        {
            int val, row = 0, col = 0, dimension;

            //preprocess dimension...
            string line;
            while(getline(f, line) && line.find("EDGE_WEIGHT_SECTION") == string::npos)
                if(line.find("DIMENSION: ") != string::npos)
                    dimension = stoi(line.substr(10,15));

//        f >> dimension;
            //initialize matrix vector
            matrix.resize(dimension);

            //input values from file
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    f >> val;
                    matrix[i].push_back(val);
                }
            }
        }else
            cout << "\n\tNie ma takiego pliku.\n";

        f.close();

        return matrix;
    }

    void print_matrix(std::vector<std::vector<int>> matrix)
    {
        cout << endl <<matrix.size() << endl;
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix[i].size(); j++) {
                cout<< setw(2) << right << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }


    ~Algorithm()
    {
        cost_matrix.clear();
        bestSolution.path.clear();
        cout <<"Algorithm deleted\n";
    }
};


#endif //CODE_ALGORITHM_CPP
