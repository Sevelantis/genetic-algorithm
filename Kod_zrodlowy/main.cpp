#include <iostream>
#include <bitset>
#include <list>
#include <cstdlib>
#include <thread>

#include "Algorithm.cpp"
#include "SW.h"


using namespace std;

//----------------------MAIN-------------------------

void runGA(vector<vector<int>> matrix, int stop, double fCrossover, double fMutation, int fPopulation,
           Algorithm::CrossoverType crossoverType,
           double fCut, double fWeakProb)
{
    Algorithm *pAlg = new Algorithm(matrix, stop, fCrossover, fMutation, fPopulation, crossoverType);
    pAlg->setParams(fCut, fWeakProb);
    pAlg->GA();
}

int main(int argc, char *argv[])
{
    srand(time_t(NULL));

    Timer timer;
    Algorithm *pAlgorithm = nullptr;
    SW *pSW = nullptr;

    vector< vector<int> > matrix;
    int stop = 180;
    double fCrossover = 0.8;
    double fMutation = 0.01;
    int fPopulation = 800;
    Algorithm::CrossoverType crossoverType = Algorithm::PMX;

    double cool_rate = 1.2;
    bool running = true;
    char choice, ch;
    string filename = "tsp/ftv47.atsp";

    matrix = Algorithm::read_from_file(filename.c_str());

    //testing
//    pAlgorithm = new Algorithm(matrix, stop, fCrossover, fMutation, fPopulation, crossoverType);
//    pAlgorithm->GA();
//    pAlgorithm->print_best_path();
//    cout << pAlgorithm->time_found << "s";
//    return 0;

    //info data
    double err, errs[3][10];
    int ans_correct[3] = {1776, 2755, 2465};
    int ans_curr;
    int best_answers[3], answers[3][10];
    double avg_answers[3];
    int anscorr = -1;

    //turbo read
    while(running)
    {
        cout <<"\n--- MENU GLOWNE ---\n";
        cout <<"1.Wczytanie danych z pliku i wyswietlenie(aktualnie: "<< filename << ")\n"
               "2.Kryterum stopu[s]\n"
               "3.Wielkosc populacji poczatkowej\n"
               "4.Wspolczynnik mutacji\n"
               "5.Wspolczynnik krzyzowania\n"
               "6.Wybor metody krzyzowania\n"
               "7.Uruchom algorytm GA\n"
               "8.Uruchom algorytm SW\n"
               "0.Wyjscie\n";
        cout <<"--- MENU GLOWNE ---\n";
        cout << "User choice:";
        cin >> choice;
        switch(choice)
        {
            case '1': //wczytanie z pliku danych
                cout << "Podaj nazwe pliku(uzupelnic pole X: X.atsp): ";
                cin >> filename;
                filename.insert(0, "tsp/");
                filename.append(".atsp");
                matrix = Algorithm::read_from_file(filename.c_str());
                break;

            case '2': // kryterium stopu
                cout <<"obecnie: " << stop <<"[s], podaj nowe [s], (0,inf): " <<endl;
                cin >> stop;
                if(stop<=0)
                    cout <<"podaj poprawna wartosc";
                break;

            case '3': //wielkosc populacji poczatkowej
                cout << "obecnie: " << fPopulation << ", podaj nowy wspolczynnik wielkosci populacji: " << endl;
                cin >> fPopulation;
                break;

            case '4': // wspolczynnik mutacji
                cout << "obecnie: " << fMutation << ", podaj nowy wspolczynnik mutacji: " << endl;
                cin >> fMutation;
                break;

            case '5': // wspolczynnik krzyzowania
                cout << "obecnie: " << fCrossover << ", podaj nowy wspolczynnik krzyzowania: " << endl;
                cin >> fCrossover;
                break;

            case '6'://wybor metody krzyzowania
                cout << "obecnie: ";
                if (crossoverType == Algorithm::OX)
                    cout << "OX";
                else
                    cout <<"PMX";
                cout << "\n1 -> OX\n2 - > PMX\n";
                cin >> ch;

                if(ch == '1')
                    crossoverType = Algorithm::OX;
                else if(ch =='2')
                    crossoverType = Algorithm::PMX;
                else
                    cout << "Nie ma takiej opcji.";

                break;

            case '7'://ga
                pAlgorithm = new Algorithm(matrix, stop, fCrossover, fMutation, fPopulation, crossoverType);
                timer.start();
                pAlgorithm->GA();
                timer.stop();

                if(filename.compare("ftv47.atsp"))
                    anscorr = 1776;
                if(filename.compare("ftv170.atsp"))
                    anscorr = 2755;
                if(filename.compare("rbg403.atsp"))
                    anscorr = 2465;


                ans_curr = pAlgorithm->cost(pAlgorithm->bestSolution.path);
                err = (double) abs(ans_curr - anscorr) *100/ anscorr;
                cout << "---------------------------\n";
                if (crossoverType == Algorithm::OX)
                    cout << "\nKRZYZOWANIE: OX\n";
                else
                    cout <<"\nKRZYZOWANIE: PMX\n";
                cout << "WYNIK ALGORYTMU: "
                     << ans_curr << "\t"
                                   "\nblad wzgledny: " << err << "%"
                     << "\nCZAS WYKRYCIA: " << pAlgorithm->time_found / 1000000000.0 << "[s]\n";
                cout <<"CZAS CALKOWIY: ";
                timer.info();
                cout << "SCIEZKA: \n";
                pAlgorithm->print_best_path();
                cout << "---------------------------\n";
                break;

            case '8'://sa
                pSW = new SW(matrix, stop, cool_rate);
                timer.start();
                pSW->find_best_path();
                timer.stop();

                if(filename.compare("ftv47.atsp"))
                    anscorr = 1776;
                if(filename.compare("ftv170.atsp"))
                    anscorr = 2755;
                if(filename.compare("rbg403.atsp"))
                    anscorr = 2465;


                ans_curr = pAlgorithm->cost(pAlgorithm->bestSolution.path);
                err = (double) abs(ans_curr - anscorr) *100/ anscorr;
                cout << "---------------------------\n";
                cout << "WYNIK ALGORYTMU: "
                     << ans_curr << "\t"
                                   "\nblad wzgledny: " << err << "%"
                     << "\nCZAS WYKRYCIA: " << pAlgorithm->time_found / 1000000000.0 << "[s]\n";
                cout <<"CZAS CALKOWIY: ";
                timer.info();
                cout << "SCIEZKA: \n";
                pAlgorithm->print_best_path();
                cout << "---------------------------\n";
                break;

            case '9'://test all
                for (int i = 0; i < 3; ++i)
                {
                    if (i==0)
                    {
                        matrix = Algorithm::read_from_file("tsp/ftv47.atsp");
                        cout << "||ftv47:";
                        stop = 60;
                    } else if (i==1)
                    {
                        matrix = Algorithm::read_from_file("tsp/ftv170.atsp");
                        cout << "||ftv170:";
                        stop = 180;
                    } else if (i==2)
                    {
                        matrix = Algorithm::read_from_file("tsp/rbg403.atsp");
                        cout << "||rbg403:";
                        stop = 180;
                    }
                    for (int j = 0; j < 10; j++)
                    {
                        cout << " population: " << fPopulation << endl;
                        pAlgorithm = new Algorithm(matrix, stop, fCrossover, fMutation, fPopulation, crossoverType);

                        timer.start();
                        pAlgorithm->GA();
                        timer.stop();

                        ans_curr = pAlgorithm->cost(pAlgorithm->bestSolution.path);
                        answers[i][j] = ans_curr;
                        err = (double) abs(ans_curr - ans_correct[i]) *100/ ans_correct[i];
                        errs[i][j] = err;
                        if (best_answers[i] > ans_curr)
                            best_answers[i] = ans_curr;

                        cout << "\n------------" << i << ","<< j << "------------\n";
                        cout << "WYNIK ALGORYTMU: "
                             << ans_curr
                             << " WYNIK OCZEKIWANY: " << ans_correct[i] << "\t";
                        cout << "BLAD WZGL.: " << err << "%"
                             << "\nCZAS WYKRYCIA: " << pAlgorithm->time_found / 1000000000.0 << "[s]\n";
                        cout <<"CZAS CALKOWIY: ";
                        timer.info();
                        cout << "SCIEZKA: \n";
                        pAlgorithm->print_best_path();
                        cout << "---------------------------\n";
//                        break;
                    }
                }
                //
                fill(avg_answers, avg_answers+10, 0);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 10; ++j)
                        avg_answers[i] += answers[i][j];

                    avg_answers[i] /= 10.0;
                }
                cout << "\n------------POSUMOWANIE------------\n";
                cout << "ftv47 :: SREDNIA WYNIKOW: "<<avg_answers[0]<< ", OCZEKIWANO: "<<ans_correct[0]<< "\n";
                cout << "ftv170 :: SREDNIA WYNIKOW: "<<avg_answers[1]<<", OCZEKIWANO: "<<ans_correct[1]<< "\n";
                cout << "rbg403 :: SREDNIA WYNIKOW: "<<avg_answers[2]<< ", OCZEKIWANO: "<<ans_correct[2]<<"\n";
                cout << "------------POSUMOWANIE------------\n";
                break;

            case '0':
                running = false;
                break;

            default: break;
        }//switch choice
    }//while running

    //free memory
    delete pAlgorithm;
    pAlgorithm = nullptr;

    system("pause");
    return 0;
}
