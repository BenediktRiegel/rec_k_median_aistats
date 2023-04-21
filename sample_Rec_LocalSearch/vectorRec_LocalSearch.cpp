#include "vectorRec_LocalSearch.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <numeric>
#include "fixedDouble.h"
#include <unordered_set>
#include "kMSolution.h"
#include <omp.h>
#include <stdlib.h>

using namespace std;


double calculate_clientcost(int j, const vector<int>& S, vector<int>* C, vector<vector<double>>* dFtoC) {
    //fixedDouble minCost = (*dFtoC).at(S.at(0)).at(j);
    //fixedDouble tempCost;
    double minCost = (*dFtoC).at(S.at(0)).at(j);
    double tempCost;
    for (int i = 0; i < S.size(); ++i) {
        tempCost = (*dFtoC).at(S.at(i)).at(j);
        if (tempCost < minCost) {
            minCost = tempCost;
        }
    }
    return minCost;
}

double calculate_servicecost(const vector<int>& S, vector<int>* C, vector<vector<double>>* dFtoC) {
    //fixedDouble cost("0");
    double cost = 0;
    for (int j : *C) {
        cost += calculate_clientcost(j, S, C, dFtoC);
    }
    return cost;
}


int maxOf(vector<int>* vec) {
    int max = (*vec).at(0);
    for (int el : (*vec)){
        if (max < el) {
            max = el;
        }
    }
    return max;
}


double calculate_single_reccost(int i, const vector<int>& S, vector<vector<double>>* dFtoF, double lam) {
    //fixedDouble cost("0");
    double cost = 0;
    for (int f : S) {
        cost += (*dFtoF).at(i).at(f);
    }
    return cost * lam;
}


double calculate_reccost(const vector<int>& S, vector<vector<double>>* dFtoF, double lam) {
    //fixedDouble cost("0");
    double cost = 0;
    for (int i1 = 0; i1 < S.size(); ++i1) {
        for (int i2 = i1; i2 < S.size(); ++i2) {
            cost += (*dFtoF).at(S.at(i1)).at(S.at(i2));
        }
    }
    return cost * lam;
}

double calculate_clientcost(int j, const vector<int>& S, vector<int>* C, vector<vector<double>>* dFtoC,
                            vector<int>* serving_f) {
    int serv_i = S.at(0);
    //fixedDouble minCost = (*dFtoC).at(S.at(0)).at(j);
    //fixedDouble tempCost;
    double minCost = (*dFtoC).at(S.at(0)).at(j);
    double tempCost;
    for (int i = 0; i < S.size(); ++i) {
        tempCost = (*dFtoC).at(S.at(i)).at(j);
        if (tempCost < minCost) {
            minCost = tempCost;
            serv_i = S.at(i);
        }
    }
    (*serving_f)[j] = serv_i;
    return minCost;
}

double calculate_servicecost(const vector<int>& S, vector<int>* C, vector<vector<double>>* dFtoC, vector<int>* serving_f) {
    //fixedDouble cost("0");
    double cost = 0;
    for (int j : *C) {
        cost += calculate_clientcost(j, S, C, dFtoC, serving_f);
    }
    return cost;
}


double calculate_cost(const vector<int>& S, vector<int>* C, vector<vector<double>>* dFtoC,
                           vector<vector<double>>* dFtoF, vector<int>* serving_f, double lam) {
    double cost = calculate_servicecost(S, C, dFtoC, serving_f) + calculate_reccost(S, dFtoF, lam);
    return cost;
}

pair<vector<int>, double> update_service_cost(int newi, int removedi, const kMSolution& S, vector<int>* C,
                                                      vector<vector<double>>* dFtoC, vector<int> serving_f) {

    //fixedDouble service_cost = S.service_cost;
    double service_cost = 0;
    for (int j = 0; j < serving_f.size(); ++j) {
        if ((*dFtoC).at(serving_f.at(j)).at(j) >= (*dFtoC).at(newi).at(j)) {
            //service_cost -= (*dFtoC).at(serving_f.at(j)).at(j);
            serving_f[j] = newi;
            service_cost += (*dFtoC).at(newi).at(j);
        }
        else if (serving_f.at(j) == removedi) {
            //service_cost -= (*dFtoC).at(removedi).at(j);
            service_cost += calculate_clientcost(j, S.solution, C, dFtoC, &serving_f);
        }
        else {
            service_cost += (*dFtoC).at(serving_f.at(j)).at(j);
        }
    }

    return pair<vector<int>, double> (serving_f, service_cost);
}

double update_rec_cost(double current_cost, int newi, int removedi, double lam,
                            vector<int> S, vector<vector<double>>* dFtoF, int del_index) {
    S.erase(S.begin() + del_index);
    return current_cost - calculate_single_reccost(removedi, S, dFtoF, lam) + calculate_single_reccost(newi, S, dFtoF, lam);
}


pair<kMSolution, vector<int>> getRandomS(vector<int>* F, int k, double lam, vector<int>* C, vector<vector<double>>* dFtoC,
                                         vector<vector<double>>* dFtoF, mt19937 random_engine, vector<int>* serving_f) {
    kMSolution S;

    cout << omp_get_thread_num() << ": getRandomS" << endl;
    //int maxC = maxOf(C) + 1;
    cout << omp_get_thread_num() << ": maxOf(C) + 1 = " << maxOf(C) + 1 << endl;
    cout << omp_get_thread_num() << ": (*C).size() = " << (*C).size() << endl;
    int maxC = maxOf(C) + 1;;
    //(*serving_f).reserve(maxC);
    cout << omp_get_thread_num() << ": \tserving_f" << endl;
    for (int j = 0; j < maxC; ++j){
        (*serving_f).push_back(0);
    }
    vector<int> solution;
    vector<int> indices;
    //indices.reserve((*F).size());
    cout << omp_get_thread_num() << ": \tindices" << endl;
    for (int i = 0; i < (*F).size(); ++i) {
        indices.push_back(i);
    }
    //iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), random_engine);
    cout << omp_get_thread_num() << ": \tsolution" << endl;
    for (int i = 0; i < k; ++i) {
        solution.push_back((*F).at(indices.at(i)));
    }
    S.solution = solution;
    S.service_cost = calculate_servicecost(S.solution, C, dFtoC, serving_f);
    S.other_cost = calculate_reccost(S.solution, dFtoF, lam);

    vector<int> notS;
    cout << omp_get_thread_num() << ": \tnotS" << endl;
    for (int i = k; i < (*F).size(); ++i) {
        notS.push_back((*F).at(indices.at(i)));
    }
    cout << omp_get_thread_num() << ": notS worked" << endl;

    return pair<kMSolution, vector<int>>(S, notS);
}


pair<kMSolution, int> localsearchRec(vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoC, int k, double lam,
                           vector<vector<double>>* dFtoF) {

    /*
    cout << "f = " << "[" << (*f).at(0);
    for (int i = 1; i < (*f).size(); ++i) {
        cout << ", " << (*f).at(i);
    }
    cout << "]\n" << endl;
    */

    if (k == (*F).size()) {
        kMSolution S;
        S.solution = (*F);
        S.service_cost = calculate_servicecost(S.solution, C, dFtoC);
        S.other_cost = calculate_reccost(S.solution, dFtoF, lam);
        return pair<kMSolution, int> (S, 1);
    }

    random_device rd;
    mt19937 random_engine(rd());

    vector<int> serving_f;
    serving_f.reserve((*C).size());
    auto SAndNotS = getRandomS(F, k, lam, C, dFtoC, dFtoF, random_engine, &serving_f);
    cout << omp_get_thread_num() << ": getRandomS worked" << endl;

    kMSolution S = SAndNotS.first;
    vector<int> notS = SAndNotS.second;
    int notSsize = notS.size();
    uniform_int_distribution<> rd_S(0, k - 1);
    uniform_int_distribution<> rd_notS(0, notSsize - 1);
    /*
    cout << "starting S: [" << S.solution.at(0);
    for (int i = 1; i < S.solution.size(); ++i) {
        cout << ", " << S.solution.at(i);
    }
    cout << "], cost: " << S.cost() << "\n" << endl;
    */
    cout << omp_get_thread_num() << ": get tempS" << endl;
    kMSolution tempS;
    tempS.solution = S.solution;
    cout << omp_get_thread_num() << ": tempS service cost" << endl;
    tempS.service_cost = S.service_cost;
    cout << omp_get_thread_num() << ": tempS other cost" << endl;
    tempS.other_cost = S.other_cost;
    cout << omp_get_thread_num() << ": tempS it worked" << endl;
    bool finished = false;

    int currenti;
    int newi;
    int num_itr = 0;
    while (!finished) {
        ++num_itr;
        finished = true;
        cout << omp_get_thread_num() << ": get offsets" << endl;
        int offset_S = rd_S(random_engine);
        //int offset_S = rand() % (k - 1);
        cout << omp_get_thread_num() << " received offset_S with: " << offset_S << endl;
        int offset_notS = rd_notS(random_engine);

        //int offset_notS = rand() % (notSsize - 1);
        cout << omp_get_thread_num() << " reced offset_notS with: " << offset_notS << endl;
        cout << omp_get_thread_num() << ": S.size()=" << tempS.solution.size() << ", notS.size()=" << notS.size() << ", total F=" << (*F).size() << endl;
        //shuffle(notS.begin(), notS.end(), random_engine);
        //shuffle(S.solution.begin(), S.solution.end(), random_engine);
        for (int Si = 0; Si < k; ++Si) {
            currenti = Si + offset_S;
            if (currenti >= k) {
                currenti -= k;
            }
            cout << omp_get_thread_num() << ": currenti: " << currenti << " smaller " << k << "?" << endl;
            for (int notSi = 0; notSi < notSsize; ++notSi) {
                newi = notSi + offset_notS;
                if (newi >= notSsize) {
                    newi -= notSsize;
                }
                cout << omp_get_thread_num() << ": newi: " << newi << " smaller " << notS.size() << "?" << endl;
                tempS.solution[currenti] = notS.at(newi);
                //newi, removedi, S, C, dFtoC, f_serves
                //cout << "within itr update_service_cost" << endl;
                auto f_servesAndsC = update_service_cost(tempS.solution.at(currenti), S.solution.at(currenti), tempS, C, dFtoC, serving_f);
                tempS.service_cost = f_servesAndsC.second;
                //update rec cost
                //tempS.other_cost += (calculate_single_reccost(tempS.solution.at(currenti), tempS.solution, dFtoF, lam) -
                //        calculate_single_reccost(S.solution.at(currenti), S.solution, dFtoF, lam));
                tempS.other_cost = calculate_reccost(tempS.solution, dFtoF, lam);
                /*
                cout << "swapping " << S.solution.at(currenti) << " with " << notS.at(newi)
                     << " results in a cost of: " << tempS.cost() << "\n" << endl;
                */
                if (tempS.cost() < S.cost()) {
                    if (tempS.cost() < 0) {
                        cout << "negative cost on iteration " << num_itr << endl;
                        cout << "tempS.cost() < 0, with " << to_string(tempS.cost()) << " = " << to_string(tempS.service_cost) << " + " << to_string(tempS.other_cost) << endl;
                    }
                    notS[newi] = S.solution.at(currenti);
                    S.solution[currenti] = tempS.solution.at(currenti);
                    S.other_cost = tempS.other_cost;
                    S.service_cost = tempS.service_cost;
                    finished = false;
                    serving_f = f_servesAndsC.first;
                    /*
                    cout << "new S: [" << S.solution.at(0);
                    for (int i = 1; i < S.solution.size(); ++i) {
                        cout << ", " << S.solution.at(i);
                    }
                    cout << "], cost: " << S.cost() << "\n" << endl;
                     */
                    break;
                } else {
                    tempS.solution[currenti] = S.solution.at(currenti);
                    tempS.service_cost = S.service_cost;
                    tempS.other_cost = S.other_cost;
                }
            }
            if (!finished) {
                break;
            }
        }
    }
    cout << "LocalSearch took " << num_itr << " iteration to converge" << endl;

    return pair<kMSolution, int> (S, num_itr);
}
