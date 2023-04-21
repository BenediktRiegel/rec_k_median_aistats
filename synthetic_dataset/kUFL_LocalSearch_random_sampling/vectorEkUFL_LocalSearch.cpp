#include "vectorEkUFL_LocalSearch.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <numeric>
#include "fixedDouble.h"
#include <unordered_set>
#include "kMSolution.h"

using namespace std;


double calculate_openingcost(const vector<int>& S, vector<double>* f) {
    //fixedDouble cost = fixedDouble("0");
    double cost = 0;
    for (int i : S) {
        cost += (*f).at(i);
    }
    return cost;
}

double calculate_clientcost(int j, const vector<int>& S, vector<int>* C, map<int, map<int, double>>* dFtoC,
                            map<int, int>* serving_f) {
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

double calculate_servicecost(const vector<int>& S, vector<int>* C, map<int, map<int, double>>* dFtoC, map<int, int>* serving_f) {
    //fixedDouble cost = fixedDouble("0");
    double cost = 0;
    for (int j : *C) {
        cost += calculate_clientcost(j, S, C, dFtoC, serving_f);
    }
    return cost;
}


double calculate_cost(const vector<int>& S, vector<int>* C, map<int, map<int, double>>* dFtoC, vector<double>* f,
                      map<int, int>* serving_f) {
    double cost = calculate_servicecost(S, C, dFtoC, serving_f) + calculate_openingcost(S, f);
    return cost;
}

pair<map<int, int>, double> update_service_cost(int newi, int removedi, const kMSolution& S, vector<int>* C,
                                              map<int, map<int, double>>* dFtoC, map<int, int> serving_f) {

    //fixedDouble service_cost = S.service_cost;
    double service_cost = 0;
    for (int j : *C) {
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

    return pair<map<int, int>, double> (serving_f, service_cost);
}

double update_opening_cost(double current_cost, int newi, int removedi, vector<double>* f) {
    return current_cost - (*f).at(removedi) + (*f).at(newi);
}


pair<kMSolution, vector<int>> getRandomS(vector<int>* F, int k, vector<int>* C, map<int, map<int, double>>* dFtoC,
                                         vector<double>* f, mt19937 random_engine, map<int, int>* serving_f) {
    kMSolution S;

    // init serving f with all 0
    for (int j : (*C)){
        (*serving_f)[j] = 0;
    }
    // init solution
    vector<int> solution;
    vector<int> indices((*F).size());
    for (int i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), random_engine);
    for (int i = 0; i < k; ++i) {
        solution.push_back((*F)[indices[i]]);
    }
    S.solution = solution;
    //init costs
    S.service_cost = calculate_servicecost(S.solution, C, dFtoC, serving_f);
    S.other_cost = calculate_openingcost(S.solution, f);

    //init notS
    vector<int> notS;
    for (int i = k; i < (*F).size(); ++i) {
        notS.push_back((*F)[indices[i]]);
    }

    return pair<kMSolution, vector<int>>(S, notS);
}


pair<kMSolution, int> localsearchkUFL(vector<int>* C, vector<int>* F, map<int, map<int, double>>* dFtoC, int k, vector<double>* f) {

    /*
    cout << "f = " << "[" << (*f).at(0);
    for (int i = 1; i < (*f).size(); ++i) {
        cout << ", " << (*f).at(i);
    }
    cout << "]\n" << endl;
    */
    /*
    cout << "C.size = " << (*C).size() << endl;
    cout << "F.size = " << (*F).size() << endl;
    cout << "dFtoC.size = (" << dFtoC->size() << ", " << dFtoC->begin()->second.size() << ")" << endl;
    cout << "k = " << k << endl;
    cout << "f.size = " << (*f).size() << endl;
     */

    random_device rd;
    mt19937 random_engine(rd());

    map<int, int> serving_f;
    auto SAndNotS = getRandomS(F, k, C, dFtoC, f, random_engine, &serving_f);

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

    kMSolution tempS;
    tempS.solution = S.solution;
    tempS.service_cost = S.service_cost;
    tempS.other_cost = S.other_cost;
    bool finished = false;

    int currenti;
    int newi;
    int num_itr = 0;
    while (!finished) {
        ++num_itr;
        finished = true;
        int offset_S = rd_S(random_engine);
        int offset_notS = rd_notS(random_engine);
        //shuffle(notS.begin(), notS.end(), random_engine);
        //shuffle(S.solution.begin(), S.solution.end(), random_engine);
        for (int Si = 0; Si < k; ++Si) {
            currenti = Si + offset_S;
            if (currenti >= k) {
                currenti -= k;
            }
            for (int notSi = 0; notSi < notSsize; ++notSi) {
                newi = notSi + offset_notS;
                if (newi >= notSsize) {
                    newi -= notSsize;
                }
                tempS.solution[currenti] = notS.at(newi);
                //newi, removedi, S, C, dFtoC, f_serves
                //auto f_servesAndsC = update_service_cost(tempS.solution.at(currenti), S.solution.at(currenti), tempS, C, dFtoC, serving_f);
                //tempS.service_cost = f_servesAndsC.second;
                //tempS.other_cost = update_opening_cost(tempS.other_cost, tempS.solution.at(currenti), S.solution.at(currenti), f);
		        auto f_servesAndsC = update_service_cost(tempS.solution.at(currenti), S.solution.at(currenti), tempS, C, dFtoC, serving_f);
                tempS.service_cost = f_servesAndsC.second;
                tempS.other_cost = calculate_openingcost(tempS.solution, f);
                /*
                cout << "swapping " << S.solution.at(currenti) << " with " << notS.at(newi)
                     << " results in a cost of: " << tempS.cost() << "\n" << endl;
                */
                if (tempS.cost() < S.cost()) {
                    if (tempS.cost() < 0) {
                        cout << "negative cost on iteration " << num_itr << endl;
                        cout << "tempS.cost() < 0, with " << stringValue(tempS.cost()) << " = " << stringValue(tempS.service_cost) << " + " << stringValue(tempS.other_cost) << endl;
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
    //cout << "LocalSearch took " << num_itr << " iteration to converge" << endl;
    S.service_cost = calculate_servicecost(S.solution, C, dFtoC, &serving_f);
    return pair<kMSolution, int> (S, num_itr);
}
