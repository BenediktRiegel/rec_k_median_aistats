#include "RandomSampling.h"
#include <iostream>
#include <vector>
#include "kMSolution.h"
#include <random>
#include <algorithm>
#include "vectorEkUFL_LocalSearch.h"
#include "vectorRec_LocalSearch.h"

using namespace std;

double calculate_RCcost(const vector<int>& S, vector<vector<double>>* dFtoF, double lam) {
    double cost = 0;
    for (int i1 = 0; i1 < S.size(); ++i1) {
        for (int i2 = i1 + 1; i2 < S.size(); ++i2) {
            cost += (*dFtoF).at(S.at(i1)).at(S.at(i2));
        }
    }
    cost = cost * lam;
    return cost;
}

pair<kMSolution, vector<int>> recsolve(vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoF,
                    map<int, map<int, double>>* dAtoC, double lam, int k, int sampleAmount) {
    cout << "start recsolve" << endl;
    kMSolution S;
    S.service_cost = -1;
    S.other_cost = -1;
    int indexM = 0;
    int maxM = (*F).size();
    vector<int> itr;
    // guessing the median in F
    random_device rd;
    mt19937 random_engine(rd());

    vector<int> indices((*F).size());
    for (int i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    shuffle(indices.begin(), indices.end(), random_engine);
    cout << "F.size(): " << maxM << endl;
    uniform_int_distribution<> rd_S(0, maxM-1);
    cout << "randomly sampling m in F" << endl;
#pragma omp declare reduction(minKMS : kMSolution : omp_out = (omp_out.service_cost == -1 || (omp_in.cost() < omp_out.cost())) ? omp_in : omp_out) initializer (omp_priv=omp_orig)
#pragma omp parallel for reduction(minKMS:S)
    for (int p = 0; p < sampleAmount; ++p) {
//        int rd_int = rd_S(random_engine);
        int rd_int = indices[p];
        int m = (*F)[rd_int];
//#pragma omp atomic
//        indices[rd_int] += 1;
        cout << "rd_int = " << rd_int << " => m = " << m << endl;
        indexM += 1;
        // print the number of the median guessed
        cout << "JKRedAlgo: " << indexM << " of " << maxM << endl;
        // time_t start = time(nullptr);
        vector<double> f;
        // calculate the opening cost of each i in F
        for (int i : *F) {
            f.push_back(((k - 1) * lam * (*dFtoF).at(i).at(m)));
        }
        pair<kMSolution, int> SAndItr = localsearchkUFL(C, F, dAtoC, k, &f);
        kMSolution tempS = SAndItr.first;
        pair<kMSolution, int> SAndItr_ = localsearchRec(C, F, dAtoC, k, lam, dFtoF, tempS.solution);
        itr.push_back(SAndItr.second + SAndItr_.second);
        tempS = SAndItr_.first;
        // tempS.other_cost = calculate_RCcost(tempS.solution, dFtoF, lam);
        // print out the facilities of the solution and what the costs are for the used median
        cout << "m: " << m << " tempCost: " << stringValue(tempS.cost()) << endl;
        // check, if this solution is better than previous solutions
        if (S.service_cost == -1 || tempS.cost() < S.cost()) {
            S = tempS;
            //cout << "newCost: " << stringValue(S.cost()) << endl;
        }
    }
//    double prob = 0;
//    for (int i = 0; i < 5; i++) {
//        prob += indices[i];
//    }
//    prob = prob / sampleAmount;
//    cout << "empirical Pr(0 <= m <= 4) = " << prob << endl;
    /*
    for (int i = 1; i < S.solution.size(); ++i) {
        cout << ", " << S.solution[i];
    }
    cout << "]" << endl;
     */
//    S.service_cost = -1;
//    S.other_cost = -1;
//    S.solution.push_back(-1);
//    itr.push_back(-1);
    return pair<kMSolution, vector<int>> (S, itr);
}
