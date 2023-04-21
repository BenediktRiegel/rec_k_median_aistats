#include "ReductionAlgo.h"
#include <iostream>
#include <vector>
#include "kMSolution.h"
#include <set>
//#include "CharLiKUFL.h"
#include "vectorEkUFL_LocalSearch.h"
#include <omp.h>
#include "fixedDouble.h"

using namespace std;

double calculate_RCcost(const vector<int>& S, vector<vector<double>>* dFtoF, double lam) {
    //fixedDouble cost = fixedDouble("0");
    double cost = 0;
    for (int i1 = 0; i1 < S.size(); ++i1) {
        for (int i2 = i1 + 1; i2 < S.size(); ++i2) {
            cost += (*dFtoF).at(S.at(i1)).at(S.at(i2));
        }
    }
    cost = cost * lam;
    return cost;
}

pair<kMSolution, pair<vector<int>, vector<int>>> recsolve(const vector<int>& sampled_C, vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoF,
                                       vector<vector<double>>* dAtoC, double lam, int k,
                                       const vector<vector<int>>& nearest_k, const vector<int>& nearest_f) {
    cout << "start recsolve" << endl;
    kMSolution S;
    S.service_cost = -1;
    S.other_cost = -1;
    S.Fused = -1;
    int indexM = 0;
    int maxM = (*F).size();
    vector<int> itr;
    vector<int> differentFSize;
    // guessing the median in F
    cout << "guessing median in F" << endl;
#pragma omp declare reduction(minKMS : kMSolution : omp_out = (omp_out.service_cost == -1 || (omp_in.cost() < omp_out.cost())) ? omp_in : omp_out) initializer (omp_priv=omp_orig)
#pragma omp parallel for reduction(minKMS:S)
    for (int p = 0; p < (*F).size(); ++p) {
        int m = (*F).at(p);
        indexM += 1;
        // print the number of the median guessed
        cout << "RedAlgo: " << indexM << " of " << maxM << endl;
        // time_t start = time(nullptr);
        vector<double> f;
        // calculate the opening cost of each i in F
        for (int i : *F) {
            f.push_back(((k - 1) * lam * (*dFtoF).at(i).at(m)));
        }
        vector<int> Fsubvec;
        // Get Fsubvec and newC
        {
            set<int> Fsubset;
            //cout << "Fsubset.insert k nearest for m" << endl;
            //cout << "Csets.insert k nearest for m" << endl;
            for (int i = 0; i < k; ++i) {
                Fsubset.insert(nearest_k.at(m).at(i));
            }
            //cout << "Fsubset.insert nearest f" << endl;
            //cout << "Csets.insert nearest f" << endl;
            for (int i = 0; i < sampled_C.size(); ++i) {
                Fsubset.insert(nearest_f.at(sampled_C.at(i)));
            }
            //cout << "Fset to Fvec" << endl;
            for (int el : Fsubset) {
                Fsubvec.push_back(el);
            }
        }
        cout << "Fsubvec.size(): " << Fsubvec.size() << endl;
        pair<kMSolution, int> SAndItr = localsearchkUFL(C, &Fsubvec, dAtoC, k, &f);

        kMSolution tempS = SAndItr.first;
#pragma omp critical
        {
            itr.push_back(SAndItr.second);
            differentFSize.push_back(Fsubvec.size());
        }
        tempS.other_cost = calculate_RCcost(tempS.solution, dFtoF, lam);
        tempS.Fused = Fsubvec.size();
        // print out the facilities of the solution and what the costs are for the used median
        cout << "m: " << m << " tempCost: " << stringValue(tempS.cost()) << endl;
        cout << "iterations to converge = " << SAndItr.second << endl;
        // check, if this solution is better than previous solutions
        if (S.service_cost == -1 || tempS.cost() < S.cost()) {
            S = tempS;
            //cout << "newCost: " << stringValue(S.cost()) << endl;
        }
    }
    /*
    for (int i = 1; i < S.solution.size(); ++i) {
        cout << ", " << S.solution[i];
    }
    cout << "]" << endl;
     */

    return pair<kMSolution, pair<vector<int>, vector<int>>> (S, pair<vector<int>, vector<int>> (itr, differentFSize));
}
