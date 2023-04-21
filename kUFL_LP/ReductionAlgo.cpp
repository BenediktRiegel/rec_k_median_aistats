#include "ReductionAlgo.h"
#include <iostream>
#include <set>
#include <vector>
#include "kMSolution.h"
#include <time.h>
#include "CharLiKUFL.h"
#include <omp.h>
#include "fixedDouble.h"

using namespace std;

double calculate_RCcost(const vector<int>& S, map<int, map<int, double>>* dFtoF, fixedDouble lam) {
    //fixedDouble cost = 0;
    double cost = 0;
    for (int i1 = 0; i1 < S.size(); ++i1) {
        for (int i2 = i1; i2 < S.size(); ++i2) {
            cost += (*dFtoF).at(S.at(i1)).at(S.at(i2));
        }
    }
    cost = cost * lam.getDouble();
    return cost;
}


pair<kMSolution, bool> recsolve(vector<int>* C, vector<int>* F, map<int, map<int, double>>* dFtoF,
    map<int, map<int, double>>* dAtoC, double lam, int k, const vector<vector<int>>& G){

    kMSolution S;
    S.service_cost = -1;
    S.other_cost = -1;
    int indexM = 0;
    int maxM = (*F).size();
    bool lp_only = true;
    pair<kMSolution, bool> result = pair<kMSolution, bool>(S, lp_only);
    //cout << "max number of threads: " << omp_get_max_threads() << endl;
// guessing the median in F
//pragma omp declare reduction(minKMS : pair<kMSolution, bool> : omp_out = ((omp_out.first.service_cost == -1 || (omp_in.first.cost() < omp_out.first.cost())) ? omp_in.first : omp_out.first), (omp_in.second && omp_out.second)) initializer (omp_priv=omp_orig)
//pragma omp parallel for reduction(minKMS:result)
//#pragma omp parallel for
#pragma omp declare reduction(minKMS : pair<kMSolution, bool> : omp_out = pair<kMSolution, bool> (((omp_out.first.service_cost == -1 || (omp_in.first.cost() < omp_out.first.cost())) ? omp_in.first : omp_out.first), (omp_in.second && omp_out.second))) initializer (omp_priv=omp_orig)
#pragma omp parallel for reduction(minKMS:result)
    for (int p = 0; p < (*F).size(); ++p) {
        int m = (*F).at(p);
        //int m = 109;
    //for (int test = 0; test < 4; ++test) {
        //int m = (*F).at(test);
        //cout << "thread id: " << omp_get_thread_num() << endl;
        indexM += 1;
        // print the number of the median guessed
        std::cout << "RedAlgo: " << indexM << " of " << maxM << endl;
        // time_t start = time(nullptr);
        map<int, double> f;
        // calculate the opening cost of each i in F
        for (int i : *F) {
            f[i] = (k - 1.0) * lam * (*dFtoF)[i][m];
        }

        // print('(k - 1) * lam * dFtoF[i, m] = ' + str(k-1) + ' * ' + str(lam) + ' * ' + str(dFtoF[i, m]) + ' = ' + str((k - 1) * lam * dFtoF[i, m]))

        // use localSearch or Charikar/Li Algo, depending on which is wanted
        /*
        if LocalSearch:
            (tempS, tempCost) = LocalSearchkUFL.solve(C, F, k, dFtoC, f)
        else:*/

        pair<kMSolution, bool> SAndLP = solvekUFL(C, F, dAtoC, k, &f, G);
        kMSolution tempS = SAndLP.first;
        bool temp_lp = SAndLP.second;
        if (tempS.solution.size() != k) {
            cout << "For m = " << m << ", |S| = " << tempS.solution.size() << " != " << k << " = k" << endl;
            exit(-1);
        }
        tempS.other_cost = calculate_RCcost(tempS.solution, dFtoF, lam);
        // print out the facilities of the solution and what the costs are for the used median
        std::cout << "m: " << m << " tempCost: " << tempS.cost() << endl;
        // check, if this solution is better than previous solutions
        if (result.first.service_cost == -1 || tempS.cost() < result.first.cost()) {
            result.first = tempS;
            //std::cout <<"newCost: " << stringValue(result.first.cost()) << endl;
        }
        cout << "temp_lp: " << temp_lp << endl;
        result.second = result.second && temp_lp;
    }

    return result;
}
