#include "evaluateTwitter.h"
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <algorithm>
#include "kMSolution.h"
#include "loadsave.h"
#include "vectorRec_LocalSearch.h"
#include "fixedDouble.h"

using namespace std;


void saveResult(string path, kMSolution S, vector<int> num_itr, double duration, double lam, int k) {
	string filepath = path + "EvaluationResults/RecLocalSearchRandomSampling/k=" + stringValue(k) + "_lam=" + stringValue(lam) + ".txt";
	cout << path << endl;
	string result;
	result += "S=";
	string solString = "[" + to_string(S.solution.at(0));
	for (int i = 1; i < S.solution.size(); ++i) {
		solString += "," + to_string(S.solution.at(i));
	}
	solString += "]_serviceCost=" + stringValue(S.service_cost) + "_otherCost=" + stringValue(S.other_cost) + "_duration=" + to_string(duration) + "ms_itr=";
    solString += to_string(num_itr.at(0));
    for (int i = 1; i < num_itr.size(); ++i) {
        solString += "," + to_string(num_itr.at(i));
    }
	result += solString + "\n";

	cout << "save" << endl;
	fstream file;
	cout << filepath << endl;
	file.open(filepath, ios_base::app);
	file << result;
    cout << "done saving" << endl;
}


vector<int> getDisjunctMerge(vector<int> M1, const vector<int>& M2){
    vector<int> merge;
    copy(M1.begin(), M1.end(), back_inserter(merge));
    for (int el : M2){
        if(find(M1.begin(), M1.end(), el) == M1.end()){
            merge.push_back(el);
        }
    }
    return merge;
}


void evalTwitter(string path, int num_runs) {
	cout << "Starting evalTwitter" << endl;
	//Remember F <= C (subset equal) and F = [0, |F|)
	cout << path << endl;
	cout << "read data" << endl;
	string kpath = path + "k.txt";
	vector<int> ks = getIntVector(kpath);
	string lampath = path + "lam.txt";
	vector<double> lams = getDoubleVector(lampath);

    cout << "read C and F" << endl;
    string Cpath = path + "C.txt";
    vector<int> C = getIntVector(Cpath);

    //Since AtoC is only FtoC
    string Fpath = path + "F.txt";
    vector<int> F = getIntVector(Fpath);
    vector<int> A = getDisjunctMerge(F, C);
    string AtoCpath = path + "dAtoC.txt";
    //vector<vector<double>> dAtoC = getDistanceVector(FtoCpath);
    map<int, map<int, double>> dAtoC = getDistanceMap(AtoCpath, A, C);

    string FtoFpath = path + "dFtoF.txt";
    vector<vector<double>> dFtoF = getDistanceVector(FtoFpath);

	//cout << "dFtoC.at(0).size(): " << dAtoC.at(0).size() << endl;
    for (int k : ks) {
        for (double lam : lams) {
            cout << "k: " << k << ", lam: " << stringValue(lam) << endl;
            auto t_start = std::chrono::high_resolution_clock::now();
            kMSolution S;
            S.service_cost = -1;
            S.other_cost = -1;
            vector<int> num_itr;
#pragma omp declare reduction(minKMS : kMSolution : omp_out = (omp_out.service_cost == -1 || (omp_in.cost() < omp_out.cost())) ? omp_in : omp_out) initializer (omp_priv=omp_orig)
#pragma omp parallel for reduction(minKMS:S)
            for (int run = 0; run < num_runs; run++) {
                cout << "run " << run << endl;
                pair<kMSolution, int> SandItr = localsearchRec(&C, &F, &dAtoC, k, lam, &dFtoF);
                if (S.service_cost == -1 || SandItr.first.cost() < S.cost()) {
                    S = SandItr.first;
                }
                num_itr.push_back(SandItr.second);
                cout << "end run " << run << "\n\n\n" << endl;
            }
            auto t_end = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            saveResult(path, S, num_itr, duration, num_runs, k);
        }
    }
    /* LP rounding
    string Gpath = path + "sortedTwitter.txt";
    vector <vector<int>> G = getGraphVector(Gpath);
    cout << "start LP rounding evalutaion" << endl;
    for (int k : ks) {
        for (double lam : lams) {
            cout << "k: " << k << ", lam: " << lam << endl;
            auto t_start = std::chrono::high_resolution_clock::now();
            kMSolution S = (recsolve(&C, &F, &dFtoF, &dAtoC, lam, k, G));
            auto t_end = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();

            saveResult(path, S, duration, lam, k);
        }
    }
    */

}

