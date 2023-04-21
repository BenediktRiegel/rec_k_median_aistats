#include "evaluateTwitter.h"
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "kMSolution.h"
#include "loadsave.h"
#include "vectorRec_LocalSearch.h"
#include "fixedDouble.h"

using namespace std;


void saveResult(string path, kMSolution S, int num_itr, double duration, double lam, int k) {
	string filepath = path + "EvaluationResults/RecLocalSearch/k=" + stringValue(k) + "_lam=" + stringValue(lam) + ".txt";
	cout << path << endl;
	string result;
	result += "S=";
	string solString = "[" + to_string(S.solution.at(0));
	for (int i = 1; i < S.solution.size(); ++i) {
		solString += "," + to_string(S.solution.at(i));
	}
	solString += "]_serviceCost=" + stringValue(S.service_cost) + "_otherCost=" + stringValue(S.other_cost) + "_duration=" + to_string(duration) + "ms_itr=" + to_string(num_itr) + "\n";
	result += solString;

	cout << "save" << endl;
	fstream file;
	cout << filepath << endl;
	file.open(filepath, ios_base::app);
	file << result;
    cout << "done saving" << endl;
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
	string FtoCpath = path + "dAtoC.txt";
	vector<vector<double>> dAtoC = getDistanceVector(FtoCpath);
	
	cout << "Read C and F" << endl;
	vector<int> C = getIntVector(path + "C.txt");

	//Since AtoC is only FtoC
	string Fpath = path + "F.txt";
	vector<int> F = getIntVector(Fpath);
	cout << "F.size(): " << F.size() << endl;
	cout << "C.size(): " << C.size() << endl;
	cout << "dAtoC.size(): " << dAtoC.size() << endl;
	//cout << "dFtoC.at(0).size(): " << dAtoC.at(0).size() << endl;
    for (int run = 0; run < num_runs; ++run) {
        cout << "run " << run << endl;
        for (int k : ks) {
            for (double lam : lams) {
                cout << "k: " << k << ", lam: " << stringValue(lam) << endl;
                auto t_start = std::chrono::high_resolution_clock::now();
                pair<kMSolution, int> SandItr = localsearchRec(&C, &F, &dAtoC, k, lam, &dAtoC);
                auto t_end = std::chrono::high_resolution_clock::now();
                double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();
                kMSolution S = SandItr.first;
                int num_itr = SandItr.second;
                saveResult(path, S, num_itr, duration, lam, k);
            }
        }
        cout << "end run " << run << "\n\n\n" << endl;
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
