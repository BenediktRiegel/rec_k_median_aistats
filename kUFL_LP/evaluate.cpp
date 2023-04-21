#include "evaluate.h"
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "kMSolution.h"
#include "loadsave.h"
#include "ReductionAlgo.h"
#include "fixedDouble.h"

using namespace std;


void saveResult(string path, kMSolution S, bool lp_only, double duration, double lam, int k) {
	string filepath = path + "EvaluationResults/LP/k=" + to_string(k) + "_lam=" + to_string(lam) + ".txt";
	cout << "saving to " << filepath << endl;
	string result;
	result += "S=";
	string solString = "[" + to_string(S.solution.at(0));
	for (int i = 1; i < S.solution.size(); ++i) {
		solString += "," + to_string(S.solution.at(i));
	}
	solString += "]_serviceCost=" + to_string(S.service_cost) + "_otherCost=" + to_string(S.other_cost) + "_duration=" + to_string(duration) + "ms_lpOnly=" + to_string(lp_only);
	result += solString + "\n";

	fstream file;
	file.open(filepath, ios_base::app);
	file << result;
	file.close();
	cout << "done saving" << endl;
}


void evalFpartC(string path, int num_runs) {
	cout << "Starting evalFpartC" << endl;
	//Remember F <= C (subset equal) and F = [0, |F|)
	cout << "read data" << endl;
	map<int, map<int, double>> dAtoC = getDistanceMap(path + "dAtoC.txt");
	vector<vector<int>> G = getGraphVector(path + "G.txt");
	//string FtoFpath = path + "dFtoF.txt";
	//map<int, map<int, double>> dFtoF = getDistanceMap(FtoFpath);

	cout << "Create C and F" << endl;
	vector<int> C = getIntVector(path + "C.txt");

	vector<int> F = getIntVector(path + "F.txt");

	vector<int> ks = getIntVector(path + "k.txt");

	vector<double> lams = getDoubleVector(path + "lam.txt");

	/*
	for (int c : C){
		for (int c2 : C) {
			if (dAtoC.at(c).at(c2) < 0) {
				cout << "dAtoC[" << c << ", " << c2 << "] = " << dAtoC.at(c).at(c2) << endl;
			}
		}
	}*/

	cout << "start evalutaion" << endl;
	vector<vector<bool>> pure_LP;
	pure_LP.reserve(ks.size() * lams.size());
	for (int k : ks) {
	    vector<bool> temp;
	    temp.reserve(lams.size());
	    for (double lam : lams) {
	        temp.push_back(false);
	    }
	    pure_LP.push_back(temp);
	}
    for (int run = 0; run < num_runs; ++run) {
        for (int i = 0; i < ks.size(); ++i) {
            int k = ks.at(i);
            for (int j = 0; j < ks.size(); ++j) {
                double lam = lams.at(j);
                if (!pure_LP.at(i).at(j)) {
                    cout << "k: " << k << ", lam: " << stringValue(lam) << endl;
                    cout << "run " << run << endl;
                    auto t_start = std::chrono::high_resolution_clock::now();
                    pair<kMSolution, bool> S = (recsolve(&C, &F, &dAtoC, &dAtoC, lam, k, G));
                    auto t_end = std::chrono::high_resolution_clock::now();
                    double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();

                    saveResult(path, S.first, S.second, duration, lam, k);
                    if (S.second) {
                        cout << "Everything was LP_only" << endl;
                        pure_LP[i][j] = true;
                    }
                }
            }
            cout << "runs finished\n\n\n" << endl;
        }
    }
	cout << "Done" << endl;


}


void evalTwitter(string path, int num_runs) {
	cout << "Starting evalTwitter" << endl;
	//Remember F <= C (subset equal) and F = [0, |F|)
	cout << "read data" << endl;
	vector<int> ks = getIntVector(path + "k.txt");
	vector<double> lams = getDoubleVector(path + "lam.txt");
	string FtoCpath = path + "dFtoC.txt";
	map<int, map<int, double>> dAtoC = getDistanceMap(FtoCpath);
	string Gpath = path + "sortedTwitter.txt";
	vector<vector<int>> G = getGraphVector(Gpath);
	string FtoFpath = path + "dFtoF.txt";
	map<int, map<int, double>> dFtoF = getDistanceMap(FtoFpath);

	cout << "Create C and F" << endl;
	vector<int> C;
	C.reserve(dAtoC[0].size());
	for (int i = 0; i < dAtoC[0].size(); ++i) {
		C.push_back(i);
	}
	
	vector<int> F = getIntVector(path + "F.txt");

	cout << "start evalutaion" << endl;
	for (int run = 0; run < num_runs; ++run) {
		cout << "run " << run << endl;
		for (int k : ks) {
			for (double lam : lams) {
				cout << "k: " << k << ", lam: " << stringValue(lam) << endl;
				auto t_start = std::chrono::high_resolution_clock::now();
				pair<kMSolution, bool> S = (recsolve(&C, &F, &dFtoF, &dAtoC, lam, k, G));
				auto t_end = std::chrono::high_resolution_clock::now();
				double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	
				saveResult(path, S.first, S.second, duration, lam, k);
			}
		}
		cout << "end run " << run << "\n\n\n" << endl;
	}


}
