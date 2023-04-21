#include "evaluateTwitter.h"
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "kMSolution.h"
#include "loadsave.h"
#include "Rec_sample.h"
#include "fixedDouble.h"
#include "SampleData.h"

using namespace std;


void saveResult(string path, kMSolution S, vector<int> itr, double duration, double lam, int k, int sample_amount, vector<int> differentFSize) {
    //Use this for D sampling
    //string filepath = path + "sampleEvaluationResults/RecLocalSearch/k=" + to_string(k) + "_lam=" + to_string(lam) + ".txt";
    //Use this for uniform sampling
    string filepath = path + "uniformSampleEvaluationResults/RecLocalSearch/k=" + to_string(k) + "_lam=" + to_string(lam) + ".txt";
	cout << path << endl;
	string result;
	result += "S=";
	string solString = "[" + to_string(S.solution.at(0));
	for (int i = 1; i < S.solution.size(); ++i) {
		solString += "," + to_string(S.solution.at(i));
	}
    solString += "]_serviceCost=" + to_string(S.service_cost) + "_otherCost=" + to_string(S.other_cost) + "_duration=" + to_string(duration) + "ms_itr=";
    solString += to_string(itr.at(0));
    for (int i = 1; i < itr.size(); ++i) {
        solString += "," + to_string(itr.at(i));
    }
    result += solString + "_sampleAmount=" + to_string(sample_amount) + "_FAmount=" + to_string(S.Fused) + "_differentFSize=" + to_string(differentFSize.at(0));
    for (int i = 1; i < differentFSize.size(); ++i) {
        result += "," + to_string(differentFSize.at(i));
    }
    result += "\n";
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
	string kpath = path + "sample_k.txt";
	vector<int> ks = getIntVector(kpath);
	string lampath = path + "sample_lam.txt";
	vector<double> lams = getDoubleVector(lampath);
	string FtoCpath = path + "dAtoC.txt";
	vector<vector<double>> dAtoC = getDistanceVector(FtoCpath);
    vector<vector<int>> nearest_k = getIntVecVec(path + "nearest_k.txt");
    vector<int> nearest_f = getIntVector(path + "nearest_f.txt");
    //vector<vector<int>> G = getGraphVector(path + "G.txt");
	
	cout << "Read C and F" << endl;
    vector<int> sample_amounts = getIntVector(path + "sample_amounts.txt");

	//Since AtoC is only FtoC
	string Fpath = path + "F.txt";
	vector<int> F = getIntVector(Fpath);
	vector<int> C = getIntVector(path + "C.txt");

	cout << "ks.size(): " << ks.size() << endl;
	cout << "lams.size()" << lams.size() << endl;
    cout << "nearest_f.size(): " << nearest_f.size() << endl;
    cout << "nearest_k.size(): " << nearest_k.size() << endl;
    cout << "F.size(): " << F.size() << endl;
    cout << "C.size(): " << C.size() << endl;
    cout << "test" << endl;
    cout << "sample_amounts.size(): " << sample_amounts.size() << endl;
    cout << "dAtoC.size(): " << dAtoC.size() << endl;
    for (int run = 0; run < num_runs; ++run) {
        cout << "run " << run << endl;
        for (int amount : sample_amounts) {
            cout << "num_to_be_sampled: " << amount << endl;
            auto t_uni_start = std::chrono::high_resolution_clock::now();
            // Use this line for D-Sampling
            //vector<int> sampled_C = fullMatrix_D_sampling(C, dAtoC, amount);
            // Use this line for uniform sampling
            vector<int> sampled_C = uniform_sampling(C, amount);
            auto t_uni_end = std::chrono::high_resolution_clock::now();
            double uni_duration = std::chrono::duration<double, std::milli>(t_uni_end - t_uni_start).count();
            cout << "seconds it took for uniform_sampling: " << uni_duration / 1000 << endl;

            for (int k : ks) {
                for (double lam : lams) {
                    cout << "k: " << k << ", lam: " << stringValue(lam) << ", sampled_C.size(): " << sampled_C.size() << endl;
                    auto t_start = std::chrono::high_resolution_clock::now();
                    pair<kMSolution, pair<vector<int>, vector<int>>> SandItr = rec_solve(sampled_C, &C, &F, &dAtoC, k, lam, &dAtoC, nearest_k, nearest_f);
                    //pair<kMSolution, int> SandItr = rec_solve(&C, &F, &dAtoC, k, lam, &dAtoC);
                    auto t_end = std::chrono::high_resolution_clock::now();
                    double duration = std::chrono::duration<double, std::milli>(t_end - t_start).count();

                    saveResult(path, SandItr.first, SandItr.second.first, duration, lam, k, amount, SandItr.second.second);
                }
            }
        }
        cout << "end run " << run << "\n\n\n" << endl;
    }
}
