#include "evaluateTwitter.h"
#include <iostream>
#include "fixedDouble.h"
#include <omp.h>

using namespace std;

int main(int argc, char *argv[]){
    if (argc != 4) {
        cout << "Exactly 4 arguments are needed, first being the default argument and the second being the path and the third is the number of runs and the fourth the sample amount!" << endl;
        exit(1);
    }
    string path = argv[1];
    int num_runs = stoi(argv[2]);
    int sampleAmount = stoi(argv[3]);
    /*
    vector<int> ks = {2, 4, 8};
    vector<fixedDouble> lams;
    lams.push_back(fixedDouble("0"));
    fixedDouble currentlam = fixedDouble("1");
    for (int i = 1; i < 7; ++i) {
        lams.push_back(currentlam / fixedDouble("10"));
        currentlam = currentlam * fixedDouble("0.2");
    }
    */
    evalTwitter(path, num_runs, sampleAmount);
    cout << "end\n\n\n" << endl;
    return 0;
}

