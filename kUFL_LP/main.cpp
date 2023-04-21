#include <iostream>
#include <vector>
#include "evaluate.h"
#include <string>
//#include "testChariLi.h"
//#include "fixedDouble.h"

using namespace std;


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Exactly 3 parameters are required one is the default and the other is a path to the data and the number of runs" << endl;
        exit(1);
    }
    /*
    string path = argv[1];
    string GOrNoG = argv[2];
    vector<int> ks = { 2, 4, 8 };
    vector<double> lams;
    lams.push_back(0.0);
    double currentlam = 0.02;
    for (int i = 0; i < 7; ++i) {
        lams.push_back(currentlam / 10);
        currentlam *= 0.2;
    }
    if (GOrNoG == "1") {
        evalTwitter(path, ks, lams);
    }
    else {
        evalFpartC(path, ks, lams);
    }
    */
    string path = argv[1];
    int num_runs = stoi(argv[2]);
    evalFpartC(path, num_runs);
    //for (int i = 0; i < 8; ++i) {
    //test_ChariLi(path, 5, fixedDouble(0, 0), 0);
    //}
    //fixedDouble sol = fixedDouble(5 - 1) * fixedDouble("1.0") * 1;
    //cout << stringValue(sol) << endl;
    cout << "end\n\n\n" << endl;
    
    return 0;
}
