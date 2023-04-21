#include "evaluateTwitter.h"
#include <iostream>
#include "fixedDouble.h"
#include <omp.h>

using namespace std;

int main(int argc, char *argv[]){
    if (argc != 3) {
        cout << "Exactly 3 arguments are needed, first being the default argument and the second being the path and the third is the number of runs" << endl;
        exit(1);
    }
    string path = argv[1];
    int num_runs = stoi(argv[2]);

    evalTwitter(path, num_runs);
    cout << "end\n\n\n" << endl;
    return 0;
}

