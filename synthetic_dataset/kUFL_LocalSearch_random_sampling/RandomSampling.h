#include "kMSolution.h"
#include <vector>
#include <map>
#include "fixedDouble.h"

using namespace std;

pair<kMSolution, vector<int>> recsolve(vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoF,
                                       map<int, map<int, double>>* dAtoC, double lam, int k, int sampleAmount);
