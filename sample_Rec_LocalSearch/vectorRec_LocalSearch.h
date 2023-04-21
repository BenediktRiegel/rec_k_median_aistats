#include <vector>
#include "kMSolution.h"
#include "fixedDouble.h"
using namespace std;

pair<kMSolution, int> localsearchRec(vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoC, int k, double lam,
                          vector<vector<double>>* dFtoF);
