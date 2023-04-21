#include <vector>
#include <map>
#include "kMSolution.h"
#include "fixedDouble.h"

using namespace std;

pair<kMSolution, int> localsearchRec(vector<int>* C, vector<int>* F, map<int, map<int, double>>* dFtoC, int k, double lam,
                          vector<vector<double>>* dFtoF, const vector<int>& startingS);
