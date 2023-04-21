#include "kMSolution.h"
#include <vector>
#include <map>
#include "fixedDouble.h"

using namespace std;

pair<kMSolution, pair<vector<int>, vector<int>>> recsolve(const vector<int>& sampled_C, vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoF,
                                       vector<vector<double>>* dAtoC, double lam, int k,
                                       const vector<vector<int>>& nearest_k, const vector<int>& nearest_f);

