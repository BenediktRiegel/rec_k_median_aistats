#include <vector>
#include "kMSolution.h"

using namespace std;

pair<kMSolution, pair<vector<int>, vector<int>>> rec_solve(const vector<int>& sampled_C, vector<int>* C, vector<int>* F, vector<vector<double>>* dFtoC, int k, double lam,
                                        vector<vector<double>>* dFtoF, const vector<vector<int>>& nearest_k, const vector<int>& nearest_f);
