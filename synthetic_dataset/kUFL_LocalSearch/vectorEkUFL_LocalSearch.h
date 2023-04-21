#include <vector>
#include <map>
#include "kMSolution.h"
#include "fixedDouble.h"
using namespace std;

pair<kMSolution, int> localsearchkUFL(vector<int>* C, vector<int>* F, map<int, map<int, double>>* dFtoC, int k, vector<double>* f);
