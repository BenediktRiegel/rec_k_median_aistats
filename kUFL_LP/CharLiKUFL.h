#ifndef BAALGO_CHARLIKUFL_H
#define BAALGO_CHARLIKUFL_H
#include <vector>
#include <map>
#include "kMSolution.h"
using namespace std;

pair<kMSolution, bool> solvekUFL(vector<int>* inC, vector<int>* inF, map<int, map<int, double>>* indAtoC,
    int ink, map<int, double>* inf, const vector<vector<int>>& G);

kMSolution testkUFLNoLP(vector<int>* inC, vector<int>* inF, map<int, map<int, double>>* indAtoC,
    int ink, map<int, double>* inf, const vector<vector<int>>& G, map<int, map<int, fixedDouble>>* inx,
    map<int, fixedDouble>* iny);
#endif //BAALGO_CHARLIKUFL_H
