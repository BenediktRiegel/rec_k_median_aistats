#ifndef BAALGO_LOADSAVE_H
#define BAALGO_LOADSAVE_H
#include <string>
#include <map>
#include <vector>
#include "fixedDouble.h"

using namespace std;

map<int, map<int, double>> getDistanceMap(string path);
vector<vector<double>> getDistanceVector(string path);
vector<int> getIntVector(string path);
vector<double> getDoubleVector(string path);
vector<vector<int>> getGraphVector(string path);
vector<vector<int>> getIntVecVec(string path);

void reverseFile(string inpath, string outpath);
void cutFileAfter(string inpath, string outpath, int cutafter);
void FCtoAC(string inpath, string outpath);


#endif //BAALGO_LOADSAVE_H
