#ifndef BAALGO_LOADSAVE_H
#define BAALGO_LOADSAVE_H
#include <string>
#include <map>
#include <vector>
#include "fixedDouble.h"

using namespace std;

map<int, map<int, double>> getDistanceMapByFollowers(string path);
map<int, map<int, double>> getDistanceMap(string path, vector<int> M1, vector<int> M2);
vector<vector<double>> getDistanceVector(string path);
vector<int> getIntVector(string path);
vector<double> getDoubleVector(string path);
vector<vector<int>> getGraphVector(string path);

void reverseFile(string inpath, string outpath);
void cutFileAfter(string inpath, string outpath, int cutafter);
void FCtoAC(string inpath, string outpath);


#endif //BAALGO_LOADSAVE_H
