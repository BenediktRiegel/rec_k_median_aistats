#ifndef BAALGO_LOADSAVE_H
#define BAALGO_LOADSAVE_H
#include <string>
#include <map>
#include <vector>
#include "fixedDouble.h"

using namespace std;

vector<vector<double>> getDistanceVector(string path);
map<int, map<int, double>> getDistanceMap(string path);
vector<int> getIntVector(string path);
vector<double> getDoubleVector(string path);
vector<vector<int>> getGraphVector(string path);

void reverseFile(string inpath, string outpath);
void cutFileAfter(string inpath, string outpath, int cutafter);
void FCtoAC(string inpath, string outpath);

string join(vector<int> input, char insert);


#endif //BAALGO_LOADSAVE_H
