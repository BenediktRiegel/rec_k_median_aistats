#ifndef BAALGO_ROUNDGRAPH_H
#define BAALGO_ROUNDGRAPH_H
#include <map>
#include <vector>
#include <set>
#include "fixedDouble.h"
#include <iostream>

using namespace std;

//all parameters for the bipartit rounding procedure
class roundgraph {
    public:
        map<int, map<int, fixedDouble>> G;
        int numE;
        int root;
        map<int, int> parentNode;
        map<int, map<int, vector<pair<int, int>>>> EdgeToCycle;
        map<int, map<int, vector<int>>> CycleToEdges;
        set<int> removedOriginalFs;

        vector<vector<int>> return_root_paths(int node){
            vector<vector<int>> paths;
            for (auto tuple : G.at(node)) {
                vector<int> path;
                path.push_back(node);
                int neighbour = tuple.first;
                path.push_back(neighbour);
                while (neighbour != root) {
                    neighbour = parentNode.at(neighbour);
                    path.push_back(neighbour);
                }
                paths.push_back(path);
            }
            return paths;
        }

        map<int, fixedDouble> get_none_zero_nodes() {
            map<int, fixedDouble> result;
            for (auto el : G) {
                fixedDouble node_sum(0,0);
                for (auto el2 : el.second) {
                    node_sum += el2.second;
                }
                if (node_sum != 0) {
                    cout << "node: " << el.first << " has a node_sum of " << stringValue(node_sum) << endl;
                    for (auto el2 : el.second) {
                        cout << "\t" << el2.first << ", " << stringValue(el2.second) << endl;
                    }
                }
                result[el.first] = node_sum;
            }
            return result;
        }
};


#endif //BAALGO_ROUNDGRAPH_H
