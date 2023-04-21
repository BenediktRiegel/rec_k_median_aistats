#include "BFS.h"
#include <utility>
#include <vector>
#include <queue>
#include <random>
#include <iostream>
#include <fstream>
#include "loadsave.h"
using namespace std;

vector<int> bfs(const vector<vector<int>>& G, int v){
    //cout << "reserve space for color" << endl;
    //bool color [(*G).size()];
    vector<bool> color;
    color.reserve(G.size());
    for (int i = 0; i < G.size(); ++i) {
        color[i] = false;
    }
    //cout << "reserve space for distance" << endl;
    vector<int> distance;
    distance.reserve(G.size());
    for (int i = 0; i < G.size(); ++i) {
        distance.push_back(0);
    }
    //cout << "start BFS computation" << endl;
    color[v] = true;
    queue<int> Q;
    Q.push(v);
    while (!Q.empty()){
        int currentNode = Q.front();
        Q.pop();
        for (int neighbour : G[currentNode]) {
            if (!color[neighbour]){
                color[neighbour] = true;
                distance[neighbour] = distance[currentNode] + 1;
                Q.push(neighbour);
            }
        }
    }

    return distance;
}


void writeToFile(string path, int v, vector<int> distance){
    string towrite = to_string(v) + "\t" + to_string(distance[0]);
    for (int i = 1; i < distance.size(); ++i) {
        towrite.append("," + to_string(distance[i]));
    }
    fstream file;
    file.open(path, ios_base::out);
    file << towrite;
}


void checkDistances(vector<int> distance, int v){
    int numZero = 0;
    for (int el : distance) {
        if (el == 0) {
            numZero += 1;
        }
    }
    if (distance[v] == 0 && numZero == 1){
        cout << "it seems legit" << endl;
    }
    cout << numZero << endl;
}


long bfstestTime(vector<vector<int>> G, int v) {
    time_t start = time(nullptr);
    checkDistances(bfs(move(G), v), v);
    long duration = time(nullptr) - start;
    return duration;
}
