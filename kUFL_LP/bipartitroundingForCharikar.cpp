#include "bipartitroundingForCharikar.h"
#include <iostream>
#include <random>
#include <map>
#include "roundgraph.h"
#include "fixedDouble.h"

using namespace std;
// G = {node:(next, weight) for every {node, next} in E}

random_device rd;
mt19937 mt(rd());
//mt19937 mt(0);
uniform_real_distribution<double> dist(0.0, 1.0);

bool checkIfInMap(const map<int, map<int, vector<int>>>& values, int firstIndex, int secondIndex) {
    if (values.find(firstIndex) != values.end()) {
        return values.at(firstIndex).find(secondIndex) != values.at(firstIndex).end();
    }

    return false;
}


bool checkIfInSet(const set<int>& values, int value) {
    return values.find(value) != values.end();
}

/*
def debugInG(G, v, u):
if v in G:
if u in G[v]:
print('(' + str(v) + ', ' + str(u) + ') is still in G')
else:
print('(' + str(v) + ', ' + str(u) + ') not in G')
else:
print('(' + str(v) + ', ' + str(u) + ') not in G')
*/

//Update the cycle structures
void updateCycles(int vertexA, int vertexB, roundgraph* rgraph) {
    for (auto cycle : (*rgraph).EdgeToCycle[vertexA][vertexB]) {
        //if (!(*rgraph).CycleToEdges[cycle.first][cycle.second].empty()) {
        if (checkIfInMap((*rgraph).CycleToEdges, cycle.first, cycle.second)) {
            //cout << "removed cycle: (" << cycle.first << ", " << cycle.second << ")" << endl;
            (*rgraph).CycleToEdges[cycle.first].erase(cycle.second);

            if ((*rgraph).CycleToEdges[cycle.first].empty()) {
                (*rgraph).CycleToEdges.erase(cycle.first);
            }
        }
    }
}

void cleanGraph(roundgraph* rgraph, map<int, map<int, int>>* result) {
    int removedEdges = 0;

    //Get all edges that are already 1 or 0, so don't need rounding
    for (const auto& secG : (*rgraph).G) {
        for (const auto& thirG : secG.second) {
            if (thirG.second == 0 or thirG.second == 1) {
               (*result)[secG.first][thirG.first] = (int) thirG.second.getDouble();
            }
        }
    }

    //Now update all the cycles that depended on those edges
    for (const auto& secRes : *result) {
        for (const auto& thirRes : secRes.second) {
            (*rgraph).G[secRes.first].erase(thirRes.first);
            updateCycles(secRes.first, thirRes.first, rgraph);
            ++removedEdges;
        }
    }

    vector<int> toDel;
    //delete vertices that have no edges anymore
    //find these vertices
    for (auto i : (*rgraph).G) {
        if (i.second.empty()) {
            toDel.push_back(i.first);
        }
    }
    //now delete them
    for (int i : toDel) {
        (*rgraph).G.erase(i);
    }
    removedEdges /= 2;
    (*rgraph).numE -= removedEdges;
}


/*
//Not needed anymore
        def copyGraph(G):
temp = {i: G[i].copy() for i in G.keys()}
return temp
*/


bool doesEdgeExist(const map<int, map<int, fixedDouble>>& G, int nodeA, int nodeB) {
    return G.at(nodeA).find(nodeB) != G.at(nodeA).end();
}


// Goes up parentNode until root and if every edge exists, then it returns the path
vector<int> getPathToRoot(const roundgraph& rgraph, int node){
    vector<int> path;
    path.push_back(node);
    int temp;
    while (node != rgraph.root) {
        temp = rgraph.parentNode.at(node);
        if (doesEdgeExist(rgraph.G, node, temp)) {
            node = temp;
            path.push_back(node);
        } else {
            path = vector<int>();
            break;
        }
    }
    return path;
}


/*
//not needed anymore
vector<int> getCycle(const roundgraph& rgraph, int node) {
    neighbours = list(G[node].keys())
    vector<int> cycle;
    vector<int> path1;
    vector<int> path2;
    int n1 = 0;
    // try and get two neighbour <-- ... --> root paths
    for i in range(0, len(neighbours)):
        path1 = getPathToRoot(G, root, neighbours[i], parentNode)
        if len(path1) != 0:
            n1 = i + 1
            break
        else:
            path1 = list()
    if n1 != 0:
        for i in range(n1, len(neighbours)):
            path2 = getPathToRoot(G, root, neighbours[i], parentNode)
            if len(path2) != 0:
                break
            else:
                path2 = list()

    //if two exist, then there is a cycle
    if len(path1) != 0 and len(path2) != 0:
        cycle.append(node)
        cycle.extend(path1)
        path2 = path2[:-1]
        path2.reverse()
        cycle.extend(path2)
        cycle.append(node)

    return cycle
}
*/

// def getSomeCycle(G, root, removedOriginalFs, parentNode):
//     cycle = list()
//     for node in removedOriginalFs:
//         if len(G[node].keys()) >= 2:
//             cycle = getCycle(G, root, node, parentNode)
//         if len(cycle) != 0:
//             break
//     return cycle


pair<int, int> getSomeCycle(const map<int, map<int, vector<int>>>& CycleToEdges) {
    pair<int, int> cycle;
    if (!CycleToEdges.empty()) {
        auto c = *CycleToEdges.begin();
        cycle.first = c.first;
        cycle.second = (*c.second.begin()).first;
    }
    return cycle;
}


// Destroys the Graph, if not copied
// so we copy it ;)
vector<int> maximalPath(map<int, map<int, fixedDouble>> G) {
    vector<int> pathFront;
    vector<int> pathBack;
    // first node
    int start = G.begin()->first;
    int end = start;
    pathFront.push_back(start);
    int nextN = start;
    while (!G[nextN].empty()) {
        start = G[nextN].begin()->first;
        G[nextN].erase(start);
        G[start].erase(nextN);
        pathFront.push_back(start);
        nextN = start;
    }

    nextN = end;
    while (!G[nextN].empty()) {
        end = G[nextN].begin()->first;
        G[nextN].erase(end);
        G[end].erase(nextN);
        pathBack.push_back(end);
        nextN = end;
    }

    vector<int> path;
    for (auto itr = pathFront.rbegin(); itr != pathFront.rend(); ++itr) {
        path.push_back(*itr);
    }
    for (int i : pathBack) {
        path.push_back(i);
    }

    return path;
}


/*
// Not needed anymore
def maximalPath_Copy(G):
    return maximalPath(copyGraph(G))
*/


vector<int> getCycleOrMaxPath(roundgraph* rgraph) {
    pair<int, int> cycle = getSomeCycle((*rgraph).CycleToEdges);
    //cout << "cycle: " << cycle.first << ", " << cycle.second << endl;
    vector<int> nodes = (*rgraph).CycleToEdges[cycle.first][cycle.second];
    //does the cycle even exist?
    if (nodes.empty()) {
        //If note get a maximalPath
        nodes = maximalPath((*rgraph).G);
        //print('path: ' + str(nodes))
    }
    //else:
        //print('cycle: ' + str(nodes))
    return nodes;
}


vector<fixedDouble> getMinMaxM1M2(const map<int, map<int, fixedDouble>>& G, const vector<int>& pathE) {
    fixedDouble minM1(0, 0);
    fixedDouble minM2(0, 0);
    fixedDouble maxM1(0, 0);
    fixedDouble maxM2(0, 0);
    int indicator = 0;
    fixedDouble temp(0, 0);
    // get max and min from M1 and M2 to get alpha and beta
    for (int i = 0; i < pathE.size() - 1; ++i) {
        temp = G.at(pathE[i]).at(pathE[i + 1]);
        // in M1
        if (indicator % 2 == 0) {
            if ((temp < minM1) || (minM1 == 0)) {
                minM1 = temp;
            }
            if (temp > maxM1) {
                maxM1 = temp;
            }
        }
            // in M2
        else {
            if ((temp < minM2) || (minM2 == 0)) {
                minM2 = temp;
            }
            if (temp > maxM2) {
                maxM2 = temp;
            }
        }

        indicator = (indicator + 1) % 2;
    }

    return {minM1, maxM1, minM2, maxM2};
}


pair<fixedDouble, fixedDouble> getAlphaBeta(const map<int, map<int, fixedDouble>>& G,
                                  const vector<int>& pathE) {

    //0: minM1, 1: maxM1, 2: minM2, 3: maxM2
    vector<fixedDouble> minsandmax = getMinMaxM1M2(G, pathE);
    // AlphaAndBeta
    //beta = min(minM1, 1 - maxM2);
    fixedDouble onefd(1, 0);
    fixedDouble beta = min(minsandmax[0], onefd - minsandmax[3]);
    fixedDouble alpha;
    // It can happen, that there is only one edge in M1.
    //if (minM2 == 0) {
    if (minsandmax[2] == 0) {
        //alpha = maxM1;
        alpha = onefd - minsandmax[1];
    }
    else {
        //alpha = min(minM2, 1 - maxM1)
        alpha = min(minsandmax[2], onefd - minsandmax[1]);
    }

    //cout << "Alpha: " << to_string(alpha) << ", Beta: " << to_string(beta) << endl;
    fixedDouble alphaProb = (beta / (alpha + beta));
    //cout << "Alpha-Prob: " << stringValue(alphaProb) << endl;
    fixedDouble outcome = stofixedd(to_string(dist(mt)));
    //cout << "outcome: " << stringValue(outcome) << endl;
    if (outcome < alphaProb) {
        // go alpha
        beta = fixedDouble(0, 0);
        //cout << "Alpha was chosen" << endl;
    }
    else {
        // go beta
        alpha = fixedDouble(0, 0);
        //cout << "Beta was chosen" << endl;
    }
    return pair<fixedDouble, fixedDouble> (alpha, beta);
}


int roundEdges(roundgraph* rgraph, const vector<int>& pathE, map<int, map<int, int>>* result) {
    int edgesRemoved = 0;
    // print('numE: ' + str(numE))
    auto alphabeta = getAlphaBeta((*rgraph).G, pathE);
    fixedDouble alpha = alphabeta.first;
    fixedDouble beta = alphabeta.second;
    // update edge weights
    int indicator = 0;
    for (int i = 0; i < pathE.size() - 1; ++i) {
        if (doesEdgeExist((*rgraph).G, pathE[i], pathE[i + 1])) {
            // in M1
            if (indicator % 2 == 0) {
                (*rgraph).G[pathE[i]][pathE[i + 1]] =
                        (*rgraph).G[pathE[i]][pathE[i + 1]] + alpha - beta;

                (*rgraph).G[pathE[i + 1]][pathE[i]] = (*rgraph).G[pathE[i]][pathE[i + 1]];

            }
            // in M2
            else {
                if (doesEdgeExist((*rgraph).G, pathE[i + 1], pathE[i])) {

                    (*rgraph).G[pathE[i]][pathE[i + 1]] =
                            (*rgraph).G[pathE[i]][pathE[i + 1]] + beta - alpha;

                    (*rgraph).G[pathE[i + 1]][pathE[i]] = (*rgraph).G[pathE[i]][pathE[i + 1]];

                }
            }

            // print('(' + str(pathE[i]) + ', ' + str(pathE[i + 1]) + ')')
            if ((*rgraph).G[pathE[i]][pathE[i + 1]] == 0 ||
                    (*rgraph).G[pathE[i]][pathE[i + 1]] == 1) {

                //cout << "edge (" << pathE[i] << ", " << pathE[i + 1] << ") was rounded to " \
                << to_string((*rgraph).G[pathE[i]][pathE[i + 1]]) << "\n" << endl;
                // update result
                (*result)[pathE[i]][pathE[i + 1]] = (int) (*rgraph).G[pathE[i]][pathE[i + 1]].getDouble();
                (*result)[pathE[i + 1]][pathE[i]] = (int) (*rgraph).G[pathE[i]][pathE[i + 1]].getDouble();
                //print('edge (' + str(pathE[i]) + ', ' + str(pathE[i + 1]) + ') removed')
                // update graph G
                edgesRemoved += 1;
                (*rgraph).G[pathE[i]].erase(pathE[i + 1]);
                (*rgraph).G[pathE[i + 1]].erase(pathE[i]);

                if ((*rgraph).G[pathE[i]].empty()) {
                    (*rgraph).G.erase(pathE[i]);
                    if (checkIfInSet((*rgraph).removedOriginalFs, pathE[i])) {
                        (*rgraph).removedOriginalFs.erase(pathE[i]);
                    }
                }
                if ((*rgraph).G[pathE[i + 1]].empty()) {
                    (*rgraph).G.erase(pathE[i + 1]);
                    if (checkIfInSet((*rgraph).removedOriginalFs, pathE[i + 1])) {
                        (*rgraph).removedOriginalFs.erase(pathE[i + 1]);
                    }
                }
                //updateCycles
                updateCycles(pathE[i], pathE[i + 1], rgraph);
                //cout << "cycles updated" << endl;
            }
        }
        indicator = (indicator + 1) % 2;
    }
    //cout << "edges removed: " << edgesRemoved << endl;
    return edgesRemoved;
}


map<int, map<int, int>> bipartitroundingForCharikar::solve(roundgraph* rgraph) {
    /*
    int testEdges = 0;
    for i in G.keys():
        testEdges += len(G[i].keys())
    testEdges /= 2 */
    // print('parentNode: ' + str(parentNode))
    // print('numE before cleaning: ' + str(numE))
    // print('G before clening: ' + str(G))
    //cout << "numE: " << (*rgraph).numE << endl;
    map<int, map<int, int>> result;
    /*
    vector<int> rneighbours;
    fixedDouble preRootSum("0");
    cout << "root neighbours: " << endl;
    for (auto el : (*rgraph).G.at((*rgraph).root)) {
        rneighbours.push_back(el.first);
        preRootSum += el.second;
        cout << "\t" << el.first << endl;
    }
    */
    cleanGraph(rgraph, &result);

    // print('G after clening: ' + str(G))
    // print('numE after cleaning: ' + str(numE))
    int index = 0;
    vector<int> pathE;
    while ((*rgraph).numE != 0) {
        index += 1;
        // print(str(index) + '. path')
        pathE = getCycleOrMaxPath(rgraph);
        /*
        cout << "\n\npath: " << endl;
        for (int i : pathE) {
            cout << i << ", ";
        }
        cout << "\n\n" << endl;
        */
        // print(pathE)
        (*rgraph).numE -= roundEdges(rgraph, pathE, &result);
        // if cycle is not None:
        //     updateCycles(pathE, EdgeToCycle, CycleToEdges)
        // print('numE: ' + str(numE))
    }
    //root test
    //cout << "root check" << endl;
    //cout << "preRoot sum: " << stringValue(preRootSum) << endl;

    return result;
}
