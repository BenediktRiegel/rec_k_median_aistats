#include "CharLiKUFL.h"
#include "kMSolution.h"
#include "roundgraph.h"
#include "BFS.h"
#include "fixedDouble.h"
#include "bipartitroundingForCharikar.h"
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include "loadsave.h"
#include <omp.h>
#include <ilcplex/ilocplex.h>

using namespace std;


int nextFreeFacility = 0;
vector<int>* C;
vector<int>* F;
map<int, map<int, double>>* dAtoC;
map<int, double>* f;
int k;

map<int, map<int, fixedDouble>>* x;
map<int, fixedDouble>* y;

vector<int>* newToOldA;

#pragma omp threadprivate(nextFreeFacility)
#pragma omp threadprivate(C)
#pragma omp threadprivate(F)
#pragma omp threadprivate(dAtoC)
#pragma omp threadprivate(f)
#pragma omp threadprivate(k)
#pragma omp threadprivate(x)
#pragma omp threadprivate(y)
#pragma omp threadprivate(newToOldA)

typedef IloArray<IloNumVarArray> NumVarMatrix;


void solveLP() {
    IloEnv env;
    IloModel Model(env);
    //Setup X and Y
    //cout << "Setup X and Y" << endl;
    NumVarMatrix cpX(env, (*F).size());
    for (int i = 0; i < (*F).size(); ++i) {
        cpX[i] = IloNumVarArray(env, (*C).size(), 0, 1, ILOFLOAT);
    }
    IloNumVarArray cpY(env, (*F).size(), 0, 1, ILOFLOAT);

    //Sum to minimize
    //cout << "sum to minimize" << endl;
    IloExpr expr0(env);
    for (int i = 0; i < (*F).size(); ++i) {
        for (int j = 0; j < (*C).size(); ++j) {
            expr0 += (cpX[i][j] * (*dAtoC)[(*F)[i]][(*C)[j]]);
        }
        //cout << i << endl;
        expr0 += ((cpY[i]) * ((*f)[(*F)[i]]));
    }
    Model.add(IloMinimize(env, expr0));

    //one facility per client
    //cout << "one facility per client" << endl;
    for (int j = 0; j < (*C).size(); ++j) {
        IloExpr exp1(env);
        for (int i = 0; i < (*F).size(); ++i) {
            exp1 += cpX[i][j];
        }
        Model.add(exp1 == 1);
    }

    //k facilities
    //cout << "k facilities" << endl;
    IloExpr exp2(env);
    for (int i = 0; i < (*F).size(); ++i) {
        exp2 += cpY[i];
    }
    Model.add(exp2 == k);

    //For a client to connect to a facility, the facility has to be open enough.
    //cout << "For a client to connect to a facility, the facility has to be open enough" << endl;
    for (int i = 0; i < (*F).size(); ++i) {
        for (int j = 0; j < (*C).size(); ++j) {
            Model.add(cpX[i][j] <= cpY[i]);
        }
    }

    IloCplex cplex(Model);
    cplex.setOut(env.getNullStream());
    //cout << "solving Model" << endl;
    if (!cplex.solve()) {
        env.error() << "Failed to solve the LP!" << endl;
        throw(-1);
    }

    //get X values
    //cout << "get X values" << endl;
    for (int i = 0; i < (*F).size(); ++i) {
        for (int j = 0; j < (*C).size(); ++j) {
            (*x)[(*F)[i]][(*C)[j]] = fixedDouble(cplex.getValue(cpX[i][j]));
        }
    }

    //get Y values
    //cout << "get Y values" << endl;
    for (int i = 0; i < (*F).size(); ++i) {
        //(*y).set((*F).at(i), cplex.getValue(cpY[i]));
        (*y)[(*F).at(i)] = fixedDouble(cplex.getValue(cpY[i]));
    }
    /*
    cout << "done" << endl;

    cout << "x: " << endl;
    for (int i : (*F)) {
        for (int j : (*C)) {
            if ((*x).at(i).at(j) != 0) {
                cout << "\tx[" << i << "][" << j << "] = " << stringValue((*x).at(i).at(j)) << endl;
            }
        }
    }
    cout << "y: " << endl;
    for (int i : (*F)) {
        if ((*y).at(i) != 0) {
            cout << "\ty[" << i << "] = " << stringValue((*y).at(i)) << endl;
        }
    }
    */

    env.end();
}


vector<int> checkIfLPHasIntegerSolution() {
    vector<int> openFacilities;
    for (int i : (*F)) {
        if (y->at(i) == 1) {
            openFacilities.push_back(i);
        }
        else if (y->at(i) != 0) {
            break;
        }
    }
    return openFacilities;
}

bool keyInMap(map<int, fixedDouble> mapp, int key) {
    if (mapp.find(key) == mapp.end()) {
        return false;
    }
    return true;
}


bool keyInMap(map<int, double> mapp, int key) {
    if (mapp.find(key) == mapp.end()) {
        return false;
    }
    return true;
}


int maxof(const vector<int>& stack){
    int temp = 0;
    for (int el : stack){
        if (el > temp){
            temp = el;
        }
    }
    return temp;
}

int minof(const vector<int>& stack){
    int temp = -1;
    for (int el : stack){
        if (el < temp || temp == -1){
            temp = el;
        }
    }
    return temp;
}


void initilizeNextFreeFacility() {
    int tempF = maxof(*F);
    int tempC = maxof(*C);
    if (tempF > tempC) {
        nextFreeFacility = tempF;
    }
    else {
        nextFreeFacility = tempC;
    }
}


int getNextFreeFacility() {
    nextFreeFacility += 1;
    //(*y).push_back(0);
    (*y)[nextFreeFacility] = 0.0;
    // First fill all x with 0
    //vector<double> temp;
    for (int j : *C) {
        //temp.push_back(0.0);
        (*x)[nextFreeFacility][j] = 0.0;
    }
    //(*x).push_back(temp);
    return nextFreeFacility;
}


map<int, double> computeDav() {
    map<int, double> dav;
    double temp;
    for (int j : *C) {
        temp = 0;
        for (int i : *F) {
            temp += dAtoC->at(i).at(j) * x->at(i).at(j).getDouble();
        }
        dav[j] = temp;
    }
    return dav;
}


double newdAtoC(int a, int j) {
    return (*dAtoC)[(*newToOldA)[a]][j];
}


set<int> splitOneFacility(int i) {
    // set of every facility, that will be newly opened
    set<int> newFacilities;
    // first we want every different x exactly once
    set<fixedDouble> differentX;
    set<int> connectedC;
    for (int j : *C) {
        if ((*x)[i][j] != 0) {
            differentX.insert((*x)[i][j]);
            connectedC.insert(j);
        }
    }
    // only continue, if there are different x's
    if (differentX.size() > 1) {
        // c++ set is already sorted
        vector<fixedDouble> differences_betw_x;
        // compute the difference between the current X and the X before in the list differentX
        {
            auto itr = differentX.begin();
            differences_betw_x.push_back(*itr);
            ++itr;
            auto itrBefore = differentX.begin();
            //itrBefore is always one element before itr, there for computing the difference between neighbouring elements
            //Since differentX is sorted from smallest to biggest, we always get something >= 0
            for (;itr != differentX.end(); ++itr) {
                differences_betw_x.push_back(*itr - *itrBefore);
                ++itrBefore;
            }
        }
        // Now create and open up the new facilities. Also establish connection cost between clients and the facility
        {
            auto itrDiffX = differentX.begin(); //to iterate over the elements
            for (int i2 = 0; i2 < differentX.size(); ++i2) {
                int facility = getNextFreeFacility();
                (*newToOldA).push_back(i);
                newFacilities.insert(facility);
                //(*y).set(facility, differences_betw_x[i2]);
                (*y)[facility] = differences_betw_x[i2];
                set<int> to_keep_C;
                // connect clients to facility, if x_ij is higher or equal, than differentX. Otherwise delete them from the set.
                for (int j : connectedC) {
                    if ((*x)[i][j] > *itrDiffX) {
                        (*x)[facility][j] = (*y)[facility];
                        to_keep_C.insert(j);
                    }
                    else if ((*x)[i][j] == *itrDiffX) {
                        (*x)[facility][j] = (*y)[facility];
                        (*x)[i][j] = 0; // Remove connections between the original facility and all clients
                    }
                    else {
                        (*x)[i][j] = 0; // Remove connections between the original facility and all clients
                    }
                }
                connectedC = to_keep_C;
                ++itrDiffX; //keep iterator updated. For-loop has exactly the size of the itr.
            }
        }
        /*
        cout << "original: " << i << endl;
        fixedDouble sum("0");
        for (auto el : newFacilities) {
            sum += (*y).at(el);
            cout << "\t" << el << ": " << stringValue((*y).at(el)) << endl;
        }
        cout << "\tsum: " << stringValue(sum) << endl;
        */
        // Remove opening of the original facility
        //(*y).set(i, 0);
        //(*y)[i] = 0;
        // The distances between the new facilities and the clients are identical to the original facility i
        // and can be accessed with the method newdAtoC, with the help of newToOldA
    }
    // since C, F, dAtoC, f, x, y are complex datatypes, there is no need to return them
    // newFacilities has to be added to F in the method up, because of a for i in F loop.
    return newFacilities;

}


void initNewToOldA(){
    newToOldA = new vector<int>();
    for (int i = 0; i <= max(maxof(*F), maxof(*C)); ++i) {
        (*newToOldA).push_back(i);
    }
}


vector<int> splitFacilities(map<int, set<int>>* originalFToSplits,
                     set<int>* removedOriginalFs) {
    initNewToOldA();
    /*
    cout << "newToOldA: ";
    for (auto el : *newToOldA) {
        cout << el << ", ";
    }
    cout << endl;*/
    vector<int> newF;
    // Split every facility, if necessary. But remember the original structure, with the above variables
    for (int i : *F) {
        set<int> newFacility = splitOneFacility(i);
        if (!newFacility.empty()) {
            (*removedOriginalFs).insert(i);
            (*originalFToSplits)[i] = newFacility;
            // add every new facility
            for (int i2 : newFacility) {
                newF.push_back(i2);
            }
        }
        else {
            //add facilities that weren't split
            newF.push_back(i);
        }
    }
    return newF;
    // since C, dAtoC, f, x, y are complex datatypes, there is no need to return them
    // OBACHT
}


int getMinDavInCTemp(const map<int, double>& dav, const map<int, bool>& flags) {
    int minJ = -1;
    //fixedDouble minValue(-1);
    //fixedDouble otherValue;
    double minValue = -1;
    double otherValue;
    for (int j : *C) {
        if (!flags.at(j)) {        // if j hasn't got "deleted" yet
            otherValue = dav.at(j);
            if (minValue == -1 or minValue > otherValue) {
                minJ = j;
                minValue = otherValue;
            }
        }
    }
    return minJ;
}


map<int, double> vectorToMap(vector<int> vec) {
    map<int, double> result;
    int index = 0;
    for (int distance : vec) {
        result[index] = distance;
        ++index;
    }
    return result;
}


void addDistanceVectorToTable(vector<double> vector, int node) {
    for (int a = 0; a < dAtoC->size(); ++a) {
        // add the distances
        (*dAtoC)[a][node] = vector[a];
        (*dAtoC)[node][a] = vector[a];
    }
    // dAtoC is a complex data structure, so there is no need to return anything!
    //OBACHT
}


bool checkForNotFlagged(const map<int, bool>& flags) {
    for (auto flag : flags) {
        if (!flag.second) {
            return true;
        }
    }
    return false;
}


set<int> filterClients(const map<int, double>& dav, const vector<vector<int>>& G) {
    map<int, bool> flags;
    for (int j : *C) {
        flags[j] = false;
    }
    set<int> CP;
    //fixedDouble four("4");
    // delete j from C'' <=> set flag to True for j
    while (checkForNotFlagged(flags)) {
        // search for min dav, add it to CP and remove it from C'' <=> flag it
        int newJ = getMinDavInCTemp(dav, flags);
        CP.insert(newJ);
        flags[newJ] = true;      // "del" newJ from C''
        // remove every j with dAtoC[j, newJ] <= 4 * dav[j]
        map<int, double> dOfNewJ = (*dAtoC)[newJ];
        //cout << newJ << ": ";
        for (int j : *C) {
            if (!flags[j]) {    // check if j is already "deleted" from C''
                if (!keyInMap(dOfNewJ, newJ)) {    // calculate the distances, if we are missing one
                    vector<int> newJDistances = bfs(G, newJ);
                    // vector<double> doubleVec(newJDistances.begin(), newJDistances.end());
                    dOfNewJ = vectorToMap(newJDistances);
                }
                if (dOfNewJ[j] <= 4 * dav.at(j)) {
                    flags[j] = true;     // "del" j from C''
                    //cout << j << ", ";
                }
            }
        }
        //cout << endl;
    }

    return CP;
}


// Changed from looking through newF, to looking through F, since fractional facilities have the same
// distance as their split facility.
double maxdistanceToF(int j) {
    double maxd = (*dAtoC)[(*F)[0]][j];
    for (int i = 1; i < (*F).size(); ++i) {
        if (maxd < newdAtoC((*F)[i], j)){
            maxd = newdAtoC((*F)[i], j);
        }
    }
    return maxd;
}


map<int, double> computeR(const set<int>& CP) {

    // R = {j: 0.5 * min(d[j, other] for other in CP) for j in CP}
    map<int, double> R;
    for (int j : CP) {
        int minJ = -1;
        for (int j2 : CP) {
            if (j != j2 && (minJ == -1 || (*dAtoC)[j][minJ] > (*dAtoC)[j][j2])) {
                minJ = j2;
            }
        }
        if (minJ == -1) {
            R[j] = maxdistanceToF(j) / 2.0;
        }
        else {
            R[j] = ((*dAtoC)[j].at(minJ) / 2.0);
        }
    }

    return R;
}


void bundlingFacilities(const set<int>& CP, const vector<int>& newF, map<int, set<int>>* bundles,
                        set<int>* unbundledFacilities) {

    map<int, double> R = computeR(CP);

    int minClient;
    //fixedDouble onepfive = fixedDouble(15, 1);
    for (int i : newF){
        if ((*y).at(i) != 0) {
            minClient = -1;
            for (int j : CP) {
                //cout << "Abstand " << i << " to " << j << " is " << stringValue(newdAtoC(i, j)) << endl;
                if ((*x).at(i).at(j) != 0 &&
                    (minClient == -1 || newdAtoC(i, j) < newdAtoC(i, minClient)) &&
                    newdAtoC(i, j) <= (1.5 * R[j])) {

                    minClient = j;
                }
            }
            if (minClient != -1) {
                (*bundles)[minClient].insert(i);
            }
            else {
                (*unbundledFacilities).insert(i);
            }
        }
    }
}


void matchingBundles(const set<int>& CP, set<int>* unmatchedBundles, set<pair<int, int>>* M) {

    (*unmatchedBundles) = CP;
    while ((*unmatchedBundles).size() > 1) {
        pair<int, int> minPair = pair<int, int>(-1, -1);
        for (int j1 : *unmatchedBundles) {
            for (int j2 : *unmatchedBundles) {
                if (j1 != j2 &&
                    (minPair.first == -1 ||
                     (*dAtoC).at(j1).at(j2) < (*dAtoC).at(minPair.first).at(minPair.second))) {
                    minPair = pair<int, int>(j1, j2);
                }
            }
        }

        (*M).insert(minPair);
        (*unmatchedBundles).erase(minPair.first);
        (*unmatchedBundles).erase(minPair.second);
    }
}


fixedDouble vol(int j, const map<int, set<int>>& bundles) {
    fixedDouble result = 0;
    for (int i : bundles.at(j)) {
        result += y->at(i);
    }

    return result;
}


int graphMatchings(map<int, map<int, fixedDouble>>* G, int root, const set<pair<int, int>>& M,
                   const map<int, set<int>>& bundles, map<int, int>* parentNode){
    int numE = 0;
    // set up matchs in graph
    int a, b;
    for (pair<int, int> p : M) {
        a = p.first;
        b = p.second;
        int gammaNode = getNextFreeFacility();
        int alphaNode = getNextFreeFacility();
        int betaNode = getNextFreeFacility();
        fixedDouble volA = vol(a, bundles);
        fixedDouble volB = vol(b, bundles);
        //cout << "vol(" << a << ") = " << stringValue(volA) << endl;
        //cout << "vol(" << b << ") = " << stringValue(volB) << endl;

        // parentNode
        (*parentNode)[gammaNode] = root;
        (*parentNode)[alphaNode] = gammaNode;
        (*parentNode)[betaNode] = gammaNode;


        // Connecting match (A, B) to root: GammaNode <--> root
        (*G)[root][gammaNode] = volA + volB - 1;
        (*G)[gammaNode][root] = volA + volB - 1;
        // Connecting match (A, B) to bundle A: GammaNode <--> AlphaNode
        fixedDouble oneFixedDouble = fixedDouble(1, 0);
        (*G)[gammaNode][alphaNode] = oneFixedDouble - volA;
        (*G)[alphaNode][gammaNode] = oneFixedDouble - volA;
        //Connecting match (A, B) to bundle B: GammaNode <--> BetaNode
        (*G)[gammaNode][betaNode] = oneFixedDouble - volB;
        (*G)[betaNode][gammaNode] = oneFixedDouble - volB;
        numE += 3;

        // Connecting facilities i in bundle A to bundle A: i <--> AlphaNode
        for (int i : bundles.at(a)) {
            (*parentNode)[i] = alphaNode;
            (*G)[alphaNode][i] = (*y)[i];
            (*G)[i][alphaNode] = (*y)[i];
            numE += 1;
        }

        // Connecting facilities i in bundle B to bundle B: i <--> BetaNode
        for (int i : bundles.at(b)) {
            (*parentNode)[i] = betaNode;
            (*G)[betaNode][i] = (*y)[i];
            (*G)[i][betaNode] = (*y)[i];
            numE += 1;
        }
    }
    // G is complex
    return numE;
}

int graphUnmatchedBundles(map<int, map<int, fixedDouble>>* G, int root, const set<int>& unmatchedBundles,
                          const map<int, set<int>>& bundles, map<int, int>* parentNode) {

    int numE = 0;
    // set up unmatched bundle
    for (int c : unmatchedBundles){
        int muNode = getNextFreeFacility();
        int phiNode = getNextFreeFacility();

        // parentNode
        (*parentNode)[muNode] = root;
        (*parentNode)[phiNode] = muNode;

        fixedDouble volC = vol(c, bundles);
        //cout << "vol(" << c << ") = " << stringValue(volC) << endl;

        // Connecting unmatched bundle C to root: GammaNode <--> root
        (*G)[root][muNode] = volC;
        (*G)[muNode][root] = volC;
        // Connecting negating root connection: GammaNode <--> AlphaNode
        fixedDouble oneFixedDouble = fixedDouble(1, 0);
        (*G)[muNode][phiNode] = oneFixedDouble - volC;
        (*G)[phiNode][muNode] = oneFixedDouble - volC;
        numE += 2;

        // Connecting facilities i in bundle C to negated bundle C
        for (int i : bundles.at(c)) {
            (*parentNode)[i] = phiNode;
            (*G)[phiNode][i] = (*y)[i];
            (*G)[i][phiNode] = (*y)[i];
            numE += 1;
        }
    }
    // G is complex
    return numE;
}


int graphUnbundledFacilities(map<int, map<int, fixedDouble>>* G, int root,
                             set<int>* unbundledFacilities, map<int, int>* parentNode) {
    int numE = 0;
    // set up unbundled facilities in graph
    for (int i : *unbundledFacilities) {
        (*parentNode)[i] = root;
        (*G)[root][i] = y->at(i);
        (*G)[i][root] = y->at(i);
        numE += 1;
    }
    return numE;
}


vector<int> getCycle(int original, int omega1, int omega2, const map<int, int>& parentNode) {
    vector<int> startPath = {original};
    vector<int> endPath = {original};
    int parent1 = omega1;
    int parent2 = omega2;
    startPath.push_back(parent1);
    endPath.push_back(parent2);
    // parent of root is root
    while (parent1 != parent2) {
        if (parent1 != parentNode.at(parent1)) {
            parent1 = parentNode.at(parent1);
            startPath.push_back(parent1);
        }
        if (parent2 != parentNode.at(parent2)) {
            parent2 = parentNode.at(parent2);
            endPath.push_back(parent2);
        }
    }
    endPath.pop_back();
    for (auto itr = endPath.rbegin(); itr != endPath.rend(); ++itr){
        startPath.push_back(*itr);
    }
    return startPath;
}


int graphSplitFacility(map<int, map<int, fixedDouble>>* G,
                       const set<int>& fractionalFs,
                       map<int, int>* parentNode, int original,
                       map<int, map<int, vector<pair<int, int>>>>* EdgeToCycle,
                       map<int, map<int, vector<int>>>* CycleToEdges){
    int numE = 0;
    //G[original] = dict()
    (*parentNode)[original] = original;
    fixedDouble oneFixedDouble(1, 0);
    vector<int> newSplits;
    for (int fraction : fractionalFs) {
        // Buffer vertices
        int omega = getNextFreeFacility();
        newSplits.push_back(omega);
        (*parentNode)[omega] = fraction;

        // edge omega <--> fractional
        (*G)[fraction][omega] = oneFixedDouble - (*y)[fraction];
        (*G)[omega][fraction] = oneFixedDouble - (*y)[fraction];
        numE += 1;
        // edge omega <--> split
        (*G)[omega][original] = (*y)[fraction];
        (*G)[original][omega] = (*y)[fraction];
        numE += 1;
    }

    //get cycles
    for (int omega1 : newSplits) {
        for (int omega2 : newSplits) {
            if (omega1 != omega2) {
                if ((*CycleToEdges)[omega1][omega2].empty()) {
                    vector<int> cycle = getCycle(original, omega1, omega2, *parentNode);
                    (*CycleToEdges)[omega1][omega2] = cycle;
                    (*CycleToEdges)[omega2][omega1] = cycle;
                    auto cyclepair1 = pair<int, int> (omega1, omega2);
                    auto cyclepair2 = pair<int, int> (omega2, omega1);
                    for (int i = 0; i < cycle.size() - 1; ++i) {
                        (*EdgeToCycle)[cycle[i]][cycle[i + 1]].push_back(cyclepair1);
                        (*EdgeToCycle)[cycle[i]][cycle[i + 1]].push_back(cyclepair2);
                        (*EdgeToCycle)[cycle[i + 1]][cycle[i]].push_back(cyclepair1);
                        (*EdgeToCycle)[cycle[i + 1]][cycle[i]].push_back(cyclepair2);
                    }
                }
            }
        }
    }

    return numE;
}


int graphSplitFacilities(map<int, map<int, fixedDouble>>* G,
                         const set<int>& removedOriginalFs,
                         const map<int, set<int>>& originalFToSplits,
                         map<int, int>* parentNode,
                         map<int, map<int, vector<pair<int, int>>>>* EdgeToCycle,
                         map<int, map<int, vector<int>>>* CycleToEdges) {
    int numE = 0;
    // set up unbundled facilities in graph
    for (int original : removedOriginalFs) {
        (*parentNode)[original] = original;
        numE += graphSplitFacility(G, originalFToSplits.at(original),
                                   parentNode, original, EdgeToCycle, CycleToEdges);
    }
    // G is complex
    return numE;
}


roundgraph createGraph(const vector<int>& newF, const set<pair<int, int>>& M,
                 const set<int>& unmatchedBundles, set<int>* unbundledFacilities,
                 const map<int, set<int>>& bundles, const set<int>& removedOriginalFs,
                 const map<int, set<int>>& originalFToSplits) {

    /*
    cout << "creating graph + y sum" << endl;
    fixedDouble sum = 0;
    for (auto el : *y) {
        sum += el.second;
    }
    cout << "All y's sum up to: " << stringValue(sum) << endl;
    */
    roundgraph result;
    result.root = getNextFreeFacility();
    result.parentNode[result.root] = result.root;
    result.numE = 0;
    // create graph part for the matchings and update the number of edges
    result.numE += graphMatchings(&result.G, result.root, M, bundles, &result.parentNode);
    // create graph part for the unmatched bundles and update the number of edges
    result.numE += graphUnmatchedBundles(&result.G, result.root, unmatchedBundles, bundles,
                                         &result.parentNode);
    // create graph part for the unbundled facilities and update the number of edges
    result.numE += graphUnbundledFacilities(&result.G, result.root, unbundledFacilities,
                                            &result.parentNode);
    // create graph part for the split facilities and update the number of edges
    // Also get the cycles
    result.numE += graphSplitFacilities(&result.G, removedOriginalFs, originalFToSplits,
                                        &result.parentNode, &result.EdgeToCycle,
                                        &result.CycleToEdges);
    return result;
}


map<int, int> cleanRoundingSolution(const map<int, map<int, int>>& roundingSolution) {
    map<int, int> result;
    for (int i : *F) {
        result[i] = 0;
        if (roundingSolution.find(i) != roundingSolution.end()) {
            for (auto secondlayer : roundingSolution.at(i)) {
                result[i] += secondlayer.second; //secondlayer.second is the edgeweight of [i][secondlayer.first]
            }
        }
    }
    return result;
}



vector<int> sample(const vector<int>& newF, const set<pair<int, int>>& M, const set<int>& unmatchedBundles,
                   set<int>* unbundledFacilities, const map<int, set<int>>& bundles,
                   const set<int>& removedOriginalFs, const map<int, set<int>>& originalFToSplits) {

    // Create the bipartit graph for the dependent rounding procedure
    cout << "create graph" << endl;
    roundgraph rgraph = createGraph(newF, M, unmatchedBundles, unbundledFacilities, bundles, \
                                    removedOriginalFs, originalFToSplits);

    //rgraph.get_none_zero_nodes();
    //cout << stringValue(none_zero_map.at(rgraph.root)) << endl;
    //cout << "y.at(251) = " << stringValue((*y).at(251)) << endl;
    //cout << "root of G: " << rgraph.root << endl;
    /*
    auto root_paths = rgraph.return_root_paths(0);
    for (const auto& el : root_paths) {
        cout << join(el, ',') << endl;
    }*/

    //DEBUG
    /*
    auto root = rgraph.root;
    fixedDouble sum("0");
    for (auto el : rgraph.G.at(root)) {
        sum += el.second;
    }

    //DEBUG
    cout << "sum at the root is " << stringValue(sum) << endl;
    sum = fixedDouble("0");
    for (int i : newF) {
        sum += (*y).at(i);
    }
    cout << "the sum of all facilities y is " << stringValue(sum) << endl;
    */
    // Round
    map<int, map<int, int>> roundingSolution = bipartitroundingForCharikar::solve(&rgraph);
    map<int, int> roundedVar = cleanRoundingSolution(roundingSolution);
    // Open the facilities
    vector<int> openedF;
    // if y[i] == 0, then i is not in the graph G and therefore not in roundingSolution.
    // To avoid a KeyError y[i] == 0 is tested

    /*
    //DEBUG
    sum = fixedDouble("0");
    fixedDouble other_sum("0");
    fixedDouble zeroy("0");
    for (auto el : roundedVar) {
        if (find((*F).begin(), (*F).end(), el.first) != (*F).end()) {
            other_sum += el.second;
            if ((*y).find(el.first) != (*y).end()) {
                if ((*y).at(el.first) != 0) {
                    zeroy += el.second;
                    cout << stringValue(el.second) << endl;
                }
            }
        }
        sum += el.second;
    }
    cout << "the sum of the roundedVar is " << stringValue(sum) << endl;
    cout << "the other_sum is " << stringValue(other_sum) << endl;
    cout << "the zeroy is " << stringValue(zeroy) << endl;
    cout << "both sums should be " << k << endl;
    */
    for (int i : *F) {
        if ((*y)[i] != 0 and roundedVar[i] == 1) {
            openedF.push_back(i);
        }
    }

    return openedF;
}


double compute_clientCost(const vector<int>& openedF, int j) {
    double mindistance(-1);
    //int mini = -1;
    for (int i : openedF) {
        if (mindistance == -1 || (mindistance > dAtoC->at(i).at(j))){
            mindistance = dAtoC->at(i).at(j);
            //mini = i;
        }
    }
    //cout << "\t" << j << " is serveded by " << mini << " with " << stringValue(mindistance) << endl;
    return mindistance;
}


double compute_openingCost(const vector<int>& openedF) {
    //fixedDouble cost(0);
    double cost = 0;
    //cout << "compute f: " << endl;
    for (int i : openedF) {
        /*
        if (cost + f->at(i) < 0) {
            cout << "opening cost becomes negative for " << cost << " + " << stringValue(f->at(i)) << " = " << cost + f->at(i);
        }*/
        cost += f->at(i);
        //cout << "\t" << stringValue(f->at(i)) << endl;
    }
    //cout << "\tsum: " << stringValue(cost) << endl;
    return cost;
}

double compute_serviceCost(const vector<int>& openedF) {
    //fixedDouble cost(0);
    double cost = 0;
    //cout << "compute service cost: " << endl;
    for (int j : *C) {
        /*
        if (cost + compute_clientCost(openedF, j) < 0) {
            cout << "opening cost becomes negative for " << cost << " + " << compute_clientCost(openedF, j) << " = " << cost + compute_clientCost(openedF, j);
        }*/
        cost += compute_clientCost(openedF, j);
        //cout << "\t" << stringValue(compute_clientCost(openedF, j)) << endl;
    }
    //cout << "\tsum: " << stringValue(cost) << endl;
    return cost;
}


pair<kMSolution, bool> solvekUFL(vector<int>* inC, vector<int>* inF, map<int, map<int, double>>* indAtoC,
    int ink, map<int, double>* inf, const vector<vector<int>>& G) {
    C = inC;
    F = inF;
    dAtoC = indAtoC;
    k = ink;
    f = inf;

    x = new map<int, map<int, fixedDouble>>();
    y = new map<int, fixedDouble>();


    // 1. solve LP
    cout << "solve LP" << endl;
    solveLP();

    /*
    cout << "x: " << endl;
    for (int i : (*F)) {
        for (int j : (*C)) {
            cout << "\tx[" << i << "][" << j << "] = " << (*x).at(i).at(j) << endl;
        }
    }
    cout << "y: " << endl;
    for (int i : (*F)) {
        cout << "\ty[" << i << "] = " << (*y).at(i) << endl;
    }
    */

    // x[i, j] = the service connection between facility i and client j
    // y[i] = how much facility i is open

    cout << "check for integer solution" << endl;
    vector<int> openedF = checkIfLPHasIntegerSolution();
    bool lp_only = true;
    // If the lp already solved it perfectly, then there is no need to continue
    if (openedF.size() != k) {
        lp_only = false;
        // compute dav, before the facilities are split. (That makes it easier)
        //Because of datatypes
        cout << "dav" << endl;
        map<int, double> dav = computeDav();
        // 2. split facilities
        initilizeNextFreeFacility();

        set<int> removedOriginalFs;
        set<int> removedFs;
        map<int, set<int>> originalFToSplits;
        //Because of datatypes
        cout << "split F" << endl;
        vector<int> newF = splitFacilities(&originalFToSplits, &removedOriginalFs);
        //DEBUG
        /*
        for (auto el : originalFToSplits) {
            cout << el.first << "got split into: ";
            for (auto el2 : el.second) {
                cout << el2 << ", ";
            }
        }
        cout << endl;
        */
        // 3. Filtering phase
        //Because of datatypes
        cout << "filterClients" << endl;
        set<int> CP = filterClients(dav, G);
        // 4. Bundling phase
        map<int, set<int>> bundles;
        set<int> unbundledFacilities;
        cout << "bundling F" << endl;
        bundlingFacilities(CP, newF, &bundles, &unbundledFacilities);

        //DEBUG
        /*
        cout << "251 is in the following bundles: " << endl;
        for (auto el : bundles) {
            if(find(el.second.begin(), el.second.end(), 251) != el.second.end()) {
                cout << "\t" << el.first << endl;
            }
        }
        cout << "bundle 251: " << endl;
        for (auto el : bundles.at(251)) {
            cout << "\t" << el << ", " << stringValue((*x).at(251).at(el)) << endl;
        }
        */
        delete(newToOldA);
        // 5. Matching phase
        set<int> unmatchedBundles;
        set<pair<int, int>> M;
        //Because of datatypes
        cout << "matching bundles" << endl;
        matchingBundles(CP, &unmatchedBundles, &M);
        /*
        for (auto el : M) {
            if (el.first == 251) {
                cout << "251 is matched with " << el.second << endl;
            }
            else if (el.second == 251) {
                cout << "251 is matched with " << el.first << endl;
            }
        }
         */
        // for j in CP:
        //     print('U_' + str(j) + ': ' + str(bundles[j]))
        // for j in unmatchedBundles:
        //     print('unmatched bundle ' + str(j))
        // print('M: ' + str(M))
        // print('unbundledFacilities: ' + str(unbundledFacilities))
        // for i in removedOriginalFs:
        //     print('splitted ' + str(i) + ': ' + str(originalFToSplits[i]))
        // 6. Sampling phase
        cout << "sample" << endl;
        openedF = sample(newF, M, unmatchedBundles, &unbundledFacilities, bundles, removedOriginalFs,
            originalFToSplits);
    }
    // 7. Compute cost and return
    delete(x);
    delete(y);
    cout << "compute costs" << endl;
    double service_cost = compute_serviceCost(openedF);
    double opening_cost = compute_openingCost(openedF);
    kMSolution solution;
    solution.set(openedF, service_cost, opening_cost);
    return pair<kMSolution, bool> (solution, lp_only);
}


kMSolution testkUFLNoLP(vector<int>* inC, vector<int>* inF, map<int, map<int, double>>* indAtoC,
    int ink, map<int, double>* inf, const vector<vector<int>>& G, map<int, map<int, fixedDouble>>* inx,
    map<int, fixedDouble>* iny) {
    C = inC;
    F = inF;
    dAtoC = indAtoC;
    k = ink;
    f = inf;

    x = inx;
    y = iny;
    // x[i, j] = the service connection between facility i and client j
    // y[i] = how much facility i is open
    cout << "check for integer solution" << endl;
    vector<int> openedF = checkIfLPHasIntegerSolution();
    // If the lp already solved it perfectly, then there is no need to continue
    if (openedF.size() != k) {
        cout << "no integer solution" << endl;
        // compute dav, before the facilities are split. (That makes it easier)
        //Because of datatypes
        cout << "compute Dav" << endl;
        map<int, double> dav = computeDav();
        // 2. split facilities
        cout << "init next free f" << endl;
        initilizeNextFreeFacility();

        set<int> removedOriginalFs;
        set<int> removedFs;
        map<int, set<int>> originalFToSplits;
        //Because of datatypes
        vector<int> newF = splitFacilities(&originalFToSplits, &removedOriginalFs);
        // 3. Filtering phase
        //Because of datatypes
        cout << "filter clients" << endl;
        set<int> CP = filterClients(dav, G);
        // 4. Bundling phase
        map<int, set<int>> bundles;
        set<int> unbundledFacilities;
        cout << "bundling facilities" << endl;
        bundlingFacilities(CP, newF, &bundles, &unbundledFacilities);
        delete(newToOldA);
        /*
        fixedDouble fsum("0");
        cout << "bundles: " << endl;
        for (auto el : bundles) {
            cout << "\t" << el.first << ": ";
            for (auto i : el.second) {
                fsum += (*y).at(i);
                cout << "(" << i << ", " << stringValue((*y).at(i)) << "), ";
            }
            cout << endl;
        }
        for (auto el : unbundledFacilities) {
            fsum += (*y).at(el);
        }
        cout << "summed unbundled f and bundled f. The resulting y sum is: " << stringValue(fsum) << endl;
         */
        // 5. Matching phase
        set<int> unmatchedBundles;
        set<pair<int, int>> M;
        //Because of datatypes
        cout << "matching bundles" << endl;
        matchingBundles(CP, &unmatchedBundles, &M);
        // for j in CP:
        //     print('U_' + str(j) + ': ' + str(bundles[j]))
        // for j in unmatchedBundles:
        //     print('unmatched bundle ' + str(j))
        // print('M: ' + str(M))
        // print('unbundledFacilities: ' + str(unbundledFacilities))
        // for i in removedOriginalFs:
        //     print('splitted ' + str(i) + ': ' + str(originalFToSplits[i]))
        // 6. Sampling phase
        cout << "sample" << endl;
        openedF = sample(newF, M, unmatchedBundles, &unbundledFacilities, bundles, removedOriginalFs, \
            originalFToSplits);
    }
    cout << "finished ChariLi" << endl;
    cout << "x and y deleted" << endl;
    // 7. Compute cost and return
    double service_cost = compute_serviceCost(openedF);
    double opening_cost = compute_openingCost(openedF);
    cout << "service cost: " << service_cost << endl;
    cout << "opening_cost: " << opening_cost << endl;
    kMSolution solution;
    solution.set(openedF, service_cost, opening_cost);
    cout << "solution set" << endl;
    return solution;
}
