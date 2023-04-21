#include "loadsave.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include "fixedDouble.h"
#include <cassert>
using namespace std;

vector<string> splitString(string s, char splitter);

map<int, map<int, double>> getDistanceMapByFollowers(string path) {
    map<int, map<int, double>> distances;
    fstream newfile;
    newfile.open(path, ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()) {   //checking whether the file is open
        string tp;
        int index = 0;
        while (getline(newfile, tp)) { //read data from file object and put it into string.
            vector<string> followers = splitString(tp, ',');
            cout << followers.size() << endl;
            map<int, double> fIDs;
            for (int i = 0; i < followers.size(); ++i) {
                fIDs[i] = stod(followers.at(i));
            }
            distances[index] = fIDs;
            ++index;
        }
        newfile.close(); //close the file object.
    }
    return distances;
}

map<int, map<int, double>> getDistanceMap(string path, vector<int> M1, vector<int> M2) {
    map<int, map<int, double>> distances;
    fstream newfile;
    newfile.open(path, ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()) {   //checking whether the file is open
        string tp;
        int i1 = 0;
        while (getline(newfile, tp)) { //read data from file object and put it into string.
            vector<string> row = splitString(tp, ',');
            assert(row.size() == M2.size());
            map<int, double> indexed_row;
            for (int i2 = 0; i2 < M2.size(); ++i2) {
                indexed_row[M2[i2]] = stod(row.at(i2));
            }
            distances[M1[i1]] = indexed_row;
            ++i1;
        }
        cout << "i1 = " << i1 << ", M1.size() = " << M1.size() << endl;
        assert(i1 == M1.size());
        newfile.close(); //close the file object.
    }
    return distances;
}

vector<vector<double>> getDistanceVector(string path) {
    vector<vector<double>> distances;
    fstream newfile;
    newfile.open(path, ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()) {   //checking whether the file is open
        string tp;
        while (getline(newfile, tp)) { //read data from file object and put it into string.
            vector<string> followers = splitString(tp, ',');
            vector<double> fIDs;
            for (int i = 0; i < followers.size(); ++i) {
                fIDs.push_back(stod(followers.at(i)));
            }
            distances.push_back(fIDs);
        }
        newfile.close(); //close the file object.
    }
    return distances;
}

vector<int> getIntVector(string path) {
    fstream newfile;
    vector<int> result;
    newfile.open(path, ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()) {   //checking whether the file is open
        string tp;
        getline(newfile, tp); //read data from file object and put it into string.
        vector<string> stringvec = splitString(tp, ',');
        for (auto& el : stringvec) {
            result.push_back(stoi(el));
        }
        
        newfile.close(); //close the file object.
    }
    return result;
}

vector<double> getDoubleVector(string path) {
    fstream newfile;
    vector<double> result;
    newfile.open(path, ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()) {   //checking whether the file is open
        string tp;
        getline(newfile, tp); //read data from file object and put it into string.
        vector<string> stringvec = splitString(tp, ',');
        for (auto& el : stringvec) {
            result.push_back(stod(el));
        }

        newfile.close(); //close the file object.
    }
    return result;
}

vector<vector<int>> getGraphVector(string path){
    //auto G = new vector<vector<int>>;
    vector<vector<int>> G;
    fstream newfile;
    newfile.open(path,ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){   //checking whether the file is open
        string tp;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            vector<string> followers = splitString(splitString(tp, '\t')[1],',');
            vector<int> fIDs;
            for (auto & follower : followers) {
                fIDs.push_back(stoi(follower));
            }
            G.push_back(fIDs);
        }
        newfile.close(); //close the file object.
    }
    return G;
}


vector<string> splitString(string s, char splitter){
    stringstream ss(s);
    vector<string> result;
    //string print = "Split \"" + s + "\" at \'" + splitter + "\' into ";
    while(ss.good()){
        string substr;
        getline(ss, substr, splitter);
        result.push_back(substr);
    }
    //print.append(vectorToString(result));
    //cout << print << endl;

    return result;
}

string join(vector<int> input, char insert){
    string result;
    result = "";
    if (!input.empty()){
        result += to_string(input[0]);
        for (int i = 1; i < input.size(); ++i) {
            result += insert + to_string(input[i]);
        }
    }
    return result;
}

string join(vector<double> input, char insert){
    string result;
    result = "";
    if (!input.empty()){
        result += to_string(input[0]);
        for (int i = 1; i < input.size(); ++i) {
            result += insert + to_string(input[i]);
        }
    }
    return result;
}


string join(vector<string> input, char insert){
    string result;
    result = "";
    if (!input.empty()){
        result += input[0];
        for (int i = 1; i < input.size(); ++i) {
            result += insert + input[i];
        }
    }
    return result;
}


void write(string path, vector<vector<int>> content){
    fstream file;
    file.open(path, ios_base::out);
    string strcontent;
    strcontent = "";
    if (!content.empty()){
        strcontent += join(content[0], ',');
        for (int i = 1; i < content.size(); ++i) {
            strcontent += '\n' + join(content[i], ',');
        }
    }
    file << strcontent;
    file.close();
}


void cutFileAfter(string inpath, string outpath, int cutafter){
    fstream infile;
    fstream outfile;
    infile.open(inpath, ios_base::in);
    outfile.open(outpath, ios_base::out);

    if (infile.is_open() && outfile.is_open()){   //checking whether the file is open
        string tp;
        int index = 0;
        vector<string> allLines;
        while(index < 500 && getline(infile, tp)){ //read data from file object and put it into string.
            vector<string> line = splitString(tp,',');
            vector<string> newLine;
            unsigned long maxEl = line.size();
            if (maxEl > cutafter) {
                maxEl = cutafter;
            }
            for (int i = 0; i < maxEl; ++i) {
                newLine.push_back(line[i]);
            }
            allLines.push_back(join(newLine, ','));
            ++index;
        }
        outfile << join(allLines, '\n');
        outfile.close(); //close outfile object
        infile.close(); //close infile object.
    }
}


void reverseFile(string inpath, string outpath) {
    fstream infile;
    fstream outfile;
    infile.open(inpath, ios_base::in);
    outfile.open(outpath, ios_base::out);

    if (infile.is_open() && outfile.is_open()){   //checking whether the file is open
        string tp;
        int index = 0;
        vector<string> allLines;
        while(index < 500 && getline(infile, tp)){ //read data from file object and put it into string.
            vector<string> line = splitString(tp,',');
            reverse(line.begin(), line.end());
            allLines.push_back(join(line, ','));
            ++index;
        }
        reverse(allLines.begin(), allLines.end());
        outfile << join(allLines, '\n');
        outfile.close(); //close outfile object
        infile.close(); //close infile object.
    }
}


void FCtoAC(string inpath, string outpath){
    cout << "starting FCtoAC, this can take a while" << endl;
    fstream infile;
    fstream outfile;
    infile.open(inpath, ios_base::in);
    vector<vector<string>>* dFtoC = new vector<vector<string>>;
    string tp;
    int index = 0;
    while(getline(infile, tp)){ //read data from file object and put it into string.
        vector<string> line = splitString(tp,',');
        dFtoC->push_back(line);
        if (index % 50 == 0) {
            cout << "reading line " << index << endl;
        }
        ++index;
    }
    cout << index << " lines read"  << endl;
    infile.close();
    int numA = dFtoC->at(0).size();
    int numF = dFtoC->size();
    int numOnlyC = numA - numF;
    //First Cs, then Fs
    vector<string> clientLine;
    cout << "starting out part" << endl;
    for (int i = 0; i < numOnlyC; ++i) {
        clientLine.emplace_back("-1");
    }
    outfile.open(outpath, ios_base::app);
    for (int a = 0; a < numA; ++a) {
        vector<string> temp = clientLine;
        for (int i = 0; i < numF; ++i) {
            temp.push_back(dFtoC->at(i).at(a));
        }
        temp[a] = "0";
        outfile << (join(temp, ',') + '\n');
        outfile.flush();
        if (a % 1000) {
            cout << a << endl;
        }
    }
    outfile.close();

}
