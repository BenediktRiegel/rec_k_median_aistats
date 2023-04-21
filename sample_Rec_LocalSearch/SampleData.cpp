#include "SampleData.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include "BFS.h"

using namespace std;


random_device rd_device;
mt19937 random_engine(rd_device());


void print_vec(const vector<unsigned long>& vec) {
	string result = "[" + to_string(vec.at(0));
	for (int i = 1; i < vec.size(); ++i) {
		result += ", " + to_string(vec.at(i));
	}
	result += "]";
	cout << result << endl;
}


bool el_in_vec(int el, vector<int> vec) {
    if(std::find(vec.begin(), vec.end(), el) != vec.end()) {
        return true;
    }
    return false;
}


vector<int> uniform_sampling(vector<int> C, int new_amount) {
    if (C.size() <= new_amount) {
        return C;
    }
    shuffle(C.begin(), C.end(), random_engine);
    vector<int> sampled_C;
    sampled_C.reserve(new_amount);
    for (int i = 0; i < new_amount; ++i) {
        sampled_C.push_back(C.at(i));
    }
    sort(sampled_C.begin(), sampled_C.end());
    return sampled_C;
}


void remove_by_value(vector<int>* vec, int el) {
    (*vec).erase(std::remove((*vec).begin(), (*vec).end(), el), (*vec).end());
}


void save_new_c(const string& path, int new_c) {
    fstream file;
    file.open((path + "new_C.txt"), ios_base::app);
    if (file.is_open()) {
        file << "," <<  to_string(new_c);
    }
    file.close();
}


double compute_score_sum(const vector<unsigned long> inv_scores, const vector<int>& old_C) {
    cout << "compute_score_sum" << endl;
	double score_sum = 0;
	for (int c : old_C) {
		score_sum += (1.0 / inv_scores.at(c));
	}
	if (score_sum < 0) {
		cout << "negative score_sum!!!" << endl;
	}
	return score_sum;
}


void init_D_sampling(vector<unsigned long>* scores, vector<int>* old_C, vector<int>* new_C,
                     const vector<vector<int>>& G) {
    //init random generator
    cout << "init D sampling" << endl;
    uniform_int_distribution<> rd_int(0, (*old_C).size() - 1);

    // Chose first new_c randomly and uniformly for all old_c
    int rd_index = rd_int(random_engine);
    (*new_C).push_back((*old_C).at(rd_index));
    (*old_C).erase((*old_C).begin() + rd_index);
    //compute distances between old_C and new_C
    cout << "do bfs for " << (*new_C).at(0) << endl;
    vector<int> distances = bfs(G, (*new_C).at(0));
    cout << "update scores" << endl;
    for (int i = 0; i < distances.size(); ++i) {
        if (i == rd_index) {
            (*scores).push_back(0);
        }
        else {
            (*scores).push_back(distances.at(i));
        }
    }
}


void sample_next_C(vector<unsigned long>* inv_scores, vector<int>* old_C, vector<int>* new_C,
               const vector<vector<int>>& G, double rd) {

    cout << "function sample next C" << endl;
    auto prob = compute_score_sum(*inv_scores, *old_C) * rd;
    double current = 0;
    int to_add_index = 0;
    cout << "chose next C" << endl;
    for (int index = 0; index < (*old_C).size(); ++index) {
        current += 1.0 / (*inv_scores).at((*old_C).at(index));
        if (prob < current) {
            to_add_index = index;
            break;
        }
    }
    int to_add_c = (*old_C).at(to_add_index);
    cout << "do bfs for " << to_add_c << endl;
    vector<int> distances = bfs(G, to_add_c);
    (*old_C).erase((*old_C).begin() + to_add_index);

    cout << "update scores" << endl;
    for (int old_c : (*old_C)) {
        (*inv_scores)[old_c] += distances.at(old_c);
        if ((*inv_scores).at(old_c) < 0) {
            cout << "negative score!!!" << endl;
            exit(-1);
        }
    }
    (*new_C).push_back(to_add_c);
}


void print_actual_score_sum(vector<unsigned long> inv_scores, const vector<int>& old_C) {
	cout << "actual_score_sum:\n\t" << compute_score_sum(inv_scores, old_C) << endl;
}


vector<int> D_sampling(vector<int> old_C, const vector<vector<int>>& G, int new_amount) {
    if (old_C.size() <= new_amount) {
        new_amount = old_C.size();
        return old_C;
    }
    vector<unsigned long> scores;
    vector<int> new_C;
    init_D_sampling(&scores, &old_C, &new_C, G);
	cout << "added " << new_C.at(0) << endl;
	//cout << "scores\n\t";
	//print_vec(scores);
    new_amount -= 1;
	print_actual_score_sum(scores, old_C);
    //init random prob generator
    uniform_real_distribution<> rd_prob(0, 1);
    for (int i = 0; i < new_amount; ++i) {
        sample_next_C(&scores, &old_C, &new_C, G, rd_prob(random_engine));
		cout << "added " << new_C.back() << endl;
		//cout << "scores\n\t";
		//print_vec(scores);
		//print_actual_score_sum(scores, old_C);
    }

    return new_C;
}



double fullMatrix_compute_score_sum(const vector<double>& scores, const vector<int>& old_C) {
    double score_sum = 0;
    for (int c : old_C) {
        score_sum += scores.at(c);
    }
    if (score_sum < 0) {
        cout << "negative score_sum!!!" << endl;
    }
    return score_sum;
}


void fullMatrix_init_D_sampling(vector<double>* scores, vector<int>* old_C, vector<int>* new_C,
                                const vector<vector<double>>& dAtoC) {
    //init random generator
    uniform_int_distribution<> rd_int(0, (*old_C).size() - 1);

    // Chose first new_c randomly and uniformly for all old_c
    int rd_index = rd_int(random_engine);
    (*new_C).push_back((*old_C).at(rd_index));
    (*old_C).erase((*old_C).begin() + rd_index);
    //compute distances between old_C and new_C
    vector<double> distances = dAtoC.at((*new_C).at(0));
    for (int i = 0; i < distances.size(); ++i) {
        if (i == rd_index) {
            (*scores).push_back(0);
        }
        else {
            (*scores).push_back(distances.at(i));
        }
    }
}


void fullMatrix_sample_next_C(vector<double>* scores, vector<int>* old_C, vector<int>* new_C,
                              const vector<vector<double>>& dAtoC, double rd) {

    auto prob = fullMatrix_compute_score_sum(*scores, *old_C) * rd;
    double current = 0;
    int to_add_index = 0;
    for (int index = 0; index < (*old_C).size(); ++index) {
        current += (*scores).at((*old_C).at(index));
        if (prob < current) {
            to_add_index = index;
            break;
        }
    }
    //cout << "dAtoC.size(): " << dAtoC.size() << endl;
    int to_add_c = (*old_C).at(to_add_index);

    vector<double> distances = dAtoC.at(to_add_c);
    (*old_C).erase((*old_C).begin() + to_add_index);

    for (int old_c : (*old_C)) {
        // distances.at(old_c) <=> distance between old_c and to_add_c
        if ((*scores).at(old_c) > distances.at(old_c)) {
            (*scores)[old_c] = distances.at(old_c);
        }
        if ((*scores).at(old_c) < 0) {
            cout << "negative score!!!" << endl;
            exit(-1);
        }
    }
    (*new_C).push_back(to_add_c);
}


vector<int> fullMatrix_D_sampling(vector<int> old_C, const vector<vector<double>>& dAtoC, int new_amount) {
    if (old_C.size() < new_amount) {
        return old_C;
    }
    vector<double> scores;
    vector<int> new_C;
    //cout << "init D_sampling" << endl;
    fullMatrix_init_D_sampling(&scores, &old_C, &new_C, dAtoC);
    //cout << "added " << new_C.at(0) << endl;
    //cout << "scores\n\t";
    //print_vec(scores);
    new_amount -= 1;
    //print_actual_score_sum(scores, old_C);
    //init random prob generator
    uniform_real_distribution<> rd_prob(0, 1);
    for (int i = 0; i < new_amount; ++i) {
        fullMatrix_sample_next_C(&scores, &old_C, &new_C, dAtoC, rd_prob(random_engine));
        //cout << "added " << new_C.back() << endl;
        //cout << "scores\n\t";
        //print_vec(scores);
        //print_actual_score_sum(scores, old_C);
    }

    return new_C;
}