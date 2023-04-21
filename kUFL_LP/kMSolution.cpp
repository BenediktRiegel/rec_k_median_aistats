#include "kMSolution.h"
#include <vector>

// Date member function
void kMSolution::set(std::vector<int> solution, double service_cost, double other_cost) {
    this->solution = solution;
    this->service_cost = service_cost;
    this->other_cost = other_cost;
}
