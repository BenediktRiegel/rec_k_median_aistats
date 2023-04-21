#ifndef BAALGO_KMSOLUTION_H
#define BAALGO_KMSOLUTION_H
#include <vector>
#include "fixedDouble.h"

class kMSolution {
public:
    std::vector<int> solution;
    fixedDouble service_cost;
    fixedDouble other_cost;

    void set(std::vector<int> solution, double service_cost, double other_cost);

    fixedDouble cost() {
        return service_cost + other_cost;
    }
};


#endif //BAALGO_KMSOLUTION_H
