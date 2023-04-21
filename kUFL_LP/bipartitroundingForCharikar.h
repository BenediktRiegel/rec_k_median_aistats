#ifndef BAALGO_BIPARTITROUNDINGFORCHARIKAR_H
#define BAALGO_BIPARTITROUNDINGFORCHARIKAR_H
#include "roundgraph.h"

class bipartitroundingForCharikar {
public:
    static map<int, map<int, int>> solve(roundgraph* rgraph);
};


#endif //BAALGO_BIPARTITROUNDINGFORCHARIKAR_H
