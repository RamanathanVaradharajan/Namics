#pragma once

#include <iostream>
#include <random>
#include "../../namics.h"

using namespace std;

class Random {
public:
    Random();

    explicit Random(int seed);

    int getInt(int min, int max);

    Real getReal(int min, int max);

    int getIntExcludeArray(int min, int max, const vector<int>& exclude_values);


private:
    mt19937 engine;
};