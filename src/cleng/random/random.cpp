#include "random.h"

Random::Random() {
    random_device rd;
    engine = mt19937(rd());
}

Random::Random(int seed) {
    engine = mt19937(seed);
}

int Random::getInt(int min, int max) {
    return uniform_int_distribution<int>(min, max)(engine);
}

Real Random::getReal(int min, int max) {
    return uniform_real_distribution<Real>(min, max)(engine);
}

int Random::getIntExcludeArray(int min, int max, const vector<int>& exclude_values) {
    int out = getInt(min, max);

    int coincidence = 0;
    for (auto exclude_value : exclude_values) if (out == exclude_value) coincidence++;

    while (coincidence) {
        out = getInt(min, max);
        coincidence = 0;
        for (auto exclude_value : exclude_values) if (out == exclude_value) coincidence++;
    }

    return out;
}