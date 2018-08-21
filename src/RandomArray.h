#ifndef RANDOMARRAY_H
#define RANDOMARRAY_H

#include "RandomNumbers.h"
#include <array>


template<uint size>
auto &generateRandomArray(const std::array<double, size> &minArray,
                          const std::array<double, size> &maxArray,) {

    std::array<double, size> randomArray;
    auto engine = getRandomEngine();
    for (int i = 0; i < size; ++i) {
      randomArray[i] = std::uniform_real_distribution<double>(minArray[i], maxArray[i])(engine);
    }
    return randomArray;
}


#endif // RANDOMARRAY_H
