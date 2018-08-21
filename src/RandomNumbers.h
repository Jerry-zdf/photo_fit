#ifndef RANDOM_NUMBERS_H_
#define RANDOM_NUMBERS_H_

#include <random>

/**
 * This function returns Mersenne Twister random number engine initialized with
 * time of first invocation.
 */
std::mt19937 &getRandomEngine();

#endif
