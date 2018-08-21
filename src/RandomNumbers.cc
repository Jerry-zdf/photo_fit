#include <ctime>
#include <random>
#include "RandomNumbers.h"

std::mt19937 &getRandomEngine() {
  static std::mt19937
    engine(static_cast<typename std::mt19937::result_type>(std::time(0)));
  return engine;
}
