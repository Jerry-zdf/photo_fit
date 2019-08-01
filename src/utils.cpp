#include "utils.h"

std::chrono::duration<double> Clock::restart() {
    const auto dur = duration();
    _start         = std::chrono::system_clock::now();
    return dur;
}

std::chrono::duration<double> Clock::duration() const {
    const auto end = std::chrono::system_clock::now();
    return end - _start;
}

std::ostream& operator<<(std::ostream& os, const Clock& rhs) {
    os << rhs.duration().count() << " s";
    return os;
}

std::mt19937& get_random_engine() noexcept {
    static std::random_device rd;
    static std::mt19937 engine{rd()};
    return engine;
}