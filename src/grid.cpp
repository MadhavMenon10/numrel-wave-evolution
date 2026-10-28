#include "grid.hpp"
#include <cmath>
#include <stdexcept>

Grid::Grid(double L, double delta_x): L_(L), delta_x_(delta_x) {
    if (L_ <= 0.0 || delta_x_ <= 0.0) {
        throw std::runtime_error("L and/or delta_x must be positive");
    }
    N_ = static_cast<int>(std::lround(L_ / delta_x_));
    if (std::abs(N_ * delta_x_ - L_) > 0.000000001) {
        throw std::runtime_error("Division wasn't accurate enough");
    }
}

double Grid::get_L() const {
    return L_;
}

double Grid::get_delta_x() const {
    return delta_x_;
}

int Grid::get_N() const {
    return N_;
}

double Grid::position(int i) const {
    return i * delta_x_;
}

int Grid::wrap(int i) const {
    return ((i % N_) + N_) % N_; // The first modulo moves i between [-N_, N_], adding N then moves that to [0, 2N_ - 1], the final modulo brings that to [0, N_ - 1]
}
