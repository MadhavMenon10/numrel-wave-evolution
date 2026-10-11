#include "grid.hpp"
#include <cmath>
#include <stdexcept>

Grid::Grid(double L, double dx): L_(L), dx_(dx) {
    N_ = std::lround(L_ / dx_);
    if ((N_ * dx_ - L_) < 0.000000001) {
        throw std::runtime_error("Division wasn't accurate enough")
    }
}

double Grid::get_L() const {
    return L_;
}

double Grid::get_dx() const {
    return dx_;
}

int Grid::get_N() const {
    return N_;
}

double Grid::position(int i) const {
    return i * dx_;
}

int Grid::wrap(int i) const {
    return ((i % N_) + N_) % N_; // The first modulo moves i between [-N_, N_], adding N then moves that to [0, 2N_ - 1], the final modulo brings that to [0, N_ - 1]
}
