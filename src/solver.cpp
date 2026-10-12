#include "solver.hpp"
#include <stdexcept>

WaveSolver::WaveSolver(Grid grid, double delta_t, int stencil_order): grid_(grid), phi_(grid_.get_N()), pi_(grid_.get_N()), t_(0), delta_t_(delta_t), stencil_order_(stencil_order) {
    if (stencil_order != 2 && stencil_order != 4) {
        throw std::invalid_argument("Invalid stencil order passed");
    }
}

const std::vector<double>& WaveSolver::get_field() const {
    return phi_;
}

const std::vector<double>& WaveSolver::get_momentum() const {
    return pi_;
}

double WaveSolver::get_t() const {
    return t_;
}

void WaveSolver::set_initial_conditions(const std::function<double(double)>& phi_0, const std::function<double(double)>& phi_t0) {
    int N = grid_.get_N();
    for (int i {}; i < N; ++i) {
        auto x_i {grid_.position(i)};
        phi_[i] = phi_0(x_i);
        pi_[i] = -phi_t0(x_i);
    }
    t_ = 0;
}


