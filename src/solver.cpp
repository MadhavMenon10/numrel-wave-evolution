#include "solver.hpp"
#include "stencils.hpp"
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
    int N {grid_.get_N()};
    for (int i {}; i < N; ++i) {
        auto x_i {grid_.position(i)};
        phi_[i] = phi_0(x_i);
        pi_[i] = -phi_t0(x_i);
    }
    t_ = 0;
}

void WaveSolver::compute_derivatives(const std::vector<double>& phi, const std::vector<double>& pi, std::vector<double>& dphi_dt, std::vector<double>& dpi_dt, int stencil_order) const {
    int N {grid_.get_N()};
    dphi_dt.resize(N);
    dpi_dt.resize(N);
    for (int i {}; i < N; ++i) {
        dphi_dt[i] = -pi[i];
    }
    if (stencil_order == 2) {
        Stencils::second_order(phi, dpi_dt, grid_);
    } else if (stencil_order == 4) {
        Stencils::fourth_order(phi, dpi_dt, grid_);
    } else {
        throw std::runtime_error("Invalid stencil order");
    }
    for (int i {}; i < N; ++i) {
        dpi_dt[i] *= -1; // The stencil call doesn't negate. We have dpi/dt = -d^2phi/dx^2
    }
}

void WaveSolver::step(int stencil_order) {
    int N {grid_.get_N()};
    compute_derivatives(phi_, pi_, k1_.phi, k1_.pi, stencil_order);
    for (int i {}; i < N; ++i) {
        temp_.phi[i] = phi_[i] + delta_t_/2 * k1_.phi[i];
        temp_.pi[i] = pi_[i] + delta_t_/2 * k1_.pi[i];
    }
    compute_derivatives(temp_.phi, temp_.pi, k2_.phi, k2_.pi, stencil_order);
    for (int i {}; i < N; ++i) {
        temp_.phi[i] = phi_[i] + delta_t_/2 * k2_.phi[i];
        temp_.pi[i] = pi_[i] + delta_t_/2 * k2_.pi[i];
    }
    compute_derivatives(temp_.phi, temp_.pi, k3_.phi, k3_.pi, stencil_order);
    for (int i {}; i < N; ++i) {
        temp_.phi[i] = phi_[i] + delta_t_ * k3_.phi[i];
        temp_.pi[i] = pi_[i] + delta_t_ * k3_.pi[i];
    }
    compute_derivatives(temp_.phi, temp_.pi, k4_.phi, k4_.pi, stencil_order);
    for (int i {}; i < N; ++i) {
        phi_[i] = phi_[i] + delta_t_/6 * (k1_.phi[i] + 2*k2_.phi[i] + 2*k3_.phi[i] + k4_.phi[i]);
        pi_[i] = pi_[i] + delta_t_/6 * (k1_.pi[i] + 2*k2_.pi[i] + 2*k3_.pi[i] + k4_.pi[i]);
    }
    t_ += delta_t_;
 }
