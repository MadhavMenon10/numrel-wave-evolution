#include "stencils.hpp"
#include <stdexcept>

void Stencils::second_order(const std::vector<double>& field, std::vector<double>& out, const Grid& grid) {
    const int N {grid.get_N()};
    if (static_cast<int>(field.size()) != N) {
        throw std::runtime_error("Field does not have the right size");
    }
    out.resize(N);
    const double delta_x {grid.get_delta_x()};
    const double delta_x_squared {delta_x * delta_x};
    for (int i {}; i < N; ++i) {
        out[i] = (field[grid.wrap(i + 1)] + field[grid.wrap(i - 1)] - 2 * field[i]) / (delta_x_squared);
    }
}

void Stencils::fourth_order(const std::vector<double>& field, std::vector<double>& out, const Grid& grid) {
    const int N = grid.get_N();
    if (static_cast<int>(field.size()) != N) {
        throw std::runtime_error("Field does not have the right size");
    }
    const double delta_x = grid.get_delta_x();
    const double delta_x_squared {delta_x * delta_x};
    out.resize(N);
    for (int i {}; i < N; ++i) {
        out[i] = (-field[grid.wrap(i + 2)] + 16.0 * field[grid.wrap(i + 1)] - 30.0 * field[i] + 16.0 * field[grid.wrap(i - 1)] - field[grid.wrap(i - 2)]) / (12.0 * delta_x_squared);
    }
}
