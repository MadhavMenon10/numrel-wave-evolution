#include "stencils.hpp"

void Stencils::second_order(const std::vector<double>& field, std::vector<double>& out, const Grid& grid) {
    int N {grid.get_N()};
    if (static_cast<int>(field.size()) != N) {
        throw std::runtime_error("Field does not have the right size");
    }
    out.resize(N);
    double delta_x {grid.get_delta_x()};
    double delta_x_squared {delta_x * delta_x};
    for (int i {}; i < N; ++i) {
        out[i] = (field[grid.wrap(i + 1)] + field[grid.wrap(i - 1)] - 2 * field[i]) / (delta_x_squared);
    }
}
