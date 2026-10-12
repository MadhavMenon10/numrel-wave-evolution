#include "grid.hpp"
#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>

int main(int argc, char* argv[]) {
    try {
        if (argc != 7) {
            std::cerr << "delta_x, delta_t, stencil_order, t_final, interval, output_path";
            return 1;
        }
        double delta_x {std::stod(argv[1])};
        double delta_t {std::stod(argv[2])};
        int stencil_order {std::stoi(argv[3])};
        double t_final {std::stod(argv[4])};
        double save_interval {std::stod(argv[5])};
        std::string output_path {argv[6]};
        if (delta_t <= 0.0 || save_interval <= 0.0) {
            throw std::invalid_argument("Delta_t or save_interval must be positive");
        }
        const double L {20.0};
        int num_steps {static_cast<int>(std::lround(t_final / delta_t))};
        int save_every {std::max(static_cast<int>(std::lround(save_interval / delta_t)), 1)};
        if (!(std::abs(num_steps * delta_t - t_final) < 0.0000000001)) {
            throw std::invalid_argument("Values don't match");
        }
        auto phi_0 = [](double x) {return std::exp(-((x-10)*(x-10)));};
        auto phi_t0 = [](double x) {return 0.0;};
        Grid grid(L, delta_x);
        WaveSolver solver(grid, delta_t, stencil_order);
        solver.set_initial_conditions(phi_0, phi_t0);
        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("Cannot open output file: " + output_path);
        }
        out << std::setprecision(17);
        auto write_row = [&]() {
            const auto& phi = solver.get_field();
            out << solver.get_t();
            for (double value : phi) {
                out << ',' << value;
            }
            out << '\n';
        };
        write_row();
        for (int i {}; i < num_steps; ++i) {
            solver.step(stencil_order);
            if (i % save_every == 0 || i == num_steps) {
                write_row();
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
