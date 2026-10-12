#ifndef WAVE_SOLVER_H
#define WAVE_SOLVER_H
#include "grid.hpp"
#include <vector>
#include <functional>


// This is used in our step function 
struct State {
    std::vector<double> phi;
    std::vector<double> pi;
    State(int N): phi(N), pi(N) {}
};

class WaveSolver {
    public:
        WaveSolver(Grid grid, double delta_t, int stencil_order);
        const std::vector<double>& get_field() const;
        const std::vector<double>& get_momentum() const;
        double get_t() const;
        void set_initial_conditions(const std::function<double(double)>& init_phi, const std::function<double(double)>& init_phit0);
        void compute_derivatives(const std::vector<double>& phi, const std::vector<double>& pi, std::vector<double>& dphi_dt, std::vector<double>& dpi_dt, int stencil_order) const;
        void step(int stencil_order);

    private:
        Grid grid_;
        std::vector<double> phi_;
        std::vector<double> pi_;
        double t_;
        double delta_t_;
        int stencil_order_;
        State k1_;
        State k2_;
        State k3_;
        State k4_;
        State temp_;
};

#endif
