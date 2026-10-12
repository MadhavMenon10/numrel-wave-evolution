#ifndef WAVE_SOLVER_H
#define WAVE_SOLVER_H
#include "grid.hpp"
#include <vector>
#include <functional>

class WaveSolver {
    public:
        WaveSolver(Grid grid, double delta_t, int stencil_order);
        const std::vector<double>& get_field() const;
        const std::vector<double>& get_momentum() const;
        double get_t() const;
        void set_initial_conditions(const std::function<double(double)>& init_phi, const std::function<double(double)>& init_phit0);

    private:
        Grid grid_;
        std::vector<double> phi_;
        std::vector<double> pi_;
        double t_;
        double delta_t_;
        int stencil_order_;
};

#endif
