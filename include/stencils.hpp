#ifndef STENCILS_H
#define STENCILS_H
#include "grid.hpp"
#include <vector>

namespace Stencils {
     // We use an output vector to avoid the overhead in returning a vector each time
    void second_order(const std::vector<double>& field, std::vector<double>& out, const Grid& grid);
    void fourth_order(const std::vector<double>& field, std::vector<double>& out, const Grid& grid);
};



#endif
