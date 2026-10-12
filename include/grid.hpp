#ifndef GRID_H
#define GRID_H

class Grid {
    public:
        Grid(double L, double delta_x);
        double get_L() const;
        double get_delta_x() const;
        int get_N() const;
        double position(int i) const;
        /* Maps any integer index onto [0, N-1] by wrapping around the grid,
         * so stencils can circle back (as we have periodic boundary conditions)
         */
        int wrap(int i) const;

    private:
        double L_;
        double delta_x_;
        int N_;
};


#endif
