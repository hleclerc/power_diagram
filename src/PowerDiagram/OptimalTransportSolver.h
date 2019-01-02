#pragma once

#include "get_der_integrals_wrt_weights.h"
#include "get_centroids.h"

namespace PowerDiagram {

/**
*/
template<class Grid,class Bounds>
class OptimalTransportSolver {
public:
    using           CP                    = typename Grid::CP; ///< convex polyhedron
    using           Pt                    = typename Grid::Pt; ///< convex polyhedron
    using           TF                    = typename Grid::TF;
    using           TI                    = typename Grid::TI;

    /* */           OptimalTransportSolver( Grid *grid, Bounds *bounds );

    template        <class VO>
    void            display               ( VO& vtk_output, const Pt *positions, const TF *weights, TI nb_diracs ); ///< result in `new_weights`
    void            solve                 ( const Pt *positions, TF *weights, TI nb_diracs ); ///< result in `new_weights`

    // input parameters
    std::size_t     max_nb_iter;

    Bounds&         bounds;
    Grid&           grid;

    // by products
    std::vector<TF> old_weights;
    std::vector<TI> m_offsets;
    std::vector<TI> m_columns;
    std::vector<TF> m_values;
    std::vector<TF> v_values;
    std::vector<TF> dw;
};

} // namespace PowerDiagram

#include "OptimalTransportSolver.tcc"

