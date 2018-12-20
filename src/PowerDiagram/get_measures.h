#pragma once

#include "SpaceFunctions/Constant.h"
#include <vector>

namespace PowerDiagram {

/**
   We assume that grid has already been initialized by diracs
*/
template<class TF,class Grid,class Bounds,class Diracs>
void get_measures( std::vector<TF> &measures, Grid &grid, Bounds &bounds, const Diracs &diracs ) {
    measures.resize( diracs.size() );
    grid.for_each_laguerre_cell( [&]( auto &lc, auto num_dirac ) {
        TF measure = 0;
        bounds.for_each_intersection( lc, [&]( auto &cp, SpaceFunctions::Constant<TF> space_func ) {
            measure += space_func.coeff * cp.measure();
        } );
        measures[ num_dirac ] = measure;
    }, bounds.englobing_convex_polyhedron(), diracs );
}

} // namespace PowerDiagram
