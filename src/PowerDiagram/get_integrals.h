#pragma once

#include "SpaceFunctions/Constant.h"
#include "FunctionEnum.h"
#include <vector>

namespace PowerDiagram {

/**
   We assume that grid has already been initialized by diracs
*/
template<class TF,class Grid,class Bounds,class Pt,class Func>
void get_integrals( TF *res, Grid &grid, Bounds &bounds, const Pt *positions, const TF *weights, std::size_t nb_diracs, const Func &func ) {
    grid.for_each_laguerre_cell( [&]( auto &lc, auto num_dirac ) {
        TF measure = 0;
        bounds.for_each_intersection( lc, [&]( auto &cp, SpaceFunctions::Constant<TF> space_func ) {
            measure += space_func.coeff * cp.measure( func );
        } );
        res[ num_dirac ] = measure;
    }, bounds.englobing_convex_polyhedron(), positions, weights, nb_diracs );
}

template<class TF,class Grid,class Bounds,class Pt>
void get_integrals( TF *res, Grid &grid, Bounds &bounds, const Pt *positions, const TF *weights, std::size_t nb_diracs ) {
    get_integrals( res, grid, bounds, positions, weights, nb_diracs, FunctionEnum::Unit() );
}

} // namespace PowerDiagram
