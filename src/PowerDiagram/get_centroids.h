#pragma once

#include "SpaceFunctions/Constant.h"
#include "system/Assert.h"
#include "FunctionEnum.h"
#include <vector>

namespace PowerDiagram {

/**
   We assume that grid has already been initialized by diracs
*/
template<class Grid,class Bounds,class Pt,class TF,class TI,class Func,class CB>
void get_centroids( Grid &grid, Bounds &bounds, const Pt *positions, const TF *weights, TI nb_diracs, const Func &radial_func, const CB &cb ) {
    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0, int ) {
        TF mass = 0;
        Pt centroid;
        for( std::size_t d = 0; d < Grid::dim; ++d )
            centroid[ d ] = 0;
        bounds.for_each_intersection( lc, [&]( auto &cp, SpaceFunctions::Constant<TF> ) {
            cp.add_centroid_contrib( centroid, mass, FunctionEnum::func_for_final_cp_integration( radial_func ) );
        } );
        cb( centroid / TF( mass + ( mass == 0 ) ), mass, num_dirac_0 );
    }, bounds.englobing_convex_polyhedron(), positions, weights, nb_diracs, false, need_ball_cut( radial_func ) );
}

} // namespace PowerDiagram
