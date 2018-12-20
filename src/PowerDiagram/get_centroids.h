#pragma once

#include "SpaceFunctions/Constant.h"
#include "system/Assert.h"
#include "FunctionEnum.h"
#include <vector>

namespace PowerDiagram {

/**
   We assume that grid has already been initialized by diracs
*/
template<class Grid,class Bounds,class Diracs,class CB>
void get_centroids( Grid &grid, Bounds &bounds, const Diracs &diracs, const CB &cb ) {
    using Pt = typename Grid::Pt;
    using TF = typename Grid::TF;

    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
        TF mass = 0;
        Pt centroid;
        for( std::size_t d = 0; d < Grid::dim; ++d )
            centroid[ d ] = 0;
        bounds.for_each_intersection( lc, [&]( auto &cp, SpaceFunctions::Constant<TF> ) {
            cp.add_centroid_contrib( centroid, mass );
        } );
        cb( centroid / TF( mass + ( mass == 0 ) ), mass, num_dirac_0 );
    }, bounds.englobing_convex_polyhedron(), diracs );
}

} // namespace PowerDiagram
