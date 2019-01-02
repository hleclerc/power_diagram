#include "get_der_integrals_wrt_weights.h"
#include "OptimalTransportSolver.h"
#include "EigenSolver.h"
#include "system/Tick.h"

namespace PowerDiagram {

template<class Grid, class Bounds>
OptimalTransportSolver<Grid, Bounds>::OptimalTransportSolver(Grid *grid, Bounds *bounds) : bounds( *bounds ), grid( *grid ) {
    max_nb_iter = 100;
}

template<class Grid, class Bounds>
void OptimalTransportSolver<Grid, Bounds>::solve( const Pt *positions, TF *weights, TI nb_diracs ) {
    using std::max;
    using std::abs;

    old_weights.resize( nb_diracs );
    for( std::size_t i = 0; i < nb_diracs; ++i )
        old_weights[ i ] = weights[ i ];

    for( std::size_t num_iter = 0; num_iter < max_nb_iter; ++num_iter ) {
        // grid
        tick << "grid";
        grid.update( positions, weights, nb_diracs, num_iter == 0, true );
        tick >> "grid";

        // der
        tick << "der";
        auto t0 = tick.get_time();
        int error = get_der_integrals_wrt_weights( m_offsets, m_columns, m_values, v_values, grid, bounds, positions, weights, nb_diracs );
        P( tick.elapsed_since( t0 ) );
        tick >> "der";
        m_values[ 0 ] *= 2;

        // go back if pb
        if ( error ) {
            TF ratio = 0.1;
            for( std::size_t i = 0; i < nb_diracs; ++i )
                weights[ i ] = ( 1 - ratio ) * old_weights[ i ] + ratio * weights[ i ];
            continue;
        }
        for( std::size_t i = 0; i < nb_diracs; ++i )
            old_weights[ i ] = weights[ i ];

        // solve
        EigenSolver es;
        tick << "solve";
        es.solve( dw, m_offsets, m_columns, m_values, v_values );
        tick >> "solve";

        TF mdw = 0;
        for( std::size_t i = 0; i < nb_diracs; ++i ) {
            mdw = max( mdw, abs( dw[ i ] ) );
            weights[ i ] -= dw[ i ];
        }

        P( mdw );
        if ( mdw < 1e-6 )
            break;
    }
}

template<class Grid, class Bounds> template<class VO>
void OptimalTransportSolver<Grid,Bounds>::display( VO &vtk_output, const Pt *positions, const TF *weights, TI nb_diracs ) {
    grid.update( positions, weights, nb_diracs );

    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
        lc.display( vtk_output, { weights[ num_dirac_0 ] } );
    }, bounds.englobing_convex_polyhedron(), positions, weights, nb_diracs );
}

//    template<class Diracs>
//    TF volume( const Diracs &diracs ) {
//        grid.init( diracs );

//        std::vector<TF> volumes( diracs.size() );
//        grid.for_each_laguerre_cell( [&]( CP &lc, std::size_t num ) {
//            bounds.for_each_intersection( lc, [&]( CP &cp, auto space_func ) {
//                volumes[ num ] += cp.integration( space_func );
//            } );
//        }, bounds.englobing_convex_polyhedron(), diracs );

//        TF res = 0;
//        for( auto v : volumes )
//            res += v;
//        return res;
//    }

} // namespace PowerDiagram
