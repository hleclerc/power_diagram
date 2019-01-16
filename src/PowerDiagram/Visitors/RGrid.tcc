#include "../system/StaticRange.h"
#include "../system/Span.h"
#include "IntSpiral.h"
#include "RGrid.h"
#include <cmath>

namespace PowerDiagram {
namespace Visitor {

template<class Pc>
RGrid<Pc>::RGrid( TF cell_size ) : cell_size( cell_size ) {
    inv_cell_size = TF( 1 ) / cell_size;
    ball_cut = allow_ball_cut;
}

template<class Pc>
typename RGrid<Pc>::TI RGrid<Pc>::index_cell( Pt pos ) const {
    TI index = 0;
    for( std::size_t d = 0; d < dim; ++d )
        index += acc_grid_lengths[ d ] * TI( inv_cell_size * ( pos[ d ] - min_point[ d ] ) );
    return index;
}

template<class Pc> template<class D>
void RGrid<Pc>::init( const D &diracs ) {
    using std::ceil;
    using std::min;
    using std::max;

    for( std::size_t d = 0; d < dim; ++d ) {
        min_point[ d ] = + std::numeric_limits<TF>::max();
        max_point[ d ] = - std::numeric_limits<TF>::max();
    }
    min_weight = + std::numeric_limits<TF>::max();
    max_weight = - std::numeric_limits<TF>::max();

    for( const auto &dirac : diracs ) {
        for( std::size_t d = 0; d < dim; ++d ) {
            min_point[ d ] = min( min_point[ d ], dirac.pos[ d ] );
            max_point[ d ] = max( max_point[ d ], dirac.pos[ d ] );
        }
        min_weight = min( min_weight, dirac.weight );
        max_weight = max( max_weight, dirac.weight );
    }
    max_weight *= 1 + std::numeric_limits<TF>::epsilon();

    //
    acc_grid_lengths[ 0 ] = 1;
    for( std::size_t d = 0; d < dim; ++d ) {
        grid_lengths[ d ] = ceil( inv_cell_size * ( max_point[ d ] - min_point[ d ] ) * ( 1 + std::numeric_limits<TF>::epsilon() ) );
        acc_grid_lengths[ d + 1 ] = acc_grid_lengths[ d ] * grid_lengths[ d ];
    }

    // set up the grid: get nb dirac per cell
    di_offsets.resize( acc_grid_lengths[ dim ] + 1 );
    for( std::size_t i = 0; i < di_offsets.size(); ++i )
        di_offsets[ i ] = 0;
    for( std::size_t i = 0; i < diracs.size(); ++i )
        ++di_offsets[ index_cell( diracs[ i ].pos ) ];

    // suffix scan
    for( std::size_t i = 0, acc = 0; i < di_offsets.size(); ++i ) {
        TI old = acc;
        acc += di_offsets[ i ];
        di_offsets[ i ] = old;
    }

    // get indices
    di_values.resize( di_offsets.back() );
    for( std::size_t i = 0; i < diracs.size(); ++i )
        di_values[ di_offsets[ index_cell( diracs[ i ].pos ) ]++ ] = i;

    // shift ng_offsets (to get the suffix scan again)
    if ( di_offsets.size() ) {
        for( TI i = di_offsets.size(); --i; )
            di_offsets[ i ] = di_offsets[ i - 1 ];
        di_offsets[ 0 ] = 0;
    }
}

template<class Pc> template<class D>
void RGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, std::size_t num )> &cb, const CP &starting_lc, const D &diracs ) {
    for_each_laguerre_cell( [&]( CP &cp, std::size_t num, int ) {
        cb( cp, num );
    }, starting_lc, diracs );
}

template<class Pc> template<class D>
void RGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, std::size_t num, int num_thread )> &cb, const CP &starting_lc, const D &diracs ) {
    using std::sqrt;
    using std::ceil;

    if ( ball_cut == 0 )
        TODO; // for now, RGrid is mainly developed for the ball_cut case

    auto plane_cut = [&]( CP &lc, TI i0, TI i1 ) {
        const auto &d0 = diracs[ i0 ];
        const auto &d1 = diracs[ i1 ];
        Pt V = d1.pos - d0.pos;
        TF d = norm_2( V );
        TF i = 1.0 / d;
        Pt N = i * V;
        TF x = 0.5 * ( d + i * ( d0.weight - d1.weight ) );
        lc.plane_cut( d0.pos + x * N, N, i1 );
    };

    // vectors for stuff that will be reused inside the execution threads
    int nb_threads = thread_pool.nb_threads(), nb_jobs = 4 * nb_threads;

    // for each item
    thread_pool.execute( nb_jobs, [&]( std::size_t num_job, int num_thread ) {
        CP lc;

        TI beg_cell = ( num_job + 0 ) * nb_cells() / nb_jobs;
        TI end_cell = ( num_job + 1 ) * nb_cells() / nb_jobs;
        for( TI num_cell = beg_cell; num_cell < end_cell; ++num_cell ) {
            Pi cell_coords;
            for( TI d = 0, c = num_cell; d < dim; ++d ) {
                cell_coords[ d ] = c % grid_lengths[ d ];
                c /= grid_lengths[ d ];
            }

            for( TI num_dirac : Span<TI>{ di_values.data() + di_offsets[ num_cell + 0 ], di_values.data() + di_offsets[ num_cell + 1 ] } ) {
                if ( ball_cut == 0 )
                    TODO;
                lc = starting_lc;

                // diracs in the current cell
                for( TI num_cr_dirac : Span<TI>{ di_values.data() + di_offsets[ num_cell + 0 ], di_values.data() + di_offsets[ num_cell + 1 ] } )
                    if ( num_cr_dirac != num_dirac )
                        plane_cut( lc, num_dirac, num_cr_dirac );

                // spiral
                IntSpiral<dim> is;
                TF max_dist = sqrt( diracs[ num_dirac ].weight ) + sqrt( max_weight );
                is.for_each_until( ceil( inv_cell_size * max_dist ), [&]( auto offset ) {
                    TI cr_cell_index = 0;
                    for( TI d = 0; d < dim; ++d ) {
                        static_assert( std::is_signed<TI>::value == false, "" );
                        TI cr_cell_coord = cell_coords[ d ] + offset[ d ];
                        if ( cr_cell_coord >= grid_lengths[ d ] )
                            return;
                        cr_cell_index += acc_grid_lengths[ d ] * cr_cell_coord;
                    }

                    // cuts
                    for( TI num_cr_dirac : Span<TI>{ di_values.data() + di_offsets[ cr_cell_index + 0 ], di_values.data() + di_offsets[ cr_cell_index + 1 ] } )
                        plane_cut( lc, num_dirac, num_cr_dirac );
                } );

                //
                if ( ball_cut )
                    lc.ball_cut( diracs[ num_dirac ].pos, sqrt( diracs[ num_dirac ].weight ), num_dirac );
                else
                    lc.sphere_center = positions[ num_dirac ];

                //
                cb( lc, num_dirac, num_thread );
            }
        }
    } );
}


template<class Pc> template<class V>
void RGrid<Pc>::display( V &vtk_output ) const {
    TODO;
    //    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
    //        const Grid &grid = grids[ num_grid ];

    //        for( TI num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
    //            Pt p;
    //            for( int d = 0; d < dim; ++d )
    //                p[ d ] = grid.cells[ num_cell ].pos[ d ];

    //            TF a = 0, b = grid.cells[ num_cell ].size;
    //            switch ( dim ) {
    //            case 2:
    //                vtk_output.add_lines( {
    //                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
    //                    Point2<TF>{ p[ 0 ] + b, p[ 1 ] + a },
    //                    Point2<TF>{ p[ 0 ] + b, p[ 1 ] + b },
    //                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + b },
    //                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
    //                }, { TF( num_grid ) } );
    //                break;
    //            case 3:
    //                TODO;
    //                break;
    //            default:
    //                TODO;
    //            }
    //        }
    //    }
}

} // namespace Visitor
} // namespace PowerDiagram
