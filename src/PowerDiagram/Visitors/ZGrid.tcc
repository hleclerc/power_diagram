#include "../system/StaticRange.h"
#include "../system/RadixSort.h"
#include "../system/Span.h"
#include "FrontZgrid.h"
#include "ZGrid.h"
#include <cmath>

// #define DISPLAY_nb_explored_cells

extern const std::uint32_t morton_256_2D_x[ 256 ];
extern const std::uint32_t morton_256_2D_y[ 256 ];
extern const std::uint32_t morton_256_3D_x[ 256 ];
extern const std::uint32_t morton_256_3D_y[ 256 ];
extern const std::uint32_t morton_256_3D_z[ 256 ];

namespace PowerDiagram {
namespace Visitor {


template<class Pc>
ZGrid<Pc>::ZGrid( std::size_t max_diracs_per_cell, TF max_delta_weight_per_grid ) : max_delta_weight_per_grid( max_delta_weight_per_grid ), max_diracs_per_cell( max_diracs_per_cell ) {
    ball_cut = allow_ball_cut;
}

template<class Pc>
void ZGrid<Pc>::update( const Pt *positions, const TF *weights, std::size_t nb_diracs, bool positions_have_changed, bool weights_have_changed ) {
    if ( positions_have_changed || weights_have_changed ) {
        update_the_limits( positions, weights, nb_diracs );
        fill_the_grids   ( positions, weights, nb_diracs );
    }
}

template<class Pc>
bool ZGrid<Pc>::may_cut( const CP &lc, TI i0, const Grid &cr_grid, const Cell &cr_cell, const Pt *positions, const TF *weights ) {
    using std::sqrt;
    using std::max;
    using std::min;
    using std::pow;
    using std::abs;

    auto c0 = positions[ i0 ];
    auto w0 = weights  [ i0 ];

    //
    if ( ball_cut ) {
        TF md = 0; // min dist pow 2
        for( size_t d = 0; d < dim; ++d ) {
            TF o = c0[ d ] - cr_cell.pos[ d ];
            if ( o > 0 )
                o = max( TF( 0 ), o - cr_cell.size );
            md += pow( o, 2 );
        }

        return md < pow( sqrt( cr_grid.max_weight ) + sqrt( w0 ), 2 );
    }

    //
    if ( w0 > cr_grid.min_weight ) {
        TF d2 = w0 - cr_grid.min_weight;
        if ( d2 > 0 ) {
            //
            TF max_cur_dist = 0;
            for( std::size_t num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point )
                max_cur_dist = max( max_cur_dist, norm_2_p2( lc.point( num_lc_point ) ) );
            if ( max_cur_dist > 4 * d2 ) {
                TF dx = cr_cell.pos[ 0 ] - c0[ 0 ];
                TF dy = cr_cell.pos[ 1 ] - c0[ 1 ];
                TF x2 = min( pow( dx, 2 ), pow( dx + cr_cell.size, 2 ) );
                TF y2 = min( pow( dy, 2 ), pow( dy + cr_cell.size, 2 ) );
                if ( x2 + y2 < d2 )
                    return true;
            }
        }
    }

    //
    TF l0 = cr_cell.size;
    auto test_line = [&]( Pt A, Pt B, Pt p ) {
        TF s2 = l0 * l0;
        TF s1 = dot( p - A, B - A );
        TF s0 = ( norm_2_p2( p - A ) - cr_grid.max_weight ) - ( norm_2_p2( c0 - p ) - w0 );
        if ( s1 > 0 && s1 < s2 && s2 * s0 - s1 * s1 < 0 )
            return true;
        return s0 < 0 || s2 - 2 * s1 + s0 < 0;
    };

    TF dx = cr_cell.pos[ 0 ] < c0[ 0 ] ? l0 : 0;
    TF dy = cr_cell.pos[ 1 ] < c0[ 1 ] ? l0 : 0;

    Pt PA { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ] + dy };
    Pt PB { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ] + dy };

    Pt PC { cr_cell.pos[ 0 ] + dx, cr_cell.pos[ 1 ]      };
    Pt PD { cr_cell.pos[ 0 ] + dx, cr_cell.pos[ 1 ] + l0 };

    for( std::size_t num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point ) {
        if ( test_line( PA, PB, lc.point( num_lc_point ) ) ) return true;
        if ( test_line( PC, PD, lc.point( num_lc_point ) ) ) return true;
    }

    //    Pt PA { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ]      };
    //    Pt PB { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ]      };
    //    Pt PC { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ] + l0 };
    //    Pt PD { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ] + l0 };

    //    for( std::size_t num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point ) {
    //        if ( test_line( PA, PB, lc.point( num_lc_point ) ) ) return true;
    //        if ( test_line( PB, PC, lc.point( num_lc_point ) ) ) return true;
    //        if ( test_line( PC, PD, lc.point( num_lc_point ) ) ) return true;
    //        if ( test_line( PD, PA, lc.point( num_lc_point ) ) ) return true;
    //    }

    return false;
}

template<class Pc>
int ZGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, std::size_t num )> &cb, const CP &starting_lc, const Pt *positions, const TF *weights, std::size_t nb_diracs, bool stop_if_void_lc ) {
    return for_each_laguerre_cell( [&]( CP &cp, std::size_t num, int ) {
        cb( cp, num );
    }, starting_lc, positions, weights, nb_diracs, stop_if_void_lc );
}

template<class Pc>
int ZGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, std::size_t num, int num_thread )> &cb, const CP &starting_lc, const Pt *positions, const TF *weights, std::size_t nb_diracs, bool stop_if_void_lc ) {
    using Front = FrontZgrid<ZGrid>;
    using std::sqrt;

    auto plane_cut = [&]( CP &lc, TI i0, TI i1 ) {
        Pt V = positions[ i1 ] - positions[ i0 ];
        TF n = norm_2_p2( V );
        TF x = TF( 1 ) + ( weights[ i0 ] - weights[ i1 ] ) / n;
        TF i = TF( 1 ) / sqrt( n );
        lc.plane_cut( positions[ i0 ] + TF( 0.5 ) * x * V, i * V, i1 );
    };

    // vectors for stuff that will be reused inside the execution threads
    int nb_threads = thread_pool.nb_threads(), nb_jobs = 4 * nb_threads;
    std::vector<std::vector<std::vector<TI>>> visited( nb_threads );
    std::vector<TI> op_counts( nb_threads, 0 );

    for( int num_thread = 0; num_thread < nb_threads; ++num_thread ) {
        visited[ num_thread ].resize( grids.size() );
        for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid )
            visited[ num_thread ][ num_grid ].resize( grids[ num_grid ].cells.size(), op_counts[ num_thread ] );
    }

    #ifdef DISPLAY_nb_explored_cells
    std::vector<std::size_t> nb_explored_cells( nb_threads, 0 );
    #endif // DISPLAY_nb_explored_cells

    // for each item
    int err = 0;
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
        Grid &grid = grids[ num_grid ];

        thread_pool.execute( nb_jobs, [&]( std::size_t num_job, int num_thread ) {
            Front front( op_counts[ num_thread ], visited[ num_thread ] );
            CP lc;

            TI beg_cell = ( num_job + 0 ) * ( grid.cells.size() - 1 ) / nb_jobs;
            TI end_cell = ( num_job + 1 ) * ( grid.cells.size() - 1 ) / nb_jobs;
            for( TI num_cell = beg_cell; num_cell < end_cell && err == 0; ++num_cell ) {
                const Cell &cell = grid.cells[ num_cell ];
                for( TI num_dirac : Span<TI>{ grid.dpc_values.data() + cell.dpc_offset, grid.dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } ) {
                    // start of lc: cut with nodes in the same cell
                    lc = starting_lc;
                    for( TI num_cr_dirac : Span<TI>{ grid.dpc_values.data() + cell.dpc_offset, grid.dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } )
                        if ( num_cr_dirac != num_dirac )
                            plane_cut( lc, num_dirac, num_cr_dirac );

                    // front
                    front.init( grids, num_grid, num_cell );

                    // => neighbors in the same grid
                    for( TI num_ng_cell : Span<TI>{ grid.ng_indices.data() + grid.ng_offsets[ num_cell + 0 ], grid.ng_indices.data() + grid.ng_offsets[ num_cell + 1 ] } )
                        front.push_without_check( num_grid, num_ng_cell, grids );

                    // => items from the grids made for != weight (containing the dirac position)
                    for( std::size_t num_pa_grid = 0; num_pa_grid < grids.size(); ++num_pa_grid ) {
                        if ( num_pa_grid != num_grid ) {
                            Grid &pa_grid = grids[ num_pa_grid ];
                            TI num_pa_cell = pa_grid.cell_index_vs_dirac_number[ num_dirac ];

                            // cut with items in pa cell
                            front.set_visited( grids, num_pa_grid, num_pa_cell );
                            for( TI num_cr_dirac : Span<TI>{ pa_grid.dpc_values.data() + pa_grid.cells[ num_pa_cell + 0 ].dpc_offset, pa_grid.dpc_values.data() + pa_grid.cells[ num_pa_cell + 1 ].dpc_offset } )
                                plane_cut( lc, num_dirac, num_cr_dirac );

                            // add neighbors in the front
                            for( TI num_ng_cell : Span<TI>{ pa_grid.ng_indices.data() + pa_grid.ng_offsets[ num_pa_cell + 0 ], pa_grid.ng_indices.data() + pa_grid.ng_offsets[ num_pa_cell + 1 ] } )
                                front.push_without_check( num_pa_grid, num_ng_cell, grids );
                        }
                    }

                    // neighbors
                    while( ! front.empty() ) {
                        typename Front::Item cr = front.pop();
                        #ifdef DISPLAY_nb_explored_cells
                        ++nb_explored_cells[ num_thread ];
                        #endif // DISPLAY_nb_explored_cells

                        const Grid &cr_grid = grids[ cr.num_grid ];
                        const Cell &cr_cell = cr_grid.cells[ cr.num_cell ];

                        // if no cut is possible, we don't go further.
                        if ( may_cut( lc, num_dirac, cr_grid, cr_cell, positions, weights ) == false )
                            continue;

                        // if we have diracs in cr, do the cuts
                        for( TI num_cr_dirac : Span<TI>{ cr_grid.dpc_values.data() + cr_cell.dpc_offset, cr_grid.dpc_values.data() + cr_grid.cells[ cr.num_cell + 1 ].dpc_offset } )
                            plane_cut( lc, num_dirac, num_cr_dirac );

                        // update the front
                        for( TI num_ng_node : Span<TI>{ cr_grid.ng_indices.data() + cr_grid.ng_offsets[ cr.num_cell + 0 ], cr_grid.ng_indices.data() + cr_grid.ng_offsets[ cr.num_cell + 1 ] } )
                            front.push( cr.num_grid, num_ng_node, grids );
                    }

                    //
                    if ( ball_cut )
                        lc.ball_cut( positions[ num_dirac ], sqrt( weights[ num_dirac ] ), num_dirac );

                    //
                    if ( stop_if_void_lc && lc.empty() ) {
                        err = 1;
                        break;
                    }

                    //
                    cb( lc, num_dirac, num_thread );
                }
            }
        } );
    }

    #ifdef DISPLAY_nb_explored_cells
    TI nen = 0, tot = 0;
    for( TI n : nb_explored_cells )
        nen += n;
    for( const Grid &grid : grids )
        tot += grid.cells.size();
    P( tot, nen / nb_diracs );
    #endif // DISPLAY_nb_explored_cells

    return err;
}

template<class Pc>
bool ZGrid<Pc>::check_sanity( const Pt *positions ) const {
    if ( grids.size() == 0 )
        return false;

    // check diracs appear only once
    std::vector<bool> c( grids[ 0 ].cell_index_vs_dirac_number.size(), false );
    for( const Grid &grid : grids ) {
        ASSERT( grid.cell_index_vs_dirac_number.size() == grids[ 0 ].cell_index_vs_dirac_number.size(), "" );
        for( std::size_t num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
            for( TI num_dirac : Span<TI>{ grid.dpc_values.data() + grid.cells[ num_cell + 0 ].dpc_offset, grid.dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } ) {
                ASSERT( c[ num_dirac ] == false, "" );
                c[ num_dirac ] = true;
            }
        }
    }

    // check cell_index_vs_dirac_number: diracs must be inside
    for( std::size_t num_dirac = 0; num_dirac < grids[ 0 ].cell_index_vs_dirac_number.size(); ++num_dirac ) {
        for( const Grid &grid : grids ) {
            const Cell &cell = grid.cells[ grid.cell_index_vs_dirac_number[ num_dirac ] ];
            for( std::size_t d = 0; d < dim; ++d ) {
                ASSERT( positions[ num_dirac ][ d ] <= cell.pos[ d ] + cell.size, "" );
                ASSERT( positions[ num_dirac ][ d ] >= cell.pos[ d ], "" );
            }
        }
    }

    return true;
}

template<class Pc>
void ZGrid<Pc>::update_the_limits( const Pt *positions, const TF *weights, std::size_t nb_diracs ) {
    using std::min;
    using std::max;

    for( std::size_t d = 0; d < dim; ++d ) {
        min_point[ d ] = + std::numeric_limits<TF>::max();
        max_point[ d ] = - std::numeric_limits<TF>::max();
    }
    min_weight = + std::numeric_limits<TF>::max();
    max_weight = - std::numeric_limits<TF>::max();

    for( std::size_t num_dirac = 0; num_dirac < nb_diracs; ++num_dirac ) {
        for( std::size_t d = 0; d < dim; ++d ) {
            min_point[ d ] = min( min_point[ d ], positions[ num_dirac ][ d ] );
            max_point[ d ] = max( max_point[ d ], positions[ num_dirac ][ d ] );
        }
        min_weight = min( min_weight, weights[ num_dirac ] );
        max_weight = max( max_weight, weights[ num_dirac ] );
    }

    int nb_grids = 1 + floor( ( max_weight - min_weight ) / max_delta_weight_per_grid );
    grids.resize( nb_grids );
    P( nb_grids );
    
    //
    grid_length = 0;
    for( std::size_t d = 0; d < dim; ++d )
        grid_length = max( grid_length, max_point[ d ] - min_point[ d ] );
    grid_length *= 1 + std::numeric_limits<TF>::epsilon();

    step_length = grid_length / ( TZ( 1 ) << nb_bits_per_axis );
}


template<class Pc>
void ZGrid<Pc>::update_neighbors( TI num_grid ) {
    Grid &grid = grids[ num_grid ];

    // make a list of requests to get the neighbors
    znodes.resize( 0 );
    znodes.reserve( 2 * dim * zcells.size() );
    for( TI num_cell = 0; num_cell < zcells.size() - 1; ++num_cell ) {
        constexpr TZ f00 = ~ ( ( TZ( 1 ) << dim * nb_bits_per_axis ) - 1 );
        StaticRange<dim>::for_each( [&]( auto d ) {
            TZ z0 = zcells[ num_cell + 0 ].zcoords;
            TZ z1 = zcells[ num_cell + 1 ].zcoords;
            TZ nz = ng_zcoord( z0, z1 - z0, d );
            if ( ( nz & f00 ) == 0 ) // test for overflow
                znodes.push_back( { nz, num_cell } );
        } );
    }

    znodes.reserve( 2 * znodes.size() );
    ZNode *sorted_znodes = radix_sort( znodes.data() + znodes.size(), znodes.data(), znodes.size(), N<sizeof_zcoords>(), rs_tmps );

    // helper function to get the neighbors for each node
    auto for_each_ng = [&]( auto cb ) {
        for( TI i = 0, j = 0; i < znodes.size(); ++i ) {
            // find first node with zcoords > sorted_znodes[ i ].zcoords
            while ( zcells[ j ].zcoords <= sorted_znodes[ i ].zcoords )
                ++j;

            // first node
            TI index_node = sorted_znodes[ i ].index;
            TI index_nbor = j - 1;
            cb( index_node, index_nbor );

            // next touching ones
            TZ off = zcells[ index_node + 1 ].zcoords - zcells[ index_node ].zcoords;
            TZ lim = zcells[ index_nbor     ].zcoords + off;
            if ( zcells[ index_nbor + 1 ].zcoords < lim ) {
                std::array<TZ,dim> tgts;
                StaticRange<dim>::for_each( [&]( auto d ) {
                    tgts[ d ] = zcells[ index_node ].zcoords & _ZcoordsOnesOnAxis<d.val>::value;
                    tgts[ d ] = ( tgts[ d ] | _ZcoordsZerosOnAxis<d.val>::value ) + off;
                    tgts[ d ] = tgts[ d ] & _ZcoordsOnesOnAxis<d.val>::value;
                } );

                ++index_nbor;
                do {
                    bool touching = false;
                    StaticRange<dim>::for_each( [&]( auto d ) {
                        TZ val = zcells[ index_nbor ].zcoords & _ZcoordsOnesOnAxis<d.val>::value;
                        touching |= val == tgts[ d ];
                    } );
                    if ( touching )
                        cb( index_node, index_nbor );
                } while ( zcells[ ++index_nbor ].zcoords < lim );
            }
        }
    };

    //
    grid.ng_offsets.resize( zcells.size() );
    for( TI i = 0; i < zcells.size(); ++i )
        grid.ng_offsets[ i ] = 0;

    // get count
    for_each_ng( [&]( std::size_t index_node, std::size_t index_nbor ) {
        ++grid.ng_offsets[ index_node ];
        ++grid.ng_offsets[ index_nbor ];
    } );

    // suffix scan
    for( TI i = 0, acc = 0; i < zcells.size(); ++i ) {
        TI v = acc;
        acc += grid.ng_offsets[ i ];
        grid.ng_offsets[ i ] = v;
    }

    // get indices
    grid.ng_indices.resize( grid.ng_offsets.back() );
    for_each_ng( [&]( std::size_t index_node, std::size_t index_nbor ) {
        grid.ng_indices[ grid.ng_offsets[ index_node ]++ ] = index_nbor;
        grid.ng_indices[ grid.ng_offsets[ index_nbor ]++ ] = index_node;
    } );

    // shift ng_offsets (to get the suffix scan again)
    if ( grid.ng_offsets.size() ) {
        for( TI i = grid.ng_offsets.size(); --i; )
            grid.ng_offsets[ i ] = grid.ng_offsets[ i - 1 ];
        grid.ng_offsets[ 0 ] = 0;
    }
}

template<class Pc>
void ZGrid<Pc>::fill_grid_using_zcoords( TI num_grid, const Pt *positions, const TF *weights, std::size_t nb_diracs ) {
    using std::round;
    using std::ceil;
    using std::pow;
    using std::min;
    using std::max;

    Grid &grid = grids[ num_grid ];

    // min and max_weights
    if ( grids.size() > 1 ) {
        grid.min_weight = + std::numeric_limits<TF>::max();
        grid.max_weight = - std::numeric_limits<TF>::max();
        for( TI index : grid.dirac_indices ) {
            grid.min_weight = min( grid.min_weight, weights[ index ] );
            grid.max_weight = max( grid.max_weight, weights[ index ] );
        }
    } else {
        grid.min_weight = min_weight;
        grid.max_weight = max_weight;
    }

    // get zcoords for each dirac
    znodes.resize( 0 );
    znodes.reserve( 2 * nb_diracs );
    if ( grids.size() > 1 ) {
        for( TI index : grids[ num_grid ].dirac_indices )
            znodes.push_back( { zcoords_for( positions[ index ] ), index } );
    } else {
        for( TI index = 0; index < nb_diracs; ++index )
            znodes.push_back( { zcoords_for( positions[ index ] ), index } );
    }

    // prepare cell_index_vs_dirac_number => we will set the values for the diracs in this grid
    grid.cell_index_vs_dirac_number.resize( nb_diracs, 666000 );

    // sorting w.r.t. zcoords
    znodes.reserve( 2 * znodes.size() );
    ZNode *sorted_znodes = radix_sort( znodes.data() + znodes.size(), znodes.data(), znodes.size(), N<sizeof_zcoords>(), rs_tmps );

    // fill `cells` with zcoords
    int level = 0;
    TZ prev_z = 0;
    zcells.resize( 0 );
    zcells.reserve( znodes.size() );
    grid.dpc_values.resize( 0 );
    grid.dpc_values.reserve( znodes.size() );
    for( TI index = max_diracs_per_cell; ; ) {
        if ( index >= znodes.size() ) {
            while ( prev_z < ( TZ( 1 ) << dim * nb_bits_per_axis ) ) {
                for( ; ; ++level ) {
                    TZ m = TZ( 1 ) << dim * ( level + 1 );
                    if ( level == nb_bits_per_axis || prev_z & ( m - 1 ) ) {
                        TZ new_prev_z = prev_z + ( TZ( 1 ) << dim * level );

                        ZNode cell;
                        cell.zcoords = prev_z;

                        cell.index = grid.dpc_values.size();
                        for( TI n = index - max_diracs_per_cell; n < znodes.size(); ++n ) {
                            if ( sorted_znodes[ n ].zcoords >= prev_z && sorted_znodes[ n ].zcoords < new_prev_z ) {
                                grid.cell_index_vs_dirac_number[ sorted_znodes[ n ].index ] = zcells.size();
                                grid.dpc_values.push_back( sorted_znodes[ n ].index );
                                ++index;
                            }
                        }

                        zcells.push_back( cell );
                        prev_z = new_prev_z;
                        break;
                    }
                }
            }
            break;
        }

        // level too high ?
        for( ; ; --level ) {
            TZ m = TZ( 1 ) << dim * ( level + 1 );
            if ( sorted_znodes[ index ].zcoords >= prev_z + m )
                break;
            ASSERT( level, "Seems not possible to have $max_diracs_per_cell considering the discretisation (some points are too close)" );
        }

        // look for a level before the one that will take the $max_diracs_per_cell next points or that will lead to an illegal cell
        for( ; ; ++level ) {
            TZ m = TZ( 1 ) << dim * ( level + 1 );
            if ( sorted_znodes[ index ].zcoords < prev_z + m || ( prev_z & ( m - 1 ) ) ) {
                TZ new_prev_z = prev_z + ( TZ( 1 ) << dim * level );

                ZNode zcell;
                zcell.zcoords = prev_z;

                zcell.index = grid.dpc_values.size();
                for( TI n = index - max_diracs_per_cell, l = index; n < l; ++n ) {
                    if ( sorted_znodes[ n ].zcoords >= prev_z && sorted_znodes[ n ].zcoords < new_prev_z ) {
                        grid.cell_index_vs_dirac_number[ sorted_znodes[ n ].index ] = zcells.size();
                        grid.dpc_values.push_back( sorted_znodes[ n ].index );
                        ++index;
                    }
                }

                zcells.push_back( zcell );
                prev_z = new_prev_z;
                break;
            }
        }
    }

    // add an ending cell
    ZNode zcell;
    zcell.index = grid.dpc_values.size();
    zcell.zcoords = TZ( 1 ) << dim * nb_bits_per_axis;
    zcells.push_back( zcell );
}

template<class Pc>
void ZGrid<Pc>::fill_the_grids( const Pt *positions, const TF *weights, std::size_t nb_diracs ) {
    static_assert( sizeof( TZ ) >= sizeof_zcoords, "zcoords types is not large enough" );

    // assign diracs to grids
    if ( grids.size() > 1 ) {
        for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid )
            grids[ num_grid ].dirac_indices.clear();
        for( std::size_t num_dirac = 0; num_dirac < nb_diracs; ++num_dirac ) {
            int num_grid = ( weights[ num_dirac ] - min_weight ) / max_delta_weight_per_grid;
            grids[ num_grid ].dirac_indices.push_back( num_dirac );
        }
    }

    // set grid content
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
        fill_grid_using_zcoords( num_grid, positions, weights, nb_diracs );
        update_neighbors       ( num_grid );
        repl_zcoords_by_ccoords( num_grid );
    }

    // cousins
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid )
        find_englobing_cousins( num_grid, positions );
}

template<class Pc> template<class C>
typename ZGrid<Pc>::TZ ZGrid<Pc>::zcoords_for( const C &pos ) {
    std::array<TZ,dim> c;
    for( int d = 0; d < dim; ++d )
        c[ d ] = TZ( TF( TZ( 1 ) << nb_bits_per_axis ) * ( pos[ d ] - min_point[ d ] ) / grid_length );

    TZ res = 0;
    switch ( dim ) {
    case 1:
        res = c[ 0 ];
        break;
    case 2:
        for( int o = 0; o < nb_bits_per_axis; o += 8 )
            res |= TZ( morton_256_2D_x[ ( c[ 0 ] >> o ) & 0xFF ] |
                       morton_256_2D_y[ ( c[ 1 ] >> o ) & 0xFF ] ) << dim *  o;
        break;
    case 3:
        for( int o = 0; o < nb_bits_per_axis; o += 8 )
            res |= TZ( morton_256_3D_x[ ( c[ 0 ] >> o ) & 0xFF ] |
                       morton_256_3D_y[ ( c[ 1 ] >> o ) & 0xFF ] |
                       morton_256_3D_z[ ( c[ 2 ] >> o ) & 0xFF ] ) << dim *  o;
        break;
    default:
        TODO;
    }

    return res;
}

template<class Pc> template<class V>
void ZGrid<Pc>::display( V &vtk_output ) const {
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
        const Grid &grid = grids[ num_grid ];

        for( TI num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
            Pt p;
            for( int d = 0; d < dim; ++d )
                p[ d ] = grid.cells[ num_cell ].pos[ d ];

            TF a = 0, b = grid.cells[ num_cell ].size;
            switch ( dim ) {
            case 2:
                vtk_output.add_lines( {
                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
                    Point2<TF>{ p[ 0 ] + b, p[ 1 ] + a },
                    Point2<TF>{ p[ 0 ] + b, p[ 1 ] + b },
                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + b },
                    Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
                }, { TF( num_grid ) } );
                break;
            case 3:
                TODO;
                break;
            default:
                TODO;
            }
        }
    }
}

template<class Pc> template<int axis>
typename ZGrid<Pc>::TZ ZGrid<Pc>::ng_zcoord( TZ zcoords, TZ off, N<axis> ) const {
    TZ ff0 = _ZcoordsZerosOnAxis<axis>::value;
    TZ res = ( ( zcoords | ff0 ) + off ) & ~ ff0;
    return res | ( zcoords & ff0 );

}

template<class Pc>
void ZGrid<Pc>::repl_zcoords_by_ccoords( TI num_grid ) {
    Grid &grid = grids[ num_grid ];

    // convert zcoords to cartesian coords
    grid.cells.resize( zcells.size() );
    for( TI num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
        const ZNode &p = zcells[ num_cell + 0 ];
        const ZNode &n = zcells[ num_cell + 1 ];

        Cell &c = grid.cells[ num_cell ];
        c.size = step_length * round( pow( n.zcoords - p.zcoords, 1.0 / dim ) );
        c.zcoords = p.zcoords;
        c.dpc_offset = p.index;

        StaticRange<dim>::for_each( [&]( auto d ) {
            c.pos[ d ] = TF( 0 );
            StaticRange<nb_bits_per_axis>::for_each( [&]( auto i ) {
                c.pos[ d ] += ( p.zcoords & ( TZ( 1 ) << ( dim * i + d ) ) ) >> ( ( dim - 1 ) * i + d );
            } );
            c.pos[ d ] = min_point[ d ] + step_length * c.pos[ d ];
        } );
    }

    Cell &c = grid.cells.back();
    c.dpc_offset = zcells.back().index;
    c.zcoords = zcells.back().zcoords;
    c.size = 0;
    c.pos = max_point;
}

template<class Pc>
void ZGrid<Pc>::find_englobing_cousins( TI num_grid, const Pt *positions ) {
    if ( grids.size() == 1 )
        return;

    // znodes for diracs of current grid
    znodes.clear();
    for( TI num_dirac : grids[ num_grid ].dirac_indices )
        znodes.push_back( { zcoords_for( positions[ num_dirac ] ), num_dirac } );

    // sort znodes
    znodes.reserve( 2 * znodes.size() );
    ZNode *out = radix_sort( znodes.data() + znodes.size(), znodes.data(), znodes.size(), N<sizeof_zcoords>(), rs_tmps );

    //
    for( std::size_t num_ot_grid = 0; num_ot_grid < grids.size(); ++num_ot_grid ) {
        if ( num_ot_grid != num_grid ) {
            for( std::size_t i = 0, j = 0; i < znodes.size(); ++i ) {
                while ( grids[ num_ot_grid ].cells[ j ].zcoords <= out[ i ].zcoords )
                    ++j;
                grids[ num_ot_grid ].cell_index_vs_dirac_number[ out[ i ].index ] = j - 1;
            }
        }
    }
}

} // namespace Visitor
} // namespace PowerDiagram
