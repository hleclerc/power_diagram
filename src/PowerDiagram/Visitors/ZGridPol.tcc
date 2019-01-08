#include "../system/StaticRange.h"
#include "../system/RadixSort.h"
#include "../system/Span.h"

#include "internal/ZCoords.h"
#include "ZGridPol.h"

#include <Eigen/Cholesky>
#include <cmath>

namespace PowerDiagram {
namespace Visitor {


template<class Pc>
ZGridPol<Pc>::ZGridPol( TI max_diracs_per_cell ) : max_diracs_per_cell( max_diracs_per_cell ) {
}

template<class Pc>
void ZGridPol<Pc>::update( const Pt *positions, const TF *weights, TI nb_diracs, bool positions_have_changed, bool weights_have_changed ) {
    if ( positions_have_changed ) {
        _update_the_limits( positions, nb_diracs );
        _fill_the_grid( positions, nb_diracs );
        _update_ngs( positions, nb_diracs );
    }

    if ( positions_have_changed || weights_have_changed ) {
        _update_the_polynomials( positions, weights, nb_diracs );
    }
}

template<class Pc>
void ZGridPol<Pc>::_update_the_limits( const Pt *positions, TI nb_diracs ) {
    using std::min;
    using std::max;

    // min and max point
    for( TI d = 0; d < dim; ++d ) {
        _min_point[ d ] = + std::numeric_limits<TF>::max();
        _max_point[ d ] = - std::numeric_limits<TF>::max();
    }
    for( TI num_dirac = 0; num_dirac < nb_diracs; ++num_dirac ) {
        for( TI d = 0; d < dim; ++d ) {
            _min_point[ d ] = min( _min_point[ d ], positions[ num_dirac ][ d ] );
            _max_point[ d ] = max( _max_point[ d ], positions[ num_dirac ][ d ] );
        }
    }

    // grid and step length
    _grid_length = 0;
    for( TI d = 0; d < dim; ++d )
        _grid_length = max( _grid_length, _max_point[ d ] - _min_point[ d ] );
    _grid_length *= 1 + std::numeric_limits<TF>::epsilon();

    _inv_step_length = ( TZ( 1 ) << nb_bits_per_axis ) / _grid_length;
    _step_length = TF( 1 ) / _inv_step_length;
}

template<class Pc>
void ZGridPol<Pc>::_fill_the_grid( const Pt *positions, TI nb_diracs ) {
    static_assert( sizeof( TZ ) >= sizeof_zcoords, "zcoords types is not large enough" );

    using std::round;
    using std::ceil;
    using std::pow;
    using std::min;
    using std::max;

    // get zcoords for each dirac
    _zcoords.resize( 0 );
    _zcoords.reserve( 2 * nb_diracs );
    for( TI index = 0; index < nb_diracs; ++index )
        _zcoords.push_back( { _zcoords_for( positions[ index ] ), index } );

    // sorting w.r.t. zcoords
    _zcoords.reserve( 2 * _zcoords.size() );
    _sorted_zcoords = radix_sort( _zcoords.data() + _zcoords.size(), _zcoords.data(), _zcoords.size(), N<sizeof_zcoords>(), _tmp_rs );

    // fill `cells` with zcoords
    int level = 0;
    TZ prev_z = 0;
    _cells.resize( 0 );
    _dpc_values.resize( 0 );
    for( TI index = max_diracs_per_cell; ; ) {
        if ( index >= _zcoords.size() ) {
            while ( prev_z < ( TZ( 1 ) << dim * nb_bits_per_axis ) ) {
                for( ; ; ++level ) {
                    TZ m = TZ( 1 ) << dim * ( level + 1 );
                    if ( level == nb_bits_per_axis || prev_z & ( m - 1 ) ) {
                        // prepare a new cell
                        Cell cell;
                        cell.dpc_offset = _dpc_values.size();
                        cell.zcoords = prev_z;

                        // register the corresponding diracs
                        TZ new_prev_z = prev_z + ( TZ( 1 ) << dim * level );
                        for( TI n = index - max_diracs_per_cell; n < _zcoords.size(); ++n ) {
                            if ( _sorted_zcoords[ n ].zcoords >= prev_z && _sorted_zcoords[ n ].zcoords < new_prev_z ) {
                                _dpc_values.push_back( _sorted_zcoords[ n ].index );
                                ++index;
                            }
                        }

                        _cells.push_back( cell );
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
            if ( _sorted_zcoords[ index ].zcoords >= prev_z + m )
                break;
            ASSERT( level, "Seems not possible to have $max_diracs_per_cell considering the discretisation (some points are too close)" );
        }

        // look for a level before the one that will take the $max_diracs_per_cell next points or that will lead to an illegal cell
        for( ; ; ++level ) {
            TZ m = TZ( 1 ) << dim * ( level + 1 );
            if ( _sorted_zcoords[ index ].zcoords < prev_z + m || ( prev_z & ( m - 1 ) ) ) {
                // prepare a new cell
                Cell cell;
                cell.dpc_offset = _dpc_values.size();
                cell.zcoords = prev_z;

                // register the corresponding diracs
                TZ new_prev_z = prev_z + ( TZ( 1 ) << dim * level );
                for( TI n = index - max_diracs_per_cell, l = index; n < l; ++n ) {
                    if ( _sorted_zcoords[ n ].zcoords >= prev_z && _sorted_zcoords[ n ].zcoords < new_prev_z ) {
                        _dpc_values.push_back( _sorted_zcoords[ n ].index );
                        ++index;
                    }
                }

                // register the cell
                _cells.push_back( cell );
                prev_z = new_prev_z;
                break;
            }
        }
    }

    // add an ending cell
    Cell cell;
    cell.zcoords = TZ( 1 ) << dim * nb_bits_per_axis;
    cell.dpc_offset = _dpc_values.size();
    _cells.push_back( cell );

    // set pos and size in each cell
    for( TI num_cell = 0; num_cell < _cells.size() - 1; ++num_cell ) {
        Cell &c = _cells[ num_cell + 0 ];
        Cell &n = _cells[ num_cell + 1 ];

        // pos
        StaticRange<dim>::for_each( [&]( auto d ) {
            c.pos[ d ] = TF( 0 );
            StaticRange<nb_bits_per_axis>::for_each( [&]( auto i ) {
                c.pos[ d ] += ( c.zcoords & ( TZ( 1 ) << ( dim * i + d ) ) ) >> ( ( dim - 1 ) * i + d );
            } );
            c.pos[ d ] = _min_point[ d ] + _step_length * c.pos[ d ];
        } );

        // size
        c.size = _step_length * round( pow( n.zcoords - c.zcoords, 1.0 / dim ) );
    }
}

template<class Pc>
void ZGridPol<Pc>::_update_ngs( const Pt *positions, TI nb_diracs ) {
    // make a list of requests to get the neighbors
    _tmp_zn.resize( 0 );
    _tmp_zn.reserve( 2 * dim * _cells.size() );
    for( TI num_cell = 0; num_cell < _cells.size() - 1; ++num_cell ) {
        constexpr TZ f00 = ~ ( ( TZ( 1 ) << dim * nb_bits_per_axis ) - 1 );
        StaticRange<dim>::for_each( [&]( auto d ) {
            TZ z0 = _cells[ num_cell + 0 ].zcoords;
            TZ z1 = _cells[ num_cell + 1 ].zcoords;
            TZ nz = _ng_zcoord( z0, z1 - z0, d );
            if ( ( nz & f00 ) == 0 ) // test for overflow
                _tmp_zn.push_back( { nz, num_cell } );
        } );
    }

    _tmp_zn.reserve( 2 * _tmp_zn.size() );
    ZNode *sorted_znodes = radix_sort( _tmp_zn.data() + _tmp_zn.size(), _tmp_zn.data(), _tmp_zn.size(), N<sizeof_zcoords>(), _tmp_rs );

    // helper function to get the neighbors for each node
    auto for_each_ng = [&]( auto cb ) {
        for( TI i = 0, j = 0; i < _tmp_zn.size(); ++i ) {
            // find first node with zcoords > sorted_znodes[ i ].zcoords
            while ( _cells[ j ].zcoords <= sorted_znodes[ i ].zcoords )
                ++j;

            // first node
            TI index_node = sorted_znodes[ i ].index;
            TI index_nbor = j - 1;
            cb( index_node, index_nbor );

            // next touching ones
            TZ off = _cells[ index_node + 1 ].zcoords - _cells[ index_node ].zcoords;
            TZ lim = _cells[ index_nbor     ].zcoords + off;
            if ( _cells[ index_nbor + 1 ].zcoords < lim ) {
                std::array<TZ,dim> tgts;
                StaticRange<dim>::for_each( [&]( auto d ) {
                    using Zooa = typename ZCoords<TZ,dim,nb_bits_per_axis,sizeof_zcoords>::template _ZcoordsOnesOnAxis<d.val>;
                    tgts[ d ] = _cells[ index_node ].zcoords & Zooa::value;
                    tgts[ d ] = ( tgts[ d ] | Zooa::value ) + off;
                    tgts[ d ] = tgts[ d ] & Zooa::value;
                } );

                ++index_nbor;
                do {
                    bool touching = false;
                    StaticRange<dim>::for_each( [&]( auto d ) {
                        using Zooa = typename ZCoords<TZ,dim,nb_bits_per_axis,sizeof_zcoords>::template _ZcoordsOnesOnAxis<d.val>;
                        TZ val = _cells[ index_nbor ].zcoords & Zooa::value;
                        touching |= val == tgts[ d ];
                    } );
                    if ( touching )
                        cb( index_node, index_nbor );
                } while ( _cells[ ++index_nbor ].zcoords < lim );
            }
        }
    };

    //
    _ng_offsets.resize( _cells.size() );
    for( TI i = 0; i < _cells.size(); ++i )
        _ng_offsets[ i ] = 0;

    // get count
    for_each_ng( [&]( TI index_node, TI index_nbor ) {
        ++_ng_offsets[ index_node ];
        ++_ng_offsets[ index_nbor ];
    } );

    // suffix scan
    for( TI i = 0, acc = 0; i < _cells.size(); ++i ) {
        TI v = acc;
        acc += _ng_offsets[ i ];
        _ng_offsets[ i ] = v;
    }

    // get indices
    _ng_values.resize( _ng_offsets.back() );
    for_each_ng( [&]( TI index_node, TI index_nbor ) {
        _ng_values[ _ng_offsets[ index_node ]++ ] = index_nbor;
        _ng_values[ _ng_offsets[ index_nbor ]++ ] = index_node;
    } );

    // shift ng_offsets (to get the suffix scan again)
    if ( _ng_offsets.size() ) {
        for( TI i = _ng_offsets.size(); --i; )
            _ng_offsets[ i ] = _ng_offsets[ i - 1 ];
        _ng_offsets[ 0 ] = 0;
    }
}

template<class Pc>
void ZGridPol<Pc>::_update_the_polynomials( const Pt *positions, const TF *weights, TI nb_diracs ) {
    std::vector<Pt> cp_positions( nb_diracs );
    std::vector<TF> cp_weights( nb_diracs );
    for( TI i = 0; i < nb_diracs; ++i ) {
        cp_positions[ i ] = positions[ i ];
        cp_weights[ i ] = weights[ i ];
    }

    _polynomials.clear();
    _subdivide_add_poly_rec( cp_positions.data(), cp_weights.data(), nb_diracs, _min_point, _max_point );
}

template<class Pc>
void ZGridPol<Pc>::_subdivide_add_poly_rec( Pt *positions, TF *weights, TI nb_diracs, Pt p0, Pt p1 ) {
    using EM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    using EV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
    using std::pow;

    constexpr int nc = dim * ( dim + 3 ) / 2 + 1;

    // try with all the position and all the weights
    EM M( nc, nc );
    for( TI r = 0; r < nc; ++r )
        for( TI c = 0; c < nc; ++c )
            M( r, c ) = 0;
    EV V( nc );
    for( TI r = 0; r < nc; ++r )
        V[ r ] = 0;

    //
    for( TI n = 0; n < nb_diracs; ++n ) {
        std::array<TF,nc> coeffs;
        coeffs[ 0 ] = 1;
        coeffs[ 1 ] = positions[ n ].x;
        coeffs[ 2 ] = positions[ n ].y;
        coeffs[ 3 ] = positions[ n ].x * positions[ n ].x;
        coeffs[ 4 ] = positions[ n ].x * positions[ n ].y;
        coeffs[ 5 ] = positions[ n ].y * positions[ n ].y;

        for( TI r = 0; r < nc; ++r ) {
            for( TI c = 0; c < nc; ++c )
                M( r, c ) += coeffs[ r ] * coeffs[ c ];
            V[ r ] += coeffs[ r ] * weights[ n ];
        }
    }

    P( M );
    P( V );

    // solve and update the weights
    Eigen::Cholesky<EM> C;
    C.compute( M );
    EV D = C.solve( V );

    P( D );
}

//template<class Pc>
//bool ZGrid<Pc>::may_cut( const CP &lc, TI i0, const Grid &cr_grid, const Cell &cr_cell, const Pt *positions, const TF *weights ) {
//    using std::sqrt;
//    using std::max;
//    using std::min;
//    using std::pow;
//    using std::abs;

//    auto c0 = positions[ i0 ];
//    auto w0 = weights  [ i0 ];

//    return true;

//    //
//    if ( ball_cut ) {
//        TF md = 0; // min dist pow 2
//        for( size_t d = 0; d < dim; ++d ) {
//            TF o = c0[ d ] - cr_cell.pos[ d ];
//            if ( o > 0 )
//                o = max( TF( 0 ), o - cr_cell.size );
//            md += pow( o, 2 );
//        }

//        return md < pow( sqrt( cr_grid.max_weight ) + sqrt( w0 ), 2 );
//    }

//    //
//    if ( w0 > cr_grid.min_weight ) {
//        TF d2 = w0 - cr_grid.min_weight;
//        if ( d2 > 0 ) {
//            //
//            TF max_cur_dist = 0;
//            for( TI num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point )
//                max_cur_dist = max( max_cur_dist, norm_2_p2( lc.point( num_lc_point ) ) );
//            if ( max_cur_dist > 4 * d2 ) {
//                TF dx = cr_cell.pos[ 0 ] - c0[ 0 ];
//                TF dy = cr_cell.pos[ 1 ] - c0[ 1 ];
//                TF x2 = min( pow( dx, 2 ), pow( dx + cr_cell.size, 2 ) );
//                TF y2 = min( pow( dy, 2 ), pow( dy + cr_cell.size, 2 ) );
//                if ( x2 + y2 < d2 )
//                    return true;
//            }
//        }
//    }

//    //
//    TF l0 = cr_cell.size;
//    auto test_line = [&]( Pt A, Pt B, Pt p ) {
//        TF s2 = l0 * l0;
//        TF s1 = dot( p - A, B - A );
//        TF s0 = ( norm_2_p2( p - A ) - cr_grid.max_weight ) - ( norm_2_p2( c0 - p ) - w0 );
//        if ( s1 > 0 && s1 < s2 && s2 * s0 - s1 * s1 < 0 )
//            return true;
//        return s0 < 0 || s2 - 2 * s1 + s0 < 0;
//    };

//    TF dx = cr_cell.pos[ 0 ] < c0[ 0 ] ? l0 : 0;
//    TF dy = cr_cell.pos[ 1 ] < c0[ 1 ] ? l0 : 0;

//    Pt PA { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ] + dy };
//    Pt PB { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ] + dy };

//    Pt PC { cr_cell.pos[ 0 ] + dx, cr_cell.pos[ 1 ]      };
//    Pt PD { cr_cell.pos[ 0 ] + dx, cr_cell.pos[ 1 ] + l0 };

//    for( TI num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point ) {
//        if ( test_line( PA, PB, lc.point( num_lc_point ) ) ) return true;
//        if ( test_line( PC, PD, lc.point( num_lc_point ) ) ) return true;
//    }

//    //    Pt PA { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ]      };
//    //    Pt PB { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ]      };
//    //    Pt PC { cr_cell.pos[ 0 ] + l0, cr_cell.pos[ 1 ] + l0 };
//    //    Pt PD { cr_cell.pos[ 0 ]     , cr_cell.pos[ 1 ] + l0 };

//    //    for( TI num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point ) {
//    //        if ( test_line( PA, PB, lc.point( num_lc_point ) ) ) return true;
//    //        if ( test_line( PB, PC, lc.point( num_lc_point ) ) ) return true;
//    //        if ( test_line( PC, PD, lc.point( num_lc_point ) ) ) return true;
//    //        if ( test_line( PD, PA, lc.point( num_lc_point ) ) ) return true;
//    //    }

//    return false;
//}

//template<class Pc>
//int ZGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, TI num )> &cb, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc ) {
//    return for_each_laguerre_cell( [&]( CP &cp, TI num, int ) {
//        cb( cp, num );
//    }, starting_lc, positions, weights, nb_diracs, stop_if_void_lc );
//}

//template<class Pc>
//int ZGrid<Pc>::for_each_laguerre_cell( const std::function<void( CP &, TI num, int num_thread )> &cb, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc ) {
//    using Front = FrontZgrid<ZGrid>;
//    using std::sqrt;

//    auto plane_cut = [&]( CP &lc, TI i0, TI i1 ) {
//        Pt V = positions[ i1 ] - positions[ i0 ];
//        TF n = norm_2_p2( V );
//        TF x = TF( 1 ) + ( weights[ i0 ] - weights[ i1 ] ) / n;
//        TF i = TF( 1 ) / sqrt( n );
//        lc.plane_cut( positions[ i0 ] + TF( 0.5 ) * x * V, i * V, i1 );
//    };

//    // vectors for stuff that will be reused inside the execution threads
//    int nb_threads = thread_pool.nb_threads(), nb_jobs = 4 * nb_threads;
//    std::vector<std::vector<std::vector<TI>>> visited( nb_threads );
//    std::vector<TI> op_counts( nb_threads, 0 );

//    for( int num_thread = 0; num_thread < nb_threads; ++num_thread ) {
//        visited[ num_thread ].resize( grids.size() );
//        for( TI num_grid = 0; num_grid < grids.size(); ++num_grid )
//            visited[ num_thread ][ num_grid ].resize( grids[ num_grid ].cells.size(), op_counts[ num_thread ] );
//    }

//    #ifdef DISPLAY_nb_explored_cells
//    std::vector<TI> nb_explored_cells( nb_threads, 0 );
//    #endif // DISPLAY_nb_explored_cells

//    // for each item
//    int err = 0;
//    for( TI num_grid = 0; num_grid < grids.size(); ++num_grid ) {
//        Grid &grid = grids[ num_grid ];

//        thread_pool.execute( nb_jobs, [&]( TI num_job, int num_thread ) {
//            Front front( op_counts[ num_thread ], visited[ num_thread ] );
//            CP lc;

//            TI beg_cell = ( num_job + 0 ) * ( grid.cells.size() - 1 ) / nb_jobs;
//            TI end_cell = ( num_job + 1 ) * ( grid.cells.size() - 1 ) / nb_jobs;
//            for( TI num_cell = beg_cell; num_cell < end_cell && err == 0; ++num_cell ) {
//                const Cell &cell = grid.cells[ num_cell ];
//                for( TI num_dirac : Span<TI>{ _dpc_values.data() + cell.dpc_offset, _dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } ) {
//                    // start of lc: cut with nodes in the same cell
//                    lc = starting_lc;
//                    for( TI num_cr_dirac : Span<TI>{ _dpc_values.data() + cell.dpc_offset, _dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } )
//                        if ( num_cr_dirac != num_dirac )
//                            plane_cut( lc, num_dirac, num_cr_dirac );

//                    // front
//                    front.init( grids, num_grid, num_cell, positions[ num_dirac ], weights[ num_dirac ] );

//                    // => neighbors in the same grid
//                    for( TI num_ng_cell : Span<TI>{ grid.ng_indices.data() + grid.ng_offsets[ num_cell + 0 ], grid.ng_indices.data() + grid.ng_offsets[ num_cell + 1 ] } )
//                        front.push_without_check( num_grid, num_ng_cell, grids );

//                    // => items from the grids made for != weight (containing the dirac position)
//                    for( TI num_pa_grid = 0; num_pa_grid < grids.size(); ++num_pa_grid ) {
//                        if ( num_pa_grid != num_grid ) {
//                            Grid &pa_grid = grids[ num_pa_grid ];
//                            TI num_pa_cell = pa_grid.cell_index_vs_dirac_number[ num_dirac ];

//                            // cut with items in pa cell
//                            front.set_visited( grids, num_pa_grid, num_pa_cell );
//                            for( TI num_cr_dirac : Span<TI>{ pa__dpc_values.data() + pa_grid.cells[ num_pa_cell + 0 ].dpc_offset, pa__dpc_values.data() + pa_grid.cells[ num_pa_cell + 1 ].dpc_offset } )
//                                plane_cut( lc, num_dirac, num_cr_dirac );

//                            // add neighbors in the front
//                            for( TI num_ng_cell : Span<TI>{ pa_grid.ng_indices.data() + pa_grid.ng_offsets[ num_pa_cell + 0 ], pa_grid.ng_indices.data() + pa_grid.ng_offsets[ num_pa_cell + 1 ] } )
//                                front.push_without_check( num_pa_grid, num_ng_cell, grids );
//                        }
//                    }

//                    // neighbors
//                    while( ! front.empty() ) {
//                        typename Front::Item cr = front.pop();
//                        #ifdef DISPLAY_nb_explored_cells
//                        ++nb_explored_cells[ num_thread ];
//                        #endif // DISPLAY_nb_explored_cells

//                        const Grid &cr_grid = grids[ cr.num_grid ];
//                        const Cell &cr_cell = cr_grid.cells[ cr.num_cell ];

//                        // if no cut is possible, we don't go further.
//                        if ( may_cut( lc, num_dirac, cr_grid, cr_cell, positions, weights ) == false )
//                            continue;

//                        // if we have diracs in cr, do the cuts
//                        for( TI num_cr_dirac : Span<TI>{ cr__dpc_values.data() + cr_cell.dpc_offset, cr__dpc_values.data() + cr_grid.cells[ cr.num_cell + 1 ].dpc_offset } )
//                            plane_cut( lc, num_dirac, num_cr_dirac );

//                        // update the front
//                        for( TI num_ng_node : Span<TI>{ cr_grid.ng_indices.data() + cr_grid.ng_offsets[ cr.num_cell + 0 ], cr_grid.ng_indices.data() + cr_grid.ng_offsets[ cr.num_cell + 1 ] } )
//                            front.push( cr.num_grid, num_ng_node, grids );
//                    }

//                    //
//                    if ( ball_cut )
//                        lc.ball_cut( positions[ num_dirac ], sqrt( weights[ num_dirac ] ), num_dirac );

//                    //
//                    if ( stop_if_void_lc && lc.empty() ) {
//                        err = 1;
//                        break;
//                    }

//                    //
//                    cb( lc, num_dirac, num_thread );
//                }
//            }
//        } );
//    }

//    #ifdef DISPLAY_nb_explored_cells
//    TI nen = 0, tot = 0;
//    for( TI n : nb_explored_cells )
//        nen += n;
//    for( const Grid &grid : grids )
//        tot += grid.cells.size();
//    P( tot, nen / nb_diracs );
//    #endif // DISPLAY_nb_explored_cells

//    return err;
//}

//template<class Pc>
//bool ZGrid<Pc>::check_sanity( const Pt *positions ) const {
//    if ( grids.size() == 0 )
//        return false;

//    // check diracs appear only once
//    std::vector<bool> c( grids[ 0 ].cell_index_vs_dirac_number.size(), false );
//    for( const Grid &grid : grids ) {
//        ASSERT( grid.cell_index_vs_dirac_number.size() == grids[ 0 ].cell_index_vs_dirac_number.size(), "" );
//        for( TI num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
//            for( TI num_dirac : Span<TI>{ _dpc_values.data() + grid.cells[ num_cell + 0 ].dpc_offset, _dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } ) {
//                ASSERT( c[ num_dirac ] == false, "" );
//                c[ num_dirac ] = true;
//            }
//        }
//    }

//    // check cell_index_vs_dirac_number: diracs must be inside
//    for( TI num_dirac = 0; num_dirac < grids[ 0 ].cell_index_vs_dirac_number.size(); ++num_dirac ) {
//        for( const Grid &grid : grids ) {
//            const Cell &cell = grid.cells[ grid.cell_index_vs_dirac_number[ num_dirac ] ];
//            for( TI d = 0; d < dim; ++d ) {
//                ASSERT( positions[ num_dirac ][ d ] <= cell.pos[ d ] + cell.size, "" );
//                ASSERT( positions[ num_dirac ][ d ] >= cell.pos[ d ], "" );
//            }
//        }
//    }

//    return true;
//}

template<class Pc> template<class C>
typename ZGridPol<Pc>::TZ ZGridPol<Pc>::_zcoords_for( const C &pos ) {
    std::array<TZ,dim> c;
    for( int d = 0; d < dim; ++d )
        c[ d ] = TZ( _inv_step_length * ( pos[ d ] - _min_point[ d ] ) );

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
void ZGridPol<Pc>::display( V &vtk_output ) const {
    for( TI num_cell = 0; num_cell < _cells.size() - 1; ++num_cell ) {
        Pt p;
        for( int d = 0; d < dim; ++d )
            p[ d ] = _cells[ num_cell ].pos[ d ];

        TF a = 0, b = _cells[ num_cell ].size;
        switch ( dim ) {
        case 2:
            vtk_output.add_lines( {
                Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
                Point2<TF>{ p[ 0 ] + b, p[ 1 ] + a },
                Point2<TF>{ p[ 0 ] + b, p[ 1 ] + b },
                Point2<TF>{ p[ 0 ] + a, p[ 1 ] + b },
                Point2<TF>{ p[ 0 ] + a, p[ 1 ] + a },
            } );
            break;
        case 3:
            TODO;
            break;
        default:
            TODO;
        }
    }
}

template<class Pc> template<int axis>
typename ZGridPol<Pc>::TZ ZGridPol<Pc>::_ng_zcoord( TZ zcoords, TZ off, N<axis> ) const {
    using Zzoa = typename ZCoords<TZ,dim,nb_bits_per_axis,sizeof_zcoords>::template _ZcoordsZerosOnAxis<axis>;
    TZ ff0 = Zzoa::value;
    TZ res = ( ( zcoords | ff0 ) + off ) & ~ ff0;
    return res | ( zcoords & ff0 );
}

//template<class Pc>
//void ZGrid<Pc>::repl_zcoords_by_ccoords( TI num_grid ) {
//    Grid &grid = grids[ num_grid ];

//    // convert zcoords to cartesian coords
//    grid.cells.resize( zcells.size() );
//    for( TI num_cell = 0; num_cell < grid.cells.size() - 1; ++num_cell ) {
//        const ZNode &p = zcells[ num_cell + 0 ];
//        const ZNode &n = zcells[ num_cell + 1 ];

//        Cell &c = grid.cells[ num_cell ];
//        c.size = step_length * round( pow( n.zcoords - p.zcoords, 1.0 / dim ) );
//        c.zcoords = p.zcoords;
//        c.dpc_offset = p.index;

//        StaticRange<dim>::for_each( [&]( auto d ) {
//            c.pos[ d ] = TF( 0 );
//            StaticRange<nb_bits_per_axis>::for_each( [&]( auto i ) {
//                c.pos[ d ] += ( p.zcoords & ( TZ( 1 ) << ( dim * i + d ) ) ) >> ( ( dim - 1 ) * i + d );
//            } );
//            c.pos[ d ] = min_point[ d ] + step_length * c.pos[ d ];
//        } );
//    }

//    Cell &c = grid.cells.back();
//    c.dpc_offset = zcells.back().index;
//    c.zcoords = zcells.back().zcoords;
//    c.size = 0;
//    c.pos = max_point;
//}

//template<class Pc>
//void ZGrid<Pc>::find_englobing_cousins( TI num_grid, const Pt *positions ) {
//    if ( grids.size() == 1 )
//        return;

//    // znodes for diracs of current grid
//    znodes.clear();
//    for( TI num_dirac : grids[ num_grid ].dirac_indices )
//        znodes.push_back( { zcoords_for( positions[ num_dirac ] ), num_dirac } );

//    // sort znodes
//    znodes.reserve( 2 * znodes.size() );
//    ZNode *out = radix_sort( znodes.data() + znodes.size(), znodes.data(), znodes.size(), N<sizeof_zcoords>(), rs_tmps );

//    //
//    for( TI num_ot_grid = 0; num_ot_grid < grids.size(); ++num_ot_grid ) {
//        if ( num_ot_grid != num_grid ) {
//            for( TI i = 0, j = 0; i < znodes.size(); ++i ) {
//                while ( grids[ num_ot_grid ].cells[ j ].zcoords <= out[ i ].zcoords )
//                    ++j;
//                grids[ num_ot_grid ].cell_index_vs_dirac_number[ out[ i ].index ] = j - 1;
//            }
//        }
//    }
//}

} // namespace Visitor
} // namespace PowerDiagram
