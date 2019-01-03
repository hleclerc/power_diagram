#include "../system/StaticRange.h"
#include "../system/RadixSort.h"
#include "../system/Span.h"
#include "ZGrid.h"
#include <queue>
#include <cmath>

#define DISPLAY_nb_explored_cells

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

    auto dp = positions[ i0 ];
    auto dw = weights  [ i0 ];

    //
    if ( ball_cut ) {
        // get min dist pow 2
        TF md = 0;
        for( size_t d = 0; d < dim; ++d ) {
            TF o = dp[ d ] - cr_cell.pos[ d ];
            if ( o > 0 )
                o = max( TF( 0 ), o - cr_cell.size );
            md += pow( o, 2 );
        }

        return md < pow( sqrt( cr_grid.max_weight ) + sqrt( dw ), 2 );
    }

    TF l0 = cr_cell.size;
    auto test_line = [&]( Pt A, Pt B, Pt p ) {
        TF c2 = l0 * l0;
        TF c1 = dot( B - A, p - A );
        TF c0 = norm_2_p2( p - A ) - norm_2_p2( dp - p ) + dw - cr_grid.max_weight;
        if ( c1 > 0 && c1 < c2 && c2 * c0 - c1 * c1 < 0 )
            return true;
        return c0 < 0 || c2 - 2 * c1 + c0 < 0;
    };

    Pt c0;
    for( int d = 0; d < dim; ++d )
        c0[ d ] = cr_cell.pos[ d ];

    TF dy = c0[ 1 ] < dp[ 1 ] ? l0 : 0;
    Pt PA { c0[ 0 ]     , c0[ 1 ] + dy };
    Pt PB { c0[ 0 ] + l0, c0[ 1 ] + dy };

    TF dx = c0[ 0 ] < dp[ 0 ] ? l0 : 0;
    Pt PC { c0[ 0 ] + dx, c0[ 1 ]      };
    Pt PD { c0[ 0 ] + dx, c0[ 1 ] + l0 };

    for( std::size_t num_lc_point = 0; num_lc_point < lc.nb_points; ++num_lc_point ) {
        if ( test_line( PA, PB, lc.point( num_lc_point ) ) ) return true;
        if ( test_line( PC, PD, lc.point( num_lc_point ) ) ) return true;
    }
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
    using std::sqrt;

    #define ZIndex_USES_HEAP
    #ifdef ZIndex_USES_HEAP
    struct Front {
        struct Item {
            bool operator<( const Item &that ) const { return dist > that.dist; }
            TI   num_grid;
            TI   num_cell;
            TF   dist;
        };

        Front( TI& op_count, std::vector<std::vector<TI>> &visited ) : op_count( op_count ), visited( visited ) {
        }

        void init( const std::vector<Grid> &grids, TI num_grid, TI num_cell ) {
            for( std::size_t n = 0; n < grids.size(); ++n )
                visited[ n ].resize( grids[ n ].cells.size(), op_count );

            orig_cell_pos = grids[ num_grid ].cells[ num_cell ].pos;
            visited[ num_grid ][ num_cell ] = ++op_count;
        }

        TF dist( const Cell &cell ) {
            TF res = 0;
            for( int d = 0; d < dim; ++d ) {
                TF v = cell.pos[ d ] - orig_cell_pos[ d ]; // we need abs to avoid the overflow
                res += v * v;
            }
            return res;
        }

        void push_without_check( TI num_grid, TI num_cell, const std::vector<Grid> &grids ) {
            items.push( Item{ num_grid, num_cell, dist( grids[ num_grid ].cells[ num_cell ] ) } );
            visited[ num_grid ][ num_cell ] = op_count;
        }

        void push( TI num_grid, TI num_cell, const std::vector<Grid> &grids ) {
            if ( visited[ num_grid ][ num_cell ] != op_count ) {
                visited[ num_grid ][ num_cell ] = op_count;
                push_without_check( num_grid, num_cell, grids );
            }
        }

        Item pop() {
            Item res = items.top();
            items.pop();
            return res;
        }

        bool empty() const {
            return items.empty();
        }

        Pt                            orig_cell_pos;
        TI&                           op_count;
        std::vector<std::vector<TI>>& visited;
        std::priority_queue<Item>     items;
    };
    #else
    struct Front {
        struct Item {
            TI   num_grid;
            TI   num_cell;
        };

        Front( TI nb_grids ) : visited( nb_grids ) {
            op_count = 0;
            rand = 0;
        }

        void init( const std::vector<Grid> &grids, TI num_grid, TI num_cell ) {
            for( std::size_t n = 0; n < grids.size(); ++n )
                visited[ n ].resize( grids[ n ].cells.size(), op_count );

            visited[ num_grid ][ num_cell ] = ++op_count;
        }

        void push_without_check( TI num_grid, TI num_cell, const std::vector<Grid> & ) {
            visited[ num_grid ][ num_cell ] = op_count;
            items.push_back( { num_grid, num_cell } );
        }

        void push( TI num_grid, TI num_cell, const std::vector<Grid> &grids ) {
            if ( visited[ num_grid ][ num_cell ] != op_count ) {
                visited[ num_grid ][ num_cell ] = op_count;
                push_without_check( num_grid, num_cell, grids );
            }
        }

        Item pop() {
            TI n = rand % items.size();
            Item res = items[ n ];
            rand += 13;
            items[ n ] = items[ items.size() - 1 ];
            items.pop_back();
            return res;
        }

        bool empty() const {
            return items.empty();
        }

        TI                            op_count;
        std::vector<std::vector<TI>>& visited;
        std::vector<Item>             items;
        TI                            rand;
    };
    #endif

    auto plane_cut = [&]( CP &lc, TI i0, TI i1 ) {
        Pt V = positions[ i1 ] - positions[ i0 ];
        TF d = norm_2( V );
        TF i = 1.0 / d;
        Pt N = i * V;
        TF x = 0.5 * ( d + i * ( weights[ i0 ] - weights[ i1 ] ) );
        lc.plane_cut( positions[ i0 ] + x * N, N, i1 );
    };

    // vectors for stuff that will be reused inside the execution threads
    int nb_threads = thread_pool.nb_threads(), nb_jobs = 4 * nb_threads;
    std::vector<std::vector<std::vector<TI>>> visited( nb_threads );
    std::vector<TI> op_counts( nb_threads, 0 );
    for( int n = 0; n < nb_threads; ++n )
        visited[ n ].resize( grids.size() );

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
                    // front
                    front.init( grids, num_grid, num_cell );

                    // => neighbors in the same grid
                    for( TI num_ng_node : Span<TI>{ grid.ng_indices.data() + grid.ng_offsets[ num_cell + 0 ], grid.ng_indices.data() + grid.ng_offsets[ num_cell + 1 ] } )
                        front.push_without_check( num_grid, num_ng_node, grids );

                    // => items of the grids for != weight containing the dirac
                    for( std::size_t num_ng_grid = 0; num_ng_grid < grids.size(); ++num_ng_grid ) {
                        if ( num_ng_grid != num_grid ) {
                            TI num_ng_node = grids[ num_ng_grid ].cell_index_vs_dirac_number[ num_dirac ];
                            front.push_without_check( num_ng_grid, num_ng_node, grids );
                        }
                    }

                    // start of lc: cut with nodes in the same cell
                    lc = starting_lc;
                    for( TI num_cr_dirac : Span<TI>{ grid.dpc_values.data() + cell.dpc_offset, grid.dpc_values.data() + grid.cells[ num_cell + 1 ].dpc_offset } )
                        if ( num_cr_dirac != num_dirac )
                            plane_cut( lc, num_dirac, num_cr_dirac );

                    // neighbors
                    while( ! front.empty() ) {
                        typename Front::Item cr = front.pop();
                        #ifdef DISPLAY_nb_explored_cells
                        ++nb_explored_cells[ num_thread ];
                        #endif // DISPLAY_nb_explored_cells

                        const Grid &cr_grid = grids[ cr.num_grid ];
                        const Cell &cr_cell = cr_grid.cells[ cr.num_cell ];

                        // if no cut is possible, we don't go further.
                        // The cost of the test being far from negligible, we don't make it for close cells
                        if ( may_cut( lc, num_dirac, cr_grid, cr_cell, positions, weights ) == false )
                            continue;

                        // if we have diracs in cr, do the cuts
                        for( TI num_cr_dirac : Span<TI>{ grid.dpc_values.data() + cr_cell.dpc_offset, grid.dpc_values.data() + cr_grid.cells[ cr.num_cell + 1 ].dpc_offset } )
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
    TI nen = 0;
    for( TI n : nb_explored_cells )
        nen += n;
    P( grids[ 0 ].cells.size(), nen / nb_diracs );
    #endif // DISPLAY_nb_explored_cells

    return err;
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

    if ( max_weight == min_weight )
        div_weight = 10 / max_delta_weight_per_grid;
    else
        div_weight = ( 1 - std::numeric_limits<TF>::epsilon() ) / ( max_weight - min_weight );

    int nb_grids = ceil( ( 1 / div_weight ) / max_delta_weight_per_grid );
    grids.resize( nb_grids );

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
    grid.cell_index_vs_dirac_number.resize( znodes.size(), 666 );

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

    using std::log2;

    // assign diracs to grids
    if ( grids.size() > 1 ) {
        for( std::size_t num_dirac = 0; num_dirac < nb_diracs; ++num_dirac ) {
            int num_grid = div_weight * ( weights[ num_dirac ] - min_weight );
            grids[ num_grid ].dirac_indices.push_back( num_dirac );
        }
    }

    // set grid content
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
        fill_grid_using_zcoords( num_grid, positions, weights, nb_diracs );
        update_neighbors       ( num_grid );
        repl_zcoords_by_ccoords( num_grid );
    }

    // find englobing cells for each dirac (and for each grid)
    for( std::size_t num_grid = 0; num_grid < grids.size(); ++num_grid ) {
        for( std::size_t num_ot_grid = 0; num_ot_grid < grids.size(); ++num_ot_grid ) {
            if ( num_ot_grid == num_grid )
                continue;
            znodes.clear();
            for( TI num_dirac : grids[ num_grid ].dirac_indices )
                znodes.push_back( { zcoords_for( positions[ num_dirac ] ), num_dirac } );
            TODO;
            //            grids[ num_ot_grid ].octree.find_nodes( pos_to_find, [&]( TI cell_index, TI num_dirac ) {
            //                grids[ num_ot_grid ].cell_index_vs_dirac_number[ num_dirac ] = cell_index;
            //            } );
        }
    }
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
                } ); // , { TF( num_grid ) }
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
    c.size = 0;
    c.pos = max_point;
}

} // namespace Visitor
} // namespace PowerDiagram
