#include "system/StaticCrossProdOfRanges.h"
#include "system/StaticRange.h"
#include "SpRegularGrid.h"
#include "system/LN.h"

template<class Pc> template<class TV>
void SpacePartioner<Pc,SpRegularGrid>::init( const TV &diracs ) {
    using std::min;
    using std::max;

    _grids.resize( 1 );
    Grid &grid = _grids[ 0 ];

    // get min/max of positions and weights
    grid.min_weight = + std::numeric_limits<TF>::max();
    grid.max_weight = - std::numeric_limits<TF>::max();
    for( int d = 0; d < dim; ++d ) {
        grid.min_pt[ d ] = + std::numeric_limits<TF>::max();
        grid.max_pt[ d ] = - std::numeric_limits<TF>::max();
    }
    for( const auto &dirac : diracs ) {
        grid.min_weight = min( grid.min_weight, dirac.weight );
        grid.max_weight = max( grid.max_weight, dirac.weight );
        grid.min_pt     = min( grid.min_pt    , dirac.pos    );
        grid.max_pt     = max( grid.max_pt    , dirac.pos    );
    }

    //

    // get cell size
    int p = 0;
    TF  m = 1;
    for( int d = 0; d < dim; ++d ) {
        if ( TF l = grid.max_pt[ d ] - grid.min_pt[ d ] ) {
            m *= l;
            ++p;
        }
    }
    grid.cell_size = _st.mul_cell_size * std::pow( m / diracs.size(), 1.0 / p ) *
                     ( 1.0 + 16 * std::numeric_limits<TF>::epsilon() );

    // _nb_cells and _nb_divs
    grid.nb_cells = 1;
    for( int d = 0; d < dim; ++d ) {
        grid.nb_divs[ d ] = 1 + ( grid.max_pt[ d ] - grid.min_pt[ d ] ) / grid.cell_size;
        grid.cp_nb_divs[ d ] = grid.nb_cells;
        grid.nb_cells *= grid.nb_divs[ d ];
    }

    // preparation of _offsets (get nb items / cell)
    grid.offsets.resize( grid.nb_cells + 1 );
    for( TI i = 0; i < grid.nb_cells; ++i )
        grid.offsets[ i ] = 0;
    for( TI i = 0; i < diracs.size(); ++i )
        ++grid.offsets[ grid.cell_index( diracs[ i ].pos ) ];

    // suffix scan
    for( ST i = 0, a = 0; i < grid.nb_cells; ++i ) {
        ST o = grid.offsets[ i ];
        grid.offsets[ i ] = a;
        a += o;
    }
    grid.offsets.back() = diracs.size();

    // fill _items (temporary modification of grid.offsets)
    grid.items.resize( diracs.size() );
    for( TI i = 0; i < diracs.size(); ++i )
        grid.items[ grid.offsets[ grid.cell_index( diracs[ i ].pos ) ]++ ] = i;

    // correction of grid.offsets
    for( ST i = grid.offsets.size(); --i; )
        grid.offsets[ i ] = grid.offsets[ i - 1 ];
    grid.offsets[ 0 ] = 0;
}

template<class Pc>
void SpacePartioner<Pc,SpRegularGrid>::for_each_cell( const std::function<void( const CellHandle &cell, Span<TI> cell_indices )> &f, size_t job_num, size_t job_len ) {
    for( const Grid &grid : _grids ) {
        // splitting of the space (for separation between several jobs)
        ST d_max_extent = 0;
        for( ST d = 1; d < dim; ++d )
            if ( grid.nb_divs[ d_max_extent ] < grid.nb_divs[ d ] )
                d_max_extent = d;

        CellPos beg, end;
        for( int d = 0; d < dim; ++d ) {
            if ( d == d_max_extent ) {
                beg[ d ] = ( job_num + 0 ) * grid.nb_divs[ d ] / job_len;
                end[ d ] = ( job_num + 1 ) * grid.nb_divs[ d ] / job_len;
                if ( beg[ d ] == end[ d ] )
                    return;
            } else {
                beg[ d ] = 0;
                end[ d ] = grid.nb_divs[ d ];
            }
        }

        CellPos cell_pos;
        for( int d = 0; d < dim; ++d )
            cell_pos[ d ] = beg[ d ];

        for( TI num = grid.cp_nb_divs[ d_max_extent ] * beg[ d_max_extent ]; num < grid.nb_cells; ++num ) {
            f( { cell_pos, num }, { grid.items.data() + grid.offsets[ num + 0 ], grid.items.data() + grid.offsets[ num + 1 ] } );

            for( int d = 0; d < dim; ++d ) {
                if ( ++cell_pos[ d ] < end[ d ] )
                    break;
                if ( d == d_max_extent )
                    num += ( ( grid.nb_divs[ d ] - end[ d ] ) + beg[ d ] ) * grid.cp_nb_divs[ d_max_extent ];
                cell_pos[ d ] = beg[ d ];
            }
        }
    }
}

template<class Pc> template<class FU>
void SpacePartioner<Pc,SpRegularGrid>::for_each_neighbor( const CellHandle &cell, const CellHandle &orig, const FU &f ) {
    const Grid &grid = _grids[ cell.grid ];

    CellHandle new_cell;
    StaticCrossProdOfRanges<dim,-1,2>::for_each( [&]( auto off ) {
        if ( off.has_only_zeros )
            return;

        bool ok = true;
        new_cell = cell;
        StaticRange<dim>::for_each( [&]( auto d ) {
            if ( constexpr int o = off.at( d ) ) {
                new_cell.pos[ d ] += o;
                new_cell.index += ( d ? o * grid.cp_nb_divs[ d ] : o );
                ok &= new_cell.pos[ d ] < grid.nb_divs[ d ];
            }
        } );

        if ( ok )
            f( new_cell );
    } );
}

template<class Pc>
typename SpacePartioner<Pc,SpRegularGrid>::TF SpacePartioner<Pc,SpRegularGrid>::min_sq_dist( const CellHandle &c0, const CellHandle &c1 ) const {
    using std::pow;
    using std::abs;
    if ( c0.grid == c1.grid ) {
        TF res = 0;
        for( int d = 0; d < dim; ++d ) {
            int o = abs( int( c1.pos[ d ] - c0.pos[ d ] ) );
            res += pow( TF( o - ( o != 0 ) ), 2 );
        }
        return pow( _grids[ c0.grid ].cell_size, 2 ) * res;
    }

    TODO;
    return {};
}

template<class Pc>
void SpacePartioner<Pc,SpRegularGrid>::init_visited_cells( VisitedCells &visited ) const {
    visited.by_grid.resize( _grids.size() );
    for( std::size_t i = 0; i < _grids.size(); ++i )
        visited.by_grid[ i ].resize( _grids[ i ].nb_cells, visited.cur_op_id );
    ++visited.cur_op_id;
}

template<class Pc>
void SpacePartioner<Pc, SpRegularGrid>::VisitedCells::append( const CellHandle &cell ) {
    by_grid[ cell.grid ][ cell.index ] = cur_op_id;
}

template<class Pc>
bool SpacePartioner<Pc, SpRegularGrid>::VisitedCells::operator()( const CellHandle &cell ) const {
    return by_grid[ cell.grid ][ cell.index ] == cur_op_id;
}

template<class Pc> template<class FU>
void SpacePartioner<Pc, SpRegularGrid>::for_each_direct_neighbor( const CellHandle &cell, const FU &f ) {
    const Grid &grid = _grids[ cell.grid ];

    auto test = [&]( auto ln ) {
        bool ok = true;
        CellHandle new_cell = cell;
        ln.for_each_with_cpt( [&]( auto off, auto d ) {
            if ( off ) {
                new_cell.pos[ d ] += int( off );
                ok &= new_cell.pos[ d ] < grid.nb_divs[ d ];
                new_cell.index += grid.cp_nb_divs[ d ] * int( off );
            }
        } );
        if ( ok )
            f( new_cell );
    };

    test( LN<-1, 0>() );
    test( LN<+1, 0>() );
    test( LN< 0,-1>() );
    test( LN< 0,+1>() );

    test( LN<-1,-1>() );
    test( LN<-1,+1>() );
    test( LN<+1,-1>() );
    test( LN<+1,+1>() );
}

template<class Pc> template<class FU>
void SpacePartioner<Pc, SpRegularGrid>::for_each_ring_2_neighbor( const CellHandle &cell, const FU &f ) {
    const Grid &grid = _grids[ cell.grid ];

    auto test = [&]( auto ln, auto next_ln_tuple ) {
        bool ok = true;
        CellHandle new_cell = cell;
        ln.for_each_with_cpt( [&]( auto off, auto d ) {
            if ( off ) {
                new_cell.pos[ d ] += int( off );
                ok &= new_cell.pos[ d ] < grid.nb_divs[ d ];
                new_cell.index += grid.cp_nb_divs[ d ] * int( off );
            }
        } );
        if ( ok )
            f( new_cell, next_ln_tuple );
    };

    // on the axes
    test( LN<-2, 0>(), std::tuple<LN<-3,-1>,LN<-3, 0>,LN<-3,+1>>() );
    test( LN<+2, 0>(), std::tuple<LN<+3,-1>,LN<+3, 0>,LN<+3,+1>>() );
    test( LN< 0,-2>(), std::tuple<LN<-1,-3>,LN< 0,-3>,LN<+1,-3>>() );
    test( LN< 0,+2>(), std::tuple<LN<-1,+3>,LN< 0,+3>,LN<+1,+3>>() );

    // off the axes
    test( LN<-2,-1>(), std::tuple<LN<-3,-1>,LN<-2,-2>>() );
    test( LN<-2,+1>(), std::tuple<LN<-3,+1>,LN<-2,+2>>() );
    test( LN<+2,-1>(), std::tuple<LN<+3,-1>,LN<+2,-2>>() );
    test( LN<+2,+1>(), std::tuple<LN<+3,+1>,LN<+2,+2>>() );
    test( LN<-1,-2>(), std::tuple<LN<-1,-3>,LN<-2,-2>>() );
    test( LN<+1,-2>(), std::tuple<LN<+1,-3>,LN<+2,-2>>() );
    test( LN<-1,+2>(), std::tuple<LN<-1,+3>,LN<-2,+2>>() );
    test( LN<+1,+2>(), std::tuple<LN<+1,+3>,LN<+2,+2>>() );

}

template<class Pc> template<class FU>
void SpacePartioner<Pc,SpRegularGrid>::for_each_cell_node( const CellHandle &cell, const FU &f ) {
    using std::pow;
    StaticRange<( 1 << dim )>::for_each( [&]( auto n ) {
        PT res;
        StaticRange<dim>::for_each( [&]( auto d ) {
            constexpr bool o = int( n ) & ( 1 << int( d ) );
            res[ d ] = _grids[ cell.grid ].min_pt[ d ] + ( cell.pos[ d ] + o ) * _grids[ cell.grid ].cell_size;
        } );
        f( res );
    } );
}

template<class Pc> template<class FU>
bool SpacePartioner<Pc,SpRegularGrid>::test_with_boundaries( const CellHandle &cell, const FU &test ) const {
    TF cs = _grids[ cell.grid ].cell_size;
    PT p0 = position( cell );
    PT p1 { p0.x + 1 * cs, p0.y + 0 * cs };
    PT p2 { p0.x + 1 * cs, p0.y + 1 * cs };
    PT p3 { p0.x + 0 * cs, p0.y + 1 * cs };
    return test( p0, p1 ) ||
           test( p1, p2 ) ||
           test( p2, p3 ) ||
           test( p3, p0 ) ;
}

template<class Pc> template<class Of>
typename SpacePartioner<Pc,SpRegularGrid>::CellHandle SpacePartioner<Pc,SpRegularGrid>::ng_cell( const CellHandle &cell, Of off ) {
    CellHandle res = cell;
    StaticRange<dim>::for_each( [&]( auto d ) {
        if ( constexpr int o = off.at( d ) ) {
            res.pos[ d ] += o;
            if ( res.pos[ d ] >= _grids[ cell.grid ].nb_divs[ d ] )
                res.grid = -1;
            res.index += o * _grids[ cell.grid ].cp_nb_divs[ d ];
        }
    } );
    return res;
}
