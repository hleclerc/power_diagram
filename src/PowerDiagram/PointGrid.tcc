#include "PointGrid.h"

template<class Item,class TF,int dim,class ST>
void PointGrid<Item,TF,dim,ST>::init( const std::function<void (const std::function<void(PT, Item)> &)> &f, std::size_t approx_nb_cells ) {
    // get min/max of positions
    ST nb_items = 0;
    for( ST d = 0; d < dim; ++d ) {
        _min[ d ] = + std::numeric_limits<TF>::max();
        _max[ d ] = - std::numeric_limits<TF>::max();
    }
    f( [&]( PT pos, Item item ) {
        _min = min( _min, pos );
        _max = max( _max, pos );
        ++nb_items;
    } );

    // correction of approx_nb_cells if necessary
    if ( approx_nb_cells == 0 )
        approx_nb_cells = nb_items;

    // get cell size (-> `_s`)
    int p = 0;
    TF  m = 1;
    for( ST d = 0; d < dim; ++d ) {
        if ( TF l = _max[ d ] - _min[ d ] ) {
            m *= l;
            ++p;
        }
    }
    _cell_size = std::pow( m / approx_nb_cells, 1.0 / p );

    // we want the last points to be included
    _cell_size *= ( 1.0 + 16 * std::numeric_limits<TF>::epsilon() );

    // _nb_cells and _nb_divs
    _nb_cells = 1;
    for( ST d = 0; d < dim; ++d ) {
        _nb_divs[ d ] = 1 + ( _max[ d ] - _min[ d ] ) / _cell_size;
        _cp_nb_divs[ d ] = _nb_cells;
        _nb_cells *= _nb_divs[ d ];
    }

    // preparation of _offsets (get nb items / cell)
    _offsets.resize( _nb_cells + 1 );
    for( ST i = 0; i < _nb_cells; ++i )
        _offsets[ i ] = 0;
    f( [&]( PT pos, Item item ) {
        ++_offsets[ cell_index( pos ) ];
    } );

    // suffix scan
    for( ST i = 0, a = 0; i < _nb_cells; ++i ) {
        ST o = _offsets[ i ];
        _offsets[ i ] = a;
        a += o;
    }
    _offsets.back() = nb_items;

    // fill _items
    _items.resize( nb_items );
    f( [&]( PT pos, Item item ) {
        _items[ _offsets[ cell_index( pos ) ]++ ] = item;
    } );

    // correction of _offsets
    for( ST i = _offsets.size(); --i; )
        _offsets[ i ] = _offsets[ i - 1 ];
    _offsets[ 0 ] = 0;
}

template<class Item,class TF,int dim,class ST>
ST PointGrid<Item,TF,dim,ST>::cell_index( CellPos p ) const{
    ST res = p[ 0 ];
    for( ST d = 1; d < dim; ++d )
        res += _cp_nb_divs[ d ] * p[ d ];
    return res;
}

template<class Item, class TF, int dim,class ST>
typename PointGrid<Item,TF,dim,ST>::CellPos PointGrid<Item,TF,dim,ST>::cell_pos( PT p ) const {
    CellPos res;
    for( ST d = 0; d < dim; ++d )
        res[ d ] = ( p[ d ] - _min[ d ] ) / _cell_size;
    return res;
}

template<class Item, class TF, int dim, class ST>
bool PointGrid<Item,TF,dim,ST>::is_outside( CellPos p ) const {
    for( ST d = 0; d < dim; ++d )
        if ( p[ d ] < 0 || p[ d ] >= _nb_divs[ d ] )
            return true;
    return false;
}

template<class Item,class TF,int dim,class ST>
void PointGrid<Item,TF,dim,ST>::for_each_cell( const std::function<void( CellPos cell_pos, const Item *items, ST nb_items )> &f, size_t ind_job, size_t nb_jobs ) const {
    ST d_max_extent = 0;
    for( ST d = 1; d < dim; ++d )
        if ( _nb_divs[ d_max_extent ] < _nb_divs[ d ] )
            d_max_extent = d;

    CellPos beg, end;
    for( ST d = 0; d < dim; ++d ) {
        if ( d == d_max_extent ) {
            beg[ d ] = ( ind_job + 0 ) * _nb_divs[ d ] / nb_jobs;
            end[ d ] = ( ind_job + 1 ) * _nb_divs[ d ] / nb_jobs;
            if ( beg[ d ] == end[ d ] )
                return;
        } else {
            beg[ d ] = 0;
            end[ d ] = _nb_divs[ d ];
        }
    }

    CellPos cell_pos;
    for( ST d = 0; d < dim; ++d )
        cell_pos[ d ] = beg[ d ];

    for( ST num = _cp_nb_divs[ d_max_extent ] * beg[ d_max_extent ]; num < _nb_cells; ++num ) {
        f( cell_pos, _items.data() + _offsets[ num + 0 ], _offsets[ num + 1 ] - _offsets[ num + 0 ] );

        for( ST d = 0; d < dim; ++d ) {
            if ( ++cell_pos[ d ] < end[ d ] )
                break;
            if ( d == d_max_extent )
                num += ( ( _nb_divs[ d ] - end[ d ] ) + beg[ d ] ) * _cp_nb_divs[ d_max_extent ];
            cell_pos[ d ] = beg[ d ];
        }
    }
}
