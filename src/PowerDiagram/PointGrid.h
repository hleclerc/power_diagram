#pragma once

#include "system/Span.h"
#include <functional>
#include "Point3.h"
#include "Point2.h"

/**
*/
template<class Item=std::size_t,class TF=double,int dim=3,class ST=std::size_t>
class PointGrid {
public:
    using             PT           = typename std::conditional<dim==3,Point3<TF>,Point2<TF>>::type;
    using             CellPos      = std::array<ST,dim>;

    /// Ex call: point_grid.init( [&]( const std::function<void( PT pos,Item item )> &cb ) { cb( pos_0, item_0 ); cb( pos_1, item_1 ); } );
    /// if approx_nb_cells == 0, we take the number of items
    void              init         ( const std::function<void( const std::function<void( PT pos, Item item )> &cb )> &f, std::size_t approx_nb_cells = 0 );

    ST                cell_index   ( PT p ) const { return cell_index( cell_pos( p ) ); }
    ST                cell_index   ( CellPos p ) const;
    TF                cell_size    () const { return _cell_size; }
    CellPos           cell_pos     ( PT p ) const;

    Span<Item>        items_in     ( ST cell_index ) const { return { _items.data() + _offsets[ cell_index + 0 ], _items.data() + _offsets[ cell_index + 1 ] }; }
    Span<Item>        items_in     ( CellPos p ) const { return is_outside( p ) ? Span<Item>{ 0, 0 } : items_in( cell_index( p ) ); }
    bool              is_outside   ( CellPos p ) const;

    void              for_each_cell( const std::function<void( CellPos cell_pos, const Item *items, ST nb_items )> &f, size_t ind_job = 0, size_t nb_jobs = 1 ) const;

private:
    CellPos           _cp_nb_divs; ///< cumulative product of nb divs
    TF                _cell_size;  ///<
    ST                _nb_cells;
    CellPos           _nb_divs;
    std::vector<ST>   _offsets;
    std::vector<Item> _items;
    PT                _min;
    PT                _max;
};


#include "PointGrid.tcc"
