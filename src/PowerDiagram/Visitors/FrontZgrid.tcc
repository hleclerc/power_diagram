#include "FrontZgrid.h"

template<class ZG>
FrontZgrid<ZG>::FrontZgrid( TI& op_count, std::vector<std::vector<TI>> &visited ) : op_count( op_count ), visited( visited ) {
}

template<class ZG> template<class Grid>
void FrontZgrid<ZG>::init( const std::vector<Grid> &grids, TI num_grid, TI num_cell ) {
    ++op_count;
    set_visited( grids, num_grid, num_cell );
    orig_cell_pos = grids[ num_grid ].cells[ num_cell ].pos;
}

template<class ZG> template<class Grid>
void FrontZgrid<ZG>::set_visited( const std::vector<Grid> &grids, TI num_grid, TI num_cell ) {
    visited[ num_grid ][ num_cell ] = op_count;
}

template<class ZG> template<class Cell>
typename FrontZgrid<ZG>::TF FrontZgrid<ZG>::dist( const Cell &cell ) {
    TF res = 0;
    for( int d = 0; d < dim; ++d ) {
        TF v = cell.pos[ d ] - orig_cell_pos[ d ]; // we need abs to avoid the overflow
        res += v * v;
    }
    return res;
}

template<class ZG> template<class Grid>
void FrontZgrid<ZG>::push_without_check( TI num_grid, TI num_cell, const std::vector<Grid> &grids ) {
    items.push( Item{ num_grid, num_cell, dist( grids[ num_grid ].cells[ num_cell ] ) } );
    set_visited( grids, num_grid, num_cell );
}

template<class ZG> template<class Grid>
void FrontZgrid<ZG>::push( TI num_grid, TI num_cell, const std::vector<Grid> &grids ) {
    if ( visited[ num_grid ][ num_cell ] != op_count )
        push_without_check( num_grid, num_cell, grids );
}

template<class ZG>
typename FrontZgrid<ZG>::Item FrontZgrid<ZG>::pop() {
    Item res = items.top();
    items.pop();
    return res;
}

template<class ZG>
bool FrontZgrid<ZG>::empty() const {
    return items.empty();
}
