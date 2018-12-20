#pragma once

#include "SpacePartioner.h"
#include "system/Assert.h"
#include "system/Span.h"
#include <functional>
#include "Point3.h"
#include "Point2.h"

/**
*/
struct SpRegularGrid {
    double max_weight_ratio = 2.0;
    double mul_cell_size    = 0.8;
};

/**
  Specialization of SpacePartioner for SpRegularGrid
*/
template<class Pc>
class SpacePartioner<Pc,SpRegularGrid> {
public:
    using                         TF                      = typename Pc::TF;
    using                         TI                      = typename Pc::TI;
    static constexpr int          dim                     = Pc::dim;
    using                         PT                      = typename std::conditional<dim==3,Point3<TF>,Point2<TF>>::type;
    using                         CellPos                 = std::array<std::size_t,dim>;
    struct                        CellHandle              { CellPos pos; TI index; int grid = 0; operator bool() const { return grid >= 0; } };
    struct                        VisitedCells            { void append( const CellHandle &cell ); bool operator()( const CellHandle &cell ) const; std::vector<std::vector<unsigned>> by_grid; unsigned cur_op_id = 0; };

    /*  */                        SpacePartioner          ( const SpRegularGrid *st ) : _st( *st ) {}

    template<class TV> void       init                    ( const TV &diracs );
    template<class Of> CellHandle ng_cell                 ( const CellHandle &cell, Of off );
    Span<TI>                      items_in                ( const CellHandle &cell ) const { return { _grids[ 0 ].items.data() + _grids[ 0 ].offsets[ cell.index + 0 ], _grids[ 0 ].items.data() + _grids[ 0 ].offsets[ cell.index + 1 ] }; }
    PT                            position                ( const CellHandle &cell ) const { PT res; for( int d = 0; d < dim; ++d ) res[ d ] = _grids[ cell.grid ].min_pt[ d ] + cell.pos[ d ] * _grids[ cell.grid ].cell_size; return res; }
    TF                            cell_size               ( const CellHandle &cell ) const { return _grids[ cell.grid ].cell_size; }
    TF                            min_weight              ( const CellHandle &cell ) const { return _grids[ cell.grid ].min_weight; }
    TF                            max_weight              ( const CellHandle &cell ) const { return _grids[ cell.grid ].max_weight; }
    PT                            cell_center             ( const CellHandle &cell ) const { PT res; for( int d = 0; d < dim; ++d ) res[ d ] = _grids[ cell.grid ].min_pt[ d ] + ( cell.pos[ d ] + TF( 0.5 ) ) * _grids[ cell.grid ].cell_size; return res; }
    TF                            min_sq_dist             ( const CellHandle &c0, const CellHandle &c1 ) const;
    void                          for_each_cell           ( const std::function<void( const CellHandle &cell, Span<TI> cell_indices )> &f, size_t job_num = 0, size_t job_len = 1 );
    template<class FU> void       for_each_neighbor       ( const CellHandle &cell, const CellHandle &orig, const FU &f ); ///< std::function<void( const CellHandle &ng_cell )>
    void                          init_visited_cells      ( VisitedCells &visited ) const;
    template<class FU> void       for_each_cell_node      ( const CellHandle &cell, const FU &f ); ///<
    template<class FU> bool       test_with_boundaries    ( const CellHandle &cell, const FU &test ) const;
    template<class FU> void       for_each_direct_neighbor( const CellHandle &cell, const FU &f ); ///< std::function<void( const CellHandle &ng_cell )>
    template<class FU> void       for_each_ring_2_neighbor( const CellHandle &cell, const FU &f ); ///< std::function<void( const CellHandle &ng_cell, tuple<LN0,...> )>


private:
    struct Grid {
        TI                        cell_index              ( CellPos p ) const { TI res = p[ 0 ]; for( int d = 1; d < dim; ++d ) res += cp_nb_divs[ d ] * p[ d ]; return res; }
        TI                        cell_index              ( PT      p ) const { return cell_index( cell_pos( p ) ); }
        CellPos                   cell_pos                ( PT      p ) const { CellPos res; for( ST d = 0; d < dim; ++d ) res[ d ] = ( p[ d ] - min_pt[ d ] ) / cell_size; return res; }

        TF                        min_weight;
        TF                        max_weight;
        CellPos                   cp_nb_divs; ///< cumulative product of nb divs. cp_nb_divs[ 0 ] always = 1
        TF                        cell_size;  ///<
        TI                        nb_cells;
        CellPos                   nb_divs;
        std::vector<TI>           offsets;
        PT                        min_pt;
        PT                        max_pt;
        std::vector<TI>           items;
    };

    std::vector<Grid>             _grids; ///< for each weight span
    const SpRegularGrid&          _st;
};

#include "SpRegularGrid.tcc"

