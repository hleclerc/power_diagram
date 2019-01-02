#pragma once

#include "../ConvexPolyhedron3.h"
#include "../ConvexPolyhedron2.h"
#include <functional>

namespace PowerDiagram {
namespace Visitor {

/**
  Pc is expected to contain
  - static int nb_bits_per_axis
  - static int dim
  - allow_ball_cut
  - TF => floating point type
  - TI => index type

  Comments:
  -
*/
template<class Pc>
class ZGrid {
public:
    // data from Pc
    static constexpr int    nb_bits_per_axis      = Pc::nb_bits_per_axis;
    static constexpr bool   allow_ball_cut        = Pc::allow_ball_cut;
    static constexpr int    dim                   = Pc::dim;
    using                   TF                    = typename Pc::TF;
    using                   TI                    = typename Pc::TI;

    // static definitions
    static constexpr int    sizeof_zcoords        = ( dim * nb_bits_per_axis + 7 ) / 8; ///< nb meaningful bytes in z-coordinates
    using                   CP2                   = ConvexPolyhedron2<Pc,TI>;
    using                   CP3                   = ConvexPolyhedron3<Pc,TI>;
    using                   CP                    = typename std::conditional<dim==3,CP3,CP2>::type;
    using                   Pt                    = typename CP::Pt;
    using                   TZ                    = std::uint64_t; ///< zcoords

    // methods
    /* ctor */              ZGrid                 ( std::size_t max_diracs_per_cell = 11, TF max_delta_weight_per_grid = 1e40 );

    void                    update                ( const Pt *positions, const TF *weights, std::size_t nb_diracs, bool positions_have_changed = true, bool weights_have_changed = true );
    int                     for_each_laguerre_cell( const std::function<void( CP &lc, std::size_t num )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, std::size_t nb_diracs, bool stop_if_void_lc = false ); ///< starting_lc can be a polygonal bound
    int                     for_each_laguerre_cell( const std::function<void( CP &lc, std::size_t num, int num_thread )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, std::size_t nb_diracs, bool stop_if_void_lc = false ); ///< version with num_thread

    template<class V> void  display               ( V &vtk_output ) const; ///< for debug purpose

    // values used by init
    TF                      max_delta_weight_per_grid;
    int                     max_diracs_per_cell;
    bool                    ball_cut;

private:
    struct                  Cell {
        TI                  dpc_offset;                           ///< offsets in grid.dpc_values
        TF                  size;                                 ///< integer size (must be before `pos`)
        Pt                  pos;                                  ///< integer coordinates
    };

    struct                  Grid {
        std::vector<TI>     cell_index_vs_dirac_number; ///< num cell for each dirac
        std::vector<TI>     dirac_indices;              ///< filled only if blocks.size() > 1
        std::vector<TI>     dpc_values;                 ///< diracs indices for each cell
        std::vector<TI>     ng_indices;                 ///< list of cell index of direct neighbors
        std::vector<TI>     ng_offsets;                 ///< offsets in ng_indices for each cell
        TF                  min_weight;
        TF                  max_weight;
        std::vector<Cell>   cells;
    };

    struct                  ZNode {
        void                write_to_stream( std::ostream &os ) const { os << zcoords << "[" << index << "]"; }
        TZ                  operator>>     ( int shift ) const { return zcoords >> shift; }
        TZ                  zcoords; ///<
        TI                  index;
    };

    template                <int num_axis,int _cur_bit=dim*nb_bits_per_axis-1>
    struct                  _ZcoordsZerosOnAxis {
        static constexpr TZ v_loc = TZ( _cur_bit % dim == num_axis ? 0 : 1 ) << _cur_bit;
        static constexpr TZ value = v_loc | _ZcoordsZerosOnAxis<num_axis,_cur_bit-1>::value;
    };

    template                <int num_axis>
    struct                  _ZcoordsZerosOnAxis<num_axis,-1> {
        static constexpr TZ value = 0;
    };

    /// Ex: axis = 0, dim = 3 (i.e. x) => 000... for level and free_bits ++ 001001001...
    template                <int num_axis,int _cur_bit = sizeof_zcoords - 1>
    struct                  _ZcoordsOnesOnAxis {
        static constexpr TZ v_loc = TZ( _cur_bit % dim == num_axis ? 1 : 0 ) << _cur_bit;
        static constexpr TZ value = v_loc | _ZcoordsOnesOnAxis<num_axis,_cur_bit-1>::value;
    };
    template                <int num_axis>
    struct                  _ZcoordsOnesOnAxis<num_axis,-1> {
        static constexpr TZ value = 0;
    };


    void                    fill_grid_using_zcoords( TI num_grid, const Pt *positions, const TF *weights, std::size_t nb_diracs );
    void                    repl_zcoords_by_ccoords( TI num_grid );
    void                    update_the_limits      ( const Pt *positions, const TF *weights, std::size_t nb_diracs );
    void                    update_neighbors       ( TI num_grid );
    void                    fill_the_grids         ( const Pt *positions, const TF *weights, std::size_t nb_diracs );
    template<class C> TZ    zcoords_for            ( const C &pos ); ///< floating point position
    template<int d>   TZ    ng_zcoord              ( TZ zcoords, TZ off, N<d> ) const;
    bool                    may_cut                ( const CP &lc, TI i0, const Grid &cr_grid, const Cell &cr_cell, const Pt *positions, const TF *weights );

    // true if a dirac in b1 (given max_weight in b1 and its neighbors) may cut lc

    // tmp
    using                   RsTmp                  = std::vector<std::array<std::size_t,256>>;
    RsTmp                   rs_tmps;
    std::vector<ZNode>      znodes;
    std::vector<ZNode>      zcells;

    //
    TF                      step_length;
    TF                      grid_length;
    TF                      min_weight;
    TF                      max_weight;
    Pt                      min_point;
    Pt                      max_point;
    std::vector<Grid>       grids;                     ///< for each weight span
};

} // namespace Visitor
} // namespace PowerDiagram

#include "ZGrid.tcc"

