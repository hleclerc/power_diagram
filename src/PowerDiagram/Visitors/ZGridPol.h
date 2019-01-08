#pragma once

#include "../ConvexPolyhedron3.h"
#include "../ConvexPolyhedron2.h"
#include <functional>

namespace PowerDiagram {
namespace Visitor {

/**
  Pc is expected to contain
  - static int dim              = space dimensionnality
  - using      TF               = floating point type
  - using      TI               = index type

  Comments:
  -
*/
template<class Pc>
class ZGridPol {
public:
    // data from Pc
    static constexpr int    dim                    = Pc::dim;
    using                   TF                     = typename Pc::TF;
    using                   TI                     = typename Pc::TI;

    // static definitions
    static constexpr int    nb_bits_per_axis       = 63 / dim;
    static constexpr int    sizeof_zcoords         = ( dim * nb_bits_per_axis + 7 ) / 8; ///< nb meaningful bytes in z-coordinates
    using                   CP2                    = ConvexPolyhedron2<Pc,TI>;
    using                   CP3                    = ConvexPolyhedron3<Pc,TI>;
    using                   CP                     = typename std::conditional<dim==3,CP3,CP2>::type;
    using                   Pt                     = typename CP::Pt;
    using                   TZ                     = std::uint64_t; ///< zcoords

    // methods
    /* ctor */              ZGridPol               ( TI max_diracs_per_cell = 11 );

    void                    update                 ( const Pt *positions, const TF *weights, TI nb_diracs, bool positions_have_changed = true, bool weights_have_changed = true );

    int                     for_each_laguerre_cell ( const std::function<void( CP &lc, TI num )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc = false ); ///< starting_lc can be a polygonal bound
    int                     for_each_laguerre_cell ( const std::function<void( CP &lc, TI num, int num_thread )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc = false ); ///< version with num_thread
    template<class V> void  display                ( V &vtk_output ) const; ///< for debug purpose


    // values used by init
    int                     max_diracs_per_cell;

private:
    struct                  Cell {
        TI                  dpc_offset;                 ///< offsets in grid.dpc_values
        TZ                  zcoords;
        TF                  size;
        Pt                  pos;
    };

    struct                  ZNode {
        void                write_to_stream        ( std::ostream &os ) const { os << zcoords << "[" << index << "]"; }
        TZ                  operator>>             ( int shift ) const { return zcoords >> shift; }

        TZ                  zcoords;
        TI                  index;
    };

    struct                  Poly {
        enum {              nc = dim * ( dim + 3 ) / 2 + 1 };
        std::array<TF,nc>   coeffs;
    };

    void                    _subdivide_add_poly_rec( Pt *positions, TF *weights, TI nb_diracs, Pt p0, Pt p1, int num );
    void                    _update_the_polynomials( const Pt *positions, const TF *weights, TI nb_diracs );
    void                    _update_the_limits     ( const Pt *positions, TI nb_diracs );
    void                    _fill_the_grid         ( const Pt *positions, TI nb_diracs );
    void                    _update_ngs            ( const Pt *positions, TI nb_diracs );

    template<class C> TZ    _zcoords_for           ( const C &pos ); ///< floating point position
    template<int axis> TZ   _ng_zcoord             ( TZ zcoords, TZ off, N<axis> ) const;

    //    void                    fill_grid_using_zcoords( TI num_grid, const Pt *positions, const TF *weights, TI nb_diracs );
    //    void                    repl_zcoords_by_ccoords( TI num_grid );
    //    void                    find_englobing_cousins ( TI num_grid, const Pt *positions ); ///< find englobing cells for each dirac (and for each grid). Must be done after repl_zcoords_by_ccoords
    //    void                    update_the_limits      ( const Pt *positions, const TF *weights, TI nb_diracs );
    //    void                    update_neighbors       ( TI num_grid );
    //    void                    fill_the_grids         ( const Pt *positions, const TF *weights, TI nb_diracs );
    //    template<int d>   TZ    ng_zcoord              ( TZ zcoords, TZ off, N<d> ) const;
    //    bool                    may_cut                ( const CP &lc, TI i0, const Grid &cr_grid, const Cell &cr_cell, const Pt *positions, const TF *weights );

    // true if a dirac in b1 (given max_weight in b1 and its neighbors) may cut lc

    // tmp
    using                   _RsTmp                 = std::vector<std::array<std::size_t,256>>;
    _RsTmp                  _tmp_rs;
    std::vector<ZNode>      _tmp_zn;               ///< tmp znodes

    // info from the diracs
    const ZNode*            _sorted_zcoords;       ///< zcoords for each dirac
    std::vector<ZNode>      _zcoords;              ///< zcoords for each dirac

    // grid
    std::vector<Poly>       _polynomials;          ///<
    std::vector<TI>         _dpc_values;           ///< diracs indices for each cell
    std::vector<TI>         _ng_offsets;           ///< offsets in ng_values for each cell
    std::vector<TI>         _ng_values;            ///< list of cell index of direct neighbors
    std::vector<Cell>       _cells;                ///< zcoords for each dirac

    //
    TF                      _inv_step_length;      ///< inv of smallest possible cell size
    TF                      _step_length;          ///< inv of smallest possible cell size
    TF                      _grid_length;          ///< width, height and depth of the grid (which is necessarily a square)
    Pt                      _min_point;
    Pt                      _max_point;
};

} // namespace Visitor
} // namespace PowerDiagram

#include "ZGridPol.tcc"

