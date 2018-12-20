#pragma once

#include "../ConvexPolyhedron3.h"
#include "../ConvexPolyhedron2.h"
#include <functional>

namespace PowerDiagram {
namespace Visitor {

/**
  Pc is expected to contain
  - static int dim
  - allow_ball_cut
  - TF => floating point type
  - TI => index type
*/
template<class Pc>
class RGrid {
public:
    // data from Pc
    static constexpr bool   allow_ball_cut        = Pc::allow_ball_cut;
    static constexpr int    dim                   = Pc::dim;
    using                   TF                    = typename Pc::TF;
    using                   TI                    = typename Pc::TI;

    // static definitions
    using                   CP2                   = ConvexPolyhedron2<Pc,TI>;
    using                   CP3                   = ConvexPolyhedron3<Pc,TI>;
    using                   CP                    = typename std::conditional<dim==3,CP3,CP2>::type;
    using                   Pi                    = std::array<TI,dim>;
    using                   Pt                    = typename CP::Pt;

    // methods
    /* ctor */              RGrid                 ( TF cell_size );

    template<class D> void  init                  ( const D &diracs );
    template<class D> void  for_each_laguerre_cell( const std::function<void( CP &lc, std::size_t num )> &f, const CP &starting_lc, const D &diracs ); ///< starting_lc can be a polygonal bound
    template<class D> void  for_each_laguerre_cell( const std::function<void( CP &lc, std::size_t num, int num_thread )> &f, const CP &starting_lc, const D &diracs ); ///< version with num_thread

    template<class V> void  display               ( V &vtk_output ) const; ///< for debug purpose

    TI                      nb_cells              () const { return acc_grid_lengths[ dim ]; }

    // values used by init
    bool                    ball_cut;

private:
    TI                      index_cell            ( Pt pos ) const;

    //
    std::array<TI,dim+1>    acc_grid_lengths;
    TF                      inv_cell_size;
    Pi                      grid_lengths;
    TF                      min_weight;
    TF                      max_weight;
    Pt                      min_point;
    Pt                      max_point;
    TF                      cell_size;

    // dirac indices
    std::vector<TI>         di_offsets;   ///< for each cell, position in di_values of indices
    std::vector<TI>         di_values;    ///<
};

} // namespace Visitor
} // namespace PowerDiagram

#include "RGrid.tcc"

