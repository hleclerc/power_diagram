#pragma once

#include "../ConvexPolyhedron2.h"
#include "../ConvexPolyhedron3.h"

namespace PowerDiagram {
namespace Bounds {

/**
  Currently, only support constant coeffs per polyhedron
*/
template<class Pc>
class ConvexPolyhedronAssembly {
public:
    static constexpr int   dim                        = Pc::dim;
    using                  TF                         = typename Pc::TF;
    using                  TI                         = typename Pc::TI;

    using                  CP2                        = ConvexPolyhedron2<Pc,TI>;
    using                  CP3                        = ConvexPolyhedron3<Pc,TI>;
    using                  CP                         = typename std::conditional<dim==3,CP3,CP2>::type;
    using                  Pt                         = typename CP::Pt;

    // modifications
    void                   add_box                    ( Pt p0, Pt p1, TF coeff = 1.0, TI cut_id = -1 );

    // info
    const CP&              englobing_convex_polyhedron() const;
    template<class F> void for_each_intersection      ( CP &cp, const F &f ) const; ///< f( ConvexPolyhedron, SpaceFunction )
    template<class V> void display_boundaries         ( V &vtk_output ) const;
    template<class V> void display_coeffs             ( V &vtk_output ) const;

    //
private:
    struct                 Item                       { CP polyhedron; TF coeff; };

    std::vector<Item>      items;
};

} // namespace Bounds
} // namespace PowerDiagram

#include "ConvexPolyhedronAssembly.tcc"