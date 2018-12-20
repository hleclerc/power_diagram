#include "../SpaceFunctions/Constant.h"
#include "ConvexPolyhedronAssembly.h"

namespace PowerDiagram {
namespace Bounds {

template<class Pc>
void ConvexPolyhedronAssembly<Pc>::add_box( Pt p0, Pt p1, TF coeff, TI cut_id ) {
    items.push_back( { { typename CP::Box{ p0, p1 }, cut_id }, coeff } );
}

template<class Pc>
const typename ConvexPolyhedronAssembly<Pc>::CP& ConvexPolyhedronAssembly<Pc>::englobing_convex_polyhedron() const {
    if ( items.size() == 1 )
        return items[ 0 ].polyhedron;
    TODO;
    static CP cp;
    return cp;
}

template<class Pc> template<class F>
void ConvexPolyhedronAssembly<Pc>::for_each_intersection( CP &cp, const F &f ) const {
    if ( items.size() == 1 )
        return f( cp, SpaceFunctions::Constant<TF>{ items[ 0 ].coeff } );
    TODO;
}

} // namespace Bounds
} // namespace PowerDiagram
