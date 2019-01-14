#include "../SpaceFunctions/Constant.h"
#include "ConvexPolyhedronAssembly.h"

namespace PowerDiagram {
namespace Bounds {

template<class Pc>
void ConvexPolyhedronAssembly<Pc>::add_convex_polyhedron( const std::vector<Pt> &positions, const std::vector<Pt> &normals, TF coeff, TI cut_id ) {
    using std::min;
    using std::max;
    Pt min_pos;
    Pt max_pos;
    Pt delta;
    for( std::size_t d = 0; d < dim; ++d ) {
        min_pos[ d ] = + std::numeric_limits<TF>::max();
        max_pos[ d ] = - std::numeric_limits<TF>::max();
        for( std::size_t i = 0; i < positions.size(); ++i ) {
            min_pos[ d ] = min( min_pos[ d ], positions[ i ][ d ] );
            max_pos[ d ] = max( max_pos[ d ], positions[ i ][ d ] );
        }
        delta[ d ] = max_pos[ d ] - min_pos[ d ];
    }

    CP cp{ typename CP::Box{ min_pos - delta, max_pos + delta }, cut_id };
    for( std::size_t i = 0; i < positions.size(); ++i )
        cp.plane_cut( positions[ i ], normalized( normals[ i ] ), cut_id );
    items.push_back( { cp, coeff } );
}

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

template<class Pc>
typename ConvexPolyhedronAssembly<Pc>::TF ConvexPolyhedronAssembly<Pc>::coeff_at( const Pt &pos ) const {
    TF res = 0;
    for( const Item &item : items )
        if ( item.polyhedron.contains( pos ) )
            res += item.coeff;
    return res;
}

template<class Pc> template<class F>
void ConvexPolyhedronAssembly<Pc>::for_each_intersection( CP &cp, const F &f ) const {
    if ( items.size() == 1 )
        return f( cp, SpaceFunctions::Constant<TF>{ items[ 0 ].coeff } );
    TODO;
}

template<class Pc> template<class V>
void ConvexPolyhedronAssembly<Pc>::display_boundaries( V &vtk_output ) const {
    for( const Item &item : items )
        item.polyhedron.display( vtk_output, { item.coeff }, false );
}

template<class Pc> template<class V>
void ConvexPolyhedronAssembly<Pc>::display_coeffs( V &vtk_output ) const {
    for( const Item &item : items )
        item.polyhedron.display( vtk_output, { item.coeff }, true );
}

} // namespace Bounds
} // namespace PowerDiagram
