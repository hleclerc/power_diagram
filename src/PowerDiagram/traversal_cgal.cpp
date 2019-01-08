//// nsmake avoid_inc CGAL/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include "traversal_cgal.h"

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -O5
//// nsmake lib_name CGAL
//// nsmake lib_name gmp

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using R = CGAL::Regular_triangulation_2<K>;

double traversal_cgal( const double *positions, const double *weights, int nb_diracs ) {
    std::vector<R::Weighted_point> diracs;
    for( int i = 0; i < nb_diracs; ++i )
        diracs.push_back( { { positions[ 2 * i + 0 ], positions[ 2 * i + 1 ] }, weights[ i ] } );

    R rt( diracs.begin(), diracs.end() );
    // rt.is_valid();

    double s = 0;
    for( auto v = rt.all_vertices_begin(); v != rt.all_vertices_end(); ++v ) {
        auto circulator = rt.incident_faces( v ), done( circulator );
        do {
            double v = circulator->vertex( 0 )->point().point().x();
            s += v;
        } while( ++circulator != done );
    }

    return s;
}
