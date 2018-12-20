#include "../src/PowerDiagram/Visitors/IntSpiral.h"
#include "catch_main.h"

//template<int dim>
//void test_dist() {
//    // check that points do not appear twice
//    using Pt = typename IntSpiral<dim>::Pt;
//    const int rd = 65;
//    std::set<PT> points;
//    for( IntSpiral<dim> is; is.rd2() < rd * rd; ++is ) {
//        EXPECT_EQ( 0, points.count( is.pos() ) );
//        points.insert( is.pos() );
//    }

//    // check that all expected points are present
//    PT point;
//    for( ST d = 0; d < dim; ++d )
//        point[ d ] = - rd;
//    for( ST num = 0; num < std::pow( 2 * rd + 1, dim ); ++num ) {
//        if ( IntSpiral<dim>::sq_dist( point ) < rd * rd )
//            EXPECT_EQ( 1, points.count( point ) );
//        for( ST d = 0; d < dim; ++d ) {
//            if ( ++point[ d ] <= rd )
//                break;
//            point[ d ] = - rd;
//        }
//    }
//}

TEST_CASE( "IntSpiral 2D" ) {
    IntSpiral<2> is;
    is.for_each_until( 5, []( auto p ) {
        P( p, std::sqrt( p[0] * p[0] + p[1] * p[1] ) );
    } );
}

