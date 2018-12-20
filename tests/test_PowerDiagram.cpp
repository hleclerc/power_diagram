#include "../src/PowerDiagram/PowerDiagram.h"
#include "../src/PowerDiagram/system/Tick.h"
#include <gtest/gtest.h>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -g3
//// nsmake cpp_flag -O6
//// nsmake lib_flag -g3
//// nsmake lib_flag -O6

template<class GridType>
void test_with_grid_type() {
    using TF = double;
    using CA = PcStd<2,TF,GridType>;
    using PT = typename PowerDiagram<CA>::PT;

    VtkOutput<1,TF> vo( { "num" } );

    // touching discs
    int n = 10; // 000;
    PowerDiagram<CA> pd;
    pd.add_box_shape( { -1.0, -1.0 }, { n + 0.0, n + 0.0 } );
    for( double i = 0; i < n; ++i ) {
        for( double j = 0; j < n; ++j ) {
            PT p { i, j };
            vo.add_point( p );
            pd.add_dirac( p, 1.0 );
        }
    }
    //    for( double i = 0.25; i < 1.0; i += 0.5 ) {
    //        for( double j = 0; j < n; ++j ) {
    //            PT p { i * n + 1.0 * rand() / RAND_MAX, j + 0.9 * rand() / RAND_MAX };
    //            pd.add_dirac( p, 1.0 );
    //            vo.add_point( p );
    //        }
    //    }
    // pd.set_ball_cut( true );
    tick.start( "measure" );
    std::vector<TF> measures( pd.nb_diracs() );
    pd.get_measures( measures.data() );
    tick.stop( "measure" );

    TF vol = 0;
    for( TF m : measures )
        vol += m;

    TF ref = pow( n + 1.0, 2 );
    P( ( vol - ref ) / ref );

    if ( n < 2000 )
        pd.display( vo );

    vo.save( "pd.vtk" );
}

TEST( PowerDiagram, RegularGrid ) {
    test_with_grid_type<SpZIndexGrid>();
}

//TEST( PowerDiagram, Grid ) {
//    using TF = double;
//    using CA = PcStd<2,TF>;
//    using PT = PowerDiagram<CA>::PT;
//    VtkOutput<1,TF> vo( { "num" } );

//    TF d_x = 0.5;
//    for( std::size_t i = 0; i < 1000000; ++i ) {
//        PT p = { 4.0 * rand() / RAND_MAX - 2.0, 2.0 * rand() / RAND_MAX - 1.0 };
//        PT d0_pos = { 0, 0 };
//        TF max_d1_weight = 2.0;
//        TF min_d0_weight = 1.0;
//        auto test_segment = [&]( PT A, PT B ) -> bool {
//            TF c2 = norm_2_p2( B - A );
//            TF c1 = dot( B - A, A - p );
//            TF c0 = norm_2_p2( A - p ) - norm_2_p2( d0_pos - p ) + min_d0_weight - max_d1_weight;
//            if ( c2 ) {
//                TF um = - c1 / c2;
//                if ( um > 0 && um < 1 && c2 * pow( um, 2 ) + 2 * c1 * um + c0 < 0 )
//                    return true;
//            }
//            return c0 < 0 || c2 + 2 * c1 + c0 < 0;
//        };
//        bool v = test_segment( { -d_x, 1.0 }, { +d_x, 1.0 } ) ||
//                 test_segment( { +d_x, 1.0 }, { +d_x, 2.0 } ) ||
//                 test_segment( { +d_x, 2.0 }, { -d_x, 2.0 } ) ||
//                 test_segment( { -d_x, 2.0 }, { -d_x, 1.0 } ) ;
//        vo.add_point( { p.x, p.y, 0.0 }, { TF( v ) } );
//    }
//    for( TF i = 0; i < 1; i += 0.01 ) {
//        PT p1{ - d_x + 2 * d_x * i, 1.0 };
//        PT pd = 0.5 * p1;
//        vo.add_lines( { pd - 2.0 * rot90( p1 ), pd + 2.0 * rot90( p1 ) }, { 2.0 } );
//    }

//    vo.add_lines( { { -d_x, 1.0, 0.0 }, { +d_x, 1.0, 0.0 }, { +d_x, 2.0, 0.0 }, { -d_x, 2.0, 0.0 }, { -d_x, 1.0, 0.0 } } );
//    vo.add_lines( { { -2.5, 0.0, 0.0 }, { +2.5, 0.0, 0.0 } } );
//    vo.save( "pd.vtk" );
//}

//template<int dim,class TF>
//void check_derivatives( PowerDiagram<dim,TF> &pd, TF epsilon = 1e-7 ) {
//    std::vector<std::vector<std::pair<size_t,TF>>> der( pd.nb_diracs() );
//    std::vector<std::vector<TF>> mod_measures( pd.nb_diracs() );
//    std::vector<TF> measures( pd.nb_diracs() );
//    pd.get_der_measures( measures.data(), der.data() );

//    for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
//        TF ow = pd.dirac( i ).weight;
//        pd.dirac( i ).weight -= epsilon;

//        mod_measures[ i ].resize( pd.nb_diracs() );
//        pd.get_measures( mod_measures[ i ].data() );

//        pd.dirac( i ).weight = ow;
//    }

//    for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
//        for( std::size_t j = 0; j < pd.nb_diracs(); ++j ) {
//            auto iter = std::find_if( der[ i ].begin(), der[ i ].end(), [&]( const std::pair<size_t,TF> &p ) { return p.first == j; } );
//            TF ex_der = iter == der[ i ].end() ? 0.0 : iter->second;
//            TF ap_der = ( measures[ i ] - mod_measures[ j ][ i ] ) / epsilon;
//            EXPECT_NEAR( ap_der, ex_der, 10 * epsilon );
//            // P( i, j, ap_der, ex_der, ap_der - ex_der );
//        }
//    }
//}

//TEST( PowerDiagram, LC2 ) {
//    using TF = long double;
//    using CA = PcStd<2,TF>;
//    VtkOutput<1,TF> vo( { "num" } );

//    // touching discs
//    PowerDiagram<CA> pd_0;
//    pd_0.add_bounding_simplex( { 0, 0 }, 10 );
//    pd_0.add_dirac( { 0, 0 }, 1 );
//    pd_0.add_dirac( { 1, 0 }, 1 );
//    pd_0.add_dirac( { 0, 1 }, 1 );
//    pd_0.add_dirac( { 1, 1 }, 1.5 );
//    pd_0.set_ball_cut( true );

//    check_derivatives( pd_0, TF( 1e-6 ) );
//    check_derivatives( pd_0, TF( 1e-7 ) );
//    check_derivatives( pd_0, TF( 1e-8 ) );
//    pd_0.display( vo );

//    // no intersection
//    PowerDiagram<CA> pd_1;
//    pd_1.add_bounding_simplex( { 0, 0 }, 10 );
//    pd_1.add_dirac( { 3, 0 }, std::pow( 0.4, 2 ) );
//    pd_1.add_dirac( { 4, 0 }, std::pow( 0.3, 2 ) );
//    pd_1.add_dirac( { 3, 1 }, std::pow( 0.2, 2 ) );
//    pd_1.set_ball_cut( true );

//    check_derivatives( pd_1 );
//    pd_1.display( vo );

//    // 1 intersection that makes "holes"
//    PowerDiagram<CA> pd_2;
//    pd_2.add_bounding_simplex( { 0, 0 }, 10 );
//    pd_2.add_dirac( { 6, 0 }, std::pow( 0.6, 2 ) );
//    pd_2.add_dirac( { 7, 0 }, std::pow( 0.6, 2 ) );
//    pd_2.add_dirac( { 6, 1 }, std::pow( 0.2, 2 ) );
//    pd_2.set_ball_cut( true );

//    check_derivatives( pd_2 );
//    pd_2.display( vo );

//    vo.save( "pd.vtk" );
//}


//TEST( PowerDiagram, LC3 ) {
//    using TF = long double;
//    using CA = PcStd<3,TF>;
//    const TF s = std::sqrt( 0.5 );
//    VtkOutput<1,TF> vo( { "num" } );

//    // touching spheres
//    PowerDiagram<CA> pd_0;
//    pd_0.add_bounding_simplex( { 0, 0, 0 }, 10 );
//    pd_0.add_dirac( { 0, 0, 0 }, 1 );
//    pd_0.add_dirac( { 1, 0, 0 }, 1 );
//    pd_0.add_dirac( { 0, 1, 0 }, 1 );
//    pd_0.add_dirac( { 1, 1, 0 }, 1.5 );
//    pd_0.set_ball_cut( true );

//    check_derivatives( pd_0, TF( 1e-6 ) );
//    check_derivatives( pd_0, TF( 1e-7 ) );
//    check_derivatives( pd_0, TF( 1e-8 ) );
//    pd_0.display( vo );

//    // no intersection
//    PowerDiagram<CA> pd_1;
//    pd_1.add_bounding_simplex( { 0, 0, 0 }, 10 );
//    pd_1.add_dirac( { 3, 0, 0 }, std::pow( 0.4, 2 ) );
//    pd_1.add_dirac( { 4, 0, 0 }, std::pow( 0.3, 2 ) );
//    pd_1.add_dirac( { 3, 1, 0 }, std::pow( 0.2, 2 ) );
//    pd_1.set_ball_cut( true );

//    check_derivatives( pd_1 );
//    pd_1.display( vo );

//    // 1 intersection that makes "holes"
//    PowerDiagram<CA> pd_2;
//    pd_2.add_bounding_simplex( { 0, 0, 0 }, 10 );
//    pd_2.add_dirac( { 6, 0, 0 }, std::pow( 0.6, 2 ) );
//    pd_2.add_dirac( { 7, 0, 0 }, std::pow( 0.6, 2 ) );
//    pd_2.add_dirac( { 6, 1, 0 }, std::pow( 0.2, 2 ) );
//    pd_2.set_ball_cut( true );

//    check_derivatives( pd_2 );
//    pd_2.display( vo );

//    vo.save( "pd.vtk" );
//}

////TEST( PowerDiagram, part_case_0 ) {
////    using TF = long double;
////    using PO = PowerDiagram<3,TF>::PT;

////    VtkOutput<1,TF> vo( { "num" } );

////    // touching spheres
////    PowerDiagram<TF> pd;
////    pd.add_dirac( { 0.861454, 0.861454, 0.490821 }, 0.321062 );
////    pd.add_dirac( { 0.43322 , 0.433221, 0.9246   }, 0.279741 );
////    pd.add_dirac( { 0.50188 , 0.50188 , 1.74263  }, 0.255672 );

////    TF l = 5, h = 5;
////    TF m = 1 * l / 4, n = 3 * l / 4;
////    PowerDiagram<TF>::LC lc( { 0, 0, 0 }, 100 );
////    lc.plane_cut( { 0, 0, 0 }, { -1, 0, 0 } );
////    lc.plane_cut( { 0, 0, 0 }, { 0, -1, 0 } );
////    lc.plane_cut( { 0, 0, 0 }, { 0, 0, -1 } );
////    lc.plane_cut( { l, l, h }, { +1, 0, 0 } );
////    lc.plane_cut( { l, l, h }, { 0, +1, 0 } );
////    lc.plane_cut( { l, l, h }, { 0, 0, +1 } );

////    lc.plane_cut( { m, m, 0 }, normalized( PO{ -0.5, 0, -1 } ) );
////    lc.plane_cut( { m, m, 0 }, normalized( PO{ 0, -0.5, -1 } ) );
////    lc.plane_cut( { n, n, 0 }, normalized( PO{ +0.5, 0, -1 } ) );
////    lc.plane_cut( { n, n, 0 }, normalized( PO{ 0, +0.5, -1 } ) );

////    pd.add_convex_bounds( lc );

////    //    check_derivatives( pd );

////    std::vector<TF> volumes_ex( pd.nb_diracs() );
////    std::vector<TF> volumes_ap( pd.nb_diracs() );
////    pd.get_volumes   ( volumes_ex.data() );
////    pd.get_volumes_mc( volumes_ap.data() );

////    P( volumes_ex );
////    P( volumes_ap );

////    pd.display( vo );
////    vo.save( "pd.vtk" );
////}

////TEST( PowerDiagram, part_case_1 ) {
////    using PO = PowerDiagram<TF>::PT;

////    VtkOutput<1,TF> vo( { "num" } );

////    // touching spheres
////    PowerDiagram<TF> pd;
////    pd.add_dirac( { 0.71082 , 0.710821,  0.575801 }, 0.300429 );
////    pd.add_dirac( { 1.62285 , 0.999315,  0.375209 }, 0.300064 );
////    pd.add_dirac( { 0.999313, 1.62285 ,  0.375211 }, 0.300064 );
////    pd.add_dirac( { 1.68829 , 1.68829 , -0.642599 }, 0.255261 );

////    TF l = 5, h = 5;
////    TF m = 1 * l / 4, n = 3 * l / 4;
////    PowerDiagram<TF>::LC lc( { 0, 0, 0 }, 100 );
////    lc.plane_cut( { 0, 0, 0 }, { -1, 0, 0 } );
////    lc.plane_cut( { 0, 0, 0 }, { 0, -1, 0 } );
////    lc.plane_cut( { 0, 0, 0 }, { 0, 0, -1 } );
////    lc.plane_cut( { l, l, h }, { +1, 0, 0 } );
////    lc.plane_cut( { l, l, h }, { 0, +1, 0 } );
////    lc.plane_cut( { l, l, h }, { 0, 0, +1 } );

////    lc.plane_cut( { m, m, 0 }, normalized( PO{ -0.5, 0, -1 } ) );
////    lc.plane_cut( { m, m, 0 }, normalized( PO{ 0, -0.5, -1 } ) );
////    lc.plane_cut( { n, n, 0 }, normalized( PO{ +0.5, 0, -1 } ) );
////    lc.plane_cut( { n, n, 0 }, normalized( PO{ 0, +0.5, -1 } ) );

////    pd.add_convex_bounds( lc );

////    pd.add_box_shape( { m, m, -2 }, { n, n,  0 } );
////    pd.add_box_shape( { 0, 0, -4 }, { l, l, -2 } );
////    //    check_derivatives( pd );

////    std::vector<TF> volumes_ex( pd.nb_diracs() );
////    pd.get_volumes   ( volumes_ex.data() );
////    P( volumes_ex );

////    std::vector<TF> volumes_ap( pd.nb_diracs() );
////    pd.get_volumes_mc( volumes_ap.data(), 1e8 );
////    P( volumes_ap );

////    pd.display( vo );
////    vo.save( "pd.vtk" );
////}

