#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/Visitors/ZGridPol.h"
#include "../src/PowerDiagram/triangulation.h"
#include "catch_main.h"

//// nsmake cpp_flag -march=native

template<class Pt,class TF>
void read_xyw( std::vector<Pt> &positions, std::vector<TF> &weights, std::string file ) {
    using std::min;

    std::ifstream f( file.c_str() );
    positions.resize( 0 );
    weights.resize( 0 );
    double x, y, w;
    while ( f >> x >> y >> w ) {
        positions.push_back( { x, y } );
        weights.push_back( w );
    }

    TF min_w = std::numeric_limits<TF>::max();
    for( TF w : weights )
        min_w = min( min_w, w );
    for( TF &w : weights )
        w -= min_w;
}

TEST_CASE( "ZGrid pol measures" ) {
    struct Pc     { enum { dim = 2 }; using TI = std::size_t; using TF = double; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    using  Grid   = PowerDiagram::Visitor::ZGridPol<Pc>;

    std::vector<Grid::Pt> positions;
    std::vector<Grid::TF> weights;
    read_xyw( positions, weights, "vtk/random_100000.xyw" );
    //    for( std::size_t i = 0; i < 100; ++i ) {
    //        double x = 1.0 * rand() / RAND_MAX;
    //        double y = 1.0 * rand() / RAND_MAX;
    //        positions.push_back( { x, y } );
    //        weights.push_back( sin( 2 * ( x + y ) ) );
    //    }

    //    Bounds bounds;
    //    bounds.add_box( { 0, 0 }, { 1, 1 }, 1.0, -1 );

    Grid grid( 40 );
    grid.update( positions.data(), weights.data(), positions.size() );

    VtkOutput<1> vo( { "num" } );
    grid.display( vo );
    for( std::size_t i = 0; i < weights.size(); ++i )
        vo.add_point( Point3<double>{ positions[ i ].x, positions[ i ].y, weights[ i ] }, { double( i ) } );
    vo.save( "vtk/grid.vtk" );

    //    VtkOutput<1> vo( { "num" } );
    //    std::mutex mutex;
    //    std::vector<double> volumes( positions.size() );
    //    std::vector<double> ext_perimeters( positions.size(), 0.0 );
    //    std::map<std::pair<std::size_t,std::size_t>,std::vector<double>> bms;
    //    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
    //        bounds.for_each_intersection( lc, [&]( auto &cp, auto space_func ) {
    //            mutex.lock();
    //            volumes[ num_dirac_0 ] += cp.integration( space_func );
    //            cp.display( vo, { 1.0 * num_dirac_0 } );

    //            cp.for_each_boundary_measure( FunctionEnum::Unit(), [&]( double boundary_measure, std::size_t num_dirac_1 ) {
    //                if ( num_dirac_1 == std::size_t( -1 ) ) {
    //                    ext_perimeters[ num_dirac_0 ] += boundary_measure;
    //                    return;
    //                }
    //                auto mi = std::min( num_dirac_0, num_dirac_1 );
    //                auto ma = std::max( num_dirac_0, num_dirac_1 );
    //                bms[ std::make_pair( mi, ma ) ].push_back( boundary_measure );
    //            } );

    //            mutex.unlock();
    //        } );
    //    }, bounds.englobing_convex_polyhedron(), positions.data(), weights.data(), positions.size() );

    //    double volume = 0;
    //    for( auto v : volumes )
    //        volume += v;
    //    CHECK_THAT( volume, WithinAbs<Grid::TF>( 1, 1e-6 ) );

    //    double ext_perimeter = 0;
    //    for( auto v : ext_perimeters )
    //        ext_perimeter += v;
    //    CHECK_THAT( ext_perimeter, WithinAbs<Grid::TF>( 4, 1e-6 ) );

    //    for( auto p : bms ) {
    //        REQUIRE( p.second.size() == 2 );
    //        CHECK_THAT( p.second[ 0 ], WithinAbs<Grid::TF>( p.second[ 1 ], 1e-6 ) );
    //    }

    //    grid.display( vo );
    //    // vo.save( "pd.vtk" );
}
