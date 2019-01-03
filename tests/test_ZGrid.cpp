#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "catch_main.h"
//// nsmake cpp_flag -march=native
using std::abs;

TEST_CASE( "ZGrid measures" ) {
    struct Pc     { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = double; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    using  Grid   = PowerDiagram::Visitor::ZGrid<Pc>;

    std::vector<Grid::Pt> positions;
    std::vector<Grid::TF> weights;
    for( std::size_t i = 0; i < 10; ++i ) {
        positions.push_back( { 1.0 * rand() / RAND_MAX, 1.0 * rand() / RAND_MAX } );
        weights.push_back( 1.0 );
    }

    Bounds bounds;
    bounds.add_box( { 0, 0 }, { 1, 1 }, 1.0, -1 );

    Grid grid( 5 );
    grid.update( positions.data(), weights.data(), positions.size() );

    VtkOutput<1> vo( { "num" } );
    std::mutex mutex;
    std::vector<double> volumes( positions.size() );
    std::vector<double> ext_perimeters( positions.size(), 0.0 );
    std::map<std::pair<std::size_t,std::size_t>,std::vector<double>> bms;
    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
        bounds.for_each_intersection( lc, [&]( auto &cp, auto space_func ) {
            mutex.lock();
            volumes[ num_dirac_0 ] += cp.integration( space_func );
            cp.display( vo, { 1.0 * num_dirac_0 } );

            cp.for_each_boundary_measure( FunctionEnum::Unit(), [&]( double boundary_measure, std::size_t num_dirac_1 ) {
                if ( num_dirac_1 == std::size_t( -1 ) ) {
                    ext_perimeters[ num_dirac_0 ] += boundary_measure;
                    return;
                }
                auto mi = std::min( num_dirac_0, num_dirac_1 );
                auto ma = std::max( num_dirac_0, num_dirac_1 );
                bms[ std::make_pair( mi, ma ) ].push_back( boundary_measure );
            } );
            
            mutex.unlock();
        } );
    }, bounds.englobing_convex_polyhedron(), positions.data(), weights.data(), positions.size() );

    double volume = 0;
    for( auto v : volumes )
        volume += v;
    CHECK_THAT( volume, WithinAbs<Grid::TF>( 1, 1e-6 ) );

    double ext_perimeter = 0;
    for( auto v : ext_perimeters )
        ext_perimeter += v;
    CHECK_THAT( ext_perimeter, WithinAbs<Grid::TF>( 4, 1e-6 ) );

    for( auto p : bms ) {
        REQUIRE( p.second.size() == 2 );
        CHECK_THAT( p.second[ 0 ], WithinAbs<Grid::TF>( p.second[ 1 ], 1e-6 ) );
    }

    grid.display( vo );
    // vo.save( "pd.vtk" );
}

TEST_CASE( "several ZGrids" ) {
    struct Pc     { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = double; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    using  Grid   = PowerDiagram::Visitor::ZGrid<Pc>;

    std::vector<Grid::Pt> positions;
    std::vector<Grid::TF> weights;
    for( double i = 0; i < 50; i += 2 ) {
        positions.push_back( { 0, i } );
        weights.push_back( 1.0 );
    }
    for( double i = 0; i < 50; i += 2 ) {
        positions.push_back( { 49, i } );
        weights.push_back( 2.0 );
    }

    Bounds bounds;
    bounds.add_box( { -1, -1 }, { 50, 50 }, 1.0, -1 );

    Grid grid( 2, 0.75 );
    grid.update( positions.data(), weights.data(), positions.size() );

    std::atomic<int> nb_cp( 0 );
    grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
        nb_cp++;
    }, bounds.englobing_convex_polyhedron(), positions.data(), weights.data(), positions.size() );

    CHECK( nb_cp == weights.size() );

    // VtkOutput<1> vo_grid( { "num" } );
    // grid.display( vo_grid );
    // vo_grid.save( "vtk/grid.vtk" );

    // VtkOutput<1> vo_pd( { "num" } );
    // grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
    //     bounds.for_each_intersection( lc, [&]( auto &cp, auto space_func ) {
    //         cp.display( vo_pd, { 1.0 * num_dirac_0 } );
    //     } );
    // }, bounds.englobing_convex_polyhedron(), positions.data(), weights.data(), positions.size() );
    // vo_pd.save( "vtk/pd.vtk" );
}
