#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/OptimalTransportSolver.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "set_up_diracs.h"
#include <cxxopts.hpp>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -O5
//// nsmake lib_flag -O5
// // nsmake lib_name gmp
// // nsmake lib_name mpfr


int main( int argc, char **argv ) {
    struct Pc { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = double; }; // boost::multiprecision::mpfr_float_100
    using  Pt = Point2<Pc::TF>;
    using  TF = Pc::TF;

    // options
    cxxopts::Options options( argv[ 0 ], "bench solve");
    options.add_options()
        ( "m,max-dirac-per-cell"    , "...", cxxopts::value<int>()->default_value( "11" ) )
        ( "r,max-delta-weight"      , "...", cxxopts::value<double>()->default_value( "1e40" ) )
        ( "eq-w-repartition"        , "..." )
        ( "d,distribution"          , "distribution name (regular, random, ...)", cxxopts::value<std::string>()->default_value( "regular" ) )
        ( "t,nb-threads"            , "...", cxxopts::value<int>()->default_value( "0" ) )
        ( "v,vtk-output"            , "", cxxopts::value<std::string>() )
        ( "n,nb-diracs"             , "...", cxxopts::value<double>()->default_value( "100" ) )
        ( "max-iter"                , "...", cxxopts::value<int>()->default_value( "100" ) )
        ( "o,output"                , "", cxxopts::value<std::string>() )
        ( "h,help"                  , "get help" )
        ;
    auto args = options.parse( argc, argv );

    //
    thread_pool.init( args[ "nb-threads" ].as<int>() );

    // diracs
    std::vector<Pt> positions;
    std::vector<TF> weights;
    set_up_diracs( positions, weights, args[ "distribution" ].as<std::string>(), args[ "nb-diracs" ].as<double>() );

    // grid
    using Grid = PowerDiagram::Visitor::ZGrid<Pc>;
    Grid grid( args[ "max-dirac-per-cell" ].as<int>(), args[ "max-delta-weight" ].as<double>() );
    grid.eq_rep_weight_split = args.count( "eq-w-repartition" );

    // Bounds
    using Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    Bounds bounds;
    bounds.add_box( { 0, 0 }, { 1, 1 } );

    // solve
    PowerDiagram::OptimalTransportSolver<Grid,Bounds> solver( &grid, &bounds );
    solver.max_nb_iter = args[ "max-iter" ].as<int>();
    solver.solve( positions.data(), weights.data(), weights.size() );
    //        P( solver.volume( diracs ), err );

    if ( args.count( "vtk-output" ) ) {
        VtkOutput<2> vtk_output( { "weight", "num" } );
        solver.display( vtk_output, positions.data(), weights.data(), weights.size() );
        vtk_output.save( args[ "vtk-output" ].as<std::string>() + ".vtk" );
    }

    //
    if ( args.count( "output" ) ) {
        std::string o = args[ "output" ].as<std::string>();
        std::ofstream f( o.c_str() );
        for( std::size_t i = 0; i < weights.size(); ++i ) {
            f << std::setprecision( 16 ) << positions[ i ].x << " ";
            f << std::setprecision( 16 ) << positions[ i ].y << " ";
            f << std::setprecision( 16 ) << weights  [ i ]   << "\n";
        }
    }

    // display weights, on a voronoi diagram
    if ( args.count( "vtk-output" ) ) {
        VtkOutput<1> vtk_output( { "weight" } );
        //solver.display_orig_pts( vtk_output, positions.data(), weights.data(), weights.size() );

        grid.display( vtk_output, 0.03 );
        vtk_output.add_lines( grid.proute_cells );

        vtk_output.save( args[ "vtk-output" ].as<std::string>() + "_orig_pts.vtk" );
    }
}
