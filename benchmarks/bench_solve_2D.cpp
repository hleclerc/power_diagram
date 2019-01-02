#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "../src/PowerDiagram/EigenSolver.h"
#include "../src/PowerDiagram/AmgclSolver.h"
#include "../src/PowerDiagram/VtkOutput.h"

#include "../src/PowerDiagram/get_der_integrals_wrt_weights_ap.h"
#include "../src/PowerDiagram/get_der_integrals_wrt_weights.h"
#include "../src/PowerDiagram/get_centroids.h"

#include "../src/PowerDiagram/system/Tick.h"
#include "set_up_diracs.h"

#include <cxxopts.hpp>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -O5
//// nsmake lib_flag -O5
//// nsmake lib_name gmp
//// nsmake lib_name mpfr

template<class Grid,class Bounds>
struct Solver {
    using CP = typename Grid::CP;
    using TF = typename Grid::TF;
    using TI = typename Grid::TI;

    Solver( Grid *grid, Bounds *bounds ) : bounds( *bounds ), grid( *grid ) {
    }


    template<class Diracs>
    TF iterate( Diracs &diracs, Diracs &old_diracs ) {
        using std::sqrt;
        using std::pow;

        TF relaxation = 1;

        // M, V
        while ( true ) {
            grid.init( diracs );
            // bool ok = PowerDiagram::get_der_measures_ap( m_offsets, m_columns, m_values, v_values, grid, bounds, diracs, TF( 1e-50 ) );
            int err = PowerDiagram::get_der_measures( m_offsets, m_columns, m_values, v_values, grid, bounds, diracs );
            P( relaxation );
            if ( err == 0 )
                break;
            if( old_diracs.size() == 0 )
                TODO;
            for( std::size_t r = 0; r < diracs.size(); ++r )
                diracs[ r ].weight = TF( 0.5 ) * diracs[ r ].weight + TF( 0.5 ) * old_diracs[ r ].weight;
            relaxation *= TF( 0.5 );
        }

        TF target_measure = 1.0 / diracs.size();
        for( std::size_t r = 0; r < diracs.size(); ++r )
            v_values[ r ] -= target_measure;
        if ( grid.ball_cut == 0 )
            m_values[ 0 ] *= 2;

        // solve
        EigenSolver solver;
        std::vector<TF> x( diracs.size(), TF( 0 ) );
        solver.solve( x, m_offsets, m_columns, m_values, v_values );

        // update weights and get error
        TF error = 0;
        old_diracs = diracs;
        for( std::size_t i = 0; i < diracs.size(); ++i ) {
            diracs[ i ].weight -= relaxation * x[ i ];
            error += pow( x[ i ], 2 );
        }

        //
        return sqrt( error );
    }

    template<class VO, class Diracs>
    void display( VO &vtk_output, Diracs &diracs ) {
        grid.init( diracs );

        grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
            lc.display( vtk_output, { diracs[ num_dirac_0 ].weight } );
        }, bounds.englobing_convex_polyhedron(), diracs );
    }

    template<class Diracs>
    TF volume( const Diracs &diracs ) {
        grid.init( diracs );

        std::vector<TF> volumes( diracs.size() );
        grid.for_each_laguerre_cell( [&]( CP &lc, std::size_t num ) {
            bounds.for_each_intersection( lc, [&]( CP &cp, auto space_func ) {
                volumes[ num ] += cp.integration( space_func );
            } );
        }, bounds.englobing_convex_polyhedron(), diracs );

        TF res = 0;
        for( auto v : volumes )
            res += v;
        return res;
    }

    std::vector<TI> m_offsets;
    std::vector<TI> m_columns;
    std::vector<TF> m_values;
    std::vector<TF> v_values;
    Bounds&         bounds;
    Grid&           grid;
};


int main( int argc, char **argv ) {
    struct Pc { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = double; }; // boost::multiprecision::mpfr_float_100
    struct Dirac { Point2<Pc::TF> pos; Pc::TF weight; void write_to_stream( std::ostream &os ) const { os << '[' << pos << ']'; } };

    // options
    cxxopts::Options options( argv[ 0 ], "bench solve");
    options.add_options()
        ( "m,max-dirac-per-cell"    , "...", cxxopts::value<int>()->default_value( "11" ) )
        ( "d,distribution"          , "distribution name (regular, random, ...)", cxxopts::value<std::string>()->default_value( "regular" ) )
        ( "t,nb-threads"            , "...", cxxopts::value<int>()->default_value( "0" ) )
        ( "v,vtk-output"            , "", cxxopts::value<std::string>() )
        ( "n,nb-diracs"             , "...", cxxopts::value<double>()->default_value( "100" ) )
        ( "h,help"                  , "get help" )
        ;
    auto args = options.parse( argc, argv );

    //
    thread_pool.init( args[ "nb-threads" ].as<int>() );

    // diracs
    std::vector<Dirac> diracs, old_diracs;
    set_up_diracs( diracs, args[ "distribution" ].as<std::string>(), args[ "nb-diracs" ].as<double>() );
    P( diracs.size() );

    // grid
    using Grid = PowerDiagram::Visitor::ZGrid<Pc>;
    Grid grid( args[ "max-dirac-per-cell" ].as<int>(), 1e6 );

    // Bounds
    using Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    Bounds bounds;

    bounds.add_box( { 0, 0 }, { 1, 1 } );

    // solve
    Solver<Grid,Bounds> solver( &grid, &bounds );
    for( std::size_t niter = 0; niter < 500; ++niter ) {
        double err = solver.iterate( diracs, old_diracs );
        P( solver.volume( diracs ), err );
        if ( err < 1e-4 )
            break;
    }
    if ( args.count( "vtk-output" ) ) {
        VtkOutput<1> vtk_output( { "weight" } );
        solver.display( vtk_output, diracs );
        vtk_output.save( args[ "vtk-output" ].as<std::string>() + ".vtk" );
    }

    // display weights, on a voronoi diagram
    if ( args.count( "vtk-output" ) ) {
        VtkOutput<1> vtk_output( { "weight" } );
        std::vector<Dirac> new_diracs = diracs;
        for( std::size_t i = 0; i < diracs.size(); ++i )
            new_diracs[ i ].weight = 1;

        grid.init( new_diracs );
        grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
            lc.display( vtk_output, { diracs[ num_dirac_0 ].weight } );
        }, bounds.englobing_convex_polyhedron(), new_diracs );

        vtk_output.save( args[ "vtk-output" ].as<std::string>() + "_orig_pts.vtk" );
    }
    // P( solver.volume( diracs ) );
}
