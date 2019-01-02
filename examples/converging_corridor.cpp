#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "../src/PowerDiagram/Visitors/RGrid.h"
#include "../src/PowerDiagram/EigenSolver.h"
#include "../src/PowerDiagram/AmgclSolver.h"
#include "../src/PowerDiagram/VtkOutput.h"

#include "../src/PowerDiagram/get_der_integrals_wrt_weights.h"
#include "../src/PowerDiagram/get_centroids.h"
#include "../src/PowerDiagram/system/Tick.h"
// #include <cxxopts.hpp>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -O6
//// nsmake lib_flag -O6

template<class Grid,class Bounds>
struct Solver {
    using CP = typename Grid::CP;
    using TF = typename Grid::TF;
    using TI = typename Grid::TI;

    Solver( Grid *grid, Bounds *bounds, TF target_radius ) : target_radius( target_radius ), bounds( *bounds ), grid( *grid ) {
    }

    template<class Diracs>
    TF iterate( Diracs &diracs ) {
        using std::sqrt;
        using std::pow;
        using std::max;
        using std::min;

        // M, V
        tick << "grid init";
        grid.init( diracs );
        tick >> "grid init";

        tick << "derivatives";
        PowerDiagram::get_der_measures( m_offsets, m_columns, m_values, v_values, grid, bounds, diracs );
        TODO;
        //        for( std::size_t r = 0; r < diracs.size(); ++r )
        //            v_values[ r ] -= M_PI * pow( target_radius, 2 );
        tick >> "derivatives";

        // solve
        tick << "solve";
        AmgclSolver solver;
        // EigenSolver solver;
        std::vector<TF> x( diracs.size() );
        solver.solve( x, m_offsets, m_columns, m_values, v_values );
        tick >> "solve";

        TF error = 0, max_x = 0;
        for( std::size_t i = 0; i < diracs.size(); ++i ) {
            max_x = max( max_x, x[ i ] );
            error += pow( x[ i ], 2 );
        }

        // update weights and get error
        TF coeff = min( 2.0 / max_x, 1.0 );
        for( std::size_t i = 0; i < diracs.size(); ++i )
            diracs[ i ].weight -= coeff * x[ i ];

        //
        return sqrt( error );
    }

    template<class Diracs>
    void centroids( Diracs &diracs ) {
        tick << "centroids";
        std::vector<typename Grid::Pt> centroids( diracs.size() );
        PowerDiagram::get_centroids( grid, bounds, diracs, [&]( auto p, auto m, std::size_t num ) {
            centroids[ num ] = p;
        } );
        for( size_t i = 0; i < diracs.size(); ++i )
            diracs[ i ].pos = centroids[ i ];
        tick >> "centroids";
    }

    template<class VO, class Diracs>
    void display( VO &vtk_output, Diracs &diracs ) {
        tick << "display";
        grid.init( diracs );

        grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0 ) {
            lc.display( vtk_output, { 1.0 * num_dirac_0 } );
        }, bounds.englobing_convex_polyhedron(), diracs );
        tick >> "display";
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

    TF              target_radius;
    std::vector<TI> m_offsets;
    std::vector<TI> m_columns;
    std::vector<TF> m_values;
    std::vector<TF> v_values;
    Bounds&         bounds;
    Grid&           grid;
};


int main( int argc, char **argv ) {
    using std::pow;

    struct Pc { enum { nb_bits_per_axis = 31, allow_ball_cut = 1, dim = 2 }; using TI = std::size_t; using TF = double; }; // boost::multiprecision::mpfr_float_100
    struct Dirac { Point2<Pc::TF> pos; Pc::TF weight; void write_to_stream( std::ostream &os ) const { os << '[' << pos << ']'; } };

    // diracs
    double h = 0.01;
    std::vector<Dirac> diracs;
    for( std::size_t y = 0; y < 2 / h; ++y )
        for( std::size_t x = 0; x < 2 / h; ++x )
            if ( pow( h * ( x + 0.5 ), 2 ) + pow( h * ( y + 0.5 ), 2 ) <= 4 )
                diracs.push_back( { { h * ( x + 0.5 ), h * ( y + 0.5 ) }, pow( 1.0 * h / 2, 2 ) } );
    P( diracs.size() );

    // grid
    using Grid = PowerDiagram::Visitor::ZGrid<Pc>;
    Grid grid( 16, 1e6 );

    // 0.5 * h => 2.05
    // 1.0 * h => 2.45
    // 2.0 * h => 3.34
    //    using Grid = PowerDiagram::Visitor::RGrid<Pc>;
    //    Grid grid( 0.75 * h );

    // Bounds
    using Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    Bounds bounds;

    bounds.add_box( { 0, 0 }, { 10, 10 } );

    // ZGrid, sans prise en compte taille => 1.819
    // ZGrid, avec prise en compte taille => mdpc = 20 => 1.8
    // en ne testant que la distance 1.44

    // 1.78

    // solve
    Solver<Grid,Bounds> solver( &grid, &bounds, 1.0 * h / 2 );
    for( std::size_t niter = 0; niter < 20; ++niter ) {
        if ( diracs.size() < 2000 ) {
            VtkOutput<1> vtk_output( { "num" } );
            solver.display( vtk_output, diracs );
            vtk_output.save( "vtk/cc_" + to_string( niter ) + ".vtk" );
        }

        for( auto &d : diracs )
            if ( norm_2( d.pos ) > h )
                d.pos -= h / 8 * normalized( d.pos );
        for( std::size_t nr = 0; nr < 8; ++nr ) {
            auto e = solver.iterate( diracs );
            //            P( e );
            if ( e < 1e-6 )
                break;
        }
        // P( solver.volume( diracs ) );
        solver.centroids( diracs );
    }
}
