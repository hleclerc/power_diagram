#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
#include "../src/PowerDiagram/get_der_integrals_wrt_weights_ap.h"
#include "../src/PowerDiagram/get_der_integrals_wrt_weights.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "../src/PowerDiagram/get_integrals.h"
#include <boost/multiprecision/mpfr.hpp>
#include "catch_main.h"

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -O3
//// nsmake lib_name gmp
//// nsmake lib_name mpfr

TEST_CASE( "get_der_measures", "ctor" ) {
    struct Pc     { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = boost::multiprecision::mpfr_float_100; };
    using  Grid   = PowerDiagram::Visitor::ZGrid<Pc>;
    struct Dirac  { Point2<Pc::TF> pos; Pc::TF weight; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;

    std::vector<Dirac> diracs;
    for( std::size_t i = 0; i < 100; ++i ) {
        double x = 0.5 * rand() / RAND_MAX;
        double y = 1.0 * rand() / RAND_MAX;
        diracs.push_back( { { x + 0.5 * ( x > 0.25 ), y }, 1.0 + 0.01 * i } );
    }

    Bounds bounds;
    bounds.add_box( { 0, 0 }, { 1, 1 } );

    Grid grid( 25 );

    grid.init( diracs );
    std::vector<Pc::TF> m_values, v_values;
    std::vector<Pc::TI> m_offsets, m_columns;
    PowerDiagram::get_der_measures( m_offsets, m_columns, m_values, v_values, grid, bounds, diracs );

    const Pc::TF eps = 1e-40;
    for( std::size_t r = 0; r < diracs.size(); ++r ) {
        std::vector<Pc::TF> ex_der( diracs.size(), 0 );
        for( std::size_t o = m_offsets[ r + 0 ]; o < m_offsets[ r + 1 ]; ++o )
            ex_der[ m_columns[ o ] ] = m_values[ o ];

        std::vector<Dirac> new_diracs( diracs.size() );
        for( std::size_t i = 0; i < diracs.size(); ++i )
            new_diracs[ i ] = { diracs[ i ].pos, diracs[ i ].weight + eps * ( i == r ) };

        std::vector<Pc::TF> n_values;
        PowerDiagram::get_measures( n_values, grid, bounds, new_diracs );

        std::vector<Pc::TF> ap_der( diracs.size(), 0 );
        for( std::size_t c = 0; c < diracs.size(); ++c )
            ap_der[ c ] = ( n_values[ c ] - v_values[ c ] ) / eps;

        for( std::size_t c = 0; c < diracs.size(); ++c )
            CHECK_THAT( ex_der[ c ], WithinAbs( ap_der[ c ], Pc::TF( 1e-16 ) ) );

        Pc::TF volume = 0;
        for( std::size_t c = 0; c < diracs.size(); ++c )
            volume += n_values[ c ];
        CHECK_THAT( volume, WithinAbs<Pc::TF>( 1, 1e-16 ) );
    }
}
