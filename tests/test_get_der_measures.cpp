#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
// #include "../src/PowerDiagram/get_der_integrals_wrt_weights_ap.h"
#include "../src/PowerDiagram/get_der_integrals_wrt_weights.h"
#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "../src/PowerDiagram/get_integrals.h"
#include <boost/multiprecision/mpfr.hpp>
#include "catch_main.h"

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -O3
//// nsmake lib_name gmp
//// nsmake lib_name mpfr

template<class FU>
void test_for_radial_func( const FU &fu ) {
    struct Pc     { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = boost::multiprecision::mpfr_float_100; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    using  Grid   = PowerDiagram::Visitor::ZGrid<Pc>;
    using  PT     = typename Grid::Pt;
    using  TF     = Pc::TF;
    using  TI     = Pc::TI;

    for( std::size_t num_real = 0; num_real < 10; ++num_real ) {
        std::vector<PT> positions;
        std::vector<TF> weights;
        for( std::size_t i = 0; i < 100; ++i ) {
            double x = 0.5 * rand() / RAND_MAX;
            double y = 1.0 * rand() / RAND_MAX;
            positions.push_back( { x + 0.5 * ( x > 0.25 ), y } );
            weights.push_back( 1.0 + 0.01 * i );
        }

        Bounds bounds;
        bounds.add_box( { 0, 0 }, { 1, 1 } );

        Grid grid( 25 );

        grid.update( positions.data(), weights.data(), weights.size() );
        std::vector<TF> m_values, v_values;
        std::vector<TI> m_offsets, m_columns;
        int err = PowerDiagram::get_der_integrals_wrt_weights( m_offsets, m_columns, m_values, v_values, grid, bounds, positions.data(), weights.data(), weights.size(), FunctionEnum::Unit(), false );

        const TF eps = 1e-40;
        for( std::size_t r = 0; r < weights.size(); ++r ) {
            std::vector<TF> ex_der( weights.size(), 0 );
            for( std::size_t o = m_offsets[ r + 0 ]; o < m_offsets[ r + 1 ]; ++o )
                ex_der[ m_columns[ o ] ] = m_values[ o ];

            std::vector<TF> new_weights( weights.size() );
            for( std::size_t i = 0; i < weights.size(); ++i )
                new_weights[ i ] = weights[ i ] + eps * ( i == r );

            std::vector<TF> n_values( weights.size() );
            PowerDiagram::get_integrals( &n_values[ 0 ], grid, bounds, positions.data(), new_weights.data(), weights.size(), FunctionEnum::Unit() );

            std::vector<Pc::TF> ap_der( weights.size(), 0 );
            for( std::size_t c = 0; c < weights.size(); ++c )
                ap_der[ c ] = ( n_values[ c ] - v_values[ c ] ) / eps;

            for( std::size_t c = 0; c < weights.size(); ++c )
                CHECK_THAT( ex_der[ c ], WithinAbs( ap_der[ c ], Pc::TF( 1e-16 ) ) );

            TF volume = 0;
            for( std::size_t c = 0; c < weights.size(); ++c )
                volume += n_values[ c ];
            CHECK_THAT( volume, WithinAbs<TF>( 1, 1e-16 ) );
        }
    }
}

TEST_CASE( "get_der_measures", "ctor" ) {
    test_for_radial_func( FunctionEnum::Unit() );
}
