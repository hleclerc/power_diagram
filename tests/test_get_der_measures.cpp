#include "../src/PowerDiagram/Bounds/ConvexPolyhedronAssembly.h"
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
void test_for_radial_func( const FU &fu, double expected_precision ) {
    struct Pc     { enum { nb_bits_per_axis = 31, allow_ball_cut = 0, dim = 2 }; using TI = std::size_t; using TF = boost::multiprecision::mpfr_float_100; };
    using  Bounds = PowerDiagram::Bounds::ConvexPolyhedronAssembly<Pc>;
    using  Grid   = PowerDiagram::Visitor::ZGrid<Pc>;
    using  PT     = typename Grid::Pt;
    using  TF     = typename Pc::TF;
    using  TI     = typename Pc::TI;

    using std::max;
    using std::abs;

    TF max_err = 0, max_val = 0;
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
        PowerDiagram::get_der_integrals_wrt_weights( m_offsets, m_columns, m_values, v_values, grid, bounds, positions.data(), weights.data(), weights.size(), fu, false );

        const TF eps = 1e-40;
        for( std::size_t r = 0; r < weights.size(); ++r ) {
            // ex_der
            std::vector<TF> ex_der( weights.size(), 0 );
            for( std::size_t o = m_offsets[ r + 0 ]; o < m_offsets[ r + 1 ]; ++o )
                ex_der[ m_columns[ o ] ] = m_values[ o ];

            // ap_der
            std::vector<TF> new_weights( weights.size() );
            for( std::size_t i = 0; i < weights.size(); ++i )
                new_weights[ i ] = weights[ i ] + eps * ( i == r );

            std::vector<TF> n_values( weights.size() );
            PowerDiagram::get_integrals( &n_values[ 0 ], grid, bounds, positions.data(), new_weights.data(), weights.size(), fu );

            std::vector<TF> ap_der( weights.size(), 0 );
            for( TI c = 0; c < weights.size(); ++c )
                ap_der[ c ] = ( n_values[ c ] - v_values[ c ] ) / eps;

            // tests
            for( TI c = 0; c < weights.size(); ++c ) {
                max_err = max( max_err, abs( ex_der[ c ] - ap_der[ c ] ) );
                max_val = max( max_val, abs( ap_der[ c ] ) );
            }
            //            for( std::size_t c = 0; c < weights.size(); ++c )
            //                CHECK_THAT( ex_der[ c ], WithinAbs( ap_der[ c ], TF( expected_precision ) ) );

            //            TF volume = 0;
            //            for( std::size_t c = 0; c < weights.size(); ++c )
            //                volume += n_values[ c ];
            //            CHECK_THAT( volume, WithinAbs<TF>( 1, 1e-16 ) );
        }
    }
    // 1e-10 => 2.23e-9
    // 1e-12 => 4.01e-10
    // 1e-14 => 5.03e-11
    P( fu.name(), max_err / max_val );
}

TEST_CASE( "get_der_measures", "ctor" ) {
    using TF = boost::multiprecision::mpfr_float_100;
    //    test_for_radial_func<TF>( FunctionEnum::Unit(), 1e-16 );
    //    test_for_radial_func( FunctionEnum::ExpWmR2db<TF>{ 1.0 }, 1e-8 );
    test_for_radial_func( FunctionEnum::ExpWmR2db<TF>{ 0.5 }, 1e-8 );
}
