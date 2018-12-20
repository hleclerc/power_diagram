#include "../src/PowerDiagram/PowerDiagram.h"
#include <matplotlibcpp.h>
#include <chrono>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -g3
//// nsmake cpp_flag -O6
//// nsmake lib_flag -g3
//// nsmake lib_flag -O6

int main( int argc, char **argv ) {
    using Pc = PcStd<2,double>;
    using std::pow;

    if ( argc >= 2 ) {
        double ext_off = 1.0;
        double i = std::sqrt( std::atof( argv[ 1 ] ) );

        SpRegularGrid sp;
        //        sp.mul_cell_size = 0.8;

        PowerDiagram<Pc> pd( sp );
        pd.add_box_shape( { - 1 - ext_off, - 1 - ext_off }, { i + ext_off, i + ext_off } );
        for( double y = 0; y < i; ++y )
            for( double x = 0; x < i; ++x )
                if ( pow( x - 0.5 * i, 2 ) + pow( y - 0.5 * i, 2 ) <= pow( 0.5 * i, 2 ) )
                    pd.add_dirac( { x, y }, 0.2 );
        pd.set_ball_cut( true );

        std::vector<double> measures( pd.nb_diracs() );
        std::vector<std::vector<std::pair<std::size_t,double>>> derivatives( pd.nb_diracs() );
        auto b = std::chrono::high_resolution_clock::now();

        // pd.get_measures( measures.data() );
        pd.get_der_measures( measures.data(), derivatives.data() );

        auto t = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::high_resolution_clock::now() - b ).count() / 1e6;
        P( pd.nb_diracs(), t );

        if ( pd.nb_diracs() < 1e4 ) {
            VtkOutput<1> vo( { "num" } );
            pd.display( vo );
            // pd.display_bounds( vo );
            vo.save( "scal.vtk" );
        }

    } else {
        for( double ext_off : { 1.0 } ) {
            std::vector<double> vec_1;
            std::vector<double> vec_2;
            std::vector<double> vec_s;
            std::vector<double> vec_t;
            for( std::size_t s : { 1e3, 1e4, 1e5, 1e6, 1e7 } ) {
                PowerDiagram<Pc> pd;

                double i = std::sqrt( s );
                pd.add_box_shape( { - 1 - ext_off, - 1 - ext_off }, { i + ext_off, i + ext_off } );
                for( double y = 0; y < i; ++y )
                    for( double x = 0; x < i; ++x )
                        pd.add_dirac( { x, y }, 1 );

                std::vector<double> measures( pd.nb_diracs() );
                auto b = std::chrono::high_resolution_clock::now();
                pd.get_measures( measures.data() );
                auto t = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::high_resolution_clock::now() - b ).count() / 1e6;
                P( s, t );

                //            if ( s == 1e3 ) {
                //                VtkOutput<0> vo;
                //                pd.display( vo );
                //                pd.display_bounds( vo );
                //                vo.save( "scal.vtk" );
                //            }

                //            vec_s.push_back( pd.nb_diracs() );
                //            vec_t.push_back( t );
                //            vec_1.push_back( vec_t[ 0 ] * pow( pd.nb_diracs() / vec_s[ 0 ], 1 ) );
                //            vec_2.push_back( vec_t[ 0 ] * pow( pd.nb_diracs() / vec_s[ 0 ], 2 ) );
            }
            //        matplotlibcpp::loglog( vec_s, vec_t );
            //        matplotlibcpp::loglog( vec_s, vec_1, "-" );
            //        matplotlibcpp::loglog( vec_s, vec_2, "-" );
        }
        //    matplotlibcpp::show();
    }
}

