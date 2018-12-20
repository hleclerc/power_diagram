#include "../src/PowerDiagram/Visitors/ZGrid.h"
#include "../src/PowerDiagram/system/Tick.h"
#include "../src/PowerDiagram/VtkOutput.h"
#include "../src/PowerDiagram/PcStd.h"
#include <matplotlibcpp.h>
#include <gtest/gtest.h>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -O5
//// nsmake lib_flag -O5

TEST( VisitorZIndex, part_case ) {
    struct Dirac { Point2<double> pos; double weight; };
    struct Pc {
        static constexpr int max_diracs_per_cell = 2;
        static constexpr int nb_bits_per_axis    = 20;
        static constexpr int dim                 = 2;
        using                TI                  = std::size_t;
        using                TF                  = double;
    };

    LaguerreCellVisitors::ZIndex<Pc> lv( 2.0 );

    using LC = typename Tr::LC;

    for( std::size_t n : { 1e5, 1e6 /*10, 100, 1000, 10000, 100000*/ } ) {
        VtkOutput<1> vo( { "num" } );
        //        for( std::size_t i = 0; i < n; ++i )
        //            diracs[ i ] = { { 1.0 * rand() / RAND_MAX, 1.0 * rand() / RAND_MAX }, 1.0 };

        std::size_t l = std::sqrt( n );
        std::vector<Dirac> diracs( l * l );
        for( std::size_t i = 0, c = 0; i < l; ++i )
            for( std::size_t j = 0; j < l; ++j )
                diracs[ c++ ] = { { ( i + 0.5 ) / l, ( j + 0.5 ) / l }, 1.0 };

        Pc pc;
        Tr tr( pc.sp );
        auto time = Tick::get_time();
        tr.init_nodes( diracs );

        LC bounds( LC::Box{ { 0.0, 0.0 }, { 1.0, 1.0 } } );

        std::vector<double> volumes( diracs.size(), 17 );
        tr.for_each_cell( [&]( LC &lc, std::size_t num ) {
            if ( n < 500 ) {
                lc.display( vo, { 1.0 * num } );
                vo.add_point( diracs[ num ].pos );
            }
            volumes[ num ] = lc.measure();
        }, bounds, diracs );

        double t = Tick::elapsed_since( time );

        double volume = 0;
        for( double v : volumes )
            volume += v;
        P( volume, t );
        if ( n < 500 )
            vo.save( "pd.vtk" );
    }

    //    std::vector<Dirac> diracs;
    //    for( std::size_t i = 0; i < n; ++i )
    //        diracs.push_back( { { 1.0 * rand() / RAND_MAX, 1.0 * rand() / RAND_MAX }, 1.0 } );

    //    Pc pc;
    //    Tr tr( pc.sp );
    //    tr.init( diracs );

    //    LC bounds( LC::Box{ { 0.0, 0.0 }, { 1.0, 1.0 } } );
    //    bounds.display( vo, { 1.0 } );
    //    tr.for_each_cell( [&]( LC &lc, std::size_t num ) {
    //        lc.display( vo, { 1.0 * num } );
    //        vo.add_point( diracs[ num ].pos );
    //    }, bounds, diracs );

    //    vo.save( "pd.vtk" );

    //    std::vector<double> nx, ny;
    //    matplotlibcpp::plot( x, y, "*-" );
    //    matplotlibcpp::show();
}



