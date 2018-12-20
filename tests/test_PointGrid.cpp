#include "../src/PowerDiagram/PointGrid.h"
#include <gtest/gtest.h>

TEST( SpiralGrid, traversal ) {
    PointGrid<std::size_t,double,3> pg;
    pg.init( [&]( const std::function<void( Point3<double> pos, std::size_t item )> &cb ) {
        cb( Point3<double>{ 0, 0, 0 }, 10 );
        cb( Point3<double>{ 1, 0, 0 }, 11 );
        cb( Point3<double>{ 0, 1, 0 }, 12 );
        cb( Point3<double>{ 1, 1, 0 }, 13 );
    }, 4 );

    int cpt = 10;
    pg.for_each_cell( [&]( auto cell_pos, const std::size_t *items, std::size_t nb_items ) {
        EXPECT_EQ( 1, nb_items );
        EXPECT_EQ( items[ 0 ], cpt++ );
    } );
}

