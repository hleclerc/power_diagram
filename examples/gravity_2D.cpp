//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -O3
#include "GravityPb.h"

int main() {
    using TF = double;
    using PO = Point2<TF>;

    GravityPb<2> pb( 0.5 );

    TF l = 5;
    TF h = 3;
    TF m = 1 * l / 4, n = 3 * l / 4;

    // boundaries
    typename PowerDiagram<2,TF>::LC lc( { 0, 0 }, 100 );
    lc.plane_cut( { 0, 0 }, { -1, 0 } );
    lc.plane_cut( { 0, 0 }, { 0, -1 } );
    lc.plane_cut( { l, h }, { +1, 0 } );
    lc.plane_cut( { l, h }, { 0, +1 } );
    lc.plane_cut( { m, 0 }, normalized( PO{ -0.5, -1 } ) );
    lc.plane_cut( { n, 0 }, normalized( PO{ +0.5, -1 } ) );
    pb.pd.add_convex_shape( lc );

    pb.pd.add_box_shape( { m, -2 }, { n,  0 } );
    pb.pd.add_box_shape( { 0, -4 }, { l, -2 } );

    // diracs
    for( TF y = 0; y < h; ++y )
        for( TF x = 0; x < l; ++x )
            pb.pd.add_dirac( { x + 0.5, y + 0.5 }, pb.r * pb.r );

    // simulation
    pb.save( 0 );
    for( std::size_t time_step = 1; time_step < 150; ++time_step ) {
        P( time_step );
        pb.change_pos();
        pb.get_op_radii();
        pb.move_to_centroids();
        pb.save( time_step );
    }
}

