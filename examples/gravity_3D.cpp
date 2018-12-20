#include "../src/PowerDiagram.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <matplotlibcpp.h>
#include <gtest/gtest.h>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -ffast-math
//// nsmake cpp_flag -O3

struct Pb {
    using TF = double;
    using PO = Point3<TF>;
    using TI = std::size_t;

    Pb( TF r = 0.5, TF l = 5, TF h = 3 ) : r( r ) {
        TF m = 1 * l / 4, n = 3 * l / 4;
        PowerDiagram<TF>::LC lc( { 0, 0, 0 }, 100 );
        lc.plane_cut( { 0, 0, 0 }, { -1, 0, 0 } );
        lc.plane_cut( { 0, 0, 0 }, { 0, -1, 0 } );
        lc.plane_cut( { 0, 0, 0 }, { 0, 0, -1 } );
        lc.plane_cut( { l, l, h }, { +1, 0, 0 } );
        lc.plane_cut( { l, l, h }, { 0, +1, 0 } );
        lc.plane_cut( { l, l, h }, { 0, 0, +1 } );

        lc.plane_cut( { m, m, 0 }, normalized( PO{ -0.5, 0, -1 } ) );
        lc.plane_cut( { m, m, 0 }, normalized( PO{ 0, -0.5, -1 } ) );
        lc.plane_cut( { n, n, 0 }, normalized( PO{ +0.5, 0, -1 } ) );
        lc.plane_cut( { n, n, 0 }, normalized( PO{ 0, +0.5, -1 } ) );

        pd.add_convex_bounds( lc );

        pd.add_box_shape( { m, m, -2 }, { n, n,  0 } );
        pd.add_box_shape( { 0, 0, -4 }, { l, l, -2 } );

        for( TF z = 0, c = 0; z < h; ++z )
            for( TF y = 0; y < l; ++y )
                for( TF x = 0; x < l; ++x, ++c )
                    pd.add_dirac( { x + 0.5, y + 0.5, z + 0.5 }, r * r );
    }

    void change_pos() {
        old_positions.resize( pd.nb_diracs() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            old_positions[ i ] = pd.dirac( i ).pos;
            pd.dirac( i ).pos.z -= 0.1;
        }
    }

    void save( int n ) {
        if ( n == 0 ) {
            VtkOutput<1,TF> vo( { "num" } );
            pd.display_bounds( vo );
            vo.save( "res/bounds.vtk" );
        }

        VtkOutput<1,TF> vo( { "num" } );
        pd.display( vo );
        vo.save( "res/pd_" + to_string( 100 + n ) + ".vtk" );

        std::string fna = "res/pos_" + to_string( 100 + n ) + ".dat";
        std::ofstream f( fna.c_str() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            f << pd.dirac( i ).pos << " " << pd.dirac( i ).weight << "\n";
    }

    void get_op_radii() {
        using EV = Eigen::VectorXd;

        EV best_W( pd.nb_diracs() );
        TF relaxation = 0.75;
        TF best_norm_V = std::numeric_limits<TF>::max();
        for( std::size_t num_iter = 0; ; ++num_iter ) {
            // if stuck, restart with the best weights, with a smaller relaxation
            if ( num_iter && num_iter % 20 == 0 ) {
                for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                    pd.dirac( i ).weight = best_W[ i ];
                relaxation /= 2;
            }

            //            static int cpt = 0;
            //            VtkOutput<1,TF> vo( { "num" } );
            //            pd.display( vo );
            //            vo.save( "res/iter_" + to_string( cpt++ ) + ".vtk" );

            // volumes and derivatives
            volumes.resize( pd.nb_diracs() );
            derivatives.resize( pd.nb_diracs() );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
                derivatives[ i ].resize( 0 );
                volumes[ i ] = 0;
            }
            pd.get_der_volumes( volumes.data(), derivatives.data() );

            // rhs
            Eigen::VectorXd V( pd.nb_diracs() );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                V[ i ] = 4.0 / 3.0 * M_PI * std::pow( r, 3 ) - volumes[ i ];

            //
            TF norm_V = std::sqrt( V.squaredNorm() / pd.nb_diracs() );
            if ( best_norm_V > norm_V ) {
                 best_norm_V = norm_V;
                 for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                     best_W[ i ] = pd.dirac( i ).weight;
            }

            // matrix
            Eigen::SparseMatrix<double> M( pd.nb_diracs(), pd.nb_diracs() );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                for( std::pair<TI,TF> p : derivatives[ i ] )
                    M.coeffRef( i, p.first ) += p.second;

            // cholesky
            Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> C;
            C.compute( M );

            // solve and update the weights
            Eigen::VectorXd D = C.solve( V );
            TF nD = D.lpNorm<Eigen::Infinity>();
            TF cr = relaxation; // * ( nD > r * 5e-1 ? r * 5e-1 / nD : 1 );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                pd.dirac( i ).weight += cr * D[ i ];

            //            TF mi = 1e40, ma = 0;
            //            for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            //                mi = std::min( mi, pd.dirac( i ).weight );
            //                ma = std::max( ma, pd.dirac( i ).weight );
            //            }
            //            P( mi, ma, norm_V );
            //            for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            //                std::cout << pd.dirac( i ).weight << " ";
            //            std::cout << "\n";

            P( norm_V );

            //            if ( norm_V > 10 ) {
            //                P( volumes );
            //                for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            //                    P( pd.dirac( i ).weight );
            //            }

            // stop criterion
            if ( num_iter == 60 ) {
                //                for( size_t i = 0; i < pd.nb_diracs(); ++i )
                //                    P( pd.dirac( i ).pos, pd.dirac( i ).weight );

                // 0.857755
                //                volumes_mc.resize( pd.nb_diracs() );
                //                pd.get_volumes_mc( volumes_mc.data() );
                pd.get_volumes( volumes.data() );
                P( volumes );
                //                P( volumes_mc );
            }
            if ( norm_V < 1e-4 || num_iter == 60 ) {
                //                Eigen::MatrixXd m( M );
                //                Eigen::EigenSolver<Eigen::MatrixXd> es{ m, false };
                //                P( max( volumes ) - min( volumes ), es.eigenvalues().real().maxCoeff(), es.eigenvalues().real().minCoeff() );
                break;
            }
        }
    }

    void move_to_centroids() {
        centroid_contribs.resize( pd.nb_diracs() );
        volumes.resize( pd.nb_diracs() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            centroid_contribs[ i ] = { 0, 0, 0 };
            volumes[ i ] = 0;
        }

        pd.get_centroid_contribs( centroid_contribs.data(), volumes.data() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            pd.dirac( i ).pos = centroid_contribs[ i ] / volumes[ i ];
    }

    std::vector<PO>                            centroids_mc;
    std::vector<TF>                            volumes_mc;

    std::vector<PO>                            centroid_contribs;
    std::vector<PO>                            old_positions;
    std::vector<std::vector<std::pair<TI,TF>>> derivatives;
    std::vector<TF>                            volumes;
    PowerDiagram<TF>                           pd;
    TF                                         r;
};
//    std::ifstream f("114.txt");
//    for( std::size_t i = 0; i < pb.pd.nb_diracs(); ++i )
//        f >> pb.pd.dirac( i ).pos.x >> pb.pd.dirac( i ).pos.y >> pb.pd.dirac( i ).pos.z >> pb.pd.dirac( i ).weight;
//    //    pb.save( 0 );
//    //    pb.get_op_radii();
//    //    pb.save( 1 );
//    pb.change_pos();
//    //    pb.get_op_radii();
//    //    pb.save( 2 );

//    std::vector<Pb::TF> W0{ 0.29994 , 0.300545, 0.325219, 0.283964, 0.251125 };
//    std::vector<Pb::TF> W1{ 0.300262, 0.300611, 0.331469, 0.285212, 0.251125 };
//    std::vector<double> plt[ 5 ];
//    for( double t = 0; t <= 1; t += 1.0/8 ) {
//        for( size_t i = 0; i < W0.size(); ++i )
//            pb.pd.dirac( i ).weight  = ( 1 - t ) * W0[ i ] + t * W1[ i ];
//        pb.volumes.resize( pb.pd.nb_diracs() );
//        pb.pd.get_volumes( pb.volumes.data() );
//        for( size_t i = 0; i < W0.size(); ++i )
//            plt[ i ].push_back( pb.volumes[ i ] );
//        pb.pd.get_volumes_mc( pb.volumes.data() );
//        P( pb.volumes[ 2 ] );
//        pb.save( 8 * t );
//    }

//    // 130
//    P( plt[ 2 ] );

//    //    for( size_t i = 0; i < W0.size(); ++i )
//    matplotlibcpp::plot( plt[ 2 ] );
//    matplotlibcpp::show();
//        if ( time_step == 114 ) {
//            std::ofstream f("114.txt");
//            for( std::size_t i = 0; i < pb.pd.nb_diracs(); ++i )
//                f << pb.pd.dirac( i ).pos << " " << pb.pd.dirac( i ).weight;

//        }

TEST( Gravity, part_case ) {

    Pb pb;
    // pb.save( 0 );
    for( std::size_t time_step = 0; time_step < 180; ++time_step ) {
        P( time_step );
        pb.change_pos();
        pb.get_op_radii();
        pb.save( time_step );
        pb.move_to_centroids();
    }
}

