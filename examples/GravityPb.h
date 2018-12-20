#pragma once

#include "../src/PowerDiagram/PowerDiagram.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>

/**
*/
template<int dim>
struct GravityPb {
    using TF = double;
    using PO = typename std::conditional<dim==3,Point3<TF>,Point2<TF>>::type;
    using TI = std::size_t;

    GravityPb( TF r = 0.5 ) : r( r ) {
    }

    void change_pos() {
        old_positions.resize( pd.nb_diracs() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            old_positions[ i ] = pd.dirac( i ).pos;
            pd.dirac( i ).pos[ dim - 1 ] -= 0.1;
        }
    }

    void save( int n ) {
        if ( n == 0 ) {
            VtkOutput<1,TF> vo( { "num" } );
            pd.display_bounds( vo );
            vo.save( "../res/bounds.vtk" );
        }

        VtkOutput<1,TF> vo( { "num" } );
        pd.display( vo );
        vo.save( "../res/pd_" + to_string( 100 + n ) + ".vtk" );

        std::string fna = "../res/pos_" + to_string( 100 + n ) + ".dat";
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
            measures.resize( pd.nb_diracs() );
            derivatives.resize( pd.nb_diracs() );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
                derivatives[ i ].resize( 0 );
                measures[ i ] = 0;
            }
            pd.get_der_measures( measures.data(), derivatives.data() );

            // rhs
            Eigen::VectorXd V( pd.nb_diracs() );
            for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
                V[ i ] = TF( 1 + dim ) / 3 * M_PI * std::pow( r, dim ) - measures[ i ];

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
                pd.get_measures( measures.data() );
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
        measures.resize( pd.nb_diracs() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            for( std::size_t d = 0; d < dim; ++d )
                centroid_contribs[ i ][ d ] = 0;
            measures[ i ] = 0;
        }

        pd.get_centroid_contribs( centroid_contribs.data(), measures.data() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            pd.dirac( i ).pos = centroid_contribs[ i ] / measures[ i ];
    }

    std::vector<PO>                            centroids_mc;
    std::vector<TF>                            volumes_mc;

    std::vector<PO>                            centroid_contribs;
    std::vector<PO>                            old_positions;
    std::vector<std::vector<std::pair<TI,TF>>> derivatives;
    std::vector<TF>                            measures;
    PowerDiagram<dim,TF>                       pd;
    TF                                         r;
};
