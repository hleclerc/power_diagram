#pragma once

#include "../src/PowerDiagram/system/Assert.h"
#include "../src/PowerDiagram/Point2.h"
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

template<class Pt,class TF>
void set_up_diracs( std::vector<Pt> &positions, std::vector<TF> &weights, std::string distribution, std::size_t nb_diracs ) {
    using TI = std::size_t;
    using std::max;
    using std::min;

    if ( distribution == "regular" ) {
        TI l = std::sqrt( nb_diracs );
        positions.resize( l * l );
        weights.resize( l * l );
        for( TI i = 0, c = 0; i < l; ++i ) {
            for( TI j = 0; j < l; ++j, ++c ) {
                positions[ c ] = { ( j + 0.45 + 0.1 * rand() / RAND_MAX ) / l, ( i + 0.45 + 0.1 * rand() / RAND_MAX ) / l };
                weights[ c ] = 1.0;
            }
        }
        return;
    }

    if ( distribution == "random" ) {
        positions.resize( nb_diracs );
        weights.resize( nb_diracs );
        for( TI i = 0; i < nb_diracs; ++i ) {
            positions[ i ] = { 1.0 * rand() / RAND_MAX, 1.0 * rand() / RAND_MAX };
            weights[ i ] = 1.0;
        }
        return;
    }

    if ( distribution == "split" ) {
        positions.resize( nb_diracs );
        weights.resize( nb_diracs );
        for( std::size_t i = 0; i < nb_diracs; ++i ) {
            double x = 0.5 * rand() / RAND_MAX;
            double y = 1.0 * rand() / RAND_MAX;
            positions[ i ] = { x + 0.5 * ( x > 0.25 ), y };
            weights[ i ] = 1.0;
        }
        return;
    }

    if ( distribution == "lines" ) {
        positions.resize( 0 );
        weights.resize( 0 );
        for( std::size_t i = 0; i < nb_diracs / 2; ++i ) {
            double y = 1.0 * rand() / RAND_MAX;
            double x = 0.1 + 0.05 * rand() / RAND_MAX;
            double w = 1; // 0.5 + 0.05 * rand() / RAND_MAX;
            positions.push_back( { x, y } );
            weights.push_back( w );
        }
        for( std::size_t i = 0; i < nb_diracs / 2; ++i ) {
            double y = 1.0 * rand() / RAND_MAX;
            double x = 0.9 - 0.5 * y + 0.05 * rand() / RAND_MAX;
            double w = 1; // 0.5 + 0.05 * rand() / RAND_MAX;
            positions.push_back( { x, y } );
            weights.push_back( w );
        }

        return;
    }

    if ( distribution == "concentration" ) {
        positions.resize( 0 );
        weights.resize( 0 );
        //
        for( std::size_t i = 0; i < nb_diracs / 20; ++i ) {
            double y = 1.0 * rand() / RAND_MAX;
            double x = 1.0 * rand() / RAND_MAX;
            double w = 1;
            positions.push_back( { x, y } );
            weights.push_back( w );
        }
        //
        for( std::size_t i = 0; i < 19 * nb_diracs / 20; ++i ) {
            double y = 0.7 + 0.1 * rand() / RAND_MAX;
            double x = 0.7 + 0.1 * rand() / RAND_MAX;
            double w = 1; // 0.5 + 0.05 * rand() / RAND_MAX;
            positions.push_back( { x, y } );
            weights.push_back( w );
        }

        return;
    }

    if ( distribution.size() > 4 && distribution.substr( distribution.size() - 4 ) == ".xyz" ) {
        std::ifstream f( distribution.c_str() );
        positions.resize( 0 );
        weights.resize( 0 );
        double x, y, z;
        while ( f >> x >> y >> z ) {
            if ( z > 0.1 && z < 0.2 ) {
                positions.push_back( { x, y } );
                weights.push_back( 1 );
            }
        }

        return;
    }

    if ( distribution.size() > 6 && distribution.substr( 0, 6 ) == "eur://" ) {
        std::string filename = distribution.substr( 6 );
        std::ifstream f( filename.c_str() );
        double v[ 4 ];
        double min_z = + std::numeric_limits<double>::max();
        double max_z = - std::numeric_limits<double>::max();
        while ( f.readsome( (char *)v, 4 * sizeof( double ) ) ) {
            min_z = min( min_z, v[ 2 ] );
            max_z = max( max_z, v[ 2 ] );
        }

        positions.resize( 0 );
        weights.resize( 0 );
        f.seekg( 0, std::ios_base::beg );
        while ( f.readsome( (char *)v, 4 * sizeof( double ) ) ) {
            if ( v[ 2 ] < min_z + 0.05 * ( max_z - min_z ) ) {
                positions.push_back( { v[ 0 ], v[ 1 ] } );
                weights.push_back( 1 );
            }
        }
        return;
    }

    TODO;
}

