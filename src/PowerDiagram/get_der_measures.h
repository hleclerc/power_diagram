#pragma once

#include "SpaceFunctions/Constant.h"
#include "system/ThreadPool.h"
#include "system/Assert.h"
#include "FunctionEnum.h"
#include <vector>

namespace PowerDiagram {

/**
   We assume that grid has already been initialized by diracs
*/
template<class TI,class TF,class Grid,class Bounds,class Diracs>
int get_der_measures( std::vector<TI> &m_offsets, std::vector<TI> &m_columns, std::vector<TF> &m_values, std::vector<TF> &v_values, Grid &grid, Bounds &bounds, const Diracs &diracs ) {
    using Pt = decltype( diracs[ 0 ].pos );

    struct DataPerThread {
        DataPerThread( std::size_t approx_nb_diracs ) {
            row_items    .reserve( 64 );
            offsets      .reserve( approx_nb_diracs );
            columns      .reserve( 5 * approx_nb_diracs );
            values       .reserve( 5 * approx_nb_diracs );
        }

        std::vector<std::pair<TI,TF>> row_items;
        std::vector<TI>               offsets;
        std::vector<TI>               columns;
        std::vector<TF>               values;
    };

    int nb_threads = thread_pool.nb_threads();
    std::vector<DataPerThread> data_per_threads( nb_threads, diracs.size() / nb_threads );
    std::vector<std::pair<int,TI>> pos_in_loc_matrices( diracs.size() ); // num dirac => num_thread, num sub row

    v_values.resize( diracs.size() );
    for( TF &v : v_values )
        v = 0;

    int err = grid.for_each_laguerre_cell( [&]( auto &lc, std::size_t num_dirac_0, int num_thread ) {
        DataPerThread &dpt = data_per_threads[ num_thread ];
        pos_in_loc_matrices[ num_dirac_0 ] = { num_thread, dpt.offsets.size() };

        // get local row_items (sorted)
        TF der_0 = 0;
        dpt.row_items.resize( 0 );
        Pt d0_center = diracs[ num_dirac_0 ].pos;
        TF d0_weight = diracs[ num_dirac_0 ].weight;
        bounds.for_each_intersection( lc, [&]( auto &cp, SpaceFunctions::Constant<TF> space_func ) {
            TF coeff = 0.5 * space_func.coeff;
            v_values[ num_dirac_0 ] += space_func.coeff * cp.measure();
            cp.for_each_boundary_measure( FunctionEnum::Unit(), [&]( TF boundary_measure, TI num_dirac_1 ) {
                if ( num_dirac_1 == TI( -1 ) )
                    return;
                if ( num_dirac_0 == num_dirac_1 ) {
                    der_0 += coeff * boundary_measure / sqrt( d0_weight );
                } else {
                    TI m_num_dirac_1 = num_dirac_1 % diracs.size();
                    Pt d1_center = diracs[ m_num_dirac_1 ].pos;
                    if ( std::size_t nu = num_dirac_1 / diracs.size() )
                        TODO; // d1_center = transformation( _tranformations[ nu - 1 ], d1_center );

                    TF dist = norm_2( d0_center - d1_center );
                    TF der_1 = coeff * boundary_measure / dist;
                    dpt.row_items.emplace_back( m_num_dirac_1, - der_1 );
                    der_0 += der_1;
                }
            } );
        } );
        dpt.row_items.emplace_back( num_dirac_0, der_0 );
        std::sort( dpt.row_items.begin(), dpt.row_items.end() );

        // save them in local sub matrix
        dpt.offsets.push_back( dpt.columns.size() );
        for( std::size_t i = 0; i < dpt.row_items.size(); ++i ) {
            if ( i + 1 < dpt.row_items.size() && dpt.row_items[ i ].first == dpt.row_items[ i + 1 ].first ) {
                dpt.row_items[ i + 1 ].second += dpt.row_items[ i ].second;
                continue;
            }
            dpt.columns.push_back( dpt.row_items[ i ].first  );
            dpt.values .push_back( dpt.row_items[ i ].second );
        }
    }, bounds.englobing_convex_polyhedron(), diracs, true /*stop if void laguerre cell*/ );
    if ( err )
        return err;

    // completion of local matrices
    std::size_t nnz = 0;
    for( DataPerThread &dpt : data_per_threads ) {
        dpt.offsets.push_back( dpt.columns.size() );
        nnz += dpt.columns.size();
    }

    // assembly
    m_offsets.resize( diracs.size() + 1 );
    m_columns.reserve( nnz );
    m_values .reserve( nnz );
    m_columns.resize( 0 );
    m_values .resize( 0 );
    for( std::size_t n = 0; n < diracs.size(); ++n ) {
        m_offsets[ n ] = m_columns.size();

        std::size_t lr = pos_in_loc_matrices[ n ].second;
        DataPerThread &dpt = data_per_threads[ pos_in_loc_matrices[ n ].first ];
        m_columns.insert( m_columns.end(), dpt.columns.data() + dpt.offsets[ lr + 0 ], dpt.columns.data() + dpt.offsets[ lr + 1 ] );
        m_values .insert( m_values .end(), dpt.values .data() + dpt.offsets[ lr + 0 ], dpt.values .data() + dpt.offsets[ lr + 1 ] );
    }
    m_offsets.back() = m_columns.size();

    return 0;
}

} // namespace PowerDiagram