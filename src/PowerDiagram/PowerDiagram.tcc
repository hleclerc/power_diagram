#include "system/ThreadPool.h"
#include "PowerDiagram.h"
#include <queue>

template<class Pc>
PowerDiagram<Pc>::PowerDiagram( const Pc &power_diagram_carac ) : _power_diagram_carac( power_diagram_carac ), _space_partitioner( &_power_diagram_carac.sp ) {
    _ball_cut = false;
}

template<class Pc>
void PowerDiagram<Pc>::add_convex_shape( const std::vector<CutInfo> &convex_bounds, TF orig_radius ) {
    PT c;
    for( int i = 0; i < dim; ++i )
        c[ i ] = 0;
    LC lc( c, orig_radius, TI( -1 ) );
    for( const CutInfo &ci : convex_bounds )
        lc.plane_cut( ci.O, normalized( ci.N ), TI( -1 ) );
    add_convex_shape( lc );
}

template<class Pc>
void PowerDiagram<Pc>::add_convex_shape( const LC &lc ) {
    _convex_bounds.push_back( lc );
    _convex_bounds.back().set_cut_ids( TI( -1 ) );
}

template<class Pc>
void PowerDiagram<Pc>::add_bounding_simplex( PT center, TF radius ) {
    LC lc( center, radius, TI( -1 ) );
    add_convex_shape( lc );
}

template<class Pc>
void PowerDiagram<Pc>::add_box_shape( PT _p0, PT _p1 ) {
    PT p0 = min( _p0, _p1 );
    PT p1 = max( _p0, _p1 );
    _add_box_shape( p0, p1, N<dim>() );
}


template<class Pc>
void PowerDiagram<Pc>::add_dirac( PT pos, TF weight ) {
    _diracs.push_back( Dirac{ pos, weight } );
}

template<class Pc>
void PowerDiagram<Pc>::set_ball_cut( bool ball_cut ) {
    _ball_cut = ball_cut;
}

template<class Pc>
typename PowerDiagram<Pc>::TI PowerDiagram<Pc>::nb_convex_bounds() const {
    return _convex_bounds.size();
}

template<class Pc>
typename PowerDiagram<Pc>::TF PowerDiagram<Pc>::domain_measure() const {
    TF res = 0;
    for( const LC &lc : _convex_bounds )
        res += lc.measure();
    return res;
}

template<class Pc>
void PowerDiagram<Pc>::write_to_stream( std::ostream &os ) const {
    os << _diracs;
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_der_measures_ap( const Fu &fu, TF *volumes, std::vector<std::pair<TI,TF>> *derivatives, TF epsilon ) const {
    get_measures( volumes );

    std::vector<TF> mod_volumes( nb_diracs() );
    for( TI i = 0; i < nb_diracs(); ++i ) {
        TF ow = dirac( i ).weight;
        dirac( i ).weight -= epsilon;
        get_measures( fu, mod_volumes.data() );
        dirac( i ).weight = ow;

        for( TI j = 0; j < nb_diracs(); ++j )
            if ( TF ap_der = ( volumes[ j ] - mod_volumes[ j ] ) / epsilon )
                derivatives[ j ].push_back( { i, ap_der } );
    }
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_der_measures( const Fu &fu, TF *measures, std::vector<std::pair<TI,TF>> *derivatives ) const {
    using std::sqrt;

    for( std::size_t n = 0; n < nb_diracs(); ++n )
        measures[ n ] = 0;

    for_each_cell( [&]( const LC &cs, TI index_0, int num_thread ) {
        measures[ index_0 ] += cs.measure( fu );

        // list of contributions
        TF der_0 = 0;
        PT d0_center = _diracs[ index_0 ].pos;
        TF d0_weight = _diracs[ index_0 ].weight;
        std::vector<std::pair<TI,TF>> &dv = derivatives[ index_0 ];
        cs.for_each_boundary_measure( fu, [&]( TF boundary_measure, TI index_1 ) {
            if ( index_1 == TI( -1 ) )
                return;
            if ( index_0 == index_1 ) {
                der_0 += 0.5 * boundary_measure / sqrt( d0_weight );
            } else {
                TI m_index_1 = index_1 % nb_diracs();
                PT d1_center = _diracs[ m_index_1 ].pos;
                if ( size_t nu = index_1 / nb_diracs() )
                    d1_center = transformation( _tranformations[ nu - 1 ], d1_center );

                TF dist = norm_2( d0_center - d1_center );
                TF der_1 = 0.5 * boundary_measure / dist;
                der_0 += der_1;
                dv.push_back( std::make_pair( m_index_1, - der_1 ) );
            }
        } );
        dv.push_back( std::make_pair( index_0, der_0 ) );
    } );

    // sort and sum the contributions
    for( std::size_t n = 0; n < nb_diracs(); ++n )
        _sort_and_sum( derivatives[ n ] );
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_centroid_contribs( const Fu &fu, PT *centroid_mul_by_volumes, TF *volumes ) const {
    for_each_cell( [&]( const LC &cs, TI num, int num_thread ) {
        cs.add_centroid_contrib( centroid_mul_by_volumes[ num ], volumes[ num ], fu );
    } );
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_centroids_ap( const Fu &fu, PT *centroids, TI n ) const {
    for_each_cell( [&]( const LC &cs, TI num, int num_thread ) {
        centroids[ num ] = cs.centroid_ap( fu, n );
    } );
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_measures( const Fu &fu, TF *measures ) const {
    for( std::size_t n = 0; n < nb_diracs(); ++n )
        measures[ n ] = 0;

    for_each_cell( [&]( const LC &cs, TI num, int num_thread ) {
        measures[ num ] += cs.measure( fu );
    } );
}

template<class Pc> template<class Fu>
void PowerDiagram<Pc>::get_measures_ap( const Fu &fu, TF *volumes, TI n ) const {
    for( std::size_t n = 0; n < nb_diracs(); ++n )
        volumes[ n ] = 0;

    for_each_cell( [&]( const LC &cs, TI num, int num_thread ) {
        volumes[ num ] += cs.volume_ap( fu, n );
    } );
}

template<class Pc> template<int c>
void PowerDiagram<Pc>::display( VtkOutput<c,TF> &vo, const FC &cell_data_func, bool filled, TF max_ratio_area_error ) const {
    std::mutex m;
    for_each_cell( [&]( const LC &cs, TI num, int num_thread ) {
        typename VtkOutput<c,TF>::CV cv;
        cell_data_func( cv.data(), cs, num );
        cs.display( vo, cv, filled, max_ratio_area_error, false, &m );
    } );
}

template<class Pc> template<class F>
void PowerDiagram<Pc>::for_each_cell( const F &f ) const {
    using VisitedCells = typename Sp::VisitedCells;
    using CellHandle   = typename Sp::CellHandle;
    using std::sqrt;
    using std::max;
    using std::min;

    struct CellHandleWithDist {
        bool operator<( const CellHandleWithDist &that ) const {
            return dist > that.dist;
        }
        CellHandle cell;
        double     dist;
    };

    auto plane_cut = [&]( LC &lc, std::size_t b1_index ) {
        const Dirac &d1 = _diracs[ b1_index ];
        PT V = d1.pos - lc.info.dirac.pos;
        TF D = norm_2( V );
        PT N = 1.0 / D * V;
        TF x = 0.5 * ( D + ( lc.info.dirac.weight - d1.weight ) / D );
        lc.plane_cut( lc.info.dirac.pos + x * N, N, b1_index );
    };

    struct TestWithBoundary {
        bool operator()( PT A, PT B ) const {
            TF c2 = norm_2_p2( B - A );
            TF c1 = dot( B - A, A - p );
            TF c0 = norm_2_p2( A - p ) - norm_2_p2( d0_pos - p ) + min_d0_weight - max_d1_weight;
            if ( c2 ) {
                TF um = - c1 / c2;
                if ( um > 0 && um < 1 && c2 * pow( um, 2 ) + 2 * c1 * um + c0 < 0 )
                    return true;
            }
            return c0 < 0 || c2 + 2 * c1 + c0 < 0;

        }
        TF max_d1_weight;
        TF min_d0_weight;
        PT d0_pos;
        PT p;
    };

    // true if a dirac in b1 (given max_weight in b1 and its neighbors) may cut lc
    auto may_cut_lc = [&]( const LC &lc, typename Sp::CellHandle b1 ) -> bool {
        TestWithBoundary test;
        test.max_d1_weight = _space_partitioner.max_weight( b1 );
        test.min_d0_weight = lc.info.dirac.weight;
        test.d0_pos        = lc.info.dirac.pos;

        if ( _ball_cut ) {
            PT V = _space_partitioner.cell_center( b1 ) - lc.info.dirac.pos;
            TF s = 1.415 * _space_partitioner.cell_size( b1 );
            TF D = norm_2( V );
            TF d = lc.info.dirac.weight - test.max_d1_weight;
            TF x = max( TF( 0 ), TF( 0.5 ) * ( D + d / D ) - s );
            if ( x * x > lc.info.dirac.weight )
                return false;
        }

        for( std::size_t num_lc_point = 0; num_lc_point < lc.nb_points(); ++num_lc_point ) {
            test.p = lc.point( num_lc_point );
            if ( _space_partitioner.test_with_boundaries( b1, test ) )
                return true;
        }
        return false;
    };

    // true if a dirac in b1 may cut lc (made from d0 so far)
    auto may_cut_cs = [&]( const std::vector<LC> &cs, std::size_t nb_cs, Span<std::size_t> b0_indices, typename Sp::CellHandle b1 ) -> bool {
        for( const LC &lc : cs )
            if ( may_cut_lc( lc, b1 ) )
                return true;
        return false;
    };


    int nb_threads = thread_pool.nb_threads();
    std::vector<std::vector<LC>> cs_vec( nb_threads );
    std::vector<VisitedCells> visited_vec( nb_threads );
    std::vector<std::priority_queue<CellHandleWithDist>> front_vec( nb_threads );

    for( std::vector<LC> &cs : cs_vec )
        cs.resize( 32 );

    if ( _convex_bounds.empty() ) {
        TODO;
    } else {
        int nb_jobs = 4 * nb_threads;
        for( const LC &convex_bounds : _convex_bounds ) {
            // initialization of _point_grid
            _space_partitioner.init( _diracs );

            // traversal
            thread_pool.execute( nb_jobs, [&]( std::size_t job_index, int num_thread ) {
                std::priority_queue<CellHandleWithDist> &front = front_vec[ num_thread ];
                VisitedCells &visited = visited_vec[ num_thread ];
                std::vector<LC> &cs = cs_vec[ num_thread ];

                _space_partitioner.for_each_cell( [&]( typename Sp::CellHandle b0, Span<std::size_t> b0_indices ) {
                    _space_partitioner.init_visited_cells( visited );
                    visited.append( b0 );

                    if ( b0_indices.size() == 1 ) {
                        std::size_t b0_index = b0_indices[ 0 ];
                        if ( b0_index < _diracs.size() ) {
                            LC &lc = cs[ 0 ];
                            lc = convex_bounds;
                            lc.info = { _diracs[ b0_index ], b0_index };

                            // direct neighbors
                            _space_partitioner.for_each_direct_neighbor( b0, [&]( typename Sp::CellHandle b1 ) __attribute__((always_inline)) {
                                visited.append( b1 );
                                for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                    plane_cut( lc, b1_index );
                            } );

                            // ring 2
                            _space_partitioner.for_each_ring_2_neighbor( b0, [&]( typename Sp::CellHandle b1, auto ln_tuple ) {
                                visited.append( b1 );

                                // if it's not possible to cut with a dirac in b1, we're not going to test farther cells
                                if ( may_cut_lc( lc, b1 ) ) {
                                    for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                        plane_cut( lc, b1_index );

                                    // update the front
                                    StaticRange<std::tuple_size<decltype(ln_tuple)>::value>::for_each( [&]( auto num_ng ) {
                                        auto off = std::get<num_ng.val>( ln_tuple );
                                        typename Sp::CellHandle ng_cell = _space_partitioner.ng_cell( b0, off );
                                        if ( ng_cell && visited( ng_cell ) == false ) {
                                            front.push( { ng_cell, _space_partitioner.min_sq_dist( b0, ng_cell ) } );
                                            visited.append( ng_cell );
                                        }
                                    } );
                                }
                            } );

                            // empty the front
                            while ( ! front.empty() ) {
                                CellHandle b1 = front.top().cell;
                                front.pop();

                                if ( may_cut_lc( lc, b1 ) ) {
                                    for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                        plane_cut( lc, b1_index );

                                    // update the front
                                    _space_partitioner.for_each_neighbor( b1, b0, [&]( const CellHandle &ng_cell ) {
                                        if ( visited( ng_cell ) == false ) {
                                            front.push( { ng_cell, _space_partitioner.min_sq_dist( b0, ng_cell ) } );
                                            visited.append( ng_cell );
                                        }
                                    } );
                                }
                            }

                            // ball_cut
                            if ( _ball_cut )
                                lc.sphere_cut( lc.info.dirac.pos, std::sqrt( lc.info.dirac.weight ), b0_index );

                            // callback
                            f( lc, b0_index, num_thread );
                        }
                    } else if ( b0_indices.size() ) {
                        // init a new Laguerre cell for each direct dirac (i.e. that does not come from a symetry) in the b0 box
                        std::size_t nb_cs = 0;
                        for( std::size_t b0_index : b0_indices ) {
                            if ( b0_index < _diracs.size() ) {
                                if ( nb_cs == cs.size() )
                                    cs.resize( 1 + 2 * nb_cs );
                                cs[ nb_cs ] = convex_bounds;
                                cs[ nb_cs ].info = { _diracs[ b0_index ], b0_index };
                                ++nb_cs;
                            }
                        }

                        // other diracs in the b0 cell
                        if ( nb_cs > 1 ) {
                            for( std::size_t num_in_cs_0 = 0; num_in_cs_0 < nb_cs; ++num_in_cs_0 )
                                for( std::size_t num_in_cs_1 = 0; num_in_cs_1 < nb_cs; ++num_in_cs_1 )
                                    if ( num_in_cs_1 != num_in_cs_0 )
                                        plane_cut( cs[ num_in_cs_0 ], cs[ num_in_cs_1 ].info.index );
                        }

                        // direct neighbors
                        _space_partitioner.for_each_direct_neighbor( b0, [&]( typename Sp::CellHandle b1 ) {
                            visited.append( b1 );
                            for( std::size_t num_in_cs = 0; num_in_cs < nb_cs; ++num_in_cs )
                                for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                    plane_cut( cs[ num_in_cs ], b1_index );
                        } );

                        // ring 2
                        _space_partitioner.for_each_ring_2_neighbor( b0, [&]( typename Sp::CellHandle b1, auto ln_tuple ) {
                            visited.append( b1 );

                            // if it's not possible to cut with a dirac in b1, we're not going to test farther cells
                            if ( may_cut_cs( cs, nb_cs, b0_indices, b1 ) ) {
                                for( std::size_t num_in_cs = 0; num_in_cs < nb_cs; ++num_in_cs )
                                    for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                        plane_cut( cs[ num_in_cs ], b1_index );

                                // update the front
                                StaticRange<std::tuple_size<decltype(ln_tuple)>::value>::for_each( [&]( auto num_ng ) {
                                    auto off = std::get<num_ng.val>( ln_tuple );
                                    typename Sp::CellHandle ng_cell = _space_partitioner.ng_cell( b0, off );
                                    if ( ng_cell && visited( ng_cell ) == false ) {
                                        front.push( { ng_cell, _space_partitioner.min_sq_dist( b0, ng_cell ) } );
                                        visited.append( ng_cell );
                                    }
                                } );
                            }
                        } );

                        // empty the front
                        while ( ! front.empty() ) {
                            CellHandle b1 = front.top().cell;
                            front.pop();

                            if ( may_cut_cs( cs, nb_cs, b0_indices, b1 ) ) {
                                for( std::size_t num_in_cs = 0; num_in_cs < nb_cs; ++num_in_cs )
                                    for( std::size_t b1_index : _space_partitioner.items_in( b1 ) )
                                        plane_cut( cs[ num_in_cs ], b1_index );

                                // update the front
                                _space_partitioner.for_each_neighbor( b1, b0, [&]( const CellHandle &ng_cell ) {
                                    if ( visited( ng_cell ) == false ) {
                                        front.push( { ng_cell, _space_partitioner.min_sq_dist( b0, ng_cell ) } );
                                        visited.append( ng_cell );
                                    }
                                } );
                            }
                        }

                        // ball_cut
                        if ( _ball_cut ) {
                            for( std::size_t num_in_cs = 0; num_in_cs < b0_indices.size(); ++num_in_cs ) {
                                TI b0_index = b0_indices[ num_in_cs ];
                                const Dirac &d0 = _diracs[ b0_index ];
                                LC &lc = cs[ num_in_cs ];
                                lc.sphere_cut( d0.pos, std::sqrt( d0.weight ), b0_index );
                            }
                        }

                        // callbacks
                        for( std::size_t num_in_cs = 0; num_in_cs < b0_indices.size(); ++num_in_cs ) {
                            TI b0_index = b0_indices[ num_in_cs ];
                            LC &lc = cs[ num_in_cs ];
                            f( lc, b0_index, num_thread );
                        }
                    }
                }, job_index, nb_jobs );
            } );
        }
    }
}

template<class Pc>
void PowerDiagram<Pc>::_sort_and_sum( std::vector<std::pair<TI,TF> > &dv ) {
    std::sort( dv.begin(), dv.end(), [&]( auto a, auto b ) {
        return a.first < b.first;
    } );

    size_t k = -1;
    for( TI j = 0; j < dv.size(); ++j ) {
        if ( j == 0 || dv[ j ].first != dv[ j - 1 ].first ) {
            dv[ ++k ].first = dv[ j ].first;
            dv[ k ].second = dv[ j ].second;
        } else {
            dv[ k ].second += dv[ j ].second;
        }
    }
    dv.resize( k + 1 );
}

template<class Pc>
void PowerDiagram<Pc>::_add_box_shape( PT p0,PT p1, N<2> ) {
    add_convex_shape( {
        { p0, { -1,  0 } },
        { p0, {  0, -1 } },
        { p1, { +1,  0 } },
        { p1, {  0, +1 } },
    }, 3 * max( p1 - p0 ) );
}

template<class Pc>
void PowerDiagram<Pc>::_add_box_shape( PT p0,PT p1, N<3> ) {
    add_convex_shape( {
        { p0, { -1,  0,  0 } },
        { p0, {  0, -1,  0 } },
        { p0, {  0,  0, -1 } },
        { p1, { +1,  0,  0 } },
        { p1, {  0, +1,  0 } },
        { p1, {  0,  0, +1 } },
    }, 3 * max( p1 - p0 ) );
}

template<class Pc> template<int c>
void PowerDiagram<Pc>::display_bounds( VtkOutput<c,TF> &vo ) const {
    for( const LC &lc : _convex_bounds )
        lc.display( vo );
}
