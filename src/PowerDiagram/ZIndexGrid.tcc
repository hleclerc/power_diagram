#include "system/RadixSort.h"
#include "system/Assert.h"
#include "ZIndexGrid.h"
#include <algorithm>

template<class Carac>
ZIndexGrid<Carac>::ZIndexGrid() {
}

template<class Carac>
void ZIndexGrid<Carac>::find_nodes( std::vector<TNode> &pos_to_find, const std::function<void(std::size_t,AdditionnalInfo)> &f ) {
    std::size_t nb_pos = pos_to_find.size();
    pos_to_find.reserve( 2 * nb_pos );
    TNode *out = radix_sort( pos_to_find.data() + nb_pos, pos_to_find.data(), nb_pos, N<(sizeof_zcoords+7)/8>(), rs_tmps );
    for( std::size_t i = 0, j = 0; i < nb_pos; ++i ) {
        while ( nodes[ j ].zcoords <= out[ i ].zcoords )
            ++j;
        f( j - 1, out[ i ].info );
    }
}

template<class Carac> template<class P>
typename ZIndexGrid<Carac>::TNode ZIndexGrid<Carac>::node_for( P coords, AdditionnalInfo info ) const {
    std::array<IZ,dim> c;
    for( std::size_t d = 0; d < dim; ++d )
        c[ d ] = ( 1 << nb_bits_per_axis ) * ( coords[ d ] - min_point[ d ] ) / max_length;
    return { TNode::zcoords_for( c ), info };
}

template<class Carac> template<class P>
typename ZIndexGrid<Carac>::TNode ZIndexGrid<Carac>::proj_node_for( P coords, AdditionnalInfo info ) const {
    for( std::size_t d = 0; d < dim; ++d ) {
        if ( coords[ d ] < min_point[ d ] )
            coords[ d ] = min_point[ d ];
        else if ( coords[ d ] - min_point[ d ] >= max_length )
            coords[ d ] = min_point[ d ] + max_length * ( 1 - std::numeric_limits<TF>::max() );
    }
    return node_for( coords, info );
}

template<class Carac>
Span<std::size_t> ZIndexGrid<Carac>::neighbors( const TNode *node ) const {
    std::size_t index = node - nodes.data();
    return { ng_indices.data() + ng_offsets[ index + 0 ], ng_indices.data() + ng_offsets[ index + 1 ] };
}

template<class Carac> template<class F>
void ZIndexGrid<Carac>::fill( std::size_t n, const F &func ) {
    using std::min;
    using std::max;

    Point min_point;
    Point max_point;
    for( std::size_t d = 0; d < dim; ++d ) {
        min_point[ d ] = + std::numeric_limits<TF>::max();
        max_point[ d ] = - std::numeric_limits<TF>::max();
    }
    func( [&]( auto coords, AdditionnalInfo ) {
        for( std::size_t d = 0; d < dim; ++d ) {
            min_point[ d ] = min( min_point[ d ], coords[ d ] );
            max_point[ d ] = max( max_point[ d ], coords[ d ] );
        }
    } );

    fill( n, func, min_point, max_point );
}

template<class Carac> template<class F>
void ZIndexGrid<Carac>::fill( std::size_t n, const F &func, Point min_point, Point max_point ) {
    // In terms of execution speed, all the steps are tied
    using std::min;
    using std::max;

    // min and max points
    max_length = 0;
    for( std::size_t d = 0; d < dim; ++d )
        max_length = max( max_length, max_point[ d ] - min_point[ d ] );
    max_length *= 1 + std::numeric_limits<TF>::epsilon();

    for( std::size_t d = 0; d < dim; ++d )
        max_point[ d ] = min_point[ d ] + max_length;

    this->min_point = min_point;
    this->max_point = max_point;

    // coordinates
    tmp_0.resize( 0 );
    tmp_0.reserve( n );
    func( [&]( auto coords, AdditionnalInfo info ) {
        tmp_0.push_back( node_for( coords, info ) );
    } );

    tmp_1.reserve( tmp_0.size() );
    TNode *tmp = radix_sort( tmp_1.data(), tmp_0.data(), tmp_0.size(), N<(sizeof_zcoords+7)/8>(), rs_tmps );
    // std::sort( tmp_zcoords.begin(), tmp_zcoords.end(), [&]( auto a, auto b ) { return a.zcoords < b.zcoords; } );

    // nodes
    nodes.resize( 0 );
    nodes.reserve( 2 * n );
    if ( tmp_0.size() == 1 ) {
        TODO;
    } else {
        int level = 0;
        for( IZ prev_z = 0, index = 1; ; ) {
            if ( index == tmp_0.size() ) {
                while ( prev_z < ( IZ( 1 ) << sizeof_zcoords ) ) {
                    for( ; ; ++level ) {
                        IZ m = IZ( 1 ) << dim * ( level + 1 );
                        if ( prev_z & ( m - 1 ) ) {
                            IZ new_prev_z = prev_z + ( IZ( 1 ) << dim * level );
                            nodes.push_back( {
                                prev_z,
                                index == tmp_0.size() &&
                                    tmp[ index - 1 ].zcoords >= prev_z &&
                                    tmp[ index - 1 ].zcoords < new_prev_z ? tmp[ index++ - 1 ].info : -1
                            } );
                            prev_z = new_prev_z;
                            break;
                        }
                    }
                }
                break;
            }

            // level too high ?
            for( ; ; --level ) {
                IZ m = IZ( 1 ) << dim * ( level + 1 );
                if ( tmp[ index ].zcoords >= prev_z + m )
                    break;
                ASSERT( level, "Seems not possible to have only one point per cell" );
            }

            // look for a level before the one that will take the two next points or that will lead to an illegal cell
            for( ; ; ++level ) {
                IZ m = IZ( 1 ) << dim * ( level + 1 );
                if ( tmp[ index ].zcoords < prev_z + m || ( prev_z & ( m - 1 ) ) ) {
                    IZ new_prev_z = prev_z + ( IZ( 1 ) << dim * level );
                    nodes.push_back( {
                        prev_z,
                        tmp[ index - 1 ].zcoords >= prev_z && tmp[ index - 1 ].zcoords < new_prev_z ? tmp[ index++ - 1 ].info : -1
                    } );
                    prev_z = new_prev_z;
                    break;
                }
            }
        }
    }

    nodes.push_back( { IZ( 1 ) << sizeof_zcoords, 0 } );

    // sorted neighbor requests
    tmp_0.reserve( dim * nb_nodes() );
    tmp_1.reserve( dim * nb_nodes() );
    tmp_0.resize( 0 );
    for( std::size_t index = 0; index < nb_nodes(); ++index ) {
        constexpr IZ f00 = ~ ( ( IZ( 1 ) << sizeof_zcoords ) - 1 );
        StaticRange<dim>::for_each( [&]( auto d ) {
            IZ nz = nodes[ index ].ng_zcoord( d );
            if ( ( nz & f00 ) == 0 )
                tmp_0.push_back( { nz, index } );
        } );
    }

    tmp = radix_sort( tmp_1.data(), tmp_0.data(), tmp_0.size(), N<(sizeof_zcoords+7)/8>(), rs_tmps );

    // helper function to get the neighbors for each node
    auto for_each_ng = [&]( auto cb ) {
        for( std::size_t i = 0, j = 0; i < tmp_0.size(); ++i ) {
            // find first node with zcoords > tmp_0[ i ].zcoords
            while ( nodes[ j ].zcoords <= tmp[ i ].zcoords )
                ++j;

            // first node
            std::size_t index_node = tmp[ i ].info;
            std::size_t index_nbor = j - 1;
            cb( index_node, index_nbor );

            // next touching ones
            IZ off = nodes[ index_node + 1 ].zcoords - nodes[ index_node ].zcoords;
            IZ lim = nodes[ index_nbor ].zcoords + off;
            if ( nodes[ index_nbor + 1 ].zcoords < lim ) {
                std::array<IZ,dim> tgts;
                StaticRange<dim>::for_each( [&]( auto d ) {
                    tgts[ d ] = nodes[ index_node ].zcoords & TNode::template _ZcoordsOnesOnAxis<d.val>::value;
                    tgts[ d ] = ( tgts[ d ] | TNode::template _ZcoordsZerosOnAxis<d.val>::value ) + off;
                    tgts[ d ] = tgts[ d ] & TNode::template _ZcoordsOnesOnAxis<d.val>::value;
                } );

                ++index_nbor;
                do {
                    bool touching = false;
                    StaticRange<dim>::for_each( [&]( auto d ) {
                        IZ val = nodes[ index_nbor ].zcoords & TNode::template _ZcoordsOnesOnAxis<d.val>::value;
                        touching |= val == tgts[ d ];
                    } );
                    if ( touching )
                        cb( index_node, index_nbor );
                } while ( nodes[ ++index_nbor ].zcoords < lim );
            }
        }
    };

    //
    ng_offsets.resize( nodes.size() );
    for( std::size_t i = 0; i < nodes.size(); ++i )
        ng_offsets[ i ] = 0;

    // get count
    for_each_ng( [&]( std::size_t index_node, std::size_t index_nbor ) {
        ++ng_offsets[ index_node ];
        ++ng_offsets[ index_nbor ];
    } );

    // suffix scan
    for( std::size_t i = 0, acc = 0; i < nodes.size(); ++i ) {
        std::size_t v = acc;
        acc += ng_offsets[ i ];
        ng_offsets[ i ] = v;
    }

    // get indices
    ng_indices.resize( ng_offsets.back() );
    for_each_ng( [&]( std::size_t index_node, std::size_t index_nbor ) {
        ng_indices[ ng_offsets[ index_node ]++ ] = index_nbor;
        ng_indices[ ng_offsets[ index_nbor ]++ ] = index_node;
    } );

    // shift ng_offsets (to get the suffix scan again)
    if ( ng_offsets.size() ) {
        for( std::size_t i = ng_offsets.size(); --i; )
            ng_offsets[ i ] = ng_offsets[ i - 1 ];
        ng_offsets[ 0 ] = 0;
    }
}

template<class Carac>
typename ZIndexGrid<Carac>::TF ZIndexGrid<Carac>::cell_size( const TNode *node ) const {
    TF step = max_length / ( IZ( 1 ) << nb_bits_per_axis );
    IZ d = node[ 1 ].zcoords - node[ 0 ].zcoords;
    return step * pow( d, TF( 1 ) / dim );;
}

template<class Carac>
typename ZIndexGrid<Carac>::Point ZIndexGrid<Carac>::cell_pos( const TNode *node ) const {
    Point pt;
    TF step = max_length / ( IZ( 1 ) << nb_bits_per_axis );
    StaticRange<dim>::for_each( [&]( auto d ) {
        pt[ d ] = min_point[ d ] + step * node->int_coord( d );
    } );
    return pt;
}
