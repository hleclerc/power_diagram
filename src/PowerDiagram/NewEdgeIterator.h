#pragma once

#include "ConvexPolyhedron3.h"

/**
*/
template<class TF>
class NewEdgeIterator {
public:
    using Edge = typename ConvexPolyhedron3<TF>::Edge;
    using TI   = typename ConvexPolyhedron3<TF>::TI;

    NewEdgeIterator( const std::vector<Edge> &old_edges, const std::vector<TI> &num_old_edges, std::vector<Edge> &new_edges )  : num_old_edges( num_old_edges ), old_edges( old_edges ), new_edges( new_edges ) {
        reset();
    }

    bool at_the_end() const {
        return num_old == num_old_edges.size();
    }

    NewEdgeIterator &operator++() {
        if ( ++num_ne >= 2 || old_edges[ num_old_edges[ num_old ] ].nedge[ num_ne ] == TI( -1 ) ) {
            ++num_old;
            _set_to_next_valid_pair();
        }
        return *this;
    }

    Edge *operator->() {
        return &new_edges[ num_edge() ];
    }

    Edge &operator*() {
        return new_edges[ num_edge() ];
    }

    TI num_edge() const {
        return old_edges[ num_old_edges[ num_old ] ].nedge[ num_ne ];
    }

    Edge &last() const {
        for( TI no = num_old_edges.size(); no--; ) {
            for( TI ne = 2; ne--; ) {
                TI tr = old_edges[ num_old_edges[ no ] ].nedge[ ne ];
                if ( tr != TI( -1 ) )
                    return new_edges[ tr ];
            }
        }
        // should not append
        return new_edges.back();
    }

    void reset() {
        num_old = 0;
        _set_to_next_valid_pair();
    }

    bool loop_incr() {
        operator++();
        if ( at_the_end() ) {
            reset();
            return true;
        }
        return false;
    }

private:
    /// try to start with a valid ( num_old, num_ne ) pair
    void _set_to_next_valid_pair() {
        for( ; num_old < num_old_edges.size(); ++num_old ) {
            if ( old_edges[ num_old_edges[ num_old ] ].nedge[ 0 ] != TI( -1 ) ) {
                num_ne = 0;
                break;
            }
        }
    }

    const std::vector<TI>   &num_old_edges;
    const std::vector<Edge> &old_edges;
    std::vector<Edge>       &new_edges;
    TI                       num_old;
    TI                       num_ne;
};
