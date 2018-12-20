#pragma once

#include "system/Span.h"
#include <functional>
#include "Point3.h"
#include "Point2.h"
#include "Node.h"

template<int _dim,class _TF=double,int _nb_bits_per_axis=30,class _AdditionnalInfo=std::size_t>
struct ZIndexGridCarac {
    static constexpr int     nb_bits_per_axis = _nb_bits_per_axis;
    using                    AdditionnalInfo  = _AdditionnalInfo;
    static constexpr int     dim              = _dim;
    using                    TF               = _TF;
};

/**
*/
template<class Carac=ZIndexGridCarac<2>>
class ZIndexGrid {
public:
    using                    AdditionnalInfo  = typename Carac::AdditionnalInfo;
    using                    TF               = typename Carac::TF;
    static constexpr int     nb_bits_per_axis = Carac::nb_bits_per_axis;
    static constexpr int     dim              = Carac::dim;

    static constexpr int     sizeof_zcoords   = dim * nb_bits_per_axis;
    using                    RadixSortTmps    = std::vector<std::array<std::size_t,256>>;
    using                    Point            = typename std::conditional<dim==3,Point3<TF>,Point2<TF>>::type;
    using                    TNode            = Node<Carac>;
    using                    IZ               = typename TNode::IZ;

    /**/                     ZIndexGrid       ();
    template<class F> void   fill             ( std::size_t n, const F &func, Point min_point, Point max_point );
    template<class F> void   fill             ( std::size_t n, const F &func );


    template<class P> TNode  proj_node_for    ( P coords, AdditionnalInfo info ) const;
    void                     find_nodes       ( std::vector<TNode> &pos_to_find, const std::function<void( std::size_t cell_index, AdditionnalInfo info )> &f );
    Span<std::size_t>        neighbors        ( const TNode *node ) const;
    TF                       cell_size        ( const TNode *node ) const;
    template<class P> TNode  node_for         ( P coords, AdditionnalInfo info ) const;
    Point                    cell_pos         ( const TNode *node ) const;
    Point                    cell_pos         ( std::size_t index ) const { return cell_pos( &node( index ) ); }
    std::size_t              nb_nodes         () const { return nodes.size() - 1; }
    const TNode&             node             ( std::size_t index ) const { return nodes[ index ]; }

    const TNode*             begin            () const { return nodes.data(); }
    const TNode*             end              () const { return nodes.data() + nb_nodes(); }

private:
    std::vector<std::size_t> ng_indices;
    std::vector<std::size_t> ng_offsets;
    TF                       max_length;
    Point                    min_point;
    Point                    max_point;
    RadixSortTmps            rs_tmps;
    std::vector<TNode>       tmp_0;
    std::vector<TNode>       tmp_1;
    std::vector<TNode>       nodes;
};


#include "ZIndexGrid.tcc"
