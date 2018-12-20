#include "Node.h"

extern const std::uint32_t morton_256_2D_x[ 256 ];
extern const std::uint32_t morton_256_2D_y[ 256 ];
extern const std::uint32_t morton_256_3D_x[ 256 ];
extern const std::uint32_t morton_256_3D_y[ 256 ];
extern const std::uint32_t morton_256_3D_z[ 256 ];

template<class Carac> template<int axis>
std::size_t Node<Carac>::int_coord( N<axis> ) const {
    std::size_t res = 0;
    static_assert( sizeof_zcoords <= 64, "TODO" );
    for( unsigned i = 0; i < nb_bits_per_axis; ++i ) // unroll rec ?
        res += ( zcoords & ( IZ( 1 ) << ( dim * i + axis ) ) ) >> ( ( dim - 1 ) * i + axis );
    return res;
}

template<class Carac>
void Node<Carac>::write_to_stream( std::ostream &os, int shift ) const {
    os << "c:" << zcoords;
    StaticRange<dim>::for_each( [&]( auto d ) {
        os << " " << "xyz"[ d.val ] << ":" << ( int_coord( N<d.val>() ) >> shift );
    } );
}

template<class Carac>
typename Node<Carac>::IZ Node<Carac>::zcoords_for( std::size_t x ) {
    static_assert( dim == 1, "done only for 1D" );
    return x;
}

template<class Carac>
typename Node<Carac>::IZ Node<Carac>::zcoords_for( std::size_t x, std::size_t y ) {
    static_assert( dim == 2, "done only for 2D" );
    IZ res = 0;
    for( std::size_t o = 0; o < nb_bits_per_axis; o += 8 )
        res |= IZ( morton_256_2D_x[ ( x >> o ) & 0xFF ] |
                   morton_256_2D_y[ ( y >> o ) & 0xFF ] ) << dim *  o;
    return res;
}

template<class Carac>
typename Node<Carac>::IZ Node<Carac>::zcoords_for( std::size_t x, std::size_t y, std::size_t z ) {
    static_assert( dim == 3, "done only for 3D" );
    IZ res = 0;
    for( std::size_t o = 0; o < nb_bits_per_axis; o += 8 )
        res |= IZ( morton_256_3D_x[ ( x >> o ) & 0xFF ] |
                   morton_256_3D_y[ ( y >> o ) & 0xFF ] |
                   morton_256_3D_z[ ( z >> o ) & 0xFF ] ) << dim *  o;
    return res;
}

template<class Carac> template<int axis>
typename Node<Carac>::IZ Node<Carac>::ng_zcoord( N<axis> ) const {
    IZ off = this[ 1 ].zcoords - zcoords;

    // keep values if d != axis
    IZ ff0 = _ZcoordsZerosOnAxis<axis>::value;
    IZ res = zcoords & ff0;
    res |= ( ( zcoords | ff0 ) + off ) & ~ ff0;
    return res;
}
