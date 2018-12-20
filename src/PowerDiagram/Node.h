#pragma once

#include "system/StaticRange.h"
#include "system/Stream.h"

/**
*/
template<class Carac>
class Node {
public:
    static constexpr int  nb_bits_per_axis = Carac::nb_bits_per_axis;
    using                 AdditionnalInfo  = typename Carac::AdditionnalInfo;
    static constexpr int  dim              = Carac::dim;

    static constexpr int  sizeof_zcoords   = dim * nb_bits_per_axis;
    using                 IZ               = std::uint64_t;

    /**/                  Node             ( IZ zcoords, AdditionnalInfo info ) : zcoords( zcoords ), info( info ) {}
    /**/                  Node             () {}

    // info
    void                  write_to_stream  ( std::ostream &os, int shift = 0 ) const;

    IZ                    operator>>       ( int shift ) const { return zcoords >> shift; }
    template<int axis> IZ int_coord        ( N<axis> ) const;
    template<int axis> IZ ng_zcoord        ( N<axis> ) const;


    // construction helpers
    static IZ             zcoords_for      ( std::size_t x );
    static IZ             zcoords_for      ( std::size_t x, std::size_t y );
    static IZ             zcoords_for      ( std::size_t x, std::size_t y, std::size_t z );

    static IZ             zcoords_for      ( const std::array<std::size_t,1> &v ) { return zcoords_for( v[ 0 ] ); }
    static IZ             zcoords_for      ( const std::array<std::size_t,2> &v ) { return zcoords_for( v[ 0 ], v[ 1 ] ); }
    static IZ             zcoords_for      ( const std::array<std::size_t,3> &v ) { return zcoords_for( v[ 0 ], v[ 1 ], v[ 2 ] ); }

    /// Ex: axis = 0, dim = 3 (i.e. x) => 000... for level and free_bits ++ 110110110...
    template<int num_axis,int _cur_bit = sizeof_zcoords - 1>
    struct _ZcoordsZerosOnAxis {
        static constexpr IZ v_loc = IZ( _cur_bit % dim == num_axis ? 0 : 1 ) << _cur_bit;
        static constexpr IZ value = v_loc | _ZcoordsZerosOnAxis<num_axis,_cur_bit-1>::value;
    };
    template<int num_axis>
    struct _ZcoordsZerosOnAxis<num_axis,-1> {
        static constexpr IZ value = 0;
    };

    /// Ex: axis = 0, dim = 3 (i.e. x) => 000... for level and free_bits ++ 001001001...
    template<int num_axis,int _cur_bit = sizeof_zcoords - 1>
    struct _ZcoordsOnesOnAxis {
        static constexpr IZ v_loc = IZ( _cur_bit % dim == num_axis ? 1 : 0 ) << _cur_bit;
        static constexpr IZ value = v_loc | _ZcoordsOnesOnAxis<num_axis,_cur_bit-1>::value;
    };
    template<int num_axis>
    struct _ZcoordsOnesOnAxis<num_axis,-1> {
        static constexpr IZ value = 0;
    };

    // attributes
    IZ                    zcoords;
    AdditionnalInfo       info;
};


#include "Node.tcc"
