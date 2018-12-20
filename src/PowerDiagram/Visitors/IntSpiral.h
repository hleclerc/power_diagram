#pragma once

#include "../system/Stream.h"
#include "../system/Assert.h"
#include "../Point3.h"
#include "../Point2.h"

extern const std::size_t    nb_values_2D;
extern const std::size_t    nb_values_3D;
extern const std::ptrdiff_t values_2D[];
extern const std::ptrdiff_t values_3D[];


/**
  Draws a spiral on a grid around (0,0,0) (no doubles, no forgotten points)

  The first points are sorted by min dist between the considered grid cell and the (0,0,0) one.

*/
template<int dim,int rd_ratio=8>
class IntSpiral {
public:
    using TI = std::ptrdiff_t;
    using Pt = std::array<ST,dim>;

    IntSpiral() {
    }

    template<class FU>
    void for_each_until( TI until, const FU &f ) {
        using std::pow;

        const std::ptrdiff_t *values = dim == 2 ? values_2D : values_3D;
        for( ; ; ) {
            TI sq_dist = *( values++ );
            if ( sq_dist > pow( until, 2 ) )
                break;

            Pt p;
            for( size_t d = 0; d < dim; ++d )
                p[ d ] = *( values++ );
            f( p );
        }
    }

    //    void operator++() {
    //        if ( ++_num < _nb_predef_values ) {
    //            const int *values = dim == 2 ? values_2D : values_3D;
    //            for( std::size_t d = 0; d < dim; ++d )
    //                _pos[ d ] = values[ ( dim + 1 ) * _num + d ];
    //            _rd2 = values[ ( dim + 1 ) * _num + dim ];
    //            return;
    //        }

    //        // transition  between precomputed and computed points
    //        if ( _num == _nb_predef_values ) {
    //            std::cerr << "TODO: optimize IntSpiral when _num > _nb_predef_values\n";

    //            ++_rd2;
    //            _mr2 = ( rd_ratio + 1 ) * _rd2 / rd_ratio;

    //            _brd = std::sqrt( double( _rd2 ) ) * ( 1 - 8 * std::numeric_limits<double>::epsilon() ) + 1;
    //            _erd = std::sqrt( double( _mr2 ) ) * ( 1 + 8 * std::numeric_limits<double>::epsilon() ) + std::sqrt( double( dim ) );
    //            for( int d = 1; d < dim; ++d )
    //                _pos[ d ] = - _erd;
    //            _pos[ 0 ] = - _erd - 1;
    //        }

    //        // change pos until valid
    //        while ( true ) {
    //            for( int d = 0;; ++d ) {
    //                if ( ++_pos[ d ] <= _erd )
    //                    break;
    //                if ( d == dim - 1 ) {
    //                    _rd2 = _mr2;
    //                    _mr2 = ( rd_ratio + 1 ) * _rd2 / rd_ratio;

    //                    _brd = std::sqrt( double( _rd2 ) ) * ( 1 - 8 * std::numeric_limits<double>::epsilon() ) + 1;
    //                    _erd = std::sqrt( double( _mr2 ) ) * ( 1 + 8 * std::numeric_limits<double>::epsilon() ) + std::sqrt( double( dim ) );
    //                    for( int d = 0; d < dim; ++d )
    //                        _pos[ d ] = - _erd;
    //                    break;
    //                }
    //                _pos[ d ] = - _erd;
    //            }

    //            ST sd = sq_dist( _pos );
    //            if ( sd >= _rd2 && sd < _mr2 )
    //                break;
    //        }
    //    }

    //    /// computation of min dist between the cell located at `pos` and the (0,0,0) one
    //    static ST sq_dist( PT pos ) {
    //        ST res = 0;
    //        for( std::size_t d = 0; d < dim; ++d ) {
    //            ST v = pos[ d ] - ( pos[ d ] > 0 ) + ( pos[ d ] < 0 );
    //            res += v * v;
    //        }
    //        return res;
    //    }


    //    PT pos() const { return _pos; } ///< position of the current cell
    //    ST rd2() const { return _rd2; } ///< min dist between the current cell and the (0,0,0) one

    //private:
    //    static constexpr int _nb_predef_values  = dim == 2 ? 8225 : 77359;

    //    ST                   _rd2; ///< square of minimum distance between current cell and the (0,0,0) one
    //    ST                   _mr2; ///<
    //    ST                   _brd; ///< square of minimum distance between current cell and the (0,0,0) one
    //    ST                   _erd; ///<
    //    PT                   _pos; ///< current cell
    //    std::size_t          _num; ///< num of points. Used at the beginning (giving an index in the table used for the first points)
};
