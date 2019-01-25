#include "ConvexPolyhedron2.h"
#include "AreaOutput.h"
#include <immintrin.h>

namespace PowerDiagram {

template<class Pc,class CI>
ConvexPolyhedron2<Pc,CI>::ConvexPolyhedron2( const EnglobingSimplex &es, CI cut_id ) {
    TODO;
}

template<class Pc,class CI>
ConvexPolyhedron2<Pc,CI>::ConvexPolyhedron2() : ConvexPolyhedron2( Box{ { -1e10, -1e10 }, { +1e10, +1e10 } } ) {
}

template<class Pc,class CI>
ConvexPolyhedron2<Pc,CI>::ConvexPolyhedron2( const Box &box, CI cut_id ) {
    sphere_center = TF( 0.5 ) * ( box.p0 + box.p1 );
    sphere_radius = -1;

    //
    nb_points = 4;

    // points
    points [ 0 ][ 0 ] = box.p0[ 0 ]; points [ 1 ][ 0 ] = box.p0[ 1 ];
    points [ 0 ][ 1 ] = box.p1[ 0 ]; points [ 1 ][ 1 ] = box.p0[ 1 ];
    points [ 0 ][ 2 ] = box.p1[ 0 ]; points [ 1 ][ 2 ] = box.p1[ 1 ];
    points [ 0 ][ 3 ] = box.p0[ 0 ]; points [ 1 ][ 3 ] = box.p1[ 1 ];

    // normals
    normals[ 0 ][ 0 ] =           0; normals[ 1 ][ 0 ] =          -1;
    normals[ 0 ][ 1 ] =          +1; normals[ 1 ][ 1 ] =           0;
    normals[ 0 ][ 2 ] =           0; normals[ 1 ][ 2 ] =          +1;
    normals[ 0 ][ 3 ] =          -1; normals[ 1 ][ 3 ] =           0;

    // cut_ids
    cut_ids     [ 0 ] =      cut_id;
    cut_ids     [ 1 ] =      cut_id;
    cut_ids     [ 2 ] =      cut_id;
    cut_ids     [ 3 ] =      cut_id;
}

template<class Pc, class CI>
ConvexPolyhedron2<Pc,CI>::ConvexPolyhedron2( const ConvexPolyhedron2 &that ) {
    nb_points = that.nb_points;

    for( size_t i = 0; i < nb_points; ++i ) points[ 0 ][ i ] = that.points[ 0 ][ i ];
    for( size_t i = 0; i < nb_points; ++i ) points[ 1 ][ i ] = that.points[ 1 ][ i ];
    for( size_t i = 0; i < nb_points; ++i ) cut_ids[ i ] = that.cut_ids[ i ];

    if ( store_the_normals ) {
        for( size_t i = 0; i < nb_points; ++i ) normals[ 0 ][ i ] = that.normals[ 0 ][ i ];
        for( size_t i = 0; i < nb_points; ++i ) normals[ 1 ][ i ] = that.normals[ 1 ][ i ];
        for( size_t i = 0; i < nb_points; ++i ) arcs.set( i, that.arcs[ i ] );
    }

    sphere_center = that.sphere_center;
    sphere_radius = that.sphere_radius;
    sphere_cut_id = that.sphere_cut_id;
}

template<class Pc,class CI>
typename ConvexPolyhedron2<Pc,CI>::TF ConvexPolyhedron2<Pc,CI>::integration_der_wrt_weight( FunctionEnum::ExpWmR2db<TF> fu, TF weight ) const {
    return integration( fu, weight ) / fu.eps;
}

template<class Pc,class CI>
template<class FU> typename ConvexPolyhedron2<Pc,CI>::TF ConvexPolyhedron2<Pc,CI>::integration_der_wrt_weight( FU, TF weight ) const {
    return 0;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_boundary_measure( FunctionEnum::ExpWmR2db<TF> e, const std::function<void(TF,CI)> &f, TF weight ) const {
    using std::sqrt;
    using std::erf;
    using std::exp;
    using std::pow;

    if ( nb_points == 0 ) {
        if ( sphere_radius >= 0 )
            f( 2 * M_PI * sphere_radius * exp( ( weight - pow( sphere_radius, 2 ) ) / e.eps ), sphere_cut_id );
        return;
    }

    for( size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ ) {
        if ( allow_ball_cut && arcs[ i0 ] ) {
            TF c = exp( ( weight - pow( sphere_radius, 2 ) ) / e.eps );
            f( _arc_length( point( i0 ), point( i1 ) ) * c, cut_ids[ i0 ] );
        } else {
            // Integrate[ Exp[ ( w - (bx+dx*u)*(bx+dx*u) - (by+dy*u)*(by+dy*u) ) / eps ], { u, 0, 1 } ]
            Pt P0 = point( i0 ) - sphere_center, P1 = point( i1 ) - sphere_center;
            if ( TF d2 = norm_2_p2( P1 - P0 ) ) {
                TF d1 = sqrt( d2 ), e5 = sqrt( e.eps );
                TF c = sqrt( M_PI ) * e5 / 2 * ( 1 ) * exp( ( weight - pow( P1.x * P0.y - P0.x * P1.y, 2 ) / d2 ) / e.eps ) / d1 * (
                            erf( ( P1.x * ( P1.x - P0.x ) + P1.y * ( P1.y - P0.y ) ) / e5 / d1 ) -
                            erf( ( P0.x * ( P1.x - P0.x ) + P0.y * ( P1.y - P0.y ) ) / e5 / d1 )
                            );
                f( d1 * c, cut_ids[ i0 ] );
            }
        }
    }
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_boundary_measure( FunctionEnum::R2, const std::function<void(TF,CI)> &f, TF weight ) const {
    TODO;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_boundary_measure( FunctionEnum::Unit, const std::function<void(TF,CI)> &f, TF weight ) const {
    if ( nb_points == 0 ) {
        if ( sphere_radius >= 0 )
            f( 2 * M_PI * sphere_radius, sphere_cut_id );
        return;
    }

    for( size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ ) {
        if ( allow_ball_cut && arcs[ i0 ] )
            f( _arc_length( point( i0 ), point( i1 ) ), cut_ids[ i0 ] );
        else
            f( norm_2( point( i1 ) - point( i0 ) ), cut_ids[ i0 ] );
    }
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::write_to_stream( std::ostream &os ) const {
    os << std::setprecision( 17 );
    os << "cuts: [";
    for( std::size_t i = 0; i < nb_points; ++i )
        os << ( i ? "," : "" ) << "(" << point( i ) << ")";
    os << "]";
    if ( store_the_normals ) {
        os << " nrms: [";
        for( std::size_t i = 0; i < nb_points; ++i )
            os << ( i ? "," : "" ) << "(" << normal( i ) << ")";
        os << "]";
    }
    os << " sphere center: " << sphere_center << " sphere radius: " << sphere_radius;
}

template<class Pc, class CI>
typename ConvexPolyhedron2<Pc,CI>::Pt ConvexPolyhedron2<Pc,CI>::min_position() const {
    Pt res{ + std::numeric_limits<TF>::max(), + std::numeric_limits<TF>::max() };
    for( std::size_t i = 0; i < nb_points; ++i ) {
        if ( arcs[ i ] )
            TODO;
        else
            res = min( res, point( i ) );
    }
    return res;
}

template<class Pc, class CI>
typename ConvexPolyhedron2<Pc,CI>::Pt ConvexPolyhedron2<Pc,CI>::max_position() const {
    Pt res{ - std::numeric_limits<TF>::max(), - std::numeric_limits<TF>::max() };
    for( std::size_t i = 0; i < nb_points; ++i ) {
        if ( arcs[ i ] )
            TODO;
        else
            res = max( res, point( i ) );
    }
    return res;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::set_cut_ids( CI cut_id ) {
    if ( nb_points > 64 )
        TODO;
    for( std::size_t i = 0; i < nb_points; ++i )
        cut_ids[ i ] = cut_id;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::ball_cut( Pt center, TF radius, CI cut_id ) {
    using std::sqrt;
    using std::pow;
    using std::max;
    using std::min;
    using std::abs;

    // store sphere info
    sphere_center = center;
    sphere_radius = radius;
    sphere_cut_id = cut_id;

    // void part ?
    if ( nb_points == 0 ) {
        sphere_radius = 0;
        return;
    }

    //
    if ( nb_points > 64 )
        TODO;

    // distances
    bool all_inside = true;
    alignas( 64 ) AF distances;
    for( std::size_t i = 0; i < nb_points; ++i ) {
        auto px = points[ 0 ][ i ];
        auto py = points[ 1 ][ i ];
        auto d = pow( px - center.x, 2 ) + pow( py - center.y, 2 ) - pow( radius, 2 );
        distances[ i ] = d;

        all_inside &= d <= 0;
    }

    // if all points (corners) are inside the sphere, the sphere is not going to cut anything
    if ( all_inside )
        return;

    auto find_unique_intersection = [&]( Pt p0, Pt p1 ) {
        // ( p0.x - sphere_center.x + ( p1.x - p0.x ) * t )² + ... = sphere_radius²
        TF a = norm_2_p2( p1 - p0 );
        TF b = dot( p0 - center, p1 - p0 );
        TF c = norm_2_p2( p0 - center ) - pow( radius, 2 );
        TF d = sqrt( max( TF( 0 ), b * b - a * c ) );
        TF u = ( - b + d ) / a;
        TF v = ( - b - d ) / a;
        TF t = abs( u - TF( 0.5 ) ) <= abs( v - TF( 0.5 ) ) ? u : v;
        return p0 + min( TF( 1 ), max( TF( 0 ), t ) ) * ( p1 - p0 );
    };

    auto find_two_cuts = [&]( Pt &pi0, Pt &pi1, const Pt &p0, const Pt &p1 ) {
        // ( p0.x - sphere_center.x + ( p1.x - p0.x ) * t )² + ... = sphere_radius²
        TF a = norm_2_p2( p1 - p0 );
        if ( a == 0 )
            return false;
        TF b = dot( p0 - center, p1 - p0 );
        TF c = norm_2_p2( p0 - center ) - pow( radius, 2 );
        TF s = b * b - a * c;
        if ( s <= 0 )
            return false;
        TF d = sqrt( s );
        TF u = ( - b - d ) / a;
        TF v = ( - b + d ) / a;
        if ( u > 0 && u < 1 )
            v = max( TF( 0 ), min( TF( 1 ), v ) );
        else if ( v > 0 && v < 1 )
            u = max( TF( 0 ), min( TF( 1 ), u ) );
        else
            return false;
        pi0 = p0 + u * ( p1 - p0 );
        pi1 = p0 + v * ( p1 - p0 );

        // check for different angles
        TF a0 = atan2( pi0.y - center.y, pi0.x - center.x );
        TF a1 = atan2( pi1.y - center.y, pi1.x - center.x );
        return a0 != a1;
    };

    auto cf = [&]( std::size_t i0, std::size_t i1 ) {
        const Pt &p0 = point( i0 );
        const Pt &p1 = point( i1 );

        if ( distances[ i0 ] <= 0 ) {
            if ( distances[ i1 ] <= 0 ) {
                _tmp_cuts.push_back( { LINE, cut_ids[ i0 ], normal( i0 )        , p0 } );
            } else {
                Pt nn = find_unique_intersection( p0, p1 );
                _tmp_cuts.push_back( { LINE, cut_ids[ i0 ], normal( i0 )        , p0 } );
                _tmp_cuts.push_back( { ARC ,        cut_id, { TF( 0 ), TF( 0 ) }, nn } );
            }
        } else {
            if ( distances[ i1 ] <= 0 ) {
                Pt nn = find_unique_intersection( p1, p0 );
                _tmp_cuts.push_back( { LINE, cut_ids[ i0 ], normal( i0 )        , nn } );
            } else {
                // 2 or 0 cuts
                Pt pi0, pi1;
                if ( find_two_cuts( pi0, pi1, p0, p1 ) ) {
                    _tmp_cuts.push_back( { LINE, cut_ids[ i0 ], normal( i0 )        , pi0 } );
                    _tmp_cuts.push_back( { ARC ,        cut_id, { TF( 0 ), TF( 0 ) }, pi1 } );
                }
            }
        }
    };
    _tmp_cuts.resize( 0 );
    for( std::size_t i = 1; i < nb_points; ++i )
        cf( i - 1, i );
    cf( nb_points - 1, 0 );

    // if no cut remaining, outside center => we don't keep anything
    if ( _tmp_cuts.empty() ) {
        for( std::size_t i0 = 0; i0 < nb_points; ++i0 ) {
            if ( dot( sphere_center - point( i0 ), normal( i0 ) ) > 0 ) {
                sphere_radius = 0;
                break;
            }
        }
    }

    // copy _tmp_cuts content
    nb_points = _tmp_cuts.size();
    if ( nb_points > 64 ) {
        TODO;
    } else {
        for( std::size_t i = 0; i < nb_points; ++i ) {
            for( std::size_t d = 0; d < 2; ++d ) {
                normals[ d ][ i ] = _tmp_cuts[ i ].normal[ d ];
                points[ d ][ i ] = _tmp_cuts[ i ].point[ d ];
            }
            cut_ids[ i ] = _tmp_cuts[ i ].cut_id;
            arcs.set( i, _tmp_cuts[ i ].cut_type );
        }
    }
}

template<class Pc,class CI>
bool ConvexPolyhedron2<Pc,CI>::plane_cut( Pt origin, Pt normal, CI cut_id ) {
    return plane_cut( origin, normal, cut_id, N<1>() );
}

template<class Pc,class CI> template<int no>
bool ConvexPolyhedron2<Pc,CI>::plane_cut( Pt origin, Pt normal, CI cut_id, N<no> normal_is_normalized ) {
    if ( nb_points > 64 )
        TODO;

    //    auto ox = BF( origin.x );
    //    auto oy = BF( origin.y );
    //    auto nx = BF( normal.x );
    //    auto ny = BF( normal.y );

    alignas( 64 ) AF distances;
    std::uint64_t outside = 0;
    const std::size_t alig_nb_points = 0; // nb_points - nb_points % simd_size;
    //    for( std::size_t i = 0; i < alig_nb_points; i += simd_size ) {
    //        auto px = xsimd::load_aligned( &points[ 0 ][ i ] );
    //        auto py = xsimd::load_aligned( &points[ 1 ][ i ] );
    //        auto d = ( ox - px ) * nx + ( oy - py ) * ny;
    //        auto n = d < BF( TF( 0 ) );
    //        d.store_aligned( &distances[ i ] );

    //        for( std::size_t j = 0; j < simd_size; ++j )
    //            outside |= std::uint64_t( n[ j ] ) << ( i + j );
    //    }

    for( std::size_t i = alig_nb_points; i < nb_points; ++i ) {
        auto px = points[ 0 ][ i ];
        auto py = points[ 1 ][ i ];
        auto d = ( origin.x - px ) * normal.x + ( origin.y - py ) * normal.y;
        outside |= std::uint64_t( d < 0 ) << i;
        distances[ i ] = d;
    }

    // all inside ?
    if ( outside == 0 )
        return false;

    // all outside ?
    std::size_t nb_outside = _mm_popcnt_u64( outside );
    if ( nb_outside == nb_points ) {
        nb_points = 0;
        return false;
    }

    //
    if ( normal_is_normalized.val == 0 ) {
        TF n = 1 / norm_2( normal );
        for( std::size_t i = 0; i < nb_points; ++i )
            distances[ i ] *= n;
        normal = n * normal;
    }

    // only 1 outside
    std::size_t i1 = _tzcnt_u64( outside );
    if ( nb_outside == 1 ) {
        if ( nb_points == 64 )
            TODO;

        // => creation of a new point
        std::size_t i0 = ( i1 + nb_points - 1 ) % nb_points;
        std::size_t i2 = ( i1 + 1 ) % nb_points;
        Pt p0 { points[ 0 ][ i0 ], points[ 1 ][ i0 ] };
        Pt p1 { points[ 0 ][ i1 ], points[ 1 ][ i1 ] };
        Pt p2 { points[ 0 ][ i2 ], points[ 1 ][ i2 ] };
        TF s0 = distances[ i0 ];
        TF s1 = distances[ i1 ];
        TF s2 = distances[ i2 ];

        TF m0 = s0 / ( s1 - s0 );
        TF m1 = s2 / ( s1 - s2 );

        // shift points
        for( std::size_t i = nb_points; i > i1 + 1; --i ) {
            points [ 0 ][ i ] = points [ 0 ][ i - 1 ];
            points [ 1 ][ i ] = points [ 1 ][ i - 1 ];
            cut_ids     [ i ] = cut_ids     [ i - 1 ];
            if ( store_the_normals ) {
                normals[ 0 ][ i ] = normals[ 0 ][ i - 1 ];
                normals[ 1 ][ i ] = normals[ 1 ][ i - 1 ];
            }
        }

        // modify the number of points
        ++nb_points;

        // modified or added points
        points[ 0 ][ i1 + 1 ] = p2.x - m1 * ( p1.x - p2.x );
        points[ 1 ][ i1 + 1 ] = p2.y - m1 * ( p1.y - p2.y );
        cut_ids    [ i1 + 1 ] = cut_ids[ i1 + 0 ];
        if ( store_the_normals ) {
            normals[ 0 ][ i1 + 1 ] = normals[ 0 ][ i1 + 0 ];
            normals[ 1 ][ i1 + 1 ] = normals[ 1 ][ i1 + 0 ];
        }

        points[ 0 ][ i1 + 0 ] = p0.x - m0 * ( p1.x - p0.x );
        points[ 1 ][ i1 + 0 ] = p0.y - m0 * ( p1.y - p0.y );
        cut_ids    [ i1 + 0 ] = cut_id;
        if ( store_the_normals ) {
            normals[ 0 ][ i1 + 0 ] = normal[ 0 ];
            normals[ 1 ][ i1 + 0 ] = normal[ 1 ];
        }

        return true;
    }

    // 2 points are outside
    if ( nb_outside == 2 ) {
        if ( i1 == 0 && ! ( outside & 2 ) )
            i1 = nb_points - 1;
        std::size_t i0 = ( i1 + nb_points - 1 ) % nb_points;
        std::size_t i2 = ( i1 + 1 ) % nb_points;
        std::size_t i3 = ( i1 + 2 ) % nb_points;
        Pt p0 { points[ 0 ][ i0 ], points[ 1 ][ i0 ] };
        Pt p1 { points[ 0 ][ i1 ], points[ 1 ][ i1 ] };
        Pt p2 { points[ 0 ][ i2 ], points[ 1 ][ i2 ] };
        Pt p3 { points[ 0 ][ i3 ], points[ 1 ][ i3 ] };
        TF s0 = distances[ i0 ];
        TF s1 = distances[ i1 ];
        TF s2 = distances[ i2 ];
        TF s3 = distances[ i3 ];

        TF m1 = s0 / ( s1 - s0 );
        TF m2 = s3 / ( s2 - s3 );

        // modified points
        points[ 0 ][ i1 ] = p0.x - m1 * ( p1.x - p0.x );
        points[ 1 ][ i1 ] = p0.y - m1 * ( p1.y - p0.y );
        cut_ids    [ i1 ] = cut_id;
        if ( store_the_normals ) {
            normals[ 0 ][ i1 ] = normal[ 0 ];
            normals[ 1 ][ i1 ] = normal[ 1 ];
        }

        points[ 0 ][ i2 ] = p3.x - m2 * ( p2.x - p3.x );
        points[ 1 ][ i2 ] = p3.y - m2 * ( p2.y - p3.y );

        return true;
    }

    // more than 2 points are outside, outside points are before and after bit 0
    if ( i1 == 0 && ( outside & ( 1ul << ( nb_points - 1 ) ) ) ) {
        std::size_t nb_inside = nb_points - nb_outside;
        std::size_t i3 = _tzcnt_u64( ~ outside );
        i1 = nb_inside + i3;
        std::size_t i0 = i1 - 1;
        std::size_t i2 = i3 - 1;

        Pt p0 { points[ 0 ][ i0 ], points[ 1 ][ i0 ] };
        Pt p1 { points[ 0 ][ i1 ], points[ 1 ][ i1 ] };
        Pt p2 { points[ 0 ][ i2 ], points[ 1 ][ i2 ] };
        Pt p3 { points[ 0 ][ i3 ], points[ 1 ][ i3 ] };
        TF s0 = distances[ i0 ];
        TF s1 = distances[ i1 ];
        TF s2 = distances[ i2 ];
        TF s3 = distances[ i3 ];

        TF m1 = s0 / ( s1 - s0 );
        TF m2 = s3 / ( s2 - s3 );

        // modified and shifted points
        points[ 0 ][ 0 ] = p3.x - m2 * ( p2.x - p3.x );
        points[ 1 ][ 0 ] = p3.y - m2 * ( p2.y - p3.y );
        cut_ids    [ 0 ] = cut_ids[ i2 ];
        if ( store_the_normals ) {
            normals[ 0 ][ 0 ] = normals[ 0 ][ i2 ];
            normals[ 1 ][ 0 ] = normals[ 1 ][ i2 ];
        }
        std::size_t o = 1;
        for( ; o <= nb_inside; ++o ) {
            points[ 0 ][ o ] = points[ 0 ][ i2 + o ];
            points[ 1 ][ o ] = points[ 1 ][ i2 + o ];
            cut_ids    [ o ] = cut_ids    [ i2 + o ];
            if ( store_the_normals ) {
                normals[ 0 ][ o ] = normals[ 0 ][ i2 + o ];
                normals[ 1 ][ o ] = normals[ 1 ][ i2 + o ];
            }
        }
        points[ 0 ][ o ] = p0.x - m1 * ( p1.x - p0.x );
        points[ 1 ][ o ] = p0.y - m1 * ( p1.y - p0.y );
        cut_ids    [ o ] = cut_id;
        if ( store_the_normals ) {
            normals[ 0 ][ o ] = normal[ 0 ];
            normals[ 1 ][ o ] = normal[ 1 ];
        }

        nb_points -= nb_outside - 2;
        return true;
    }

    // more than 2 points are outside, outside points do not cross `nb_points`
    std::size_t i0 = ( i1 + nb_points  - 1 ) % nb_points;
    std::size_t i2 = ( i1 + nb_outside - 1 ) % nb_points;
    std::size_t i3 = ( i1 + nb_outside     ) % nb_points;
    Pt p0 { points[ 0 ][ i0 ], points[ 1 ][ i0 ] };
    Pt p1 { points[ 0 ][ i1 ], points[ 1 ][ i1 ] };
    Pt p2 { points[ 0 ][ i2 ], points[ 1 ][ i2 ] };
    Pt p3 { points[ 0 ][ i3 ], points[ 1 ][ i3 ] };
    TF s0 = distances[ i0 ];
    TF s1 = distances[ i1 ];
    TF s2 = distances[ i2 ];
    TF s3 = distances[ i3 ];

    TF m1 = s0 / ( s1 - s0 );
    TF m2 = s3 / ( s2 - s3 );

    // modified and deleted points
    points[ 0 ][ i1 + 0 ] = p0.x - m1 * ( p1.x - p0.x );
    points[ 1 ][ i1 + 0 ] = p0.y - m1 * ( p1.y - p0.y );
    cut_ids    [ i1 + 0 ] = cut_id;
    if ( store_the_normals ) {
        normals[ 0 ][ i1 + 0 ] = normal[ 0 ];
        normals[ 1 ][ i1 + 0 ] = normal[ 1 ];
    }

    points[ 0 ][ i1 + 1 ] = p3.x - m2 * ( p2.x - p3.x );
    points[ 1 ][ i1 + 1 ] = p3.y - m2 * ( p2.y - p3.y );
    cut_ids    [ i1 + 1 ] = cut_ids[ i2 ];
    if ( store_the_normals ) {
        normals[ 0 ][ i1 + 1 ] = normals[ 0 ][ i2 ];
        normals[ 1 ][ i1 + 1 ] = normals[ 1 ][ i2 ];
    }

    std::size_t nb_to_rem = nb_outside - 2;
    for( std::size_t i = i2 + 1; i < nb_points; ++i ) {
        points[ 0 ][ i - nb_to_rem ] = points[ 0 ][ i ];
        points[ 1 ][ i - nb_to_rem ] = points[ 1 ][ i ];
        cut_ids    [ i - nb_to_rem ] = cut_ids    [ i ];
        if ( store_the_normals ) {
            normals[ 0 ][ i - nb_to_rem ] = normals[ 0 ][ i ];
            normals[ 1 ][ i - nb_to_rem ] = normals[ 1 ][ i ];
        }
    }

    // modification of the number of points
    nb_points -= nb_to_rem;

    //    //    bool all_ko = true;
    //    //    bool all_ok = true;
    //    //    for( std::size_t i = 0; i < _cuts.size(); ++i ) {
    //    //        _distances[ i ] = dot( _cuts[ i ].point - origin, normal );
    //    //        bool inside = _distances[ i ] <= 0;
    //    //        all_ko &= ! inside;
    //    //        all_ok &= inside;
    //    //    }
    //    bool all_ko_0 = true, all_ko_1 = true, all_ko_2 = true, all_ko_3 = true;
    //    bool all_ok_0 = true, all_ok_1 = true, all_ok_2 = true, all_ok_3 = true;
    //    for( int i = 0; i < (int)_cuts.size() - 3; i += 4 ) {
    //        _distances[ i + 0 ] = dot( _cuts[ i + 0 ].point - origin, normal );
    //        _distances[ i + 1 ] = dot( _cuts[ i + 1 ].point - origin, normal );
    //        _distances[ i + 2 ] = dot( _cuts[ i + 2 ].point - origin, normal );
    //        _distances[ i + 3 ] = dot( _cuts[ i + 3 ].point - origin, normal );
    //        bool inside_0 = _distances[ i + 0 ] <= 0;
    //        bool inside_1 = _distances[ i + 1 ] <= 0;
    //        bool inside_2 = _distances[ i + 2 ] <= 0;
    //        bool inside_3 = _distances[ i + 3 ] <= 0;
    //        all_ko_0 &= ! inside_0;
    //        all_ko_1 &= ! inside_1;
    //        all_ko_2 &= ! inside_2;
    //        all_ko_3 &= ! inside_3;
    //        all_ok_0 &= inside_0;
    //        all_ok_1 &= inside_1;
    //        all_ok_2 &= inside_2;
    //        all_ok_3 &= inside_3;
    //    }
    //    for( std::size_t i = _cuts.size() / 4 * 4; i < _cuts.size(); ++i ) {
    //        _distances[ i ] = dot( _cuts[ i ].point - origin, normal );
    //        bool inside = _distances[ i ] <= 0;
    //        all_ko_0 &= ! inside;
    //        all_ok_0 &= inside;
    //    }
    //    all_ko_0 &= all_ko_1;
    //    all_ko_0 &= all_ko_2;
    //    all_ko_0 &= all_ko_3;
    //    all_ok_0 &= all_ok_1;
    //    all_ok_0 &= all_ok_2;
    //    all_ok_0 &= all_ok_3;

    //    if ( all_ok_0 )
    //        return;

    //    if ( all_ko_0 ) {
    //        _cuts.resize( 0 );
    //        return;
    //    }

    //    auto cf = [&]( std::size_t i0, std::size_t i1 ) {
    //        const PT &p0 = _cuts[ i0 ].point;
    //        const PT &p1 = _cuts[ i1 ].point;
    //        TF s0 = _distances[ i0 ];
    //        TF s1 = _distances[ i1 ];
    //        if ( s0 <= 0 ) { // => p0 is inside
    //            if ( s1 <= 0 ) { // => p0 and p1 are inside
    //                _old_cuts.push_back( { SegType::line, _cuts[ i0 ].cut_id, _cuts[ i0 ].normal, p0 } );
    //            } else { // => p0 is inside, p1 is outside
    //                PT it = p0 - ( s0 / ( s1 - s0 ) ) * ( p1 - p0 );
    //                _old_cuts.push_back( { SegType::line, _cuts[ i0 ].cut_id, _cuts[ i0 ].normal, p0 } );
    //                _old_cuts.push_back( { SegType::line,             cut_id,             normal, it } );
    //            }
    //        } else if ( s1 <= 0 ) { // => p0 is outside, p1 is inside
    //            PT it = p0 + ( s0 / ( s0 - s1 ) ) * ( p1 - p0 );
    //            _old_cuts.push_back( { SegType::line, _cuts[ i0 ].cut_id, _cuts[ i0 ].normal, it } );
    //        }
    //    };
    //    _old_cuts.resize( 0 );
    //    for( size_t i = 1; i < _cuts.size(); ++i )
    //        cf( i - 1, i );
    //    cf( _cuts.size() - 1, 0 );

    //    std::swap( _old_cuts, _cuts );

    return true;
}


template<class Pc,class CI>
bool ConvexPolyhedron2<Pc,CI>::is_a_cutting_plane( Pt origin, Pt normal ) const {
    TODO;
    //    for( std::size_t i = 0; i < _cuts.size(); ++i ) {
    //        bool inside = dot( _cuts[ i ].point - origin, normal ) <= 0;
    //        if ( ! inside )
    //            return true;
    //    }
    return false;
}

template<class Pc, class CI>
bool ConvexPolyhedron2<Pc,CI>::contains( const Pt &pos ) const {
    for( TI i = 0; i < nb_points; ++i )
        if ( dot( pos - point( i ), normal( i ) ) > 0 )
            return false;
    return true;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::clear( Pt englobing_center, TF englobing_radius, CI englobing_cut_id ) {
    //    _sphere_center = englobing_center;
    //    _sphere_radius = -1;

    //    // common cut attributes
    //    _cuts.resize( 3 );
    //    for( Cut &c : _cuts ) {
    //        c.seg_type = SegType::line;
    //        c.cut_id   = englobing_cut_id;
    //    }

    //    // points
    //    const TF s3 = std::sqrt( TF( 3 ) / 4 );
    //    _cuts[ 0 ].point = { englobing_center + 4 * englobing_radius * PT{    1,   0 } };
    //    _cuts[ 1 ].point = { englobing_center + 4 * englobing_radius * PT{ -0.5,  s3 } };
    //    _cuts[ 2 ].point = { englobing_center + 4 * englobing_radius * PT{ -0.5, -s3 } };
    TODO;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::add_centroid_contrib( Pt &ctd, TF &mea ) const {
    add_centroid_contrib( ctd, mea, FunctionEnum::Unit(), SpaceFunctions::Constant<TF>{ 1.0 } );
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::add_centroid_contrib( Pt &ctd, TF &mea, FunctionEnum::Unit, SpaceFunctions::Constant<TF> sf, TF w ) const {
    using std::pow;
    //    auto arc_val = [&]( PT P0, PT P1 ) {
    //        using std::atan2;
    //        using std::pow;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        return ( a1 - a0 ) * 0.5 * pow( _sphere_radius, 2 );
    //    };

    //    auto seg_val = []( PT P0, PT P1 ) {
    //        return -0.25 * ( ( P0.x + P1.x ) * ( P0.y - P1.y ) - ( P0.x - P1.x ) * ( P0.y + P1.y ) );
    //    };

    //    TF lea = 0;
    //    if ( _cuts.empty() ) {
    //        if ( _sphere_radius > 0 )
    //            lea = 2 * M_PI * 0.5 * pow( _sphere_radius, 2 );
    //    } else {
    //        for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //            if ( _cuts[ i0 ].seg_type == SegType::arc )
    //                lea += arc_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //            else
    //                lea += seg_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        }
    //    }

    //    // centroid part
    //    ctd += lea * _sphere_center;
    //    mea += lea;

    //    if ( _cuts.empty() )
    //        return;

    //    auto arc_val_centroid = [&]( PT P0, PT P1 ) -> PT {
    //        using std::atan2;
    //        using std::pow;
    //        using std::sin;
    //        using std::cos;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        TF c = ( 0.5 ) * pow( _sphere_radius, 3 );
    //        return { c * ( sin( a1 ) - sin( a0 ) ) ,
    //                 c * ( cos( a0 ) - cos( a1 ) ) };
    //    };

    //    auto seg_val_centroid = []( PT P0, PT P1 ) -> PT {
    //        TF c = ( 0.125 ) * ( pow( P0.x + P1.x, 2 ) + pow( P0.y + P1.y, 2 ) ) + ( 0.04166666666666667 ) * ( pow( P1.x - P0.x, 2 ) + pow( P1.y - P0.y, 2 ) );
    //        return { ( P1.y - P0.y ) * c, ( P0.x - P1.x ) * c };
    //    };

    //    for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //        if ( _cuts[ i0 ].seg_type == SegType::arc )
    //            ctd += arc_val_centroid( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        else
    //            ctd += seg_val_centroid( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //    }

    // hand coded version
    if ( nb_points == 0 ) {
        if ( sphere_radius >= 0 ) {
            TF lea = M_PI * pow( sphere_radius, 2 ) * sf.coeff;
            ctd += lea * sphere_center;
            mea += lea;
        }
        return;
    }

    // triangles
    Pt A = point( 0 );
    for( size_t i = 2; i < nb_points; ++i ) {
        Pt B = point( i - 1 );
        Pt C = point( i - 0 );
        TF lea = 0.5 * ( A.x * ( B.y - C.y ) + B.x * ( C.y - A.y ) + C.x * ( A.y - B.y ) ) * sf.coeff;
        ctd += TF( TF( 1 ) / 3 * lea ) * ( A + B + C );
        mea += lea;
    }

    // arcs
    if ( allow_ball_cut )
        for( std::size_t i0 = nb_points - 1, i1 = 0; i1 < nb_points; i0 = i1++ )
            if ( arcs[ i0 ] )
                _centroid_arc( ctd, mea, point( i0 ), point( i1 ), sf.coeff );
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::add_centroid_contrib( Pt &ctd, TF &mea, FunctionEnum::ExpWmR2db<TF> func, SpaceFunctions::Constant<TF> sf, TF w ) const {
    using std::sqrt;

    // generated using nsmake run -g3 src/PowerDiagram/offline_integration/gen_approx_integration.cpp --function Gaussian --end-log-scale 10 --precision 1e-10 --centroid
    static const std::vector<std::pair<TF,std::array<TF,8>>> coeffs_centroid = {
        {             0.81, {  3.436754914e-11,     0.9999999936,     -0.499999808,     0.1666644663,   -0.04165429959,   0.008295504317,  -0.001324187054,  0.0001391024954 } },
        {           1.7956, {  3.626873334e-05,     0.9997460436,    -0.4992250072,     0.1653164327,   -0.04018994106,   0.007284377258,  -0.000905296835,  5.768165611e-05 } },
        {             2.89, {   0.002150113161,     0.9919358659,    -0.4866516465,     0.1538545005,   -0.03378503863,   0.005086512678, -0.0004760121107,  2.086282073e-05 } },
        {           4.2436, {    0.02578955805,     0.9363593867,    -0.4300017846,     0.1213726428,   -0.02246285614,   0.002686320036, -0.0001894501756,  6.001223393e-06 } },
        {           5.8564, {     0.1198526666,     0.7808092003,    -0.3189603936,    0.07700560426,   -0.01174582303,   0.001121234454, -6.151360776e-05,  1.486138102e-06 } },
        {           7.9524, {     0.3626890675,     0.4883881014,    -0.1669509705,     0.0327852988,  -0.003971431676,  0.0002952961456, -1.242844854e-05,  2.276530582e-07 } },
        {          10.3684, {     0.6658026798,     0.2144627846,   -0.06029387531,   0.009592922335, -0.0009300842507,  5.481648971e-05, -1.814679753e-06,  2.598625563e-08 } },
        {          13.9876, {     0.9052560174,    0.04867379694,   -0.01083135964,   0.001351194596, -0.0001019133051,  4.642173912e-06, -1.181253207e-07,  1.294296733e-09 } },
        {      20.16554128, {     0.9926751318,   0.002820878153, -0.0004672646601,  4.312942691e-05, -2.394479485e-06,  7.992405755e-08, -1.484493634e-09,  1.183209069e-11 } },
        {      48.13399542, {     0.9999975027,  5.298270207e-07,  -4.77053781e-08,  2.360949481e-09, -6.936850624e-11,   1.21024003e-12, -1.161161656e-14,  4.727626768e-17 } },
        {            1e+40, {     0.9999999977,  2.275142176e-18, -9.619261251e-20,  2.236212338e-21, -3.087317629e-23,  2.531671801e-25, -1.141960256e-27,  2.186301388e-30 } },
    };

    Pt tmp{ 0, 0 };
    TF eco = sqrt( func.eps );
    _r_centroid_integration( tmp.x, tmp.y, coeffs_centroid, TF( 1 ) / eco );
    ctd += pow( eco, 3 ) * sf.coeff * tmp;

    static const std::vector<std::pair<TF,std::array<TF,8>>> coeffs = {
        {           1.2996, {     0.4999999997,    -0.2499999827,    0.08333305386,   -0.02083158856,    0.00416121191, -0.0006849636014,  8.978963538e-05, -7.240855485e-06 } },
        {           2.7556, {     0.4999495695,    -0.2497751362,     0.0828948459,   -0.02034278382,   0.003819615811, -0.0005330234995,  4.938355801e-05, -2.249060705e-06 } },
        {             4.41, {     0.4981089084,    -0.2453203749,     0.0781960602,   -0.01753650964,   0.002792798345, -0.0003025475296,  1.998609405e-05, -6.059544385e-07 } },
        {           6.4516, {     0.4846982526,    -0.2243732282,    0.06401642919,   -0.01214114936,   0.001546076695, -0.0001275781801,  6.178598559e-06, -1.335000282e-07 } },
        {           8.8804, {     0.4432202492,    -0.1791255579,    0.04268794966,  -0.006509163724,  0.0006463252278, -4.062413504e-05,  1.472747438e-06, -2.351860255e-08 } },
        {          11.9716, {     0.3718824889,    -0.1222831443,     0.0231493835,  -0.002753606604,   0.000210412058, -1.007502188e-05,  2.761476809e-07, -3.315113427e-09 } },
        {          15.8404, {     0.2911023832,   -0.07417541322,    0.01080055983, -0.0009826795077,  5.719712218e-05, -2.079557094e-06,  4.317540828e-08, -3.918779009e-10 } },
        {      20.92609863, {     0.2218236671,    -0.0430074657,   0.004759441426, -0.0003288198859,  1.452268876e-05, -4.004246974e-07,  6.301737011e-09, -4.333919868e-11 } },
        {      27.62333749, {     0.1680993965,   -0.02469616128,    0.00207082746, -0.0001083996797,  3.627290314e-06, -7.577183884e-08,  9.034157059e-10, -4.706946611e-12 } },
        {      36.46397674, {     0.1273452673,   -0.01417301708,  0.0009003110811, -3.570197005e-05,  9.050274063e-07, -1.432195307e-08,  1.293590653e-10, -5.105792434e-13 } },
        {      49.03329792, {    0.09558605859,  -0.007983779363,  0.0003805362908, -1.132076813e-05,  2.152528359e-07, -2.554561589e-09,    1.7300631e-11, -5.119218133e-14 } },
        {      65.93532905, {    0.07108264816,   -0.00441516283,  0.0001564961074, -3.462204518e-06,  4.895478776e-08, -4.320482338e-10,  2.175942741e-12, -4.788052604e-15 } },
        {      88.66357763, {    0.05286114006,  -0.002441701544,   6.43609312e-05, -1.058874342e-06,  1.113422016e-08, -7.307516305e-11,  2.736893816e-13, -4.478600509e-16 } },
        {            1e+40, {    0.04232218964,  -0.001566869156,  3.313963314e-05, -4.379566085e-07,  3.703248443e-09, -1.956604449e-11,  5.905729983e-14, -7.796696243e-17 } },
    };

    TF lea = func.eps * sf.coeff * _r_polynomials_integration( coeffs, TF( 1 ) / eco );
    ctd += lea * sphere_center;
    mea += lea;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::add_centroid_contrib( Pt &ctd, TF &mea, FunctionEnum::R2, SpaceFunctions::Constant<TF> sf, TF w ) const {
    TODO;
    //    // generated using nsmake run -g3 src/PowerDiagram/offline_integration/gen_approx_integration.cpp --function R2 --end-log-scale 1000 --precision 1e-10 --centroid
    //    auto arc_val = [&]( PT P0, PT P1 ) {
    //        using std::atan2;
    //        using std::pow;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        return ( a1 - a0 ) * 0.25 * pow( _sphere_radius, 4 );
    //    };

    //    auto seg_val = []( PT P0, PT P1 ) {
    //        TF result;
    //        TF R_0 = P1.y; TF R_1 = (-1.0)*R_0; TF R_2 = P0.y; TF R_3 = (-1.0)*R_2;
    //        R_3 = R_0+R_3; TF R_4 = pow(R_3,2); R_1 = R_2+R_1; R_0 = R_2+R_0;
    //        R_2 = R_3*R_0; TF R_5 = pow(R_0,2); TF R_6 = P0.x; TF R_7 = (-1.0)*R_6;
    //        TF R_8 = P1.x; R_7 = R_8+R_7; TF R_9 = R_1*R_7; R_9 = (-1.0)*R_9;
    //        TF R_10 = pow(R_7,2); R_4 = R_10+R_4; R_10 = (-1.0)*R_8; R_10 = R_6+R_10;
    //        R_3 = R_10*R_3; R_3 = R_9+R_3; R_0 = R_10*R_0; R_0 = (-1.0)*R_0;
    //        R_6 = R_8+R_6; R_7 = R_7*R_6; R_7 = R_2+R_7; R_7 = R_3*R_7;
    //        R_7 = 2.0*R_7; R_1 = R_1*R_6; R_0 = R_1+R_0; R_4 = R_4*R_0;
    //        R_4 = (-1.0)*R_4; R_7 = R_4+R_7; R_6 = pow(R_6,2); R_5 = R_6+R_5;
    //        R_6 = 0.25; R_7 = R_6*R_7; R_7 = (1.0/24.0)*R_7; R_5 = R_6*R_5;
    //        R_0 = R_5*R_0; R_0 = -0.125*R_0; R_7 = R_0+R_7; result = R_7;
    //        return result;
    //    };

    //    TF lea = 0;
    //    if ( _cuts.empty() ) {
    //        if ( _sphere_radius > 0 )
    //            lea = 2 * M_PI * 0.25 * pow( _sphere_radius, 4 );
    //    } else {
    //        for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //            if ( _cuts[ i0 ].seg_type == SegType::arc )
    //                lea += arc_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //            else
    //                lea += seg_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        }
    //    }

    //    // centroid part
    //    ctd += lea * _sphere_center;
    //    mea += lea;

    //    if ( _cuts.empty() )
    //        return;

    //    auto arc_val_centroid = [&]( TF &r_x, TF &r_y, PT P0, PT P1 ) {
    //        using std::atan2;
    //        using std::pow;
    //        using std::sin;
    //        using std::cos;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        TF c = ( 0.25 ) * pow( _sphere_radius, 5 );
    //        r_x += c * ( sin( a1 ) - sin( a0 ) );
    //        r_y += c * ( cos( a0 ) - cos( a1 ) );
    //    };

    //    auto seg_val_centroid = []( TF &r_x, TF &r_y, PT P0, PT P1 ) {
    //        TF R_0 = P0.x; TF R_1 = (-1.0)*R_0; TF R_2 = P1.x; TF R_3 = (-1.0)*R_2;
    //        R_3 = R_0+R_3; R_0 = R_2+R_0; TF R_4 = pow(R_0,2); R_1 = R_2+R_1;
    //        R_0 = R_1*R_0; R_1 = pow(R_1,2); R_2 = 0.5; TF R_5 = P1.y;
    //        TF R_6 = (-1.0)*R_5; TF R_7 = P0.y; TF R_8 = R_5+R_7; TF R_9 = pow(R_8,2);
    //        R_9 = R_4+R_9; R_4 = pow(R_9,2); R_4 = R_2*R_4; TF R_10 = R_3*R_4;
    //        R_10 = (1.0/32.0)*R_10; TF R_11 = (-1.0)*R_7; R_11 = R_5+R_11; R_8 = R_11*R_8;
    //        R_0 = R_8+R_0; R_0 = pow(R_0,2); R_0 = 0.5*R_0; R_11 = pow(R_11,2);
    //        R_11 = R_1+R_11; R_9 = R_11*R_9; R_9 = 0.25*R_9; R_9 = R_0+R_9;
    //        R_9 = R_2*R_9; R_0 = R_3*R_9; R_0 = (1.0/12.0)*R_0; R_0 = R_10+R_0;
    //        R_11 = pow(R_11,2); R_11 = R_2*R_11; R_3 = R_3*R_11; R_3 = (1.0/160.0)*R_3;
    //        R_0 = R_3+R_0; r_y += R_0; R_6 = R_7+R_6; R_9 = R_6*R_9;
    //        R_9 = (-1.0/12.0)*R_9; R_4 = R_6*R_4; R_4 = (-1.0/32.0)*R_4; R_9 = R_4+R_9;
    //        R_11 = R_6*R_11; R_11 = (-1.0/160.0)*R_11; R_9 = R_11+R_9; r_x += R_9;
    //    };

    //    for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //        if ( _cuts[ i0 ].seg_type == SegType::arc )
    //            arc_val_centroid( ctd.x, ctd.y, _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        else
    //            seg_val_centroid( ctd.x, ctd.y, _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //    }
}

template<class Pc,class CI> template<class FU>
typename ConvexPolyhedron2<Pc,CI>::Pt ConvexPolyhedron2<Pc,CI>::centroid( const FU &f, TF w ) const {
    TF mea = 0;
    Pt ctd = { 0, 0 };
    add_centroid_contrib( ctd, mea, f, SpaceFunctions::Constant<TF>{ 1.0 }, w );
    return mea ? ctd / mea : ctd;
}

template<class Pc,class CI>
typename ConvexPolyhedron2<Pc,CI>::Pt ConvexPolyhedron2<Pc,CI>::centroid() const {
    return centroid( FunctionEnum::Unit() );
}

template<class Pc,class CI> template<class FU>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::measure( const FU &f, TF w ) const {
    return integration( f, w );
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::measure() const {
    return measure( FunctionEnum::Unit() );
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::integration( FunctionEnum::ExpWmR2db<TF> func, TF w ) const {
    using std::sqrt;
    using std::exp;

    // nsmake run -g3 src/PowerDiagram/offline_integration/gen_approx_integration.cpp --function Gaussian --end-log-scale 100 --precision 1e-10 -r 5000 -l 5000
    //    static const std::vector<std::pair<TF,std::array<TF,8>>> coeffs = {
    //        {       1.05925264, {     0.4999999999,    -0.2499999926,     0.0833332093,   -0.02083246334,   0.004163530444, -0.0006881183554,  9.192314036e-05, -7.809062357e-06 } },
    //        {       2.20344336, {     0.4999850102,    -0.2499200703,    0.08314729201,   -0.02058558165,   0.003958891852, -0.0005806804802,  5.839184767e-05, -2.974827042e-06 } },
    //        {       3.47598736, {     0.4994170562,    -0.2482442213,    0.08099722387,   -0.01902771108,   0.003269219526, -0.0003938994023,  2.972032793e-05, -1.050631798e-06 } },
    //        {         4.919524, {     0.4949170569,    -0.2394717823,     0.0736012051,   -0.01553017726,   0.002266951228, -0.0002198238933,  1.275382813e-05,  -3.34890575e-07 } },
    //        {       6.58846224, {     0.4785550305,    -0.2165548465,    0.05975708912,    -0.0108535712,   0.001312821299, -0.0001022552306,  4.653156374e-06, -9.416582146e-08 } },
    //        {         8.561476, {     0.4421210697,    -0.1781251216,    0.04229847109,  -0.006425093074,  0.0006354588594, -3.978314345e-05,  1.436662198e-06, -2.285641597e-08 } },
    //        {          10.9561, {     0.3852698768,    -0.1317380752,    0.02600847209,  -0.003233461183,  0.0002586879387,  -1.29863525e-05,  3.735955403e-07, -4.711709831e-09 } },
    //        {      13.93678224, {     0.3181270588,    -0.0887678095,    0.01417652918,  -0.001416472283,  9.063229761e-05, -3.625395206e-06,  8.287090771e-08, -8.286274423e-10 } },
    //        {      17.61111805, {     0.2554920574,   -0.05708876822,   0.007285129951, -0.0005806676738,  2.960058955e-05, -9.424080005e-07,  1.713217525e-08, -1.361529348e-10 } },
    //        {      22.52220514, {     0.2009844712,   -0.03531047147,   0.003541405849, -0.0002217667135,  8.879021478e-06, -2.219637233e-07,  3.167592646e-09, -1.975712889e-11 } },
    //        {       28.9143048, {     0.1568793508,   -0.02151233702,    0.00168391238, -8.229622634e-05,  2.571398844e-06, -5.016362552e-08,  5.586261996e-10, -2.718841774e-12 } },
    //        {      37.26425675, {     0.1219648659,   -0.01300203181,  0.0007911949871, -3.005869573e-05,  7.300820005e-07, -1.107105401e-08,  9.583080705e-11, -3.625246907e-13 } },
    //        {      48.27355954, {     0.0943931325,   -0.00778759521,  0.0003667271048, -1.078146646e-05,  2.026322207e-07, -2.377584721e-09,  1.592368025e-11, -4.660670231e-14 } },
    //        {       62.7775044, {    0.07272565045,  -0.004622566535,  0.0001677029801, -3.798210111e-06,  5.499186965e-08, -4.970502826e-10,  2.564293299e-12, -5.781202995e-15 } },
    //        {      81.95522263, {    0.05581578163,  -0.002722745924,  7.580614992e-05,  -1.31755136e-06,  1.463852101e-08, -1.015300233e-10,  4.019228325e-13, -6.952795953e-16 } },
    //        {      107.5440464, {    0.04264518344,  -0.001589325796,  3.380515859e-05, -4.488474221e-07,  3.809438887e-09, -2.018231566e-11,  6.102560026e-14, -8.063099891e-17 } },
    //        {      141.6687217, {    0.03243578947, -0.0009194043362,   1.48730617e-05, -1.501843271e-07,   9.69348859e-10, -3.905423247e-12,  8.979912011e-15, -9.022148232e-18 } },
    //        {      187.5852684, {     0.0245596726, -0.0005270865483,   6.45553521e-06,  -4.93505761e-08,  2.411364009e-10, -7.354366136e-13,   1.28003905e-15, -9.734524908e-19 } },
    //        {      249.6667076, {    0.01850050576, -0.0002990768417,  2.759000312e-06, -1.588580508e-08,  5.845965154e-11, -1.342746956e-13,  1.759976176e-16, -1.007888008e-19 } },
    //        {      334.0101702, {     0.0138646068, -0.0001679612533,  1.161071234e-06, -5.009283606e-09,  1.381212328e-11, -2.376921801e-14,  2.334117913e-17, -1.001388297e-20 } },
    //        {      449.1546486, {    0.01033699823, -9.336001117e-05,  4.811202109e-07, -1.547360118e-09,  3.180353863e-12, -4.079504339e-15,  2.985872281e-18, -9.547407599e-22 } },
    //        {      607.1126168, {   0.007667331741, -5.136155349e-05,   1.96307077e-07, -4.682271201e-10,   7.13676206e-13,   -6.7884667e-16,  3.684275478e-19, -8.734964342e-23 } },
    //        {      824.8590986, {   0.005657922505, -2.796667751e-05,  7.886892674e-08, -1.387941357e-10,  1.560765886e-13, -1.095236153e-16,   4.38495759e-20, -7.668837426e-24 } },
    //        {      1127.941873, {   0.004151014217, -1.505242016e-05,  3.113943016e-08, -4.019631457e-11,  3.315389133e-14, -1.706308445e-17,  5.010023785e-21, -6.425403112e-25 } },
    //        {      1552.351655, {   0.003025903906, -7.997947561e-06,   1.20593701e-08,  -1.13452079e-11,   6.81936027e-15, -2.557527916e-18,  5.471764213e-22, -5.113096525e-26 } },
    //        {       2150.25486, {    0.00219159211, -4.195246986e-06,  4.580868439e-09, -3.120683919e-12,  1.358204418e-15, -3.688045732e-19,  5.712525593e-23, -3.864384856e-27 } },
    //        {      2997.686458, {   0.001577131876, -2.172417423e-06,  1.706790614e-09, -8.365621434e-13,  2.619383088e-16, -5.116644973e-20,  5.700883164e-24, -2.773891843e-28 } },
    //        {      4206.093411, {   0.001127664676,  -1.11054403e-06,  6.237662941e-10, -2.185534827e-13,  4.891536228e-17, -6.829475317e-21,  5.438377121e-25, -1.891086748e-29 } },
    //        {      5947.402733, {   0.000800603188, -5.597214373e-07,  2.231610618e-10, -5.549786375e-14,  8.815514924e-18, -8.734453624e-22,  4.935443197e-26, -1.217697878e-30 } },
    //        {      8474.840042, {   0.000564025977, -2.777768998e-07,  7.800917092e-11, -1.366366998e-14,  1.528493657e-18,  -1.06644347e-22,  4.243029226e-27, -7.370538573e-32 } },
    //        {            1e+40, {  0.0004340438709, -1.647683473e-07,  3.572519365e-11, -4.838985987e-15,   4.19288882e-19,  -2.26961347e-23,  7.017013478e-28, -9.487047369e-33 } },
    //    };
    //    static const std::vector<std::pair<TF,std::array<TF,8>>> coeffs = {
    //        {       0.60155536, {              0.5,    -0.2499999998,    0.08333332835,   -0.02083327212,   0.004166280921, -0.0006930896085,  9.651114937e-05, -9.518601186e-06 } },
    //        {       1.21969936, {      0.499999749,    -0.2499976359,    0.08332363958,   -0.02081068869,   0.004133595891, -0.0006633222961,  8.050107042e-05, -5.559214563e-06 } },
    //        {         1.865956, {     0.4999879863,    -0.2499349876,    0.08317889967,   -0.02062231949,   0.003984212121, -0.0005910365678,    6.0720882e-05, -3.197148714e-06 } },
    //        {       2.55104784, {     0.4998680401,    -0.2495006372,    0.08249991871,   -0.02002810706,   0.003669670061, -0.0004903004491,  4.264492052e-05, -1.795206534e-06 } },
    //        {       3.27899664, {     0.4992759999,    -0.2479071398,     0.0806530737,   -0.01883312402,   0.003203418101, -0.0003805920641,  2.823015302e-05, -9.793549765e-07 } },
    //        {       4.05941904, {     0.4973729173,    -0.2438922166,    0.07701066726,   -0.01699097653,   0.002642461946, -0.0002777401628,  1.771655497e-05, -5.171523822e-07 } },
    //        {         4.901796, {     0.4927642042,    -0.2360053272,    0.07121125787,     -0.014615604,   0.002057157635, -0.0001909767755,   1.05523268e-05,  -2.62956983e-07 } },
    //        {       5.81581456, {      0.483695125,    -0.2231178933,    0.06334607779,   -0.01194320712,   0.001511185718,  -0.000123908626,  5.965511622e-06, -1.282335356e-07 } },
    //        {           6.8121, {      0.468542602,     -0.204935728,       0.05397902,  -0.009257459019,   0.001048318647, -7.596019975e-05,  3.201172453e-06, -5.981149575e-08 } },
    //        {       7.90959376, {     0.4463567784,    -0.1821810885,    0.04396135554,  -0.006803480729,  0.0006870703768, -4.400297459e-05,  1.628159529e-06, -2.657733041e-08 } },
    //        {       9.12764944, {     0.4172802317,    -0.1564714338,     0.0342049628,  -0.004743674949,  0.0004257763662, -2.408726943e-05,  7.836656263e-07, -1.120918546e-08 } },
    //        {      10.48982544, {     0.3827333897,    -0.1299798757,    0.02548721197,  -0.003147785924,  0.0002502575104, -1.248977523e-05,  3.573854005e-07, -4.485515696e-09 } },
    //        {      12.02979856, {     0.3450115879,    -0.1047957214,    0.01827223249,  -0.001997996098,  0.0001401800218, -6.158789834e-06,  1.548467307e-07,   -1.7051761e-09 } },
    //        {      13.77597456, {     0.3067614128,   -0.08251284335,    0.01270207258,  -0.001223489988,  7.548611198e-05, -2.912562315e-06,  6.424379464e-08, -6.201442241e-10 } },
    //        {      15.77360656, {     0.2702193106,   -0.06391748409,   0.008641681718, -0.0007303329548,  3.950480341e-05, -1.335538199e-06,  2.579879556e-08, -2.180105627e-10 } },
    //        {      18.11720762, {     0.2366997336,   -0.04901554503,   0.005798941141, -0.0004286965952,   2.02782116e-05, -5.993508952e-07,  1.012001991e-08, -7.473813834e-11 } },
    //        {      20.84748347, {     0.2059226936,   -0.03709244778,   0.003816769126, -0.0002453872538,   1.00937518e-05, -2.594160753e-07,  3.808605658e-09, -2.445544358e-11 } },
    //        {      24.02012719, {     0.1788729209,   -0.02798681236,   0.002501383316, -0.0001396828202,  4.990471392e-06, -1.113976424e-07,  1.420457476e-09, -7.921627976e-12 } },
    //        {      27.71125887, {     0.1551508658,   -0.02105565791,   0.001632291594, -7.906041436e-05,  2.449925374e-06, -4.743284537e-08,  5.245901004e-10, -2.537416911e-12 } },
    //        {      32.01079763, {     0.1343989775,    -0.0157997263,   0.001060998364, -4.451532257e-05,  1.194906762e-06, -2.003956686e-08,  1.919793459e-10, -8.043569066e-13 } },
    //        {      37.02508158, {     0.1162724094,    -0.0118251939,  0.0006869873085, -2.493538698e-05,  5.790425692e-07, -8.401019486e-09,  6.962459206e-11, -2.523589761e-13 } },
    //        {      42.88000628, {     0.1004611319,  -0.008827715631,  0.0004431030314,  -1.38958571e-05,  2.787983617e-07, -3.494777227e-09,  2.502391379e-11,   -7.8363475e-14 } },
    //        {      49.72478847, {    0.08668818832,  -0.006573085909,  0.0002846964583,  -7.70400061e-06,   1.33374275e-07, -1.442611862e-09,  8.913138413e-12, -2.408417574e-14 } },
    //        {      57.73648443, {    0.07470729156,  -0.004881717632,  0.0001822144906, -4.249246006e-06,  6.339566769e-08, -5.909165798e-10,  3.146253608e-12, -7.326195418e-15 } },
    //        {      67.12541998, {    0.06429934251,  -0.003616237423,   0.000116172981, -2.331681578e-06,  2.993982769e-08, -2.401847616e-10,  1.100626015e-12, -2.205710638e-15 } },
    //        {      78.14172309, {    0.05527016862,  -0.002671914933,  7.378183208e-05, -1.272885479e-06,  1.404890118e-08, -9.687432009e-11,  3.815665864e-13, -6.572697458e-16 } },
    //        {      91.08319137, {    0.04744776543,  -0.001969108869,  4.667838433e-05, -6.913088976e-07,  6.549961749e-09, -3.877178775e-11,  1.310947848e-13, -1.938491102e-16 } },
    //        {      106.3047774, {     0.0406800479,  -0.001447432232,  2.941739926e-05,  -3.73523289e-07,  3.034160253e-09, -1.539809263e-11,  4.463597328e-14, -5.658600671e-17 } },
    //        {      124.3901242, {    0.03481038457,  -0.001059856712,  1.843184988e-05, -2.002590994e-07,  1.391929399e-09, -6.044271214e-12,  1.499182698e-14, -1.626167819e-17 } },
    //        {      145.7398467, {    0.02973012014, -0.0007730719146,  1.148216071e-05,  -1.06543257e-07,  6.324511602e-10, -2.345453736e-12,  4.968298739e-15, -4.602415179e-18 } },
    //        {      170.9739727, {    0.02535859713, -0.0005624372127,  7.125254757e-06,   -5.6392637e-08,  2.855223598e-10, -9.031372192e-13,  1.631718951e-15, -1.289232633e-18 } },
    //        {       200.835733, {    0.02160202683, -0.0004081404218,  4.404526745e-06,  -2.96948791e-08,  1.280727103e-10, -3.450837842e-13,  5.310892721e-16, -3.574389303e-19 } },
    //        {       236.217058, {     0.0183782644, -0.0002954109441,  2.712190466e-06, -1.555619733e-08,  5.707900853e-11, -1.308396294e-13,  1.713066474e-16, -9.808384426e-20 } },
    //        {      278.5480371, {    0.01560544629, -0.0002129922842,  1.660416209e-06, -8.086370077e-09,  2.519265205e-11, -4.903181716e-14,  5.450646845e-17, -2.649722521e-20 } },
    //        {       328.888158, {    0.01322538546, -0.0001529766242,  1.010657807e-06, -4.171225483e-09,  1.101295193e-11, -1.816454687e-14,  1.711229576e-17, -7.049709491e-21 } },
    //        {       388.826323, {    0.01119389428, -0.0001095891819,  6.127924203e-07, -2.140604553e-09,  4.783402695e-12, -6.677518353e-15,  5.324183641e-18, -1.856380979e-21 } },
    //        {       460.873414, {   0.009456181234, -7.820428666e-05,  3.694011563e-07, -1.090026387e-09,  2.057533561e-12, -2.426209052e-15,   1.63404134e-18, -4.812466988e-22 } },
    //        {      546.9743257, {   0.007972798735, -5.559265316e-05,  2.213980224e-07, -5.508051103e-10,  8.765771842e-13,  -8.71466238e-16,  4.948360541e-19,  -1.22868017e-22 } },
    //        {      649.9972436, {    0.00671346195, -3.941718935e-05,  1.321816908e-07, -2.768991575e-10,  3.710536791e-13, -3.106110175e-16,  1.485060498e-19, -3.104809801e-23 } },
    //        {      774.4166184, {   0.005642141921, -2.784031284e-05,  7.845926854e-08, -1.381250975e-10,  1.555461414e-13, -1.094219448e-16,  4.396330597e-20,  -7.72384216e-24 } },
    //        {      923.8407207, {   0.004732622762, -1.958784207e-05,  4.630288565e-08,  -6.83729696e-11,  6.458264924e-14, -3.810680567e-17,  1.284182362e-20, -1.892364265e-24 } },
    //        {      1104.938519, {   0.003962064384, -1.372837832e-05,  2.716733727e-08, -3.358328721e-11,  2.655508788e-14, -1.311659025e-17,  3.700198616e-21, -4.564323633e-25 } },
    //        {      1323.239371, {   0.003310560752, -9.584643227e-06,  1.584809675e-08, -1.636904039e-11,  1.081469788e-14, -4.463246711e-18,  1.051998064e-21, -1.084233884e-25 } },
    //        {      1588.756318, {   0.002760853541, -6.665807756e-06,  9.191411185e-09, -7.916801913e-12,  4.361698734e-15, -1.501067194e-18,  2.950294803e-22, -2.535525428e-26 } },
    //        {       1910.00918, {   0.002297976342,  -4.61799786e-06,  5.300031753e-09, -3.799600767e-12,  1.742337757e-15,  -4.99070502e-19,  8.164110328e-23, -5.839688918e-27 } },
    //        {      2302.142458, {   0.001909015813, -3.186946314e-06,  3.038432671e-09, -1.809473896e-12,   6.89259893e-16, -1.639992925e-19,  2.228494041e-23, -1.324061051e-27 } },
    //        {      2781.938405, {    0.00158181228, -2.188056638e-06,   1.72848178e-09, -8.528877487e-13,  2.691782257e-16, -5.306513455e-20,  5.974233345e-24,   -2.9408656e-28 } },
    //        {      3366.062019, {    0.00130815954, -1.496464963e-06,  9.776231501e-10, -3.989272316e-13,  1.041198018e-16, -1.697425598e-20,  1.580331548e-24, -6.433128829e-29 } },
    //        {      4083.337646, {   0.001079762638, -1.019517189e-06,  5.497341672e-10, -1.851488565e-13,  3.988408536e-17, -5.366472097e-21,  4.123547986e-25,  -1.38536024e-29 } },
    //        {      4966.232498, {  0.0008889499764, -6.910116092e-07,   3.06745377e-10, -8.504985605e-14,  1.508246057e-17, -1.670606728e-21,  1.056725014e-25, -2.922483439e-30 } },
    //        {      6055.602758, {  0.0007299746758, -4.659495821e-07,  1.698428876e-10, -3.866796615e-14,  5.630552573e-18, -5.120905934e-22,  2.659628727e-26, -6.039350958e-31 } },
    //        {      7402.975029, {   0.000597887824, -3.125757638e-07,  9.331699531e-11, -1.740020137e-14,  2.075085971e-18, -1.545633751e-22,  6.574279063e-27, -1.222580169e-31 } },
    //        {      9073.477497, {   0.000488442076,  -2.08609685e-07,  5.087661061e-11, -7.749655856e-15,  7.549660601e-19, -4.593604586e-23,   1.59603844e-27, -2.424450354e-32 } },
    //        {            1e+40, {  0.0004198069977, -1.541829815e-07,  3.235297339e-11, -4.242296633e-15,  3.559564823e-19, -1.866387797e-23,  5.591108644e-28, -7.326578023e-33 } },
    //    };

    // 1e-14
    static const std::vector<std::pair<TF,std::array<TF,8>>> coeffs = {
        {       0.34621456, {              0.5,            -0.25,    0.08333333313,   -0.02083332907,   0.004166620169, -0.0006941623621,  9.824114146e-05, -1.064353408e-05 } },
        {       0.69155856, {     0.4999999962,    -0.2499999372,    0.08333288543,   -0.02083151641,   0.004162075168, -0.0006870156231,   9.16288835e-05, -7.838844456e-06 } },
        {       1.04284944, {     0.4999997964,    -0.2499980616,    0.08332525319,   -0.02081404338,   0.004137731231, -0.0006663474343,  8.171796908e-05, -5.767003935e-06 } },
        {       1.40280336, {     0.4999974983,    -0.2499831779,    0.08328367575,   -0.02074907314,   0.004076374375, -0.0006313184474,  7.052186675e-05, -4.221346776e-06 } },
        {       1.77315856, {     0.4999845213,    -0.2499197246,    0.08315015904,   -0.02059233036,   0.003965485678, -0.0005840388239,  5.927198714e-05, -3.068922585e-06 } },
        {       2.15384976, {     0.4999361595,     -0.249731191,    0.08283427733,   -0.02029745276,   0.003799837664, -0.0005280407445,  4.872357556e-05, -2.214790639e-06 } },
        {       2.54849296, {     0.4997984099,    -0.2492873908,    0.08222021709,   -0.01982443301,   0.003580745569,  -0.000467022411,  3.926207983e-05, -1.584674738e-06 } },
        {       2.95702416, {     0.4994732446,     -0.248399657,    0.08117985183,    -0.0191459712,   0.003314837927, -0.0004043890773,  3.105238402e-05, -1.122729774e-06 } },
        {       3.38118544, {     0.4988078752,    -0.2468315726,    0.07959399685,   -0.01825378867,   0.003013282551, -0.0003431532831,  2.413494741e-05, -7.873929131e-07 } },
        {       3.82437136, {     0.4975849993,    -0.2443086561,    0.07736085204,   -0.01715445003,    0.00268821593, -0.0002854181058,  1.843183252e-05, -5.456899436e-07 } },
        {           4.2849, {      0.495531194,    -0.2405586217,    0.07442364044,   -0.01587517087,   0.002353597084, -0.0002328538068,  1.384024627e-05, -3.736375296e-07 } },
        {       4.76810896, {     0.4923281959,    -0.2353367466,    0.07077216633,   -0.01445550082,   0.002022154331,  -0.000186388187,  1.021837551e-05, -2.525478759e-07 } },
        {       5.27345296, {     0.4876311013,    -0.2284509652,    0.06644294217,   -0.01294227088,   0.001704567895, -0.0001463678155,  7.414632126e-06, -1.683059531e-07 } },
        {       5.80424464, {     0.4811201303,    -0.2198180265,    0.06153416645,   -0.01139061937,   0.001410094803, -0.0001128151239,  5.289356646e-06, -1.105756074e-07 } },
        {         6.360484, {     0.4725326767,     -0.209468984,    0.05618590743,  -0.009854214061,    0.00114512029, -8.537998667e-05,  3.710325603e-06, -7.160387854e-08 } },
        {       6.94638736, {     0.4616998355,    -0.1975534438,    0.05056584391,  -0.008380784642,  0.0009132197217, -6.346917805e-05,  2.559591765e-06, -4.568917392e-08 } },
        {       7.56690064, {     0.4485229486,    -0.1842801809,    0.04483282094,  -0.007004414715,  0.0007148589203, -4.630806551e-05,   1.73435286e-06,  -2.86733235e-08 } },
        {       8.22312976, {     0.4330662472,    -0.1699835453,    0.03916286303,  -0.005754561246,  0.0005494742489, -3.317131423e-05,  1.154372867e-06, -1.769425732e-08 } },
        {       8.92097424, {     0.4155366572,    -0.1550620554,    0.03371693683,   -0.00464983491,   0.000414954939, -2.333887707e-05,  7.549250122e-07, -1.073638543e-08 } },
        {       9.66215056, {     0.3962674328,    -0.1399396606,    0.02862845566,  -0.003698198088,  0.0003081249213, -1.614017861e-05,  4.853193512e-07, -6.407127874e-09 } },
        {      10.45358224, {     0.3756936579,     -0.125030988,    0.02399645147,  -0.002898349005,  0.0002252200971, -1.098214928e-05,  3.069600438e-07, -3.762827615e-09 } },
        {        11.303044, {     0.3542383559,    -0.1106596848,    0.01986922009,  -0.002239589707,  0.0001621064248, -7.352635244e-06,  1.909548197e-07, -2.173162317e-09 } },
        {      12.21362704, {     0.3324064053,   -0.09713285151,    0.01627588236,  -0.001709069003,  0.0001150918541, -4.851783161e-06,  1.170207519e-07, -1.236037835e-09 } },
        {      13.19432976, {     0.3106763697,   -0.08467243218,    0.01321250023,   -0.00129049712,  8.076273101e-05, -3.161814378e-06,  7.078332826e-08, -6.936592254e-10 } },
        {      14.24911504, {     0.2894552592,   -0.07340681819,    0.01064839043, -0.0009661445047,   5.61352151e-05, -2.039419329e-06,  4.235387957e-08, -3.849257784e-10 } },
        {      15.38835984, {     0.2690565422,    -0.0633792946,    0.00853500448,  -0.000718594896,  3.873048269e-05, -1.304916233e-06,  2.512659996e-08, -2.116924131e-10 } },
        {      16.64095562, {     0.2498001687,    -0.0546112801,   0.006823337373, -0.0005328860258,   2.66366306e-05, -8.321844405e-07,  1.485687674e-08, -1.160409942e-10 } },
        {       18.0473482, {     0.2307716902,   -0.04659908121,    0.00537686288, -0.0003877466565,  1.789488026e-05, -5.161388742e-07,  8.506278172e-09, -6.132829672e-11 } },
        {      19.57260055, {     0.2129013767,   -0.03965842947,    0.00422106465, -0.0002807718785,  1.195166324e-05, -3.179396583e-07,  4.832627624e-09, -3.213359332e-11 } },
        {      21.22675797, {     0.1963528171,   -0.03373187163,    0.00331102871, -0.0002031056624,  7.972914098e-06, -1.955904695e-07,  2.741545123e-09, -1.681032097e-11 } },
        {       23.0503801, {     0.1809489842,   -0.02864658723,   0.002591217754, -0.0001464762135,  5.298614093e-06, -1.197814038e-07,  1.547142343e-09, -8.741801603e-12 } },
        {      25.03067229, {      0.166637471,   -0.02429432183,   0.002023718949, -0.0001053481594,   3.50941699e-06, -7.305909537e-08,  8.690141008e-10, -4.521772148e-12 } },
        {      27.18109431, {     0.1534550468,    -0.0206025758,   0.001580428766, -7.576339994e-05,  2.324208407e-06, -4.455760197e-08,  4.880694489e-10, -2.338679054e-12 } },
        {      29.55429803, {     0.1412238028,   -0.01744911359,   0.001231828658, -5.434464703e-05,  1.534238016e-06, -2.706813031e-08,    2.7285711e-10, -1.203205781e-12 } },
        {      32.13470811, {      0.129883591,   -0.01475931262,  0.0009582733998, -3.888143292e-05,  1.009542444e-06, -1.638085056e-08,  1.518657385e-10,  -6.15900814e-13 } },
        {      34.94041592, {     0.1194541228,   -0.01248417484,  0.0007454697023, -2.781824188e-05,  6.642917386e-07, -9.913284453e-09,  8.452550532e-11, -3.152721411e-13 } },
        {      37.99109239, {     0.1098618441,   -0.01055969379,  0.0005799188432, -1.990272961e-05,  4.371067936e-07, -5.999181681e-09,  4.704440324e-11, -1.613806745e-13 } },
        {      41.36135713, {     0.1009749981,  -0.008920386862,  0.0004502600427, -1.420270548e-05,  2.866861238e-07,  -3.61634818e-09,  2.606415939e-11, -8.217566964e-14 } },
        {       45.0306047, {    0.09274721012,  -0.007525886975,  0.0003489188334, -1.010925347e-05,  1.874311106e-07, -2.171662062e-09,  1.437648551e-11, -4.163314904e-14 } },
        {      49.02535846, {    0.08518985046,  -0.006349385474,  0.0002703867559, -7.195601281e-06,  1.225396638e-07, -1.304110079e-09,  7.929790968e-12, -2.109285031e-14 } },
        {      53.44327509, {    0.07819798132,  -0.005349897594,   0.000209123906, -5.108439481e-06,  7.985445425e-08, -7.800751024e-10,  4.353935728e-12, -1.063046007e-14 } },
        {      58.25931194, {     0.0717337038,  -0.004501952367,  0.0001614309572, -3.617421421e-06,  5.187256368e-08, -4.648391481e-10,  2.379994793e-12, -5.330571429e-15 } },
        {      63.59118614, {    0.06576149205,  -0.003783520406,   0.000124373303, -2.554953453e-06,  3.358640489e-08, -2.759108659e-10,  1.295032886e-12, -2.658983376e-15 } },
        {      69.41103181, {    0.06024764613,  -0.003175652196,  9.563843139e-05,  -1.79993447e-06,  2.167732288e-08, -1.631470994e-10,  7.015518513e-13, -1.319662955e-15 } },
        {      75.76350795, {    0.05519611479,  -0.002665445349,  7.354238682e-05, -1.268032533e-06,  1.399096833e-08, -9.646947379e-11,  3.800482641e-13, -6.549534426e-16 } },
        {      82.80392766, {    0.05053562362,  -0.002234325169,  5.644170779e-05,  -8.90998569e-07,  9.000713242e-09, -5.681983384e-11,   2.04940626e-13, -3.233533382e-16 } },
        {      90.49858727, {     0.0462388226,  -0.001870529892,  4.323421309e-05, -6.244728435e-07,  5.771950707e-09,  -3.33391701e-11,  1.100251633e-13, -1.588365515e-16 } },
        {      98.90828284, {    0.04230735791,  -0.001565968162,  3.311730376e-05, -4.376733541e-07,  3.701419439e-09, -1.956183586e-11,  5.906850576e-14, -7.802316264e-17 } },
        {       108.238761, {    0.03868527975,  -0.001309305322,  2.531858766e-05, -3.059563294e-07,  2.365924048e-09, -1.143308948e-11,  3.156674526e-14, -3.812549883e-17 } },
        {      118.4494266, {    0.03535050249,  -0.001093303234,  1.931919495e-05, -2.133333931e-07,  1.507474752e-09, -6.656756958e-12,  1.679495149e-14, -1.853592794e-17 } },
        {      129.6233117, {    0.03230319219, -0.0009129359983,  1.474139469e-05, -1.487504334e-07,  9.605042609e-10, -3.875803935e-12,  8.935681936e-15, -9.011832898e-18 } },
        {      142.0340762, {    0.02949959106, -0.0007613419282,  1.122652009e-05, -1.034499533e-07,  6.100074486e-10, -2.247810697e-12,  4.732450041e-15, -4.358431044e-18 } },
        {      155.6331076, {    0.02692195273, -0.0006341045529,  8.533296645e-06, -7.176167155e-08,   3.86178372e-10, -1.298682848e-12,  2.495283938e-15, -2.097271592e-18 } },
        {      170.7539328, {    0.02455374995, -0.0005274504862,  6.473594708e-06, -4.965092334e-08,  2.436840673e-10, -7.473860397e-13,  1.309672844e-15, -1.003914641e-18 } },
        {      187.3438501, {    0.02237943422, -0.0004381716198,  4.901616487e-06, -3.426513328e-08,  1.532792988e-10, -4.284820425e-13,  6.843554666e-16, -4.781310855e-19 } },
        {      205.5455918, {    0.02039766133, -0.0003640045339,  3.711360576e-06, -2.364708004e-08,  9.641394978e-11, -2.456519803e-13,  3.576025926e-16, -2.277179011e-19 } },
        {      225.8063657, {    0.01857942911, -0.0003020014504,  2.804683709e-06, -1.627702493e-08,   6.04480656e-11, -1.402833371e-13,  1.860064593e-16,  -1.07885816e-19 } },
        {      248.0642585, {    0.01691236533, -0.0002502378506,  2.115435803e-06, -1.117539734e-08,  3.777829367e-11, -7.980645232e-14,  9.632342773e-17, -5.085577065e-20 } },
        {      272.8673019, {    0.01538498429, -0.0002070792803,  1.592475177e-06, -7.652850064e-09,   2.35335918e-11, -4.522387442e-14,  4.965280769e-17, -2.384693035e-20 } },
        {      300.1503114, {    0.01398652273,  -0.000171144144,  1.196494664e-06, -5.227258744e-09,  1.461341465e-11, -2.552959682e-14,  2.548195235e-17, -1.112587226e-20 } },
        {      330.1612499, {    0.01271517828, -0.0001414449479,  8.989775805e-07, -3.570465088e-09,  9.074343164e-12, -1.441186369e-14,  1.307740541e-17, -5.190816252e-21 } },
        {      363.6408712, {    0.01155196395, -0.0001167487912,  6.741299251e-07, -2.432470262e-09,  5.616479307e-12, -8.103890018e-15,  6.680633121e-18, -2.409087359e-21 } },
        {       400.515455, {    0.01048840108, -9.624082049e-05,  5.045497467e-07, -1.652956361e-09,   3.46522528e-12, -4.539565761e-15,  3.397753419e-18, -1.112449352e-21 } },
        {      441.6977191, {   0.009516635302, -7.923294218e-05,   3.76895517e-07, -1.120332844e-09,  2.131002908e-12, -2.532981147e-15,  1.720175554e-18, -5.110020237e-22 } },
        {      487.1144738, {   0.008629339371, -6.514694101e-05,  2.809980716e-07, -7.573969609e-10,  1.306335277e-12, -1.407980671e-15,  8.670250186e-19, -2.335477032e-22 } },
        {        537.89338, {   0.007819741667, -5.349607319e-05,  2.090943982e-07, -5.107071781e-10,  7.981978772e-13, -7.795756773e-16,  4.350078141e-19, -1.061800776e-22 } },
        {      593.9656978, {   0.007081532304, -4.387241829e-05,    1.5529122e-07, -3.434877877e-10,  4.861661777e-13, -4.299988448e-16,   2.17290581e-19, -4.803100964e-23 } },
        {      655.8832349, {   0.006413012336, -3.598000697e-05,  1.153324203e-07, -2.310205639e-10,  2.961139826e-13, -2.371790346e-16,  1.085387321e-19, -2.172703147e-23 } },
        {      725.1886247, {    0.00580386964, -2.946934815e-05,  8.548945405e-08, -1.549750438e-10,  1.797703955e-13, -1.303112249e-16,  5.396779677e-20, -9.776713664e-24 } },
        {      801.8173257, {   0.005249200918, -2.410580412e-05,  6.324689132e-08, -1.036964658e-10,  1.087917363e-13, -7.132390424e-17,  2.671551289e-20, -4.377208871e-24 } },
        {      887.6856051, {    0.00474448947,  -1.96930199e-05,   4.67005936e-08, -6.920505325e-11,  6.562336737e-14, -3.888522181e-17,  1.316430622e-20, -1.949462287e-24 } },
        {      982.7496964, {   0.004285541905, -1.606737159e-05,  3.441685466e-08, -4.606836967e-11,   3.94584305e-14, -2.111942717e-17,  6.458204439e-21, -8.638630728e-25 } },
        {      1089.396429, {   0.003868501363, -1.309233003e-05,  2.531491501e-08, -3.058717863e-11,  2.364861743e-14,  -1.14255088e-17,  3.153775419e-21, -3.807917599e-25 } },
        {      1207.616327, {   0.003489793467, -1.065444634e-05,  1.858435047e-08, -2.025663505e-11,  1.412832287e-14, -6.157683686e-18,  1.533308565e-21,  -1.67010321e-25 } },
        {      1340.390351, {   0.003146135698, -8.659332248e-06,   1.36167907e-08, -1.338032321e-11,  8.413174661e-15, -3.305633741e-18,  7.420497999e-22, -7.286366426e-26 } },
        {      1487.762506, {   0.002834491335, -7.028777081e-06,  9.957901172e-09, -8.815710854e-12,  4.993997509e-15, -1.767831354e-18,   3.57533578e-22, -3.162951072e-26 } },
        {      1653.465791, {   0.002552075839, -5.697898344e-06,  7.268039282e-09, -5.793211203e-12,  2.954751492e-15, -9.417207316e-19,   1.71476591e-22, -1.365794111e-26 } },
        {      1837.624696, {    0.00229631769, -4.613085325e-06,  5.294591427e-09, -3.797282812e-12,  1.742661103e-15, -4.997496559e-19,  8.187920342e-23, -5.868030991e-27 } },
        {      2044.926539, {    0.00206486247, -3.729991349e-06,  3.849496998e-09, -2.482548216e-12,  1.024445648e-15, -2.641665018e-19,   3.89176072e-23, -2.507900058e-27 } },
        {      2275.614035, {   0.001855539647,  -3.01207679e-06,  2.793452609e-09, -1.618878004e-12,  6.003222768e-16, -1.391080908e-19,  1.841619408e-23, -1.066456353e-27 } },
        {      2535.588562, {   0.001666364922, -2.429202877e-06,   2.02318185e-09,  -1.05293508e-12,  3.506418822e-16, -7.296617388e-20,   8.67474905e-24, -4.511136593e-28 } },
        {      2825.263536, {   0.001495512111, -1.956606203e-06,  1.462494964e-09, -6.830946341e-13,  2.041563475e-16, -3.812768982e-20,  4.068137822e-24, -1.898648343e-28 } },
        {      3152.088718, {   0.001341314217, -1.573919744e-06,   1.05513996e-09, -4.420089577e-13,  1.184801409e-16, -1.984511489e-20,  1.899053545e-24, -7.948999892e-29 } },
        {      3521.252665, {   0.001201466977, -1.262824941e-06,   7.58311848e-10, -2.845404355e-13,  6.831743961e-17, -1.024968894e-20,  8.785435161e-25, -3.293862853e-29 } },
        {      3933.652075, {   0.001075506607,  -1.01191874e-06,   5.43940722e-10, -1.827043837e-13,   3.92679113e-17, -5.273733132e-21,  4.046429568e-25, -1.358049075e-29 } },
        {      4400.013293, {   0.000962132999, -8.098183552e-07,  3.894138764e-10, -1.170104401e-13,   2.24971325e-17, -2.702832669e-21,  1.855165114e-25, -5.569722401e-30 } },
        {      4921.664807, {  0.0008601556894,  -6.47249209e-07,  2.782513573e-10, -7.474677116e-14,  1.284804197e-17, -1.379973975e-21,  8.467910543e-26, -2.272841753e-30 } },
        {      5512.255815, {  0.0007684928144, -5.166481945e-07,  1.984354317e-10, -4.762454138e-14,  7.313573941e-18, -7.018040151e-22,  3.847435608e-26, -9.225976627e-31 } },
        {      6173.716691, {  0.0006861553901, -4.118700879e-07,  1.412430948e-10, -3.026644962e-14,  4.149951275e-18, -3.555591488e-22,  1.740403859e-26, -3.726264644e-31 } },
        {      6923.462044, {  0.0006122459785, -3.279178442e-07,  1.003393384e-10, -1.918503398e-14,  2.347140313e-18, -1.794323531e-22,  7.836635584e-27, -1.497069325e-31 } },
        {      7774.262904, {  0.0005455946715, -2.604061628e-07,  7.100637148e-11, -1.209835685e-14,  1.318980737e-18, -8.985320845e-23,  3.496985097e-27, -5.952990546e-32 } },
        {      8729.615806, {  0.0004858858063, -2.065282396e-07,  5.015215484e-11, -7.609966768e-15,  7.388544876e-19,  -4.48247777e-23,  1.553611506e-27, -2.355309224e-32 } },
        {      9815.000499, {  0.0004324332789, -1.635863495e-07,  3.535393539e-11, -4.774295568e-15,  4.125358171e-19, -2.227383117e-23,   6.87054319e-28, -9.269712801e-33 } },
        {            1e+40, {   0.000403624362, -1.425476136e-07,  2.876746788e-11, -3.628444588e-15,  2.928980513e-19, -1.477711462e-23,  4.260120738e-28, -5.373155215e-33 } },
    };

    return exp( w / func.eps ) * func.eps * _r_polynomials_integration( coeffs, TF( 1 ) / sqrt( func.eps ) );
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::boundary_measure() const {
    using std::pow;
    if ( nb_points == 0 )
        return sphere_radius >= 0 ? 2 * M_PI * sphere_radius : TF( 0 );

    TF res = 0;
    for( std::size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ ) {
        if ( allow_ball_cut && arcs[ i0 ] )
            res += _arc_length( point( i0 ), point( i1 ) );
        else
            res += norm_2( point( i1 ) - point( i0 ) );
    }
    return res;
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::integration( FunctionEnum::Unit, TF w ) const {
    //    // nsmake run -g3 src/PowerDiagram/offline_integration/gen_approx_integration.cpp --function Unit --end-log-scale 100 --precision 1e-10 -r 100 -l 100
    //    if ( _cuts.empty() )
    //        return _sphere_radius > 0 ? 2 * M_PI * 0.5 * pow( _sphere_radius, 2 ) : 0;

    //    auto arc_val = []( PT P0, PT P1 ) {
    //        using std::atan2;
    //        using std::pow;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        return ( a1 - a0 ) * 0.5 * pow( dot( P0, P0 ), 1 );
    //    };

    //    auto seg_val = []( PT P0, PT P1 ) {
    //        return -0.25 * ( ( P0.x + P1.x ) * ( P0.y - P1.y ) - ( P0.x - P1.x ) * ( P0.y + P1.y ) );
    //    };

    //    TF res = 0;
    //    for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //        if ( _cuts[ i0 ].seg_type == SegType::arc )
    //            res += arc_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        else
    //            res += seg_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //    }
    //    return res;

    // hand coded version:
    if ( nb_points == 0 )
        return sphere_radius > 0 ? TF( M_PI ) * pow( sphere_radius, 2 ) : TF( 0 );

    // triangles
    TF res = 0;
    Pt A = point( 0 );
    for( size_t i = 2; i < nb_points; ++i ) {
        Pt B = point( i - 1 );
        Pt C = point( i - 0 );
        TF tr2_area = A.x * ( B.y - C.y ) + B.x * ( C.y - A.y ) + C.x * ( A.y - B.y );
        res += tr2_area;
    }
    res *= 0.5;

    // arcs
    if ( allow_ball_cut ) {
        for( size_t i0 = nb_points - 1, i1 = 0; i1 < nb_points; i0 = i1++ )
            if ( arcs[ i0 ] )
                res += _arc_area( point( i0 ), point( i1 ) );
    }

    return res;
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::integration( FunctionEnum::R2, TF w ) const {
    TODO;
    return 0;
    //    using std::pow;

    //    // generated using nsmake run -g3 src/PowerDiagram/offline_integration/gen_approx_integration.cpp --function R2 --end-log-scale 100 --precision 1e-10 -r 100 -l 100
    //    if ( _cuts.empty() )
    //        return _sphere_radius > 0 ? 2 * M_PI * 0.25 * pow( _sphere_radius, 4 ) : 0;

    //    auto arc_val = []( PT P0, PT P1 ) {
    //        using std::atan2;
    //        using std::pow;
    //        TF a0 = atan2( P0.y, P0.x );
    //        TF a1 = atan2( P1.y, P1.x );
    //        if ( a1 < a0 )
    //            a1 += 2 * M_PI;
    //        return ( a1 - a0 ) * 0.25 * pow( dot( P0, P0 ), 2 );
    //    };

    //    // hand coded version
    //    auto seg_val = []( PT P0, PT P1 ) {
    //        TF result;
    //        TF R_0 = P1.y; TF R_1 = (-1.0)*R_0; TF R_2 = P0.y; TF R_3 = (-1.0)*R_2;
    //        R_3 = R_0+R_3; TF R_4 = pow(R_3,2); R_1 = R_2+R_1; R_0 = R_2+R_0;
    //        R_2 = R_3*R_0; TF R_5 = pow(R_0,2); TF R_6 = P0.x; TF R_7 = (-1.0)*R_6;
    //        TF R_8 = P1.x; R_7 = R_8+R_7; TF R_9 = R_1*R_7; R_9 = (-1.0)*R_9;
    //        TF R_10 = pow(R_7,2); R_4 = R_10+R_4; R_10 = (-1.0)*R_8; R_10 = R_6+R_10;
    //        R_3 = R_10*R_3; R_3 = R_9+R_3; R_0 = R_10*R_0; R_0 = (-1.0)*R_0;
    //        R_6 = R_8+R_6; R_7 = R_7*R_6; R_7 = R_2+R_7; R_7 = R_3*R_7;
    //        R_7 = 2.0*R_7; R_1 = R_1*R_6; R_0 = R_1+R_0; R_4 = R_4*R_0;
    //        R_4 = (-1.0)*R_4; R_7 = R_4+R_7; R_6 = pow(R_6,2); R_5 = R_6+R_5;
    //        R_6 = 0.25; R_7 = R_6*R_7; R_7 = (1.0/24.0)*R_7; R_5 = R_6*R_5;
    //        R_0 = R_5*R_0; R_0 = -0.125*R_0; R_7 = R_0+R_7; result = R_7;
    //        return result;
    //    };

    //    TF res = 0;
    //    for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ ) {
    //        if ( _cuts[ i0 ].seg_type == SegType::arc )
    //            res += arc_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //        else
    //            res += seg_val( _cuts[ i0 ].point - _sphere_center, _cuts[ i1 ].point - _sphere_center );
    //    }
    //    return res;
}

template<class Pc, class CI>
typename ConvexPolyhedron2<Pc,CI>::TF ConvexPolyhedron2<Pc,CI>::integration( SpaceFunctions::Constant<TF> cst ) const {
    return integration( FunctionEnum::Unit() ) * cst.coeff;
}

template<class Pc,class CI> template<class FU>
typename ConvexPolyhedron2<Pc,CI>::Pt ConvexPolyhedron2<Pc,CI>::centroid_ap( const FU &func, std::size_t n ) const {
    TODO;
    return {};
    //    auto rd01 = []() {
    //        return rand() / ( RAND_MAX - 1.0 );
    //    };

    //    auto rdm1 = []() {
    //        return 2.0 * rand() / ( RAND_MAX - 1.0 ) - 1.0;
    //    };

    //    auto inside_semi_planes = [&]( PT p ) {
    //        for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ )
    //            if ( _cuts[ i0 ].seg_type == SegType::line && dot( p - _cuts[ i0 ].point, _cuts[ i0 ].normal[ 0 ] ) > 0 )
    //                return false;
    //        return true;
    //    };

    //    TF count{ 0 };
    //    PT centroid{ 0, 0 };
    //    if ( _sphere_radius < 0 ) {
    //        PT mi = { + std::numeric_limits<TF>::max(), + std::numeric_limits<TF>::max() };
    //        PT ma = { - std::numeric_limits<TF>::max(), - std::numeric_limits<TF>::max() };
    //        for( size_t i1 = 0; i1 < _cuts.size(); ++i1 ) {
    //            mi = min( mi, _cuts[ i1 ].point );
    //            ma = max( ma, _cuts[ i1 ].point );
    //        }
    //        for( TI i = 0; i < n; ++i ) {
    //            PT p{ mi.x + ( ma.x - mi.x ) * rd01(), mi.y + ( ma.y - mi.y ) * rd01() };
    //            if ( inside_semi_planes( p ) ) {
    //                TF v = func( p, _sphere_center );
    //                centroid += v * p;
    //                count += v;
    //            }
    //        }
    //    } else {
    //        for( TI i = 0; i < n; ++i ) {
    //            PT p{ _sphere_center.x + _sphere_radius * rdm1(), _sphere_center.y + _sphere_radius * rdm1() };
    //            if ( norm_2( p ) <= _sphere_radius && inside_semi_planes( p ) ) {
    //                TF v = func( p, _sphere_center );
    //                centroid += v * p;
    //                count += v;
    //            }
    //        }
    //    }

    //    return centroid / ( count + ( count == 0 ) );
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::measure_ap( TF max_ratio_area_error ) const {
    TODO;
    return 0;
    //    std::vector<PT> points;
    //    for_each_approx_seg( [&]( PT p ) {
    //        points.push_back( p );
    //    }, max_ratio_area_error );

    //    TF res = 0;
    //    const PT &A = points[ 0 ];
    //    for( size_t i = 2; i < points.size(); ++i ) {
    //        const PT &B = points[ i - 1 ];
    //        const PT &C = points[ i - 0 ];
    //        res += A.x * ( B.y - C.y ) + B.x * ( C.y - A.y ) + C.x * ( A.y - B.y );
    //    }

    //    return 0.5 * res;
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::_centroid_arc( Pt &ctd, TF &mea, Pt p0, Pt p1, TF coeff ) const {
    using std::atan2;
    using std::sqrt;
    using std::max;
    using std::pow;
    using std::sin;
    using std::cos;

    p0 -= sphere_center;
    p1 -= sphere_center;

    TF a0 = atan2( p0.y, p0.x );
    TF a1 = atan2( p1.y, p1.x );
    if ( a1 < a0 )
        a1 += 2 * M_PI;

    TF d0 = norm_2( p0 - p1 ) / 2;
    TF d1 = sqrt( max( TF( 0 ), pow( sphere_radius, 2 ) - d0 * d0 ) );

    TF tri_area = d1 * d0;
    Pt tri_bary = TF( 1 ) / 3 * ( p0 + p1 );

    TF pie_area = pow( sphere_radius, 2 ) * ( a1 - a0 ) / 2;
    Pt pie_mbar = TF( 2 ) / 3 * pow( sphere_radius, 3 ) * sin( ( a1 - a0 ) / 2 ) *
            Pt{ cos( ( a1 + a0 ) / 2 ), sin( ( a1 + a0 ) / 2 ) };

    TF lea;
    if ( a1 - a0 > M_PI ) {
        ctd += coeff * ( pie_mbar + tri_area * tri_bary );
        lea = ( pie_area + tri_area ) * coeff;
    } else {
        ctd += coeff * ( pie_mbar - tri_area * tri_bary );
        lea = ( pie_area - tri_area ) * coeff;
    }

    ctd += lea * sphere_center;
    mea += lea;
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::boundary_measure_ap( TF max_ratio_area_error ) const {
    TODO;
    return 0;
    //    std::vector<PT> lines;
    //    for_each_approx_seg( [&]( PT p ) {
    //        lines.push_back( p );
    //    }, max_ratio_area_error );

    //    TF res = 0;
    //    for( size_t i = 1; i < lines.size(); ++i )
    //        res += norm_2( lines[ i ] - lines[ i - 1 ] );
    //    return res;
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::_arc_length( Pt p0, Pt p1 ) const {
    using std::atan2;
    TF a0 = atan2( p0.y - sphere_center.y, p0.x - sphere_center.x );
    TF a1 = atan2( p1.y - sphere_center.y, p1.x - sphere_center.x );
    if ( a1 < a0 )
        a1 += 2 * M_PI;
    return ( a1 - a0 ) * sphere_radius;
}

template<class Pc,class CI>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::_arc_area( Pt p0, Pt p1 ) const {
    using std::atan2;
    using std::sqrt;
    using std::pow;
    using std::max;

    TF a0 = atan2( p0.y - sphere_center.y, p0.x - sphere_center.x );
    TF a1 = atan2( p1.y - sphere_center.y, p1.x - sphere_center.x );
    if ( a1 < a0 )
        a1 += 2 * M_PI;
    TF d0 = norm_2( p0 - p1 ) / 2;
    TF d1 = sqrt( max( TF( 0 ), pow( sphere_radius, 2 ) - d0 * d0 ) );
    TF rs = pow( sphere_radius, 2 ) * ( a1 - a0 ) / 2;

    if ( a1 - a0 > M_PI )
        rs += d1 * d0;
    else
        rs -= d1 * d0;
    return rs;
}

template<class Pc, class CI>
void ConvexPolyhedron2<Pc,CI>::display_asy( std::ostream &os, const std::string &draw_info, const std::string &fill_info, bool want_fill, bool avoid_bounds, bool want_line ) const {
    using std::atan2;
    using std::sin;
    using std::cos;

    for( int fill = want_fill; fill >= 0; --fill ) {
        if ( fill == false && want_line == false )
            continue;

        bool has_avoided_line = false;
        const auto &info = fill ? fill_info : draw_info;
        if ( nb_points ) {
            os << ( fill ? "fill" : "draw" ) << "(";
            for( TI i = 0; i < nb_points; ++i ) {
                Pt p0 = point( i );
                if ( arcs[ i ] ) {
                    Pt p1 = point( ( i + 1 ) % nb_points );
                    TF a0 = atan2( p0.y - sphere_center.y, p0.x - sphere_center.x );
                    TF a1 = atan2( p1.y - sphere_center.y, p1.x - sphere_center.x );
                    if ( a1 < a0 )
                        a1 += 2 * M_PI;

                    size_t n = 10;
                    for( size_t i = 0; i < n; ++i ) {
                        TF ai = a0 + ( a1 - a0 ) * i / n;
                        os.precision( 16 );
                        os << "(" << sphere_center.x + sphere_radius * cos( ai ) << "," << sphere_center.y + sphere_radius * sin( ai ) << ")..";
                    }
                } else {
                    os.precision( 16 );
                    os << "(" << point( i )[ 0 ] << "," << point( i )[ 1 ] << ")";

                    bool last_line_avoided = avoid_bounds && cut_ids[ i ] == TI( -1 ) && fill == false;
                    has_avoided_line |= last_line_avoided;
                    if ( last_line_avoided )
                        os << "^^";
                    else {
                        os << "--";
                    }
                }
            }

            // closing
            if ( fill || has_avoided_line == false )
                os << "cycle";
            else
                os << "(" << point( 0 )[ 0 ] << "," << point( 0 )[ 1 ] << ")";

            // info
            os << ( info.empty() ? "" : "," ) << info << ");\n";
        } else if ( sphere_radius > 0 )
            os << ( fill ? "fill" : "draw" ) << "(circle((" << sphere_center.x << "," << sphere_center.y << ")," << sphere_radius << ")" << ( info.empty() ? "" : "," ) << info << ");\n";
    }
}

template<class Pc, class CI>
void ConvexPolyhedron2<Pc, CI>::intersect_with( const ConvexPolyhedron2 &cp ) {
    ASSERT( sphere_radius <= 0, "TODO: intersect ball cutted with ball cutted convex polyhedron" );
    if ( cp.nb_points ) {
        bool has_sphere_cut = false;
        for( TI i = 0; i < cp.nb_points; ++i ) {
            if ( cp.arcs[ i ] )
                has_sphere_cut = true;
            else
                plane_cut( cp.point( i ), cp.normal( i ), cp.cut_ids[ i ] );
        }

        if ( has_sphere_cut )
            ball_cut( cp.sphere_center, cp.sphere_radius, cp.sphere_cut_id );
    } else {
        if ( cp.sphere_radius > 0 ) {
            ball_cut( cp.sphere_center, cp.sphere_radius, cp.sphere_cut_id );
        } else {
            sphere_radius = -1;
            nb_points = 0;
        }
    }
}

template<class Pc,class CI> template<class V>
void ConvexPolyhedron2<Pc,CI>::display( V &vo, const typename V::CV &cell_data, bool filled, TF max_ratio_area_error, bool display_tangents ) const {
    using VTF = typename V::TF;

    std::vector<typename V::PT> lines;
    for_each_approx_seg( [&]( Pt p ) {
        if ( lines.empty() || lines.back().x != p.x || lines.back().y != p.y )
            lines.push_back( { VTF( p.x ), VTF( p.y ), VTF( 0 ) } );
    }, max_ratio_area_error );

    vo.mutex.lock();

    if ( filled )
        vo.add_polygon( lines, cell_data );
    else
        vo.add_lines( lines, cell_data );

    if ( display_tangents ) {
        //        for( TI i = 0; i < edges.size() / 2; ++i ) {
        //            const Edge &edge = edges[ 2 * i ];
        //            vo.add_arrow( nodes[ edge.n0 ].pos, edge.tangent_0, cell_data );
        //            vo.add_arrow( nodes[ edge.n1 ].pos, edge.tangent_1, cell_data );
        //        }
        TODO;
    }

    vo.mutex.unlock();
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_approx_seg( const std::function<void( Pt )> &f, TF max_ratio_area_error ) const {
    using std::atan2;

    if ( nb_points == 0 ) {
        if ( sphere_radius > 0 ) {
            size_t n = 20;
            for( size_t i = 0; i <= n; ++i ) {
                f( {
                    sphere_center.x + sphere_radius * cos( 2 * M_PI * i / n ),
                    sphere_center.y + sphere_radius * sin( 2 * M_PI * i / n )
                } );
            }
        }
        return;
    }

    for( std::size_t i = 0; i < nb_points; ++i ) {
        std::size_t j = ( i + 1 ) % nb_points;
        switch ( arcs[ i ] ) {
        case 0:
            f( point( i ) );
            break;
        case 1: {
            Pt p0 = point( i );
            Pt p1 = point( j );
            TF a0 = atan2( p0.y - sphere_center.y, p0.x - sphere_center.x );
            TF a1 = atan2( p1.y - sphere_center.y, p1.x - sphere_center.x );
            if ( a1 < a0 )
                a1 += 2 * M_PI;
            size_t n = 20; // TODO
            for( size_t i = 0; i < n; ++i ) {
                TF ai = a0 + ( a1 - a0 ) * i / n;
                f( {
                    sphere_center.x + sphere_radius * cos( ai ),
                    sphere_center.y + sphere_radius * sin( ai )
                } );
            }
            break;
        }
        }
    }
    f( point( 0 ) );
}

template<class Pc, class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_simplex( const std::function<void( CI, CI )> &f ) const {
    for( std::size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ )
        if ( allow_ball_cut == false || ( arcs[ i0 ] == false && arcs[ i1 ] == false ) )
            f( cut_ids[ i0 ], cut_ids[ i1 ] );
}

template<class Pc, class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_bound( const std::function<void(Pt,Pt,CI)> &f ) const {
    for( std::size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ )
        if ( allow_ball_cut == false || arcs[ i0 ] == false )
            f( point( i0 ), point( i1 ), cut_ids[ i0 ] );
}

template<class Pc,class CI>
void ConvexPolyhedron2<Pc,CI>::for_each_node( const std::function<void( Pt )> &f ) const {
    for( std::size_t i = 0; i < nb_points; ++i )
        f( point( i ) );
}

template<class Pc,class CI> template<class Fu>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::integration_ap( const Fu &func, std::size_t n ) const {
    TODO;
    return 0;
    //    auto rd01 = []() {
    //        return rand() / ( RAND_MAX - 1.0 );
    //    };

    //    auto inside_semi_planes = [&]( PT p ) {
    //        for( size_t i1 = 0, i0 = _cuts.size() - 1; i1 < _cuts.size(); i0 = i1++ )
    //            if ( _cuts[ i0 ].seg_type == SegType::line && dot( p - _cuts[ i0 ].point, _cuts[ i1 ].normal[ 0 ] ) > 0 )
    //                return false;
    //        return true;
    //    };

    //    PT mi, ma;
    //    TF value{ 0 };
    //    if ( _sphere_radius >= 0 ) {
    //        mi = { _sphere_center.x - _sphere_radius, _sphere_center.y - _sphere_radius };
    //        ma = { _sphere_center.x + _sphere_radius, _sphere_center.y + _sphere_radius };

    //        for( TI i = 0; i < n; ++i ) {
    //            PT p{ mi.x + ( ma.x - mi.x ) * rd01(), mi.y + ( ma.y - mi.y ) * rd01() };
    //            if ( norm_2( p - _sphere_center ) <= _sphere_radius && inside_semi_planes( p ) )
    //                value += func( p, _sphere_center );
    //        }
    //    } else {
    //        if ( _cuts.empty() )
    //            return 0;

    //        mi = { + std::numeric_limits<TF>::max(), + std::numeric_limits<TF>::max() };
    //        ma = { - std::numeric_limits<TF>::max(), - std::numeric_limits<TF>::max() };
    //        for( size_t i = 0; i < _cuts.size(); ++i ) {
    //            mi = min( mi, _cuts[ i ].point );
    //            ma = max( ma, _cuts[ i ].point );
    //        }

    //        for( TI i = 0; i < n; ++i ) {
    //            PT p{ mi.x + ( ma.x - mi.x ) * rd01(), mi.y + ( ma.y - mi.y ) * rd01() };
    //            if ( inside_semi_planes( p ) )
    //                value += func( p, _sphere_center );
    //        }
    //    }

    //    return ( ma.x - mi.x ) * ( ma.y - mi.y ) / n * value;
}

template<class Pc,class CI> template<class Coeffs>
void ConvexPolyhedron2<Pc,CI>::_r_centroid_integration( TF &r_x, TF &r_y, const Coeffs &coeffs, TF scaling ) const {
    using std::pow;

    auto cut_index_r = [&]( TF r2 ) {
        std::size_t beg = 0;
        std::size_t end = coeffs.size() - 1;
        while ( beg < end ) {
            std::size_t mid = beg + ( end - beg ) / 2;
            if ( coeffs[ mid ].first < r2 )
                beg = mid + 1;
            else
                end = mid;
        }
        return beg;
    };

    auto part_int = []( TF &r_x, TF &r_y, auto P0, auto P1, auto u0, auto u1, const auto &poly_coeffs ) {
        // generated using metil src/PowerDiagram/offline_integration/lib/gen_centroid.met
        if ( poly_coeffs.size() == 1 ) {
            TF R_0 = P0.x; TF R_1 = P1.x; R_1 = (-1.0)*R_1; R_0 = R_1+R_0;
            R_1 = u0; R_1 = (-1.0)*R_1; TF R_2 = u1; R_1 = R_2+R_1;
            R_2 = P0.y; TF R_3 = P1.y; R_3 = (-1.0)*R_3; R_2 = R_3+R_2;
            R_3 = poly_coeffs[ 0 ]; R_0 = R_3*R_0; R_0 = R_1*R_0; R_0 = 0.5*R_0;
            r_y += R_0; R_2 = R_3*R_2; R_1 = R_2*R_1; R_1 = -0.5*R_1;
            r_x += R_1;
        }
        if ( poly_coeffs.size() == 2 ) {
            TF R_0 = P1.x; TF R_1 = (-1.0)*R_0; TF R_2 = u0; TF R_3 = (-1.0)*R_2;
            TF R_4 = u1; R_3 = R_4+R_3; TF R_5 = pow(R_3,3); R_2 = R_4+R_2;
            R_4 = R_0*R_2; R_4 = 0.5*R_4; TF R_6 = -0.5*R_2; R_6 = 1.0+R_6;
            TF R_7 = P0.x; R_1 = R_7+R_1; TF R_8 = (-1.0)*R_7; R_8 = R_0+R_8;
            R_8 = pow(R_8,2); R_7 = R_7*R_6; R_4 = R_7+R_4; R_4 = pow(R_4,2);
            R_7 = poly_coeffs[ 1 ]; R_0 = poly_coeffs[ 0 ]; TF R_9 = P1.y; R_2 = R_9*R_2;
            R_2 = 0.5*R_2; TF R_10 = (-1.0)*R_9; TF R_11 = P0.y; TF R_12 = (-1.0)*R_11;
            R_12 = R_9+R_12; R_12 = pow(R_12,2); R_12 = R_8+R_12; R_12 = R_7*R_12;
            R_8 = R_1*R_12; R_8 = R_5*R_8; R_8 = (1.0/24.0)*R_8; R_6 = R_11*R_6;
            R_6 = R_2+R_6; R_6 = pow(R_6,2); R_6 = R_4+R_6; R_6 = R_7*R_6;
            R_6 = R_0+R_6; R_1 = R_1*R_6; R_1 = R_3*R_1; R_1 = 0.5*R_1;
            R_8 = R_1+R_8; r_y += R_8; R_10 = R_11+R_10; R_12 = R_10*R_12;
            R_5 = R_12*R_5; R_5 = (-1.0/24.0)*R_5; R_6 = R_10*R_6; R_3 = R_6*R_3;
            R_3 = -0.5*R_3; R_5 = R_3+R_5; r_x += R_5;
        }
        if ( poly_coeffs.size() == 3 ) {
            TF R_0 = poly_coeffs[ 2 ]; TF R_1 = poly_coeffs[ 1 ]; TF R_2 = P1.x; TF R_3 = (-1.0)*R_2;
            TF R_4 = u0; TF R_5 = (-1.0)*R_4; TF R_6 = u1; R_5 = R_6+R_5;
            TF R_7 = pow(R_5,5); TF R_8 = pow(R_5,3); R_4 = R_6+R_4; R_6 = R_2*R_4;
            R_6 = 0.5*R_6; TF R_9 = -0.5*R_4; R_9 = 1.0+R_9; TF R_10 = P0.x;
            R_3 = R_10+R_3; TF R_11 = (-1.0)*R_10; R_11 = R_2+R_11; R_2 = pow(R_11,2);
            R_10 = R_10*R_9; R_6 = R_10+R_6; R_11 = R_11*R_6; R_6 = pow(R_6,2);
            R_10 = poly_coeffs[ 0 ]; TF R_12 = P1.y; R_4 = R_12*R_4; R_4 = 0.5*R_4;
            TF R_13 = (-1.0)*R_12; TF R_14 = P0.y; TF R_15 = (-1.0)*R_14; R_15 = R_12+R_15;
            R_12 = pow(R_15,2); R_12 = R_2+R_12; R_2 = pow(R_12,2); R_2 = R_0*R_2;
            TF R_16 = R_3*R_2; R_16 = R_7*R_16; R_16 = (1.0/160.0)*R_16; R_9 = R_14*R_9;
            R_4 = R_9+R_4; R_15 = R_15*R_4; R_15 = R_11+R_15; R_15 = pow(R_15,2);
            R_15 = R_0*R_15; R_15 = 4.0*R_15; R_4 = pow(R_4,2); R_4 = R_6+R_4;
            R_0 = R_0*R_4; R_6 = 2.0*R_0; R_6 = R_1+R_6; R_6 = R_12*R_6;
            R_6 = R_15+R_6; R_15 = R_3*R_6; R_15 = R_8*R_15; R_15 = (1.0/24.0)*R_15;
            R_0 = R_1+R_0; R_0 = R_4*R_0; R_0 = R_10+R_0; R_3 = R_3*R_0;
            R_3 = R_5*R_3; R_3 = 0.5*R_3; R_15 = R_3+R_15; R_16 = R_15+R_16;
            r_y += R_16; R_13 = R_14+R_13; R_2 = R_13*R_2; R_7 = R_2*R_7;
            R_7 = (-1.0/160.0)*R_7; R_6 = R_13*R_6; R_8 = R_6*R_8; R_8 = (-1.0/24.0)*R_8;
            R_0 = R_13*R_0; R_5 = R_0*R_5; R_5 = -0.5*R_5; R_8 = R_5+R_8;
            R_7 = R_8+R_7; r_x += R_7;
        }
        if ( poly_coeffs.size() == 4 ) {
            TF R_0 = poly_coeffs[ 3 ]; TF R_1 = poly_coeffs[ 2 ]; TF R_2 = 12.0*R_1; TF R_3 = 2.0*R_1;
            TF R_4 = 4.0*R_1; TF R_5 = poly_coeffs[ 1 ]; TF R_6 = P1.x; TF R_7 = (-1.0)*R_6;
            TF R_8 = u1; TF R_9 = u0; TF R_10 = (-1.0)*R_9; R_10 = R_8+R_10;
            TF R_11 = pow(R_10,7); TF R_12 = pow(R_10,5); TF R_13 = pow(R_10,3); R_8 = R_9+R_8;
            R_9 = R_6*R_8; R_9 = 0.5*R_9; TF R_14 = -0.5*R_8; R_14 = 1.0+R_14;
            TF R_15 = P0.x; R_7 = R_7+R_15; TF R_16 = (-1.0)*R_15; R_16 = R_6+R_16;
            R_6 = pow(R_16,2); R_15 = R_15*R_14; R_9 = R_15+R_9; R_16 = R_16*R_9;
            R_9 = pow(R_9,2); R_15 = poly_coeffs[ 0 ]; TF R_17 = P0.y; TF R_18 = (-1.0)*R_17;
            R_14 = R_17*R_14; TF R_19 = P1.y; R_18 = R_19+R_18; TF R_20 = pow(R_18,2);
            R_20 = R_6+R_20; R_6 = pow(R_20,3); R_6 = R_0*R_6; TF R_21 = R_7*R_6;
            R_21 = R_11*R_21; R_21 = (1.0/896.0)*R_21; R_8 = R_19*R_8; R_8 = 0.5*R_8;
            R_14 = R_8+R_14; R_18 = R_18*R_14; R_18 = R_16+R_18; R_18 = pow(R_18,2);
            R_16 = R_0*R_18; R_16 = 144.0*R_16; R_14 = pow(R_14,2); R_14 = R_9+R_14;
            R_0 = R_0*R_14; R_9 = 36.0*R_0; R_9 = R_2+R_9; R_9 = R_20*R_9;
            R_16 = R_9+R_16; R_16 = R_20*R_16; R_9 = R_7*R_16; R_9 = R_12*R_9;
            R_9 = (1.0/1920.0)*R_9; R_2 = 3.0*R_0; R_2 = R_3+R_2; R_2 = R_14*R_2;
            R_2 = R_5+R_2; R_2 = R_20*R_2; R_20 = 12.0*R_0; R_20 = R_4+R_20;
            R_18 = R_20*R_18; R_2 = R_18+R_2; R_18 = R_7*R_2; R_18 = R_13*R_18;
            R_18 = (1.0/24.0)*R_18; R_0 = R_1+R_0; R_0 = R_14*R_0; R_0 = R_5+R_0;
            R_0 = R_14*R_0; R_0 = R_15+R_0; R_7 = R_7*R_0; R_7 = R_10*R_7;
            R_7 = 0.5*R_7; R_18 = R_7+R_18; R_9 = R_18+R_9; R_21 = R_9+R_21;
            r_y += R_21; R_19 = (-1.0)*R_19; R_17 = R_19+R_17; R_6 = R_17*R_6;
            R_11 = R_6*R_11; R_11 = (-1.0/896.0)*R_11; R_16 = R_17*R_16; R_12 = R_16*R_12;
            R_12 = (-1.0/1920.0)*R_12; R_2 = R_17*R_2; R_13 = R_2*R_13; R_13 = (-1.0/24.0)*R_13;
            R_0 = R_17*R_0; R_10 = R_0*R_10; R_10 = -0.5*R_10; R_13 = R_10+R_13;
            R_12 = R_13+R_12; R_11 = R_12+R_11; r_x += R_11;
        }
        if ( poly_coeffs.size() == 5 ) {
            TF R_0 = poly_coeffs[ 4 ]; TF R_1 = poly_coeffs[ 3 ]; TF R_2 = 24.0*R_1; TF R_3 = 3.0*R_1;
            TF R_4 = 6.0*R_1; TF R_5 = poly_coeffs[ 2 ]; TF R_6 = 2.0*R_5; TF R_7 = poly_coeffs[ 1 ];
            TF R_8 = P1.x; TF R_9 = (-1.0)*R_8; TF R_10 = u1; TF R_11 = u0;
            TF R_12 = (-1.0)*R_11; R_12 = R_10+R_12; TF R_13 = pow(R_12,9); TF R_14 = pow(R_12,7);
            TF R_15 = pow(R_12,5); TF R_16 = pow(R_12,3); R_10 = R_11+R_10; R_11 = R_8*R_10;
            R_11 = 0.5*R_11; TF R_17 = -0.5*R_10; R_17 = 1.0+R_17; TF R_18 = P0.x;
            R_9 = R_18+R_9; TF R_19 = (-1.0)*R_18; R_19 = R_8+R_19; R_8 = pow(R_19,2);
            R_18 = R_18*R_17; R_11 = R_18+R_11; R_19 = R_19*R_11; R_11 = pow(R_11,2);
            R_18 = poly_coeffs[ 0 ]; TF R_20 = P0.y; TF R_21 = (-1.0)*R_20; R_17 = R_20*R_17;
            TF R_22 = P1.y; R_21 = R_22+R_21; TF R_23 = pow(R_21,2); R_23 = R_8+R_23;
            R_8 = pow(R_23,4); R_8 = R_0*R_8; TF R_24 = R_9*R_8; R_24 = R_13*R_24;
            R_24 = (1.0/4608.0)*R_24; TF R_25 = pow(R_23,2); R_10 = R_22*R_10; R_10 = 0.5*R_10;
            R_17 = R_10+R_17; R_21 = R_21*R_17; R_21 = R_19+R_21; R_21 = pow(R_21,2);
            R_19 = R_0*R_21; R_10 = 8640.0*R_19; R_19 = 192.0*R_19; R_17 = pow(R_17,2);
            R_17 = R_11+R_17; R_0 = R_0*R_17; R_11 = 96.0*R_0; R_11 = R_2+R_11;
            R_2 = R_21*R_11; R_11 = R_23*R_11; TF R_26 = 15.0*R_11; R_10 = R_26+R_10;
            R_10 = R_25*R_10; R_25 = R_9*R_10; R_25 = R_14*R_25; R_25 = (1.0/322560.0)*R_25;
            R_11 = 5.0*R_11; R_11 = R_19+R_11; R_11 = R_21*R_11; R_19 = 4.0*R_0;
            R_19 = R_3+R_19; R_19 = R_17*R_19; R_19 = R_6+R_19; R_6 = R_17*R_19;
            R_6 = R_7+R_6; R_6 = R_23*R_6; R_19 = 2.0*R_19; R_3 = 16.0*R_0;
            R_3 = R_4+R_3; R_3 = R_17*R_3; R_19 = R_3+R_19; R_3 = R_23*R_19;
            R_3 = 3.0*R_3; R_3 = R_2+R_3; R_3 = R_23*R_3; R_3 = R_11+R_3;
            R_11 = R_9*R_3; R_11 = R_15*R_11; R_11 = (1.0/1920.0)*R_11; R_21 = R_19*R_21;
            R_6 = R_21+R_6; R_21 = R_9*R_6; R_21 = R_16*R_21; R_21 = (1.0/24.0)*R_21;
            R_0 = R_1+R_0; R_0 = R_17*R_0; R_0 = R_5+R_0; R_0 = R_17*R_0;
            R_0 = R_7+R_0; R_0 = R_17*R_0; R_0 = R_18+R_0; R_9 = R_9*R_0;
            R_9 = R_12*R_9; R_9 = 0.5*R_9; R_21 = R_9+R_21; R_11 = R_21+R_11;
            R_25 = R_11+R_25; R_24 = R_25+R_24; r_y += R_24; R_22 = (-1.0)*R_22;
            R_20 = R_22+R_20; R_8 = R_20*R_8; R_13 = R_8*R_13; R_13 = (-1.0/4608.0)*R_13;
            R_10 = R_20*R_10; R_14 = R_10*R_14; R_14 = (-1.0/322560.0)*R_14; R_3 = R_20*R_3;
            R_15 = R_3*R_15; R_15 = (-1.0/1920.0)*R_15; R_6 = R_20*R_6; R_16 = R_6*R_16;
            R_16 = (-1.0/24.0)*R_16; R_0 = R_20*R_0; R_12 = R_0*R_12; R_12 = -0.5*R_12;
            R_16 = R_12+R_16; R_15 = R_16+R_15; R_14 = R_15+R_14; R_13 = R_14+R_13;
            r_x += R_13;
        }
        if ( poly_coeffs.size() == 6 ) {
            TF R_0 = poly_coeffs[ 5 ]; TF R_1 = poly_coeffs[ 4 ]; TF R_2 = 32.0*R_1; TF R_3 = 192.0*R_1;
            TF R_4 = 4.0*R_1; TF R_5 = 8.0*R_1; TF R_6 = poly_coeffs[ 3 ]; TF R_7 = 3.0*R_6;
            TF R_8 = poly_coeffs[ 2 ]; TF R_9 = 2.0*R_8; TF R_10 = poly_coeffs[ 1 ]; TF R_11 = P1.x;
            TF R_12 = (-1.0)*R_11; TF R_13 = u1; TF R_14 = u0; TF R_15 = (-1.0)*R_14;
            R_15 = R_13+R_15; TF R_16 = pow(R_15,11); TF R_17 = pow(R_15,9); TF R_18 = pow(R_15,7);
            TF R_19 = pow(R_15,5); TF R_20 = pow(R_15,3); R_13 = R_14+R_13; R_14 = R_11*R_13;
            R_14 = 0.5*R_14; TF R_21 = -0.5*R_13; R_21 = 1.0+R_21; TF R_22 = P0.x;
            R_12 = R_22+R_12; TF R_23 = (-1.0)*R_22; R_23 = R_11+R_23; R_11 = pow(R_23,2);
            R_22 = R_22*R_21; R_14 = R_22+R_14; R_23 = R_23*R_14; R_14 = pow(R_14,2);
            R_22 = poly_coeffs[ 0 ]; TF R_24 = P1.y; R_13 = R_24*R_13; R_13 = 0.5*R_13;
            TF R_25 = (-1.0)*R_24; TF R_26 = P0.y; TF R_27 = (-1.0)*R_26; R_27 = R_24+R_27;
            R_24 = pow(R_27,2); R_11 = R_24+R_11; R_24 = pow(R_11,5); R_24 = R_0*R_24;
            TF R_28 = R_12*R_24; R_28 = R_16*R_28; R_28 = (1.0/22528.0)*R_28; TF R_29 = pow(R_11,3);
            R_21 = R_26*R_21; R_13 = R_21+R_13; R_27 = R_27*R_13; R_23 = R_27+R_23;
            R_23 = pow(R_23,2); R_27 = R_0*R_23; R_21 = 806400.0*R_27; R_27 = 28800.0*R_27;
            R_13 = pow(R_13,2); R_13 = R_14+R_13; R_0 = R_0*R_13; R_14 = 120.0*R_0;
            R_14 = R_2+R_14; R_14 = R_13*R_14; R_2 = 960.0*R_0; R_2 = R_3+R_2;
            R_3 = R_11*R_2; TF R_30 = 105.0*R_3; R_21 = R_30+R_21; R_21 = R_29*R_21;
            R_29 = R_12*R_21; R_29 = R_17*R_29; R_29 = (1.0/92897280.0)*R_29; R_3 = 42.0*R_3;
            R_27 = R_3+R_27; R_27 = R_23*R_27; R_2 = R_23*R_2; R_3 = 5.0*R_0;
            R_3 = R_4+R_3; R_3 = R_13*R_3; R_3 = R_7+R_3; R_7 = R_13*R_3;
            R_7 = R_9+R_7; R_9 = R_13*R_7; R_9 = R_10+R_9; R_9 = R_11*R_9;
            R_7 = 2.0*R_7; R_3 = 2.0*R_3; R_4 = 20.0*R_0; R_4 = R_5+R_4;
            R_4 = R_13*R_4; R_3 = R_4+R_3; R_4 = 4.0*R_3; R_14 = R_4+R_14;
            R_4 = R_23*R_14; R_14 = R_11*R_14; R_14 = 5.0*R_14; R_14 = R_2+R_14;
            R_2 = R_11*R_14; R_2 = 3.0*R_2; R_27 = R_2+R_27; R_27 = R_11*R_27;
            R_2 = R_12*R_27; R_2 = R_18*R_2; R_2 = (1.0/322560.0)*R_2; R_14 = R_23*R_14;
            R_3 = R_13*R_3; R_7 = R_3+R_7; R_3 = R_11*R_7; R_3 = 3.0*R_3;
            R_3 = R_4+R_3; R_3 = R_11*R_3; R_3 = R_14+R_3; R_14 = R_12*R_3;
            R_14 = R_19*R_14; R_14 = (1.0/1920.0)*R_14; R_23 = R_7*R_23; R_9 = R_23+R_9;
            R_23 = R_12*R_9; R_23 = R_20*R_23; R_23 = (1.0/24.0)*R_23; R_0 = R_1+R_0;
            R_0 = R_13*R_0; R_0 = R_6+R_0; R_0 = R_13*R_0; R_0 = R_8+R_0;
            R_0 = R_13*R_0; R_0 = R_10+R_0; R_0 = R_13*R_0; R_0 = R_22+R_0;
            R_12 = R_12*R_0; R_12 = R_15*R_12; R_12 = 0.5*R_12; R_23 = R_12+R_23;
            R_14 = R_23+R_14; R_2 = R_14+R_2; R_29 = R_2+R_29; R_28 = R_29+R_28;
            r_y += R_28; R_25 = R_26+R_25; R_24 = R_25*R_24; R_16 = R_24*R_16;
            R_16 = (-1.0/22528.0)*R_16; R_21 = R_25*R_21; R_17 = R_21*R_17; R_17 = (-1.0/92897280.0)*R_17;
            R_27 = R_25*R_27; R_18 = R_27*R_18; R_18 = (-1.0/322560.0)*R_18; R_3 = R_25*R_3;
            R_19 = R_3*R_19; R_19 = (-1.0/1920.0)*R_19; R_9 = R_25*R_9; R_20 = R_9*R_20;
            R_20 = (-1.0/24.0)*R_20; R_0 = R_25*R_0; R_15 = R_0*R_15; R_15 = -0.5*R_15;
            R_20 = R_15+R_20; R_19 = R_20+R_19; R_18 = R_19+R_18; R_17 = R_18+R_17;
            R_16 = R_17+R_16; r_x += R_16;
        }
        if ( poly_coeffs.size() == 7 ) {
            TF R_0 = poly_coeffs[ 6 ]; TF R_1 = poly_coeffs[ 5 ]; TF R_2 = 1920.0*R_1; TF R_3 = 40.0*R_1;
            TF R_4 = 240.0*R_1; TF R_5 = 5.0*R_1; TF R_6 = 10.0*R_1; TF R_7 = poly_coeffs[ 4 ];
            TF R_8 = 4.0*R_7; TF R_9 = poly_coeffs[ 3 ]; TF R_10 = 3.0*R_9; TF R_11 = poly_coeffs[ 2 ];
            TF R_12 = 2.0*R_11; TF R_13 = poly_coeffs[ 1 ]; TF R_14 = P1.x; TF R_15 = (-1.0)*R_14;
            TF R_16 = u0; TF R_17 = (-1.0)*R_16; TF R_18 = u1; R_17 = R_18+R_17;
            TF R_19 = pow(R_17,13); TF R_20 = pow(R_17,11); TF R_21 = pow(R_17,9); TF R_22 = pow(R_17,7);
            TF R_23 = pow(R_17,5); TF R_24 = pow(R_17,3); R_16 = R_18+R_16; R_18 = R_14*R_16;
            R_18 = 0.5*R_18; TF R_25 = -0.5*R_16; R_25 = 1.0+R_25; TF R_26 = P0.x;
            R_15 = R_26+R_15; TF R_27 = (-1.0)*R_26; R_27 = R_14+R_27; R_14 = pow(R_27,2);
            R_26 = R_26*R_25; R_18 = R_26+R_18; R_27 = R_27*R_18; R_18 = pow(R_18,2);
            R_26 = poly_coeffs[ 0 ]; TF R_28 = P1.y; R_16 = R_28*R_16; R_16 = 0.5*R_16;
            TF R_29 = (-1.0)*R_28; TF R_30 = P0.y; TF R_31 = (-1.0)*R_30; R_31 = R_28+R_31;
            R_28 = pow(R_31,2); R_28 = R_14+R_28; R_14 = pow(R_28,6); R_14 = R_0*R_14;
            TF R_32 = R_15*R_14; R_32 = R_19*R_32; R_32 = (1.0/106496.0)*R_32; TF R_33 = pow(R_28,4);
            TF R_34 = pow(R_28,2); R_25 = R_30*R_25; R_16 = R_25+R_16; R_31 = R_31*R_16;
            R_31 = R_27+R_31; R_31 = pow(R_31,2); R_27 = R_0*R_31; R_25 = 108864000.0*R_27;
            TF R_35 = 4769280.0*R_27; R_27 = 23040.0*R_27; R_16 = pow(R_16,2); R_16 = R_18+R_16;
            R_0 = R_0*R_16; R_18 = 11520.0*R_0; R_18 = R_2+R_18; R_2 = R_31*R_18;
            R_18 = R_28*R_18; TF R_36 = 945.0*R_18; R_25 = R_36+R_25; R_25 = R_33*R_25;
            R_33 = R_15*R_25; R_33 = R_20*R_33; R_33 = (1.0/40874803200.0)*R_33; R_36 = 378.0*R_18;
            R_35 = R_36+R_35; R_35 = R_31*R_35; R_18 = 9.0*R_18; R_18 = R_27+R_18;
            R_18 = R_31*R_18; R_27 = 144.0*R_0; R_27 = R_3+R_27; R_27 = R_16*R_27;
            R_3 = 1152.0*R_0; R_3 = R_4+R_3; R_3 = R_16*R_3; R_4 = 6.0*R_0;
            R_4 = R_5+R_4; R_4 = R_16*R_4; R_4 = R_8+R_4; R_8 = R_16*R_4;
            R_8 = R_10+R_8; R_10 = R_16*R_8; R_10 = R_12+R_10; R_12 = R_16*R_10;
            R_12 = R_13+R_12; R_12 = R_28*R_12; R_10 = 2.0*R_10; R_8 = 2.0*R_8;
            R_4 = 2.0*R_4; R_5 = 24.0*R_0; R_5 = R_6+R_5; R_5 = R_16*R_5;
            R_4 = R_5+R_4; R_5 = 4.0*R_4; R_5 = R_27+R_5; R_27 = R_16*R_5;
            R_5 = 6.0*R_5; R_5 = R_3+R_5; R_3 = R_28*R_5; R_3 = 7.0*R_3;
            R_3 = R_2+R_3; R_2 = R_31*R_3; R_3 = R_28*R_3; R_3 = 5.0*R_3;
            R_3 = R_18+R_3; R_18 = 3.0*R_3; R_35 = R_18+R_35; R_35 = R_34*R_35;
            R_34 = R_15*R_35; R_34 = R_21*R_34; R_34 = (1.0/92897280.0)*R_34; R_3 = R_31*R_3;
            R_5 = R_31*R_5; R_4 = R_16*R_4; R_8 = R_4+R_8; R_4 = 4.0*R_8;
            R_4 = R_27+R_4; R_27 = R_31*R_4; R_4 = R_28*R_4; R_4 = 5.0*R_4;
            R_4 = R_5+R_4; R_5 = R_28*R_4; R_5 = 3.0*R_5; R_5 = R_2+R_5;
            R_5 = R_28*R_5; R_5 = R_3+R_5; R_3 = R_15*R_5; R_3 = R_22*R_3;
            R_3 = (1.0/322560.0)*R_3; R_4 = R_31*R_4; R_8 = R_16*R_8; R_10 = R_8+R_10;
            R_8 = R_28*R_10; R_8 = 3.0*R_8; R_8 = R_27+R_8; R_8 = R_28*R_8;
            R_8 = R_4+R_8; R_4 = R_15*R_8; R_4 = R_23*R_4; R_4 = (1.0/1920.0)*R_4;
            R_31 = R_10*R_31; R_12 = R_31+R_12; R_31 = R_15*R_12; R_31 = R_24*R_31;
            R_31 = (1.0/24.0)*R_31; R_0 = R_1+R_0; R_0 = R_16*R_0; R_0 = R_7+R_0;
            R_0 = R_16*R_0; R_0 = R_9+R_0; R_0 = R_16*R_0; R_0 = R_11+R_0;
            R_0 = R_16*R_0; R_0 = R_13+R_0; R_0 = R_16*R_0; R_0 = R_26+R_0;
            R_15 = R_15*R_0; R_15 = R_17*R_15; R_15 = 0.5*R_15; R_31 = R_15+R_31;
            R_4 = R_31+R_4; R_3 = R_4+R_3; R_34 = R_3+R_34; R_33 = R_34+R_33;
            R_32 = R_33+R_32; r_y += R_32; R_29 = R_30+R_29; R_14 = R_29*R_14;
            R_19 = R_14*R_19; R_19 = (-1.0/106496.0)*R_19; R_25 = R_29*R_25; R_20 = R_25*R_20;
            R_20 = (-1.0/40874803200.0)*R_20; R_35 = R_29*R_35; R_21 = R_35*R_21; R_21 = (-1.0/92897280.0)*R_21;
            R_5 = R_29*R_5; R_22 = R_5*R_22; R_22 = (-1.0/322560.0)*R_22; R_8 = R_29*R_8;
            R_23 = R_8*R_23; R_23 = (-1.0/1920.0)*R_23; R_12 = R_29*R_12; R_24 = R_12*R_24;
            R_24 = (-1.0/24.0)*R_24; R_0 = R_29*R_0; R_17 = R_0*R_17; R_17 = -0.5*R_17;
            R_24 = R_17+R_24; R_23 = R_24+R_23; R_22 = R_23+R_22; R_21 = R_22+R_21;
            R_20 = R_21+R_20; R_19 = R_20+R_19; r_x += R_19;
        }
        if ( poly_coeffs.size() == 8 ) {
            TF R_0 = poly_coeffs[ 7 ]; TF R_1 = poly_coeffs[ 6 ]; TF R_2 = 2304.0*R_1; TF R_3 = 23040.0*R_1;
            TF R_4 = 48.0*R_1; TF R_5 = 288.0*R_1; TF R_6 = 6.0*R_1; TF R_7 = 12.0*R_1;
            TF R_8 = poly_coeffs[ 5 ]; TF R_9 = 5.0*R_8; TF R_10 = poly_coeffs[ 4 ]; TF R_11 = 4.0*R_10;
            TF R_12 = poly_coeffs[ 3 ]; TF R_13 = 3.0*R_12; TF R_14 = poly_coeffs[ 2 ]; TF R_15 = 2.0*R_14;
            TF R_16 = poly_coeffs[ 1 ]; TF R_17 = P1.x; TF R_18 = (-1.0)*R_17; TF R_19 = u1;
            TF R_20 = u0; TF R_21 = (-1.0)*R_20; R_21 = R_19+R_21; TF R_22 = pow(R_21,15);
            TF R_23 = pow(R_21,13); TF R_24 = pow(R_21,11); TF R_25 = pow(R_21,9); TF R_26 = pow(R_21,7);
            TF R_27 = pow(R_21,5); TF R_28 = pow(R_21,3); R_19 = R_20+R_19; R_20 = R_17*R_19;
            R_20 = 0.5*R_20; TF R_29 = -0.5*R_19; R_29 = 1.0+R_29; TF R_30 = P0.x;
            R_18 = R_30+R_18; TF R_31 = (-1.0)*R_30; R_31 = R_17+R_31; R_17 = pow(R_31,2);
            R_30 = R_30*R_29; R_20 = R_30+R_20; R_31 = R_31*R_20; R_20 = pow(R_20,2);
            R_30 = poly_coeffs[ 0 ]; TF R_32 = P1.y; R_19 = R_32*R_19; R_19 = 0.5*R_19;
            TF R_33 = (-1.0)*R_32; TF R_34 = P0.y; TF R_35 = (-1.0)*R_34; R_35 = R_32+R_35;
            R_32 = pow(R_35,2); R_17 = R_32+R_17; R_32 = pow(R_17,7); R_32 = R_0*R_32;
            TF R_36 = R_18*R_32; R_36 = R_22*R_36; R_36 = (1.0/491520.0)*R_36; TF R_37 = pow(R_17,5);
            TF R_38 = pow(R_17,3); R_29 = R_34*R_29; R_19 = R_29+R_19; R_35 = R_35*R_19;
            R_31 = R_35+R_31; R_31 = pow(R_31,2); R_35 = R_0*R_31; R_29 = 20118067200.0*R_35;
            TF R_39 = 1016064000.0*R_35; R_35 = 9031680.0*R_35; R_19 = pow(R_19,2); R_19 = R_20+R_19;
            R_0 = R_0*R_19; R_20 = 13440.0*R_0; R_20 = R_2+R_20; R_20 = R_19*R_20;
            R_2 = 161280.0*R_0; R_2 = R_3+R_2; R_3 = R_17*R_2; TF R_40 = 10395.0*R_3;
            R_29 = R_40+R_29; R_29 = R_37*R_29; R_37 = R_18*R_29; R_37 = R_23*R_37;
            R_37 = (1.0/25505877196800.0)*R_37; R_40 = 4620.0*R_3; R_40 = R_39+R_40; R_40 = R_31*R_40;
            R_3 = 165.0*R_3; R_3 = R_35+R_3; R_3 = R_31*R_3; R_2 = R_31*R_2;
            R_35 = 168.0*R_0; R_35 = R_4+R_35; R_35 = R_19*R_35; R_4 = 1344.0*R_0;
            R_4 = R_5+R_4; R_4 = R_19*R_4; R_5 = 7.0*R_0; R_5 = R_6+R_5;
            R_5 = R_19*R_5; R_5 = R_9+R_5; R_9 = R_19*R_5; R_9 = R_11+R_9;
            R_11 = R_19*R_9; R_11 = R_13+R_11; R_13 = R_19*R_11; R_13 = R_15+R_13;
            R_15 = R_19*R_13; R_15 = R_16+R_15; R_15 = R_17*R_15; R_13 = 2.0*R_13;
            R_11 = 2.0*R_11; R_9 = 2.0*R_9; R_5 = 2.0*R_5; R_6 = 28.0*R_0;
            R_6 = R_7+R_6; R_6 = R_19*R_6; R_5 = R_6+R_5; R_6 = 4.0*R_5;
            R_6 = R_35+R_6; R_35 = R_19*R_6; R_6 = 6.0*R_6; R_6 = R_4+R_6;
            R_4 = 8.0*R_6; R_4 = R_20+R_4; R_20 = R_31*R_4; R_4 = R_17*R_4;
            R_4 = 9.0*R_4; R_4 = R_2+R_4; R_2 = R_17*R_4; R_7 = 105.0*R_2;
            R_40 = R_7+R_40; R_40 = R_38*R_40; R_38 = R_18*R_40; R_38 = R_24*R_38;
            R_38 = (1.0/40874803200.0)*R_38; R_2 = 42.0*R_2; R_3 = R_2+R_3; R_3 = R_31*R_3;
            R_4 = R_31*R_4; R_6 = R_19*R_6; R_5 = R_19*R_5; R_9 = R_5+R_9;
            R_5 = 4.0*R_9; R_5 = R_35+R_5; R_35 = R_19*R_5; R_5 = 6.0*R_5;
            R_5 = R_6+R_5; R_6 = R_17*R_5; R_6 = 7.0*R_6; R_6 = R_20+R_6;
            R_20 = R_31*R_6; R_6 = R_17*R_6; R_6 = 5.0*R_6; R_6 = R_4+R_6;
            R_4 = R_17*R_6; R_4 = 3.0*R_4; R_3 = R_4+R_3; R_3 = R_17*R_3;
            R_4 = R_18*R_3; R_4 = R_25*R_4; R_4 = (1.0/92897280.0)*R_4; R_6 = R_31*R_6;
            R_5 = R_31*R_5; R_9 = R_19*R_9; R_11 = R_9+R_11; R_9 = 4.0*R_11;
            R_9 = R_35+R_9; R_35 = R_31*R_9; R_9 = R_17*R_9; R_9 = 5.0*R_9;
            R_9 = R_5+R_9; R_5 = R_17*R_9; R_5 = 3.0*R_5; R_5 = R_20+R_5;
            R_5 = R_17*R_5; R_5 = R_6+R_5; R_6 = R_18*R_5; R_6 = R_26*R_6;
            R_6 = (1.0/322560.0)*R_6; R_9 = R_31*R_9; R_11 = R_19*R_11; R_13 = R_11+R_13;
            R_11 = R_17*R_13; R_11 = 3.0*R_11; R_11 = R_35+R_11; R_11 = R_17*R_11;
            R_11 = R_9+R_11; R_9 = R_18*R_11; R_9 = R_27*R_9; R_9 = (1.0/1920.0)*R_9;
            R_31 = R_13*R_31; R_15 = R_31+R_15; R_31 = R_18*R_15; R_31 = R_28*R_31;
            R_31 = (1.0/24.0)*R_31; R_0 = R_1+R_0; R_0 = R_19*R_0; R_0 = R_8+R_0;
            R_0 = R_19*R_0; R_0 = R_10+R_0; R_0 = R_19*R_0; R_0 = R_12+R_0;
            R_0 = R_19*R_0; R_0 = R_14+R_0; R_0 = R_19*R_0; R_0 = R_16+R_0;
            R_0 = R_19*R_0; R_0 = R_30+R_0; R_18 = R_18*R_0; R_18 = R_21*R_18;
            R_18 = 0.5*R_18; R_31 = R_18+R_31; R_9 = R_31+R_9; R_6 = R_9+R_6;
            R_4 = R_6+R_4; R_38 = R_4+R_38; R_37 = R_38+R_37; R_36 = R_37+R_36;
            r_y += R_36; R_33 = R_34+R_33; R_32 = R_33*R_32; R_22 = R_32*R_22;
            R_22 = (-1.0/491520.0)*R_22; R_29 = R_33*R_29; R_23 = R_29*R_23; R_23 = (-1.0/25505877196800.0)*R_23;
            R_40 = R_33*R_40; R_24 = R_40*R_24; R_24 = (-1.0/40874803200.0)*R_24; R_3 = R_33*R_3;
            R_25 = R_3*R_25; R_25 = (-1.0/92897280.0)*R_25; R_5 = R_33*R_5; R_26 = R_5*R_26;
            R_26 = (-1.0/322560.0)*R_26; R_11 = R_33*R_11; R_27 = R_11*R_27; R_27 = (-1.0/1920.0)*R_27;
            R_15 = R_33*R_15; R_28 = R_15*R_28; R_28 = (-1.0/24.0)*R_28; R_0 = R_33*R_0;
            R_21 = R_0*R_21; R_21 = -0.5*R_21; R_28 = R_21+R_28; R_27 = R_28+R_27;
            R_26 = R_27+R_26; R_25 = R_26+R_25; R_24 = R_25+R_24; R_23 = R_24+R_23;
            R_22 = R_23+R_22; r_x += R_22;
        }
    };

    auto int_seg = [&]( TF &r_x, TF &r_y, Pt P0, Pt P1 ) {
        std::size_t index = cut_index_r( pow( P0.x, 2 ) + pow( P0.y, 2 ) );
        for( TF u = 0; ; ) {
            std::size_t new_index = index;
            TF new_u = 1;

            // test if line is going to cut a circle at a lower index
            if ( index ) {
                TF d = coeffs[ index - 1 ].first; d = - d; TF b = P0.y; TF R_2 = pow( b, 2 );
                TF a = (-1.0)*b; TF R_4 = P1.y; a += R_4; b *= a;
                a *= a; R_4 = P0.x; TF R_5 = pow(R_4,2); R_2 += R_5;
                d = R_2+d; R_2 = (-1.0)*R_4; R_5 = P1.x;
                R_2 = R_5+R_2; R_4 = R_4*R_2; b += R_4; R_4 = pow(b,2);
                b = (-1.0)*b; R_2 = pow(R_2,2); a = R_2+a;
                d = R_4 - a * d;
                if ( d > 0 ) {
                    TF prop_u = ( b - sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index - 1;
                        new_u = prop_u;
                    }
                    prop_u = ( b + sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index - 1;
                        new_u = prop_u;
                    }
                }
            }

            // test if line is going to cut a circle at an higher index
            if ( index < coeffs.size() - 1 ) {
                TF a; TF b; TF d;

                TF R_0 = coeffs[ index ].first; R_0 = (-1.0)*R_0; TF R_1 = P0.y; TF R_2 = pow(R_1,2);
                TF R_3 = (-1.0)*R_1; TF R_4 = P1.y; R_3 = R_4+R_3; R_1 = R_1*R_3;
                R_3 = pow(R_3,2); R_4 = P0.x; TF R_5 = pow(R_4,2); R_2 = R_5+R_2;
                R_0 = R_2+R_0; R_2 = (-1.0)*R_4; R_5 = P1.x;
                R_2 = R_5+R_2; R_4 = R_4*R_2; R_1 = R_4+R_1; R_4 = pow(R_1,2);
                R_1 = (-1.0)*R_1; b = R_1; R_2 = pow(R_2,2); R_3 = R_2+R_3;
                a = R_3; R_0 = R_3*R_0; R_0 = (-1.0)*R_0; R_0 = R_4+R_0;
                d = R_0;
                if ( d > 0 ) {
                    TF prop_u = ( b - sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index + 1;
                        new_u = prop_u;
                    }
                    prop_u = ( b + sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index + 1;
                        new_u = prop_u;
                    }
                }
            }

            // integration on sub part of the line
            part_int( r_x, r_y, P0, P1, u, new_u, coeffs[ index ].second );

            // next disc
            if ( new_u >= 1 )
                break;
            index = new_index;
            u = new_u;
        }
    };

    auto int_arc = [&]( TF &r_x, TF &r_y, Pt P0, Pt P1 ) {
        TF c = 0;
        const auto &poly_coeffs = coeffs[ cut_index_r( pow( scaling * sphere_radius, 2 ) ) ].second;
        for( std::size_t d = 0; d < poly_coeffs.size(); ++d )
            c += poly_coeffs[ d ] * pow( scaling * sphere_radius, 2 * d + 1 );
        c *= 0.5;

        TF a0 = atan2( P0.y, P0.x );
        TF a1 = atan2( P1.y, P1.x );
        if ( a1 < a0 )
            a1 += 2 * M_PI;
        r_x += c * ( sin( a1 ) - sin( a0 ) );
        r_y += c * ( cos( a0 ) - cos( a1 ) );
    };


    for( std::size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ ) {
        if ( arcs[ i0 ] )
            int_arc( r_x, r_y, scaling * ( point( i0 ) - sphere_center ), scaling * ( point( i1 ) - sphere_center ) );
        else
            int_seg( r_x, r_y, scaling * ( point( i0 ) - sphere_center ), scaling * ( point( i1 ) - sphere_center ) );
    }
}

template<class Pc,class CI> template<class Coeffs>
typename Pc::TF ConvexPolyhedron2<Pc,CI>::_r_polynomials_integration( const Coeffs &coeffs, TF scaling ) const {
    using std::atan2;
    using std::pow;

    auto cut_index_r = [&]( TF r2 ) {
        std::size_t beg = 0;
        std::size_t end = coeffs.size() - 1;
        while ( beg < end ) {
            std::size_t mid = beg + ( end - beg ) / 2;
            if ( coeffs[ mid ].first < r2 )
                beg = mid + 1;
            else
                end = mid;
        }
        return beg;
    };

    auto part_int = []( auto P0, auto P1, auto u0, auto u1, const auto &poly_coeffs ) {
        // generated using metil src/PowerDiagram/offline_integration/lib/gen.met
        if ( poly_coeffs.size() == 4 ) {
            TF result;
            TF R_0 = poly_coeffs[ 1 ]; TF R_1 = poly_coeffs[ 0 ]; TF R_2 = poly_coeffs[ 2 ]; TF R_3 = 2.0*R_2;
            TF R_4 = 4.0*R_2; TF R_5 = 12.0*R_2; TF R_6 = poly_coeffs[ 3 ]; TF R_7 = P0.x;
            TF R_8 = (-1.0)*R_7; TF R_9 = P1.x; TF R_10 = (-1.0)*R_9; R_10 = R_7+R_10;
            R_8 = R_9+R_8; TF R_11 = pow(R_8,2); TF R_12 = P0.y; TF R_13 = (-1.0)*R_12;
            TF R_14 = P1.y; R_13 = R_14+R_13; TF R_15 = pow(R_13,2); R_15 = R_11+R_15;
            R_11 = pow(R_15,2); TF R_16 = R_10*R_13; TF R_17 = (-1.0)*R_14; R_17 = R_17+R_12;
            TF R_18 = R_17*R_8; R_18 = (-1.0)*R_18; R_16 = R_18+R_16; R_18 = u0;
            TF R_19 = (-1.0)*R_18; TF R_20 = u1; R_18 = R_18+R_20; R_14 = R_14*R_18;
            R_14 = 0.5*R_14; R_9 = R_9*R_18; R_9 = 0.5*R_9; R_18 = -0.5*R_18;
            R_18 = 1.0+R_18; R_12 = R_12*R_18; R_12 = R_14+R_12; R_10 = R_10*R_12;
            R_14 = pow(R_12,2); R_12 = R_13*R_12; R_18 = R_7*R_18; R_9 = R_18+R_9;
            R_17 = R_17*R_9; R_17 = (-1.0)*R_17; R_10 = R_17+R_10; R_17 = R_15*R_10;
            R_17 = 720.0*R_17; R_18 = pow(R_9,2); R_14 = R_18+R_14; R_18 = R_6*R_14;
            R_7 = 3.0*R_18; R_7 = R_3+R_7; R_7 = R_14*R_7; R_7 = R_0+R_7;
            R_3 = R_15*R_7; R_13 = 12.0*R_18; R_13 = R_4+R_13; R_2 = R_2+R_18;
            R_2 = R_14*R_2; R_2 = R_0+R_2; R_2 = R_14*R_2; R_2 = R_1+R_2;
            R_2 = R_2*R_10; R_18 = 36.0*R_18; R_18 = R_5+R_18; R_18 = R_15*R_18;
            R_9 = R_8*R_9; R_12 = R_9+R_12; R_9 = R_16*R_12; R_9 = 4320.0*R_9;
            R_9 = R_17+R_9; R_9 = R_11*R_9; R_9 = R_6*R_9; R_7 = R_12*R_7;
            R_7 = R_16*R_7; R_7 = 2.0*R_7; R_11 = pow(R_12,2); R_13 = R_13*R_11;
            R_3 = R_13+R_3; R_3 = R_10*R_3; R_7 = R_3+R_7; R_11 = R_6*R_11;
            R_6 = 144.0*R_11; R_6 = R_18+R_6; R_6 = R_15*R_6; R_6 = R_10*R_6;
            R_6 = (1.0/30.0)*R_6; R_11 = 24.0*R_11; R_18 = R_11+R_18; R_18 = R_12*R_18;
            R_18 = R_16*R_18; R_18 = (2.0/15.0)*R_18; R_6 = R_18+R_6; R_19 = R_20+R_19;
            R_20 = pow(R_19,7); R_20 = R_9*R_20; R_20 = (1.0/322560.0)*R_20; R_9 = pow(R_19,3);
            R_9 = R_7*R_9; R_9 = (1.0/12.0)*R_9; R_2 = R_2*R_19; R_9 = R_2+R_9;
            R_19 = pow(R_19,5); R_6 = R_19*R_6; R_6 = (1.0/32.0)*R_6; R_9 = R_6+R_9;
            R_20 = R_9+R_20; result = R_20;
            return result;
        }
        if ( poly_coeffs.size() == 5 ) {
            TF result;
            TF R_0 = poly_coeffs[ 1 ]; TF R_1 = poly_coeffs[ 0 ]; TF R_2 = poly_coeffs[ 2 ]; TF R_3 = 2.0*R_2;
            TF R_4 = poly_coeffs[ 4 ]; TF R_5 = poly_coeffs[ 3 ]; TF R_6 = 3.0*R_5; TF R_7 = 6.0*R_5;
            TF R_8 = 24.0*R_5; TF R_9 = P0.x; TF R_10 = (-1.0)*R_9; TF R_11 = P1.x;
            TF R_12 = (-1.0)*R_11; R_12 = R_12+R_9; R_10 = R_11+R_10; TF R_13 = pow(R_10,2);
            TF R_14 = P0.y; TF R_15 = (-1.0)*R_14; TF R_16 = P1.y; R_15 = R_16+R_15;
            TF R_17 = pow(R_15,2); R_17 = R_13+R_17; R_13 = pow(R_17,3); TF R_18 = R_12*R_15;
            TF R_19 = (-1.0)*R_16; R_19 = R_19+R_14; TF R_20 = R_19*R_10; R_20 = (-1.0)*R_20;
            R_18 = R_20+R_18; R_20 = u0; TF R_21 = (-1.0)*R_20; TF R_22 = u1;
            R_20 = R_20+R_22; R_16 = R_16*R_20; R_16 = 0.5*R_16; R_11 = R_11*R_20;
            R_11 = 0.5*R_11; R_20 = -0.5*R_20; R_20 = 1.0+R_20; R_14 = R_14*R_20;
            R_14 = R_16+R_14; R_12 = R_12*R_14; R_16 = pow(R_14,2); R_14 = R_15*R_14;
            R_20 = R_9*R_20; R_11 = R_20+R_11; R_19 = R_19*R_11; R_19 = (-1.0)*R_19;
            R_12 = R_19+R_12; R_19 = R_17*R_12; R_20 = 40320.0*R_19; R_9 = pow(R_11,2);
            R_16 = R_9+R_16; R_9 = R_4*R_16; R_5 = R_5+R_9; R_5 = R_16*R_5;
            R_5 = R_2+R_5; R_5 = R_16*R_5; R_5 = R_0+R_5; R_5 = R_16*R_5;
            R_5 = R_1+R_5; R_5 = R_12*R_5; R_1 = 4.0*R_9; R_1 = R_6+R_1;
            R_1 = R_16*R_1; R_1 = R_3+R_1; R_3 = R_16*R_1; R_3 = R_0+R_3;
            R_0 = R_17*R_3; R_1 = 2.0*R_1; R_6 = 16.0*R_9; R_6 = R_7+R_6;
            R_6 = R_16*R_6; R_1 = R_6+R_1; R_6 = R_17*R_1; R_6 = 3.0*R_6;
            R_9 = 96.0*R_9; R_9 = R_8+R_9; R_8 = R_17*R_9; R_16 = 15.0*R_8;
            R_8 = 5.0*R_8; R_11 = R_10*R_11; R_14 = R_11+R_14; R_11 = R_18*R_14;
            R_10 = 322560.0*R_11; R_20 = R_10+R_20; R_20 = R_13*R_20; R_20 = R_4*R_20;
            R_3 = R_14*R_3; R_3 = R_18*R_3; R_3 = 2.0*R_3; R_13 = pow(R_14,2);
            R_1 = R_1*R_13; R_0 = R_1+R_0; R_0 = R_12*R_0; R_3 = R_0+R_3;
            R_4 = R_4*R_13; R_0 = 8640.0*R_4; R_0 = R_16+R_0; R_19 = R_0*R_19;
            R_19 = 2.0*R_19; R_0 = 1920.0*R_4; R_0 = R_16+R_0; R_11 = R_0*R_11;
            R_11 = 12.0*R_11; R_19 = R_11+R_19; R_19 = R_17*R_19; R_4 = 192.0*R_4;
            R_8 = R_4+R_8; R_8 = R_13*R_8; R_9 = R_13*R_9; R_6 = R_9+R_6;
            R_17 = R_17*R_6; R_17 = R_8+R_17; R_17 = R_12*R_17; R_17 = (1.0/30.0)*R_17;
            R_6 = R_14*R_6; R_6 = R_18*R_6; R_6 = (2.0/15.0)*R_6; R_17 = R_6+R_17;
            R_21 = R_22+R_21; R_22 = pow(R_21,9); R_22 = R_20*R_22; R_22 = (1.0/92897280.0)*R_22;
            R_20 = pow(R_21,7); R_20 = R_19*R_20; R_20 = (1.0/322560.0)*R_20; R_19 = pow(R_21,3);
            R_19 = R_3*R_19; R_19 = (1.0/12.0)*R_19; R_5 = R_5*R_21; R_19 = R_5+R_19;
            R_21 = pow(R_21,5); R_17 = R_21*R_17; R_17 = (1.0/32.0)*R_17; R_19 = R_17+R_19;
            R_20 = R_19+R_20; R_22 = R_20+R_22; result = R_22;
            return result;
        }
        if ( poly_coeffs.size() == 6 ) {
            TF result;
            TF R_0 = poly_coeffs[ 1 ]; TF R_1 = poly_coeffs[ 0 ]; TF R_2 = poly_coeffs[ 2 ]; TF R_3 = 2.0*R_2;
            TF R_4 = poly_coeffs[ 3 ]; TF R_5 = 3.0*R_4; TF R_6 = poly_coeffs[ 5 ]; TF R_7 = poly_coeffs[ 4 ];
            TF R_8 = 192.0*R_7; TF R_9 = 4.0*R_7; TF R_10 = 8.0*R_7; TF R_11 = 32.0*R_7;
            TF R_12 = P1.y; TF R_13 = (-1.0)*R_12; TF R_14 = P0.y; TF R_15 = (-1.0)*R_14;
            R_15 = R_12+R_15; TF R_16 = pow(R_15,2); R_13 = R_14+R_13; TF R_17 = P0.x;
            TF R_18 = (-1.0)*R_17; TF R_19 = P1.x; TF R_20 = (-1.0)*R_19; R_20 = R_17+R_20;
            TF R_21 = R_20*R_15; R_18 = R_19+R_18; TF R_22 = pow(R_18,2); R_16 = R_22+R_16;
            R_22 = pow(R_16,4); TF R_23 = pow(R_16,2); TF R_24 = R_18*R_13; R_24 = (-1.0)*R_24;
            R_21 = R_24+R_21; R_24 = u0; TF R_25 = (-1.0)*R_24; TF R_26 = u1;
            R_24 = R_24+R_26; R_12 = R_12*R_24; R_12 = 0.5*R_12; R_19 = R_19*R_24;
            R_19 = 0.5*R_19; R_24 = -0.5*R_24; R_24 = 1.0+R_24; R_14 = R_14*R_24;
            R_14 = R_12+R_14; R_20 = R_20*R_14; R_12 = pow(R_14,2); R_14 = R_15*R_14;
            R_24 = R_17*R_24; R_19 = R_24+R_19; R_13 = R_13*R_19; R_13 = (-1.0)*R_13;
            R_20 = R_13+R_20; R_13 = R_16*R_20; R_24 = 3628800.0*R_13; R_17 = pow(R_19,2);
            R_12 = R_17+R_12; R_17 = R_6*R_12; R_7 = R_7+R_17; R_7 = R_12*R_7;
            R_7 = R_4+R_7; R_7 = R_12*R_7; R_7 = R_2+R_7; R_7 = R_12*R_7;
            R_7 = R_0+R_7; R_7 = R_12*R_7; R_7 = R_1+R_7; R_7 = R_20*R_7;
            R_1 = 960.0*R_17; R_1 = R_8+R_1; R_8 = R_16*R_1; R_2 = 105.0*R_8;
            R_4 = 42.0*R_8; R_8 = 7.0*R_8; R_15 = 5.0*R_17; R_15 = R_9+R_15;
            R_15 = R_12*R_15; R_15 = R_5+R_15; R_5 = R_12*R_15; R_5 = R_3+R_5;
            R_3 = R_12*R_5; R_3 = R_0+R_3; R_0 = R_16*R_3; R_5 = 2.0*R_5;
            R_15 = 2.0*R_15; R_9 = 20.0*R_17; R_9 = R_10+R_9; R_9 = R_12*R_9;
            R_15 = R_9+R_15; R_9 = R_12*R_15; R_5 = R_9+R_5; R_9 = R_16*R_5;
            R_9 = 3.0*R_9; R_15 = 4.0*R_15; R_17 = 120.0*R_17; R_17 = R_11+R_17;
            R_17 = R_12*R_17; R_15 = R_17+R_15; R_17 = R_16*R_15; R_17 = 5.0*R_17;
            R_19 = R_18*R_19; R_14 = R_19+R_14; R_19 = R_21*R_14; R_18 = 36288000.0*R_19;
            R_18 = R_24+R_18; R_18 = R_22*R_18; R_18 = R_6*R_18; R_19 = 2.0*R_19;
            R_19 = R_13+R_19; R_3 = R_14*R_3; R_3 = R_21*R_3; R_3 = 2.0*R_3;
            R_13 = pow(R_14,2); R_6 = R_6*R_13; R_22 = 201600.0*R_6; R_22 = R_2+R_22;
            R_2 = R_21*R_22; R_2 = R_14*R_2; R_2 = 12.0*R_2; R_19 = R_22*R_19;
            R_22 = R_20*R_6; R_22 = R_16*R_22; R_22 = 604800.0*R_22; R_19 = R_22+R_19;
            R_19 = 2.0*R_19; R_2 = R_19+R_2; R_2 = R_23*R_2; R_23 = 28800.0*R_6;
            R_23 = R_4+R_23; R_23 = R_13*R_23; R_6 = 1920.0*R_6; R_8 = R_6+R_8;
            R_8 = R_13*R_8; R_5 = R_5*R_13; R_0 = R_5+R_0; R_0 = R_20*R_0;
            R_3 = R_0+R_3; R_1 = R_13*R_1; R_17 = R_1+R_17; R_1 = R_16*R_17;
            R_1 = 3.0*R_1; R_23 = R_1+R_23; R_23 = R_16*R_23; R_23 = R_20*R_23;
            R_23 = (1.0/1260.0)*R_23; R_1 = R_8+R_1; R_1 = R_14*R_1; R_1 = R_21*R_1;
            R_1 = (1.0/210.0)*R_1; R_23 = R_1+R_23; R_17 = R_13*R_17; R_15 = R_13*R_15;
            R_9 = R_15+R_9; R_16 = R_16*R_9; R_16 = R_17+R_16; R_16 = R_20*R_16;
            R_16 = (1.0/30.0)*R_16; R_9 = R_14*R_9; R_9 = R_21*R_9; R_9 = (2.0/15.0)*R_9;
            R_16 = R_9+R_16; R_25 = R_26+R_25; R_26 = pow(R_25,11); R_26 = R_18*R_26;
            R_26 = (1.0/40874803200.0)*R_26; R_18 = pow(R_25,9); R_18 = R_2*R_18; R_18 = (1.0/92897280.0)*R_18;
            R_2 = pow(R_25,7); R_23 = R_2*R_23; R_23 = (1.0/128.0)*R_23; R_2 = pow(R_25,3);
            R_2 = R_3*R_2; R_2 = (1.0/12.0)*R_2; R_7 = R_7*R_25; R_2 = R_7+R_2;
            R_25 = pow(R_25,5); R_16 = R_25*R_16; R_16 = (1.0/32.0)*R_16; R_2 = R_16+R_2;
            R_23 = R_2+R_23; R_18 = R_23+R_18; R_26 = R_18+R_26; result = R_26;
            return result;
        }
        if ( poly_coeffs.size() == 7 ) {
            TF result;
            TF R_0 = poly_coeffs[ 1 ]; TF R_1 = poly_coeffs[ 0 ]; TF R_2 = poly_coeffs[ 2 ]; TF R_3 = 2.0*R_2;
            TF R_4 = poly_coeffs[ 3 ]; TF R_5 = 3.0*R_4; TF R_6 = poly_coeffs[ 4 ]; TF R_7 = 4.0*R_6;
            TF R_8 = poly_coeffs[ 6 ]; TF R_9 = poly_coeffs[ 5 ]; TF R_10 = 1920.0*R_9; TF R_11 = 240.0*R_9;
            TF R_12 = 40.0*R_9; TF R_13 = 5.0*R_9; TF R_14 = 10.0*R_9; TF R_15 = P1.x;
            TF R_16 = (-1.0)*R_15; TF R_17 = P0.x; TF R_18 = (-1.0)*R_17; R_18 = R_18+R_15;
            TF R_19 = pow(R_18,2); R_16 = R_17+R_16; TF R_20 = P0.y; TF R_21 = (-1.0)*R_20;
            TF R_22 = P1.y; TF R_23 = (-1.0)*R_22; R_23 = R_20+R_23; TF R_24 = R_23*R_18;
            R_24 = (-1.0)*R_24; R_21 = R_22+R_21; TF R_25 = pow(R_21,2); R_19 = R_25+R_19;
            R_25 = pow(R_19,5); TF R_26 = pow(R_19,3); TF R_27 = R_21*R_16; R_24 = R_27+R_24;
            R_27 = u0; TF R_28 = (-1.0)*R_27; TF R_29 = u1; R_27 = R_27+R_29;
            R_15 = R_27*R_15; R_15 = 0.5*R_15; R_22 = R_22*R_27; R_22 = 0.5*R_22;
            R_27 = -0.5*R_27; R_27 = 1.0+R_27; R_17 = R_17*R_27; R_15 = R_17+R_15;
            R_23 = R_23*R_15; R_23 = (-1.0)*R_23; R_17 = pow(R_15,2); R_15 = R_18*R_15;
            R_27 = R_20*R_27; R_22 = R_27+R_22; R_16 = R_16*R_22; R_16 = R_23+R_16;
            R_23 = R_19*R_16; R_27 = 479001600.0*R_23; R_20 = pow(R_22,2); R_17 = R_20+R_17;
            R_20 = R_8*R_17; R_18 = 11520.0*R_20; R_18 = R_10+R_18; R_10 = R_19*R_18;
            TF R_30 = 1890.0*R_10; TF R_31 = 945.0*R_10; TF R_32 = 378.0*R_10; TF R_33 = 63.0*R_10;
            R_10 = 9.0*R_10; R_9 = R_9+R_20; R_9 = R_17*R_9; R_9 = R_6+R_9;
            R_9 = R_17*R_9; R_9 = R_4+R_9; R_9 = R_17*R_9; R_9 = R_2+R_9;
            R_9 = R_17*R_9; R_9 = R_0+R_9; R_9 = R_17*R_9; R_9 = R_1+R_9;
            R_9 = R_16*R_9; R_1 = 1152.0*R_20; R_1 = R_11+R_1; R_1 = R_17*R_1;
            R_11 = 144.0*R_20; R_11 = R_12+R_11; R_11 = R_17*R_11; R_12 = 6.0*R_20;
            R_12 = R_13+R_12; R_12 = R_17*R_12; R_12 = R_7+R_12; R_7 = R_17*R_12;
            R_7 = R_5+R_7; R_5 = R_17*R_7; R_5 = R_3+R_5; R_3 = R_17*R_5;
            R_3 = R_0+R_3; R_0 = R_19*R_3; R_5 = 2.0*R_5; R_7 = 2.0*R_7;
            R_12 = 2.0*R_12; R_20 = 24.0*R_20; R_20 = R_14+R_20; R_20 = R_17*R_20;
            R_12 = R_20+R_12; R_20 = R_17*R_12; R_7 = R_20+R_7; R_20 = R_17*R_7;
            R_5 = R_20+R_5; R_20 = R_19*R_5; R_20 = 3.0*R_20; R_7 = 4.0*R_7;
            R_12 = 4.0*R_12; R_11 = R_12+R_11; R_12 = 6.0*R_11; R_12 = R_1+R_12;
            R_1 = R_19*R_12; R_1 = 7.0*R_1; R_11 = R_17*R_11; R_7 = R_11+R_7;
            R_11 = R_19*R_7; R_11 = 5.0*R_11; R_22 = R_21*R_22; R_15 = R_22+R_15;
            R_22 = R_19*R_15; R_22 = R_8*R_22; R_22 = R_16*R_22; R_22 = 79833600.0*R_22;
            R_21 = R_24*R_15; R_17 = 5748019200.0*R_21; R_17 = R_27+R_17; R_17 = R_25*R_17;
            R_17 = R_8*R_17; R_25 = R_21+R_23; R_3 = R_15*R_3; R_3 = R_24*R_3;
            R_3 = 2.0*R_3; R_27 = pow(R_15,2); R_8 = R_8*R_27; R_14 = 58060800.0*R_8;
            R_14 = R_30+R_14; R_14 = R_24*R_14; R_14 = R_22+R_14; R_14 = R_15*R_14;
            R_22 = 29030400.0*R_8; R_22 = R_31+R_22; R_31 = R_24*R_22; R_31 = R_15*R_31;
            R_31 = 14.0*R_31; R_25 = R_22*R_25; R_14 = R_25+R_14; R_14 = 2.0*R_14;
            R_31 = R_14+R_31; R_31 = R_26*R_31; R_26 = 4769280.0*R_8; R_26 = R_32+R_26;
            R_26 = R_27*R_26; R_32 = 414720.0*R_8; R_33 = R_32+R_33; R_33 = R_27*R_33;
            R_8 = 23040.0*R_8; R_10 = R_8+R_10; R_10 = R_27*R_10; R_18 = R_27*R_18;
            R_1 = R_18+R_1; R_18 = R_19*R_1; R_18 = 5.0*R_18; R_18 = R_10+R_18;
            R_10 = 3.0*R_18; R_26 = R_10+R_26; R_23 = R_26*R_23; R_23 = 2.0*R_23;
            R_33 = R_10+R_33; R_21 = R_33*R_21; R_21 = 16.0*R_21; R_23 = R_21+R_23;
            R_23 = R_19*R_23; R_18 = R_27*R_18; R_1 = R_27*R_1; R_5 = R_5*R_27;
            R_0 = R_5+R_0; R_0 = R_16*R_0; R_3 = R_0+R_3; R_12 = R_27*R_12;
            R_11 = R_12+R_11; R_12 = R_19*R_11; R_12 = 3.0*R_12; R_12 = R_1+R_12;
            R_1 = R_19*R_12; R_1 = R_18+R_1; R_1 = R_16*R_1; R_1 = (1.0/1260.0)*R_1;
            R_12 = R_15*R_12; R_12 = R_24*R_12; R_12 = (1.0/210.0)*R_12; R_1 = R_12+R_1;
            R_11 = R_27*R_11; R_7 = R_27*R_7; R_20 = R_7+R_20; R_19 = R_19*R_20;
            R_19 = R_11+R_19; R_19 = R_16*R_19; R_19 = (1.0/30.0)*R_19; R_20 = R_15*R_20;
            R_20 = R_24*R_20; R_20 = (2.0/15.0)*R_20; R_19 = R_20+R_19; R_28 = R_29+R_28;
            R_29 = pow(R_28,13); R_29 = R_17*R_29; R_29 = (1.0/25505877196800.0)*R_29; R_17 = pow(R_28,11);
            R_17 = R_31*R_17; R_17 = (1.0/40874803200.0)*R_17; R_31 = pow(R_28,9); R_31 = R_23*R_31;
            R_31 = (1.0/92897280.0)*R_31; R_23 = pow(R_28,7); R_1 = R_23*R_1; R_1 = (1.0/128.0)*R_1;
            R_23 = pow(R_28,3); R_23 = R_3*R_23; R_23 = (1.0/12.0)*R_23; R_9 = R_9*R_28;
            R_23 = R_9+R_23; R_28 = pow(R_28,5); R_19 = R_28*R_19; R_19 = (1.0/32.0)*R_19;
            R_23 = R_19+R_23; R_1 = R_23+R_1; R_31 = R_1+R_31; R_17 = R_31+R_17;
            R_29 = R_17+R_29; result = R_29;
            return result;
        }
        if ( poly_coeffs.size() == 8 ) {
            TF result;
            TF R_0 = poly_coeffs[ 1 ]; TF R_1 = poly_coeffs[ 0 ]; TF R_2 = poly_coeffs[ 2 ]; TF R_3 = 2.0*R_2;
            TF R_4 = poly_coeffs[ 3 ]; TF R_5 = 3.0*R_4; TF R_6 = poly_coeffs[ 4 ]; TF R_7 = 4.0*R_6;
            TF R_8 = poly_coeffs[ 5 ]; TF R_9 = 5.0*R_8; TF R_10 = poly_coeffs[ 7 ]; TF R_11 = poly_coeffs[ 6 ];
            TF R_12 = 23040.0*R_11; TF R_13 = 2304.0*R_11; TF R_14 = 288.0*R_11; TF R_15 = 6.0*R_11;
            TF R_16 = 12.0*R_11; TF R_17 = 48.0*R_11; TF R_18 = P0.x; TF R_19 = (-1.0)*R_18;
            TF R_20 = P1.x; TF R_21 = (-1.0)*R_20; R_21 = R_18+R_21; R_19 = R_20+R_19;
            TF R_22 = pow(R_19,2); TF R_23 = P1.y; TF R_24 = (-1.0)*R_23; TF R_25 = P0.y;
            TF R_26 = (-1.0)*R_25; R_26 = R_26+R_23; TF R_27 = pow(R_26,2); R_27 = R_22+R_27;
            R_22 = pow(R_27,6); TF R_28 = pow(R_27,4); TF R_29 = pow(R_27,2); TF R_30 = R_21*R_26;
            R_24 = R_25+R_24; TF R_31 = R_24*R_19; R_31 = (-1.0)*R_31; R_30 = R_31+R_30;
            R_31 = u0; TF R_32 = (-1.0)*R_31; TF R_33 = u1; R_31 = R_31+R_33;
            R_23 = R_23*R_31; R_23 = 0.5*R_23; R_20 = R_20*R_31; R_20 = 0.5*R_20;
            R_31 = -0.5*R_31; R_31 = 1.0+R_31; R_25 = R_25*R_31; R_23 = R_25+R_23;
            R_21 = R_21*R_23; R_25 = pow(R_23,2); R_23 = R_26*R_23; R_31 = R_18*R_31;
            R_20 = R_31+R_20; R_24 = R_24*R_20; R_24 = (-1.0)*R_24; R_21 = R_24+R_21;
            R_24 = R_27*R_21; R_31 = 87178291200.0*R_24; R_18 = 3.0*R_24; R_26 = pow(R_20,2);
            R_25 = R_26+R_25; R_26 = R_10*R_25; TF R_34 = 161280.0*R_26; R_34 = R_12+R_34;
            R_12 = R_27*R_34; TF R_35 = 10395.0*R_12; TF R_36 = 3465.0*R_12; TF R_37 = 1155.0*R_12;
            TF R_38 = 165.0*R_12; R_12 = 11.0*R_12; TF R_39 = 13440.0*R_26; R_39 = R_13+R_39;
            R_39 = R_25*R_39; R_11 = R_11+R_26; R_11 = R_25*R_11; R_11 = R_8+R_11;
            R_11 = R_25*R_11; R_11 = R_6+R_11; R_11 = R_25*R_11; R_11 = R_4+R_11;
            R_11 = R_25*R_11; R_11 = R_2+R_11; R_11 = R_25*R_11; R_11 = R_0+R_11;
            R_11 = R_25*R_11; R_11 = R_1+R_11; R_11 = R_21*R_11; R_1 = 1344.0*R_26;
            R_1 = R_14+R_1; R_1 = R_25*R_1; R_14 = 7.0*R_26; R_14 = R_15+R_14;
            R_14 = R_25*R_14; R_14 = R_9+R_14; R_9 = R_25*R_14; R_9 = R_7+R_9;
            R_7 = R_25*R_9; R_7 = R_5+R_7; R_5 = R_25*R_7; R_5 = R_3+R_5;
            R_3 = R_25*R_5; R_3 = R_0+R_3; R_0 = R_27*R_3; R_5 = 2.0*R_5;
            R_7 = 2.0*R_7; R_9 = 2.0*R_9; R_14 = 2.0*R_14; R_15 = 28.0*R_26;
            R_15 = R_16+R_15; R_15 = R_25*R_15; R_14 = R_15+R_14; R_15 = R_25*R_14;
            R_9 = R_15+R_9; R_15 = R_25*R_9; R_7 = R_15+R_7; R_15 = R_25*R_7;
            R_5 = R_15+R_5; R_15 = R_27*R_5; R_15 = 3.0*R_15; R_7 = 4.0*R_7;
            R_9 = 4.0*R_9; R_14 = 4.0*R_14; R_26 = 168.0*R_26; R_26 = R_17+R_26;
            R_26 = R_25*R_26; R_14 = R_26+R_14; R_26 = 6.0*R_14; R_26 = R_1+R_26;
            R_1 = 8.0*R_26; R_1 = R_39+R_1; R_39 = R_27*R_1; R_39 = 9.0*R_39;
            R_26 = R_25*R_26; R_14 = R_25*R_14; R_9 = R_14+R_9; R_14 = 6.0*R_9;
            R_14 = R_26+R_14; R_26 = R_27*R_14; R_26 = 7.0*R_26; R_9 = R_25*R_9;
            R_7 = R_9+R_7; R_9 = R_27*R_7; R_9 = 5.0*R_9; R_20 = R_19*R_20;
            R_23 = R_20+R_23; R_20 = R_10*R_23; R_19 = R_30*R_23; R_25 = 1220496076800.0*R_19;
            R_25 = R_31+R_25; R_25 = R_22*R_25; R_25 = R_10*R_25; R_18 = R_19+R_18;
            R_18 = R_23*R_18; R_18 = R_10*R_18; R_18 = 2905943040.0*R_18; R_22 = R_19+R_24;
            R_20 = R_22*R_20; R_20 = 2905943040.0*R_20; R_31 = 2.0*R_19; R_24 = R_31+R_24;
            R_3 = R_3*R_23; R_3 = R_30*R_3; R_3 = 2.0*R_3; R_31 = pow(R_23,2);
            R_17 = R_21*R_31; R_17 = R_27*R_17; R_10 = R_10*R_31; R_16 = 5588352000.0*R_10;
            R_16 = R_35+R_16; R_16 = R_19*R_16; R_16 = 16.0*R_16; R_19 = 2682408960.0*R_10;
            R_19 = R_35+R_19; R_19 = R_30*R_19; R_19 = R_20+R_19; R_19 = 3.0*R_19;
            R_18 = R_19+R_18; R_18 = R_23*R_18; R_19 = 894136320.0*R_10; R_19 = R_36+R_19;
            R_22 = R_19*R_22; R_22 = 3.0*R_22; R_18 = R_22+R_18; R_18 = 2.0*R_18;
            R_16 = R_18+R_16; R_16 = R_28*R_16; R_19 = R_17*R_19; R_17 = 121927680.0*R_10;
            R_17 = R_37+R_17; R_17 = R_31*R_17; R_37 = 9031680.0*R_10; R_38 = R_37+R_38;
            R_38 = R_31*R_38; R_10 = 322560.0*R_10; R_12 = R_10+R_12; R_12 = R_31*R_12;
            R_34 = R_31*R_34; R_39 = R_34+R_39; R_34 = R_27*R_39; R_10 = 105.0*R_34;
            R_17 = R_10+R_17; R_10 = R_30*R_17; R_10 = R_23*R_10; R_10 = 16.0*R_10;
            R_24 = R_17*R_24; R_19 = R_24+R_19; R_19 = 2.0*R_19; R_10 = R_19+R_10;
            R_10 = R_29*R_10; R_29 = 42.0*R_34; R_38 = R_29+R_38; R_38 = R_31*R_38;
            R_34 = 7.0*R_34; R_34 = R_12+R_34; R_34 = R_31*R_34; R_39 = R_31*R_39;
            R_1 = R_31*R_1; R_26 = R_1+R_26; R_1 = R_27*R_26; R_1 = 5.0*R_1;
            R_1 = R_39+R_1; R_39 = R_27*R_1; R_39 = 3.0*R_39; R_38 = R_39+R_38;
            R_38 = R_27*R_38; R_38 = R_21*R_38; R_38 = (1.0/90720.0)*R_38; R_39 = R_34+R_39;
            R_39 = R_23*R_39; R_39 = R_30*R_39; R_39 = (1.0/11340.0)*R_39; R_38 = R_39+R_38;
            R_1 = R_31*R_1; R_26 = R_31*R_26; R_5 = R_5*R_31; R_0 = R_5+R_0;
            R_0 = R_21*R_0; R_3 = R_0+R_3; R_14 = R_31*R_14; R_9 = R_14+R_9;
            R_14 = R_27*R_9; R_14 = 3.0*R_14; R_14 = R_26+R_14; R_26 = R_27*R_14;
            R_26 = R_1+R_26; R_26 = R_21*R_26; R_26 = (1.0/1260.0)*R_26; R_14 = R_23*R_14;
            R_14 = R_30*R_14; R_14 = (1.0/210.0)*R_14; R_26 = R_14+R_26; R_9 = R_31*R_9;
            R_7 = R_31*R_7; R_15 = R_7+R_15; R_27 = R_27*R_15; R_27 = R_9+R_27;
            R_27 = R_21*R_27; R_27 = (1.0/30.0)*R_27; R_15 = R_23*R_15; R_15 = R_30*R_15;
            R_15 = (2.0/15.0)*R_15; R_27 = R_15+R_27; R_32 = R_33+R_32; R_33 = pow(R_32,15);
            R_33 = R_25*R_33; R_33 = (1.0/21424936845312000.0)*R_33; R_25 = pow(R_32,13); R_25 = R_16*R_25;
            R_25 = (1.0/25505877196800.0)*R_25; R_16 = pow(R_32,11); R_16 = R_10*R_16; R_16 = (1.0/40874803200.0)*R_16;
            R_10 = pow(R_32,9); R_38 = R_10*R_38; R_38 = (1.0/512.0)*R_38; R_10 = pow(R_32,7);
            R_26 = R_10*R_26; R_26 = (1.0/128.0)*R_26; R_10 = pow(R_32,3); R_10 = R_3*R_10;
            R_10 = (1.0/12.0)*R_10; R_11 = R_11*R_32; R_10 = R_11+R_10; R_32 = pow(R_32,5);
            R_27 = R_32*R_27; R_27 = (1.0/32.0)*R_27; R_10 = R_27+R_10; R_26 = R_10+R_26;
            R_38 = R_26+R_38; R_16 = R_38+R_16; R_25 = R_16+R_25; R_33 = R_25+R_33;
            result = R_33;
            return result;
        }
        TODO;
        return TF( 0 );
    };

    // TODO: start from the old u
    auto int_seg = [&]( Pt P0, Pt P1 ) {
        std::size_t index = cut_index_r( pow( P0.x, 2 ) + pow( P0.y, 2 ) );
        TF result = 0;
        for( TF u = 0; ; ) {
            std::size_t new_index = index;
            TF new_u = 1;

            // test if line is going to cut a circle at a lower index
            if ( index ) {
                TF d = coeffs[ index - 1 ].first; d = - d; TF b = P0.y; TF R_2 = pow( b, 2 );
                TF a = (-1.0)*b; TF R_4 = P1.y; a += R_4; b *= a;
                a *= a; R_4 = P0.x; TF R_5 = pow(R_4,2); R_2 += R_5;
                d = R_2+d; R_2 = (-1.0)*R_4; R_5 = P1.x;
                R_2 = R_5+R_2; R_4 = R_4*R_2; b += R_4; R_4 = pow(b,2);
                b = (-1.0)*b; R_2 = pow(R_2,2); a = R_2+a;
                d = R_4 - a * d;
                if ( d > 0 ) {
                    TF prop_u = ( b - sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index - 1;
                        new_u = prop_u;
                    }
                    prop_u = ( b + sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index - 1;
                        new_u = prop_u;
                    }
                }
            }

            // test if line is going to cut a circle at an higher index
            if ( index < coeffs.size() - 1 ) {
                TF a; TF b; TF d;

                TF R_0 = coeffs[ index ].first; R_0 = (-1.0)*R_0; TF R_1 = P0.y; TF R_2 = pow(R_1,2);
                TF R_3 = (-1.0)*R_1; TF R_4 = P1.y; R_3 = R_4+R_3; R_1 = R_1*R_3;
                R_3 = pow(R_3,2); R_4 = P0.x; TF R_5 = pow(R_4,2); R_2 = R_5+R_2;
                R_0 = R_2+R_0; R_2 = (-1.0)*R_4; R_5 = P1.x;
                R_2 = R_5+R_2; R_4 = R_4*R_2; R_1 = R_4+R_1; R_4 = pow(R_1,2);
                R_1 = (-1.0)*R_1; b = R_1; R_2 = pow(R_2,2); R_3 = R_2+R_3;
                a = R_3; R_0 = R_3*R_0; R_0 = (-1.0)*R_0; R_0 = R_4+R_0;
                d = R_0;
                if ( d > 0 ) {
                    TF prop_u = ( b - sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index + 1;
                        new_u = prop_u;
                    }
                    prop_u = ( b + sqrt( d ) ) / a;
                    if ( prop_u > u && prop_u < new_u ) {
                        new_index = index + 1;
                        new_u = prop_u;
                    }
                }
            }

            // integration on sub part of the line
            result += part_int( P0, P1, u, new_u, coeffs[ index ].second );

            // next disc
            if ( new_u >= 1 )
                break;
            index = new_index;
            u = new_u;
        }
        return result;
    };

    auto int_arc = [&]( Pt P0, Pt P1 ) {
        TF res = 0, r2 = pow( scaling * sphere_radius, 2 );
        std::size_t index = cut_index_r( r2 );
        const auto &poly_coeffs = coeffs[ index ].second;
        for( std::size_t d = 0; d < poly_coeffs.size(); ++d )
            res += poly_coeffs[ d ] * pow( r2, d + 1 );

        TF a0 = atan2( P0.y, P0.x );
        TF a1 = atan2( P1.y, P1.x );
        if ( a1 < a0 )
            a1 += 2 * M_PI;
        return ( a1 - a0 ) * res;
    };

    if ( nb_points == 0 ) {
        TF res = 0, r2 = pow( scaling * sphere_radius, 2 );
        const auto &poly_coeffs = coeffs[ cut_index_r( r2 ) ].second;
        for( std::size_t d = 0; d < poly_coeffs.size(); ++d )
            res += poly_coeffs[ d ] * pow( r2, d + 1 );
        return 2 * M_PI * res;
    }

    TF res = 0;
    for( std::size_t i1 = 0, i0 = nb_points - 1; i1 < nb_points; i0 = i1++ )
        if ( arcs[ i0 ] )
            res += int_arc( scaling * ( point( i0 ) - sphere_center ), scaling * ( point( i1 ) - sphere_center ) );
        else
            res += int_seg( scaling * ( point( i0 ) - sphere_center ), scaling * ( point( i1 ) - sphere_center ) );
    return res;
}

} // namespace PowerDiagram
