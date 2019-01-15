#pragma once

#include "Point2.h"
#include "Point3.h"

/**
*/
namespace FunctionEnum {

struct InBallW05 { template<class PT,class TF> auto operator()( PT p, PT c, TF w ) const { return norm_2_p2( p - c ) <= w; } };
struct Gaussian  { template<class PT,class TF> auto operator()( PT p, PT c, TF w ) const { using std::exp; return exp( - norm_2_p2( p - c ) ); } };
struct Unit      { template<class PT,class TF> auto operator()( PT p, PT c, TF w ) const { return 1; } };
struct R2        { template<class PT,class TF> auto operator()( PT p, PT c, TF w ) const { return norm_2_p2( p - c ); } };

//
template<class T> T    func_for_final_cp_integration( T t       ) { return t ; }
inline            Unit func_for_final_cp_integration( InBallW05 ) { return {}; }

//
template<class TV,class Pt,class TF,class TI,class Func>
void set_target_measure( TV &v_values, const Pt *positions, const TF *weights, TI nb_diracs, Func func ) {
    for( TF &v : v_values )
        v = - TF( 1 ) / nb_diracs;
}

template<class TV,class Pt,class TF,class TI>
void set_target_measure( TV &v_values, const Pt *positions, const TF *weights, TI nb_diracs, InBallW05 ) {
    for( TF &v : v_values )
        v = - TF( 1 );
}

//
template<class T> N<0> need_ball_cut( T ) { return {}; }
inline N<1> need_ball_cut( InBallW05 ) { return {}; }

}
