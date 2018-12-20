#pragma once

#include "Point2.h"
#include "Point3.h"

/**
*/
namespace FunctionEnum {

struct Gaussian  { template<class PT> auto operator()( PT p, PT c ) const { using std::exp; return exp( - norm_2_p2( p - c ) ); } };
struct Unit      { template<class PT> auto operator()( PT p, PT c ) const { return 1; } };
struct R2        { template<class PT> auto operator()( PT p, PT c ) const { return norm_2_p2( p - c ); } };

}
