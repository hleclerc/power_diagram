#pragma once

#include "Point2.h"
#include "Point3.h"

/**
*/
namespace FunctionEnum {

template<class TS>
struct ExpWmR2db {
    template<class PT,class TF>
    auto operator()( PT p, PT c, TF w ) const {
        using std::exp; return exp( ( w - norm_2_p2( p - c ) ) / eps );
    }

    const char *name() const {
        return "ExpWmR2db";
    }

    auto func_for_final_cp_integration() const {
        return *this;
    }

    N<0> need_ball_cut() const {
        return {};
    }

    TS eps;
};

struct Unit {
    template<class PT,class TF>
    auto operator()( PT p, PT c, TF w ) const {
        return 1;
    }

    const char *name() const {
        return "1";
    }

    auto func_for_final_cp_integration() const {
        return *this;
    }

    N<0> need_ball_cut() const {
        return {};
    }
};

struct R2 {
    template<class PT,class TF>
    auto operator()( PT p, PT c, TF w ) const {
        return norm_2_p2( p - c );
    }

    const char *name() const {
        return "R2";
    }

    auto func_for_final_cp_integration() const {
        return *this;
    }

    N<0> need_ball_cut() const {
        return {};
    }
};

struct InBallW05 {
    template<class PT,class TF>
    auto operator()( PT p, PT c, TF w ) const {
        return norm_2_p2( p - c ) <= w;
    }

    const char *name() const {
        return "InBallW05";
    }

    auto func_for_final_cp_integration() const {
        return Unit{};
    }

    N<1> need_ball_cut() const {
        return {};
    }
};

}
