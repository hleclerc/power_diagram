#pragma once

#include <vector>
#include <array>

/**
*/
template<class TF,int degp>
struct PolyDivIntegrator {
    struct Cut {
        TF                    pos;
        std::array<TF,degp+1> coeffs;
    };

    std::size_t               cut_index( TF pos ) const;

    std::vector<Cut>          cuts;
};

#include "PolyDivIntegrator.tcc"
