#pragma once

//#include "SpRegularGrid.h"
#include "TraversalByZIndexGrid.h"

/**
  Example of a class that can be used as a parameter for PowerDiagram
*/
template<int _dim,class _TF=double,class _Sp=TraversalByZIndexGrid>
struct PcStd {
    static constexpr int dim  = _dim;        ///< dimension of the space
    using                TF   = _TF;         ///< floating point type
    using                Sp   = _Sp;         ///< Space partitioner
    using                TI   = std::size_t; ///< integer used for the indices

    /* ctor */           PcStd( Sp sp = {} ) : sp( sp ) {}

    Sp                   sp;                 ///< instance of the space partitioner
};
