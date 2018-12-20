#include "PolyDivIntegrator.h"

template<class TF, int degp>
std::size_t PolyDivIntegrator<TF,degp>::cut_index( TF pos ) const {
    std::size_t beg = 0;
    std::size_t end = cuts.size() - 1;
    while ( beg < end ) {
        std::size_t mid = beg + ( end - beg ) / 2;
        if ( cuts[ mid ].pos < pos )
            beg = mid + 1;
        else
            end = mid;
    }
    return end;
}
