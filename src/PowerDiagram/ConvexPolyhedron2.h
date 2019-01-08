#pragma once

#include "SpaceFunctions/Constant.h"
#include "FunctionEnum.h"
#include "VtkOutput.h"
#include <functional>
#include <algorithm>
#include <bitset>

#define XSIMD_ENABLE_FALLBACK
#include <xsimd/xsimd.hpp>

namespace PowerDiagram {

/**
  Pc must contain
    - dim (2, 3, ...)
    - TI (std::size_t, ...) => index type
    - TF (double, ...) => floating point type
    - allow_ball_cut (true, false)
  CI = Cut id type

  Beware: ball_cuts must be done AFTER the plane_cuts.
*/
template<class Pc,class CI=std::size_t>
class ConvexPolyhedron2 {
public:
    static constexpr bool     allow_ball_cut            = Pc::allow_ball_cut;
    using                     TI                        = typename Pc::TI; ///< index type
    using                     TF                        = typename Pc::TF; ///< index type
    using                     Pt                        = Point2<TF>;  ///< 3D point

    // types for simd
    static constexpr TI       simd_size                 = xsimd::simd_traits<TF>::size;
    using                     VF                        = std::vector<double,xsimd::aligned_allocator<double,XSIMD_DEFAULT_ALIGNMENT>>;
    using                     BF                        = xsimd::batch<TF,simd_size>;
    using                     AF                        = std::array<TF,64>;
    using                     AC                        = std::array<CI,64>;
    using                     AB                        = std::bitset<64>;

    // types for the ctor
    struct                    EnglobingSimplex          { Pt p; TF r; };
    struct                    Box                       { Pt p0, p1; };

    /// we start from a triangle that includes the circle defined by englobing_center and englobing_radius (but this sphere is not used, it's just here to construct a triangle)
    /**/                      ConvexPolyhedron2         ( const EnglobingSimplex &es, CI cut_id = {} );
    /**/                      ConvexPolyhedron2         ( const Box &box, CI cut_id = {} );
    /**/                      ConvexPolyhedron2         ( const ConvexPolyhedron2 &that );
    /**/                      ConvexPolyhedron2         ();

    // information
    void                      for_each_boundary_measure ( FunctionEnum::Unit, const std::function<void( TF boundary_measure, CI id )> &f ) const;
    void                      for_each_approx_seg       ( const std::function<void( Pt )> &f, TF max_ratio_area_error = 1e-1 ) const; ///<
    void                      for_each_simplex          ( const std::function<void( CI num_0, CI num_1 )> &f ) const;
    void                      for_each_bound            ( const std::function<void( Pt p0, Pt p1, CI id )> &f ) const;
    void                      for_each_node             ( const std::function<void( Pt v )> &f ) const;

    void                      write_to_stream           ( std::ostream &os ) const;
    template<class V> void    display                   ( V &vo, const typename V::CV &cell_data = {}, bool filled = true, TF max_ratio_area_error = 1e-1, bool display_tangents = false ) const;
    Pt                        normal                    ( std::size_t n ) const { return { normals[ 0 ][ n ], normals[ 1 ][ n ] }; }
    Pt                        point                     ( std::size_t n ) const { return { points[ 0 ][ n ], points[ 1 ][ n ] }; }
    bool                      empty                     () const { return nb_points == 0 && sphere_radius <= 0; }

    // modifications
    void                      set_cut_ids               ( CI cut_id ); ///< replace all the cut_ids
    void                      plane_cut                 ( Pt origin, Pt normal, CI cut_id = {} );
    void                      ball_cut                  ( Pt center, TF radius, CI cut_id = {} ); ///< beware: only one sphere cut is authorized, and it must be done after all the plane cuts.
    void                      clear                     ( Pt englobing_center, TF englobing_radius, CI englobing_cut_id = {} );

    // tests
    bool                      is_a_cutting_plane        ( Pt origin, Pt normal ) const;

    // computations
    void                      add_centroid_contrib      ( Pt &ctd, TF &vol, FunctionEnum::Gaussian ) const;
    void                      add_centroid_contrib      ( Pt &ctd, TF &vol, FunctionEnum::Unit     ) const;
    void                      add_centroid_contrib      ( Pt &ctd, TF &vol, FunctionEnum::R2       ) const;
    void                      add_centroid_contrib      ( Pt &ctd, TF &vol ) const;

    TF                        boundary_measure          ( FunctionEnum::Gaussian ) const;
    TF                        boundary_measure          ( FunctionEnum::Unit     ) const;
    TF                        boundary_measure          () const;

    template<class FU> Pt     centroid                  ( const FU &f ) const;
    Pt                        centroid                  () const;

    template<class FU> TF     measure                   ( const FU &f ) const;
    TF                        measure                   () const;


    TF                        integration               ( FunctionEnum::Gaussian ) const;
    TF                        integration               ( FunctionEnum::Unit     ) const;
    TF                        integration               ( FunctionEnum::R2       ) const;

    TF                        integration               ( SpaceFunctions::Constant<TF> cst ) const;

    // approximate computations
    TF                        boundary_measure_ap       ( TF max_ratio_area_error = 1e-4 ) const; ///<
    template<class Fu> TF     integration_ap            ( const Fu &func, std::size_t n = 1e6 ) const;
    template<class Fu> Pt     centroid_ap               ( const Fu &func, std::size_t n = 1e6 ) const; ///<
    TF                        measure_ap                ( TF max_ratio_area_error = 1e-4 ) const; ///<

    // attributes
    alignas(64) AF            normals[ 2 ];
    AF                        points[ 2 ];
    AC                        cut_ids;
    AB                        arcs;

    Pt                        sphere_center;
    TF                        sphere_radius;
    CI                        sphere_cut_id;

    std::size_t               nb_points;

private:
    enum                      CutType                   { LINE = 0, ARC = 1 };
    struct                    Cut                       { int cut_type; CI cut_id; Pt normal; Pt point; };
    template<class Coeffs> TF _r_polynomials_integration( const Coeffs &coeffs ) const;
    template<class Coef> void _r_centroid_integration   ( TF &r_x, TF &r_y, const Coef &coeffs ) const;
    void                      _centroid_arc             ( Pt &ctd, TF &mea, Pt p0, Pt p1 ) const;
    TF                        _arc_length               ( Pt p0, Pt p1 ) const;
    TF                        _arc_area                 ( Pt p0, Pt p1 ) const;

    std::vector<Cut>          _tmp_cuts;
};

} // namespace PowerDiagram

#include "ConvexPolyhedron2.tcc"

