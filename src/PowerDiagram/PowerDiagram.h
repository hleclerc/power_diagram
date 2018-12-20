#pragma once

#include "ConvexPolyhedron3.h"
#include "ConvexPolyhedron2.h"
#include <functional>
#include "PcStd.h"

/**
  Pc = Carac class to define the space dimension, base scalar type, ...

  @see PowerDiagramCaracStd for an example
*/
template<class Pc = PcStd<2> >
class PowerDiagram {
public:
    using                           TF                   = typename Pc::TF;
    using                           TI                   = typename Pc::TI;
    using                           Sp                   = SpacePartioner<Pc,typename Pc::Sp>;
    static constexpr int            dim                  = Pc::dim;

    using                           PT                   = typename std::conditional<dim==3,Point3<TF>,Point2<TF>>::type;
    struct                          Dirac                { PT pos; TF weight; /* kantorovitch potential*/ void write_to_stream( std::ostream &os ) const { os << "[" << pos << "," << weight << "]"; } };
    struct                          LC_info              { Dirac dirac; std::size_t index; };
    using                           LC                   = typename std::conditional<dim==3,ConvexPolyhedron3<TF,TI,LC_info>,ConvexPolyhedron2<TF,TI,LC_info>>::type;
    using                           FC                   = std::function<void( TF *, const LC &, TI )>; ///< used in display
    struct                          CutInfo              { PT O, N; };

    /* */                           PowerDiagram         ( const Pc &power_diagram_carac = {} );

    // modifications
    void                            add_convex_shape     ( const std::vector<CutInfo> &convex_bounds, TF orig_radius = 1e40 );
    void                            add_convex_shape     ( const LC &lc );

    void                            add_bounding_simplex ( PT center, TF radius );
    void                            add_box_shape        ( PT p0, PT p1 );

    void                            add_dirac            ( PT pos, TF weight );

    void                            set_ball_cut         ( bool ball_cut ); ///< true to cut laguerre cells with disc of radii |kantorovith_potential|^0.5 during the construction

    // information
    void                            write_to_stream      ( std::ostream &os ) const;
    template<class F> void          for_each_cell        ( const F &f ) const; ///< std::function<void( const LC &cs, TI num, int num_thread )>
    // void                         for_each_cell        ( const std::function<void( const LC &cs, TI num, int num_thread )> &f ) const; ///<

    template<int c> void            display_bounds       ( VtkOutput<c,TF> &vo ) const;
    template<int c> void            display              ( VtkOutput<c,TF> &vo, const FC &cell_data_func = []( TF *v, const LC &, TI ind ){ if ( c ) v[ 0 ] = ind; }, bool filled = true, TF max_ratio_area_error = 1e-1 ) const;

    TI                              nb_convex_bounds     () const;
    TF                              domain_measure       () const;
    TI                              nb_diracs            () const { return _diracs.size(); }

    const Dirac&                    dirac                ( TI i ) const { return _diracs[ i ]; }
    Dirac&                          dirac                ( TI i ) { return _diracs[ i ]; }

    // computations
    template<class Fu> void         get_centroid_contribs( const Fu &fu, PT *centroid_mul_by_volumes, TF *volumes ) const;
    template<class Fu> void         get_der_wrt_dir_pos  ( const Fu &fu, std::vector<std::pair<TI,TF>> *derivatives ) const; ///< get derivatives / initial position
    template<class Fu> void         get_der_measures     ( const Fu &fu, TF *measures, std::vector<std::pair<TI,TF>> *derivatives ) const;
    template<class Fu> void         get_measures         ( const Fu &fu, TF *measures ) const;

    void                            get_centroid_contribs( PT *centroid_mul_by_volumes, TF *volumes )                 const { get_centroid_contribs( FunctionEnum::Unit(), centroid_mul_by_volumes, volumes ); }
    void                            get_der_wrt_dir_pos  ( std::vector<std::pair<TI,TF>> *derivatives )               const { get_der_wrt_dir_pos  ( FunctionEnum::Unit(), derivatives                      ); } ///< get derivatives / initial position
    void                            get_der_measures     ( TF *measures, std::vector<std::pair<TI,TF>> *derivatives ) const { get_der_measures     ( FunctionEnum::Unit(), measures, derivatives            ); }
    void                            get_measures         ( TF *measures )                                             const { get_measures         ( FunctionEnum::Unit(), measures                         ); }

    // approximate computations
    template<class Fu> void         get_der_measures_ap  ( const Fu &fu, TF *volumes, std::vector<std::pair<TI,TF>> *derivatives, TF epsilon = 1e-6 ) const;
    template<class Fu> void         get_centroids_ap     ( const Fu &fu, PT *centroids, TI n = 1e7 )                                                  const;
    template<class Fu> void         get_measures_ap      ( const Fu &fu, TF *volumes, TI n = 1e7 )                                                    const;

    void                            get_der_measures_ap  ( TF *volumes, std::vector<std::pair<TI,TF>> *derivatives, TF epsilon = 1e-6 ) const { get_der_measures_ap( FunctionEnum::Unit(), volumes, derivatives, epsilon ); }
    void                            get_centroids_ap     ( PT *centroids, TI n = 1e7 )                                                  const { get_centroids_ap   ( FunctionEnum::Unit(), centroids           , n       ); }
    void                            get_measures_ap      ( TF *volumes, TI n = 1e7 )                                                    const { get_measures_ap     ( FunctionEnum::Unit(), volumes             , n       ); }

private:
    using                           Trans                = std::array<TF,dim*dim>;

    void                            _add_box_shape       ( PT p0, PT p1, N<2> );
    void                            _add_box_shape       ( PT p0, PT p1, N<3> );
    static void                     _sort_and_sum        ( std::vector<std::pair<TI,TF>> &dv );

    Pc                              _power_diagram_carac;
    mutable Sp                      _space_partitioner;
    std::vector<Trans>              _tranformations;
    std::vector<LC>                 _convex_bounds;
    bool                            _ball_cut;
    std::vector<Dirac>              _diracs;
};

#include "PowerDiagram.tcc"
