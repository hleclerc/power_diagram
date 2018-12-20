#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "PowerDiagram.h"
#include <memory>

//// nsmake cpp_flag -march=native
//// nsmake cpp_flag -Wattributes
//// nsmake cpp_flag -O3
//// nsmake cpp_flag -g3

namespace py = pybind11;

template<int dim,class TF>
struct PyMeasuresAndDerivatives {
    py::array_t<TF>     measures;
    // CSR type information
    py::array_t<TF>     der_data;
    py::array_t<size_t> der_indptr;
    py::array_t<size_t> der_indices;
};

template<int dim,class TF>
struct PyDirac {
    TF x;
    TF y;
    TF z;
    TF r;
};

//
template<int dim,class TF>
class PyPowerDiagram {
public:
    using Pc = PcStd<dim,TF>;
    using PO = typename PowerDiagram<Pc>::PT;

    std::string __repr__() const {
        return ___repr__( N<dim>() );
    }

    std::string ___repr__( N<2> ) const {
        std::ostringstream ss;
        for( size_t i = 0; i < pd.nb_diracs(); ++i )
            ss << "[x=" << pd.dirac( i ).pos.x << ",y=" << pd.dirac( i ).pos.y << ",r=" << pd.dirac( i ).weight << "]";
        return ss.str();
    }

    std::string ___repr__( N<3> ) const {
        std::ostringstream ss;
        for( size_t i = 0; i < pd.nb_diracs(); ++i )
            ss << "[x=" << pd.dirac( i ).pos.x << ",y=" << pd.dirac( i ).pos.y << ",z=" << pd.dirac( i ).pos.z << ",r=" << pd.dirac( i ).weight << "]";
        return ss.str();
    }

    void add_dirac( py::array_t<TF> &pos, TF weight ) {
        auto buf_pos = pos.request();
        auto ptr_pos = (const TF *)buf_pos.ptr;
        ASSERT( pos.size() == dim, "" );

        PO p;
        for( std::size_t d = 0; d < dim; ++d )
            p[ d ] = ptr_pos[ d ];
        pd.add_dirac( p, weight );
    }

    size_t nb_diracs() const {
        return pd.nb_diracs();
    }

    PyDirac<dim,TF> dirac( size_t n ) const {
        return _dirac( n, N<dim>() );
    }

    PyDirac<dim,TF> _dirac( size_t n, N<2> ) const {
        return { pd.dirac( n ).pos.x, pd.dirac( n ).pos.y, 0, pd.dirac( n ).weight };
    }

    PyDirac<dim,TF> _dirac( size_t n, N<3> ) const {
        return { pd.dirac( n ).pos.x, pd.dirac( n ).pos.y, pd.dirac( n ).pos.z, pd.dirac( n ).weight };
    }

    void add_to_weights( py::array_t<TF> &values_to_add ) {
        auto buf_values = values_to_add.request();
        auto ptr_values = (const TF *)buf_values.ptr;
        for( size_t i = 0; i < pd.nb_diracs(); ++i )
            pd.dirac( i ).weight += ptr_values[ i ];
    }

    void add_to_positions( py::array_t<TF> &values_to_add ) {
        PO pt;
        auto buf_values = values_to_add.request();
        auto ptr_values = (const TF *)buf_values.ptr;
        for( size_t i = 0, o = 0; i < pd.nb_diracs(); ++i ) {
            for( std::size_t d = 0; d < dim; ++d )
                pt[ d ] = ptr_values[ o++ ];
            pd.dirac( i ).pos += pt;
        }
    }

    void set_positions( py::array_t<TF> &values ) {
        PO pt;
        auto buf_values = values.request();
        auto ptr_values = (const TF *)buf_values.ptr;
        for( size_t i = 0, o = 0; i < pd.nb_diracs(); ++i ) {
            for( std::size_t d = 0; d < dim; ++d )
                pt[ d ] = ptr_values[ o++ ];
            pd.dirac( i ).pos = pt;
        }
    }

    py::array_t<TF> centroids() {
        // compute
        std::vector<PO> g_res( pd.nb_diracs() );
        std::vector<TF> g_vol( pd.nb_diracs() );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i ) {
            for( std::size_t d = 0; d < dim; ++d )
                g_res[ i ][ d ] = 0;
            g_vol[ i ] = 0;
        }

        pd.get_centroid_contribs( g_res.data(), g_vol.data() );

        // save
        py::array_t<TF> res;
        res.resize( { pd.nb_diracs(), (unsigned long) dim } );
        auto buf_res = res.request();
        auto ptr_res = (TF *)buf_res.ptr;
        for( size_t i = 0, o = 0; i < pd.nb_diracs(); ++i ) {
            for( std::size_t d = 0; d < dim; ++d )
                ptr_res[ o++ ] = g_res[ i ][ d ] / g_vol[ i ];
        }

        return res;
    }

    py::array_t<TF> positions() {
        py::array_t<TF> res;
        res.resize( { (unsigned long)pd.nb_diracs(), (unsigned long)dim } );
        auto buf_res = res.request();
        auto ptr_res = (TF *)buf_res.ptr;
        for( size_t i = 0, o = 0; i < pd.nb_diracs(); ++i )
            for( std::size_t d = 0; d < dim; ++d )
                ptr_res[ o++ ] = pd.dirac( i ).pos[ d ];
        return res;
    }

    py::array_t<TF> weights() {
        py::array_t<TF> res;
        res.resize( { pd.nb_diracs() } );
        auto buf_res = res.request();
        auto ptr_res = (TF *)buf_res.ptr;
        for( size_t i = 0; i < pd.nb_diracs(); ++i )
            ptr_res[ i ] = pd.dirac( i ).weight;
        return res;
    }

    void add_convex_shape( py::array_t<TF> &pb, TF radius ) {
        auto buf_pb = pb.request();
        auto ptr_pb = (TF *)buf_pb.ptr;

        std::vector<typename PowerDiagram<Pc>::CutInfo> poly_bound( pb.shape()[ 0 ] / ( 2 * dim ) );
        for( size_t i = 0, o = 0; i < poly_bound.size(); ++i ) {
            PO cut_O;
            for( std::size_t d = 0; d < dim; ++d )
                cut_O[ d ] = ptr_pb[ o++ ];
            PO cut_N;
            for( std::size_t d = 0; d < dim; ++d )
                cut_N[ d ] = ptr_pb[ o++ ];
            poly_bound[ i ] = { cut_O, normalized( cut_N ) };
        }
        pd.add_convex_shape( poly_bound, radius );
    }

    void add_box_shape( py::array_t<TF> &pb ) {
        auto buf_pb = pb.request();
        auto ptr_pb = (TF *)buf_pb.ptr;

        PO p0;
        std::size_t o = 0;
        for( std::size_t d = 0; d < dim; ++d )
            p0[ d ] = ptr_pb[ o++ ];

        PO p1;
        for( std::size_t d = 0; d < dim; ++d )
            p1[ d ] = ptr_pb[ o++ ];
        pd.add_box_shape( p0, p1 );
    }

    void add_transformation( py::array_t<TF> &pb ) {
        //        auto buf_pb = pb.request();
        //        auto ptr_pb = (double *)buf_pb.ptr;
        //        if ( pb.shape()[ 0 ] != 3 || pb.shape()[ 1 ] != 3 )
        //            ERROR( "a tranformation should be a 3x3 matrix" );

        //        pd.add_transformation( ptr_pb[ 0 ], ptr_pb[ 1 ], ptr_pb[ 2 ], ptr_pb[ 3 ], ptr_pb[ 4 ]... );
        TODO;
    }

    size_t nb_polygonal_shapes() {
        TODO;
        return 0; // pd.nb_polygonal_shapes();
    }

    void vtk_output( const std::string &filename ) {
        VtkOutput<1,TF> vo( { "num" } );
        pd.display( vo );
        vo.save( filename );
    }

    void vtk_output_points( const std::string &filename ) {
        VtkOutput<1,TF> vo( { "num" } );
        for( std::size_t i = 0; i < pd.nb_diracs(); ++i )
            vo.add_point( pd.dirac( i ).pos, { double( i ) } );
        vo.save( filename );
    }

    void vtk_output_shapes( const std::string &filename ) {
        VtkOutput<0,TF> vo;
        pd.display_bounds( vo );
        vo.save( filename );
    }

    void set_ball_cut( bool ball_cut ) {
        pd.set_ball_cut( ball_cut );
    }

    double domain_measure() {
        return pd.domain_measure();
    }

    py::array_t<double> polygonal_shape( size_t num ) {
        //        const std::vector<PO> &pb = pd.polygonal_shape( num );

        py::array_t<double> res;
        //        res.resize( { pb.size(), 3ul } );
        //        auto buf_res = res.request();
        //        auto ptr_res = (double *)buf_res.ptr;
        //        for( size_t i = 0, o = 0; i < pb.size(); ++i ) {
        //            ptr_res[ o++ ] = pb[ i ].x;
        //            ptr_res[ o++ ] = pb[ i ].y;
        //        }
        return res;
    }

    py::array_t<double> measures() {
        py::array_t<double> res;
        res.resize( { nb_diracs() } );
        auto buf_res = res.request();
        pd.get_measures( (double *)buf_res.ptr );
        return res;
    }
  
  py::array_t<double> second_moments() {
        py::array_t<double> res;
        res.resize( { nb_diracs() } );
        auto buf_res = res.request();
        pd.get_measures( FunctionEnum::R2(), (double *)buf_res.ptr );
        return res;
    }

    PyMeasuresAndDerivatives<dim,TF> der_measures() {
        PyMeasuresAndDerivatives<dim,TF> res;

        // reservation of areas
        res.measures    .resize({ pd.nb_diracs()     }); auto buf_measures = res.measures    .request(); auto ptr_measures = (double *)buf_measures.ptr;

        // computation
        std::vector<std::vector<std::pair<size_t,TF>>> der_values( pd.nb_diracs() );
        pd.get_der_measures( ptr_measures, der_values.data() );

        // reservation of derivatives data
        size_t nz = 0;
        for( size_t i = 0; i < pd.nb_diracs(); ++i )
            nz += der_values[ i ].size();
        res.der_data   .resize({ nz                 }); auto buf_data    = res.der_data   .request(); auto ptr_data    = (double *)buf_data    .ptr;
        res.der_indptr .resize({ pd.nb_diracs() + 1 }); auto buf_indptr  = res.der_indptr .request(); auto ptr_indptr  = (size_t *)buf_indptr  .ptr;
        res.der_indices.resize({ nz                 }); auto buf_indices = res.der_indices.request(); auto ptr_indices = (size_t *)buf_indices.ptr;

        for( size_t i = 0, o = 0; ; ++i ) {
            ptr_indptr[ i ] = o;
            if ( i == pd.nb_diracs() )
                break;
            for( auto p : der_values[ i ] ) {
                ptr_indices[ o ] = p.first;
                ptr_data[ o ] = p.second;
                ++o;
            }
        }

        return res;
    }

    PowerDiagram<Pc> pd;
};

template<int dim,class M>
void reg_power_module_diagram( N<dim>, M &m ) {
    std::string pref = to_string( dim + 0 ) + "d";
    using TF = double;

    std::string p = "PowerDiagram" + pref;
    py::class_<PyPowerDiagram<dim,TF>>( m, p.c_str() )
        .def          ( py::init<>()                                                       , "" )
        .def          ( "__repr__"           , &PyPowerDiagram<dim,TF>::__repr__           , "" )
        .def          ( "add_dirac"          , &PyPowerDiagram<dim,TF>::add_dirac          , "add a dirac (a particle). Ex: p.add_dirac( [ x, y ], weight )  " )
        .def          ( "dirac"              , &PyPowerDiagram<dim,TF>::dirac              , "information on dirac by index. Ex: p.add_dirac( num_dirac )" )
        .def          ( "add_to_positions"   , &PyPowerDiagram<dim,TF>::add_to_positions   , "change the dirac positions. Ex for 2 dirac in 3D: p.add_to_positions( [ [ dx0, dy0, dz0 ], [ dx1, dy1, dz1 ] ] )  " )
        .def          ( "set_positions"      , &PyPowerDiagram<dim,TF>::set_positions      , "set the dirac positions. Ex for 2 dirac in 3D: p.set_positions( [ [ x0, y0, z0 ], [ x1, y1, z1 ] ] )" )
        .def          ( "add_to_weights"     , &PyPowerDiagram<dim,TF>::add_to_weights     , "change the dirac weight. Ex for 2 dirac in 3D: p.add_to_weights( [ dw0, dw1, ...] ] )" )
        .def          ( "nb_diracs"          , &PyPowerDiagram<dim,TF>::nb_diracs          , "number of diracs (return an integer)" )
        .def          ( "measures"           , &PyPowerDiagram<dim,TF>::measures           , "return a vector containing areas/volumes in 2D resp 3D for each power cell." )
        .def          ( "second_moments"     , &PyPowerDiagram<dim,TF>::second_moments     , "returns the second moment of the Laguerre cell wrt the site." )
        .def          ( "der_measures"       , &PyPowerDiagram<dim,TF>::der_measures       , "" )
        .def          ( "centroids"          , &PyPowerDiagram<dim,TF>::centroids          , "" )
        .def          ( "positions"          , &PyPowerDiagram<dim,TF>::positions          , "" )
        .def          ( "weights"            , &PyPowerDiagram<dim,TF>::weights            , "" )
        .def          ( "add_convex_shape"   , &PyPowerDiagram<dim,TF>::add_convex_shape   , "" )
        .def          ( "add_box_shape"      , &PyPowerDiagram<dim,TF>::add_box_shape      , "" )
        .def          ( "add_transformation" , &PyPowerDiagram<dim,TF>::add_transformation , "" )
        .def          ( "nb_polygonal_shapes", &PyPowerDiagram<dim,TF>::nb_polygonal_shapes, "" )
        .def          ( "polygonal_shape"    , &PyPowerDiagram<dim,TF>::polygonal_shape    , "" )
        .def          ( "vtk_output"         , &PyPowerDiagram<dim,TF>::vtk_output         , "" )
        .def          ( "vtk_output_points"  , &PyPowerDiagram<dim,TF>::vtk_output_points  , "" )
        .def          ( "vtk_output_shapes"  , &PyPowerDiagram<dim,TF>::vtk_output_shapes  , "" )
        .def          ( "set_ball_cut"       , &PyPowerDiagram<dim,TF>::set_ball_cut       , "p.set_ball_cut( True or False ) to cut each cell with a disc of radii kantorovitch_potential^0.5" )
        .def          ( "domain_measure"     , &PyPowerDiagram<dim,TF>::domain_measure     , "measure of the domain (defined using add_convex_shape, add_box_shape, ...)" )
    ;

    std::string v = "MeasureAndDerivatives" + pref;
    py::class_<PyMeasuresAndDerivatives<dim,TF>>( m, v.c_str() )
        .def_readwrite( "measures"           , &PyMeasuresAndDerivatives<dim,TF>::measures      )
        .def_readwrite( "der_data"           , &PyMeasuresAndDerivatives<dim,TF>::der_data      )
        .def_readwrite( "der_indptr"         , &PyMeasuresAndDerivatives<dim,TF>::der_indptr    )
        .def_readwrite( "der_indices"        , &PyMeasuresAndDerivatives<dim,TF>::der_indices   )
    ;

    std::string d = "Dirac" + pref;
    py::class_<PyDirac<dim,TF>>( m, d.c_str() )
        .def_readwrite( "x"                  , &PyDirac<dim,TF>::x                              )
        .def_readwrite( "y"                  , &PyDirac<dim,TF>::y                              )
        .def_readwrite( "z"                  , &PyDirac<dim,TF>::z                              )
        .def_readwrite( "r"                  , &PyDirac<dim,TF>::r                              )
    ;
}

PYBIND11_MODULE( power_diagram, m ) {
    m.doc() = "Power diagram tools";

    reg_power_module_diagram( N<2>(), m );
}

