#include "system/Assert.h"
#include "system/Stream.h"
#include "system/S.h"
#include "MpiInst.h"

#include <parmetis.h>
#include <metis.h>
#include <mpi.h>


namespace {

inline MPI_Datatype mpi_datatype( S<std::int8_t  > ) { return MPI_INT8_T  ; }
inline MPI_Datatype mpi_datatype( S<std::int16_t > ) { return MPI_INT16_T ; }
inline MPI_Datatype mpi_datatype( S<std::int32_t > ) { return MPI_INT32_T ; }
inline MPI_Datatype mpi_datatype( S<std::int64_t > ) { return MPI_INT64_T ; }

inline MPI_Datatype mpi_datatype( S<std::uint8_t > ) { return MPI_UINT8_T ; }
inline MPI_Datatype mpi_datatype( S<std::uint16_t> ) { return MPI_UINT16_T; }
inline MPI_Datatype mpi_datatype( S<std::uint32_t> ) { return MPI_UINT32_T; }
inline MPI_Datatype mpi_datatype( S<std::uint64_t> ) { return MPI_UINT64_T; }

}

MpiInst mpi_inst;

MpiInst::MpiInst() {
    _init = false;
}

MpiInst::~MpiInst() {
    if ( _init )
        MPI_Finalize();
}

void MpiInst::init( int argc, char **argv ) {
    if ( _init )
        return;
    _init = true;

    // init and dind out rank and size
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &_size );

    //
    mpi = this;
}

void MpiInst::send( const std::  int8_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_INT8_T  , destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const std:: uint8_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_UINT8_T , destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const std:: int32_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_INT32_T , destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const std::uint32_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_UINT32_T, destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const std:: int64_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_INT64_T , destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const std::uint64_t *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_UINT64_T, destination, tag, MPI_COMM_WORLD ); }
void MpiInst::send( const double        *data, std::size_t count, int destination, int tag ) { MPI_Send( data, count, MPI_DOUBLE  , destination, tag, MPI_COMM_WORLD ); }

void MpiInst::recv( std::  int8_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_INT8_T  , source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( std:: uint8_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_UINT8_T , source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( std:: int32_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_INT32_T , source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( std::uint32_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_UINT32_T, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( std:: int64_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_INT64_T , source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( std::uint64_t *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_UINT64_T, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }
void MpiInst::recv( double        *data, std::size_t count, int source, int tag ) { MPI_Recv( data, count, MPI_DOUBLE  , source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); }

void MpiInst::gather( std:: int32_t *dst, const std:: int32_t *src, std::size_t count, int root ) { MPI_Gather( src, count, MPI_INT32_T , dst, count, MPI_INT32_T , root, MPI_COMM_WORLD ); }
void MpiInst::gather( std:: int64_t *dst, const std:: int64_t *src, std::size_t count, int root ) { MPI_Gather( src, count, MPI_INT64_T , dst, count, MPI_INT64_T , root, MPI_COMM_WORLD ); }
void MpiInst::gather( std::uint32_t *dst, const std::uint32_t *src, std::size_t count, int root ) { MPI_Gather( src, count, MPI_UINT32_T, dst, count, MPI_UINT32_T, root, MPI_COMM_WORLD ); }
void MpiInst::gather( std::uint64_t *dst, const std::uint64_t *src, std::size_t count, int root ) { MPI_Gather( src, count, MPI_UINT64_T, dst, count, MPI_UINT64_T, root, MPI_COMM_WORLD ); }
void MpiInst::gather( double        *dst, const double        *src, std::size_t count, int root ) { MPI_Gather( src, count, MPI_DOUBLE  , dst, count, MPI_DOUBLE  , root, MPI_COMM_WORLD ); }

void MpiInst::bcast( std::  int8_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_INT8_T  , root, MPI_COMM_WORLD ); }
void MpiInst::bcast( std:: int32_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_INT32_T , root, MPI_COMM_WORLD ); }
void MpiInst::bcast( std:: int64_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_INT64_T , root, MPI_COMM_WORLD ); }
void MpiInst::bcast( std:: uint8_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_UINT8_T , root, MPI_COMM_WORLD ); }
void MpiInst::bcast( std::uint32_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_UINT32_T, root, MPI_COMM_WORLD ); }
void MpiInst::bcast( std::uint64_t *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_UINT64_T, root, MPI_COMM_WORLD ); }
void MpiInst::bcast( double        *vec, std::size_t count, int root ) { MPI_Bcast( vec, count, MPI_DOUBLE  , root, MPI_COMM_WORLD ); }

std:: int32_t MpiInst::reduction( std:: int32_t value, const std::function<std:: int32_t(std:: int32_t,std:: int32_t)> &f ) { std::vector<std:: int32_t> values( size() ); gather( values.data(), &value, 1, 0 ); bcast( values.data(), size(), 0 ); value = values[ 0 ]; for( int i = 1; i < size(); ++i ) value = f( value, values[ i ] ); return value; }
std:: int64_t MpiInst::reduction( std:: int64_t value, const std::function<std:: int64_t(std:: int64_t,std:: int64_t)> &f ) { std::vector<std:: int64_t> values( size() ); gather( values.data(), &value, 1, 0 ); bcast( values.data(), size(), 0 ); value = values[ 0 ]; for( int i = 1; i < size(); ++i ) value = f( value, values[ i ] ); return value; }
std::uint32_t MpiInst::reduction( std::uint32_t value, const std::function<std::uint32_t(std::uint32_t,std::uint32_t)> &f ) { std::vector<std::uint32_t> values( size() ); gather( values.data(), &value, 1, 0 ); bcast( values.data(), size(), 0 ); value = values[ 0 ]; for( int i = 1; i < size(); ++i ) value = f( value, values[ i ] ); return value; }
std::uint64_t MpiInst::reduction( std::uint64_t value, const std::function<std::uint64_t(std::uint64_t,std::uint64_t)> &f ) { std::vector<std::uint64_t> values( size() ); gather( values.data(), &value, 1, 0 ); bcast( values.data(), size(), 0 ); value = values[ 0 ]; for( int i = 1; i < size(); ++i ) value = f( value, values[ i ] ); return value; }
double        MpiInst::reduction( double        value, const std::function<double       (double       ,double       )> &f ) { std::vector<double       > values( size() ); gather( values.data(), &value, 1, 0 ); bcast( values.data(), size(), 0 ); value = values[ 0 ]; for( int i = 1; i < size(); ++i ) value = f( value, values[ i ] ); return value; }

std::size_t MpiInst::probe_size( int source, int tag ) {
    MPI_Status status;
    MPI_Probe( source, tag, MPI_COMM_WORLD, &status);

    int number_amount;
    MPI_Get_count( &status, MPI_CHAR, &number_amount );
    return number_amount;
}

void MpiInst::partition( std::vector<int> &partition, const std::vector<std::size_t> &node_off, const std::vector<std::size_t> &edge_indices, const std::vector<std::size_t> &edge_values, const std::vector<int> &edge_costs, const std::vector<double> &xyz, int dim, bool full_redistribution ) {
    //    if ( true ) {
    //        // std::vector<int> &partition,
    //        std::vector<std::size_t> new_node_off;
    //        std::vector<std::size_t> new_edge_indices;
    //        std::vector<std::size_t> new_edge_values ;



    //                // const std::vector<int> &edge_costs,
    //        // const std::vector<double> &xyz,
    //        // int dim, bool full_redistribution )

    //        //        std::vector<idx_t>  _vtxdist     ( node_off.begin(), node_off.end() );
    //        //        std::vector<idx_t>  _xadj        ( edge_indices.begin(), edge_indices.end() );
    //        //        std::vector<idx_t>  _adjncy      ( edge_values.begin(), edge_values.end() );
    //        //        std::vector<idx_t>  _adjwgt      ( edge_costs.begin(), edge_costs.end() );
    //        //        idx_t              *_vwgt        = NULL;
    //        //        //idx_t              *_adjwgt      = NULL;
    //        //        idx_t               _wgtflag     = 1;      // => weight on edges
    //        //        idx_t               _numflag     = 0;      // => C style numbering
    //        //        idx_t               _ndims       = dim;    //
    //        //        std::vector<real_t> _xyz         ( xyz.begin(), xyz.end() );
    //        //        idx_t               _ncon        = 1;      //
    //        //        idx_t               _nparts      = size();
    //        //        std::vector<real_t> _tpwgts      ( size(), 1.0 / size() );
    //        //        real_t              _ubvec[]     = { 1.05 };
    //        //        idx_t               _options[10] = { 0 };
    //        //        idx_t               _edgecut     = 0;
    //        //        std::vector<idx_t>  _part        ( edge_indices.size() - 1, 17 );
    //        //        MPI_Comm            _mpi_comm    = MPI_COMM_WORLD;
    //        //        real_t              _itr         = 1.0;

    //        //        if ( _adjncy.empty() ) _adjncy.push_back( 0 );
    //        //        if ( _part.empty() ) _part.push_back( 0 );

    //        //        METIS_PartGraphKway(
    //        //            &nvtxz, &ncon, _xadj.data(), _adjncy.data(), _vwgt, _adjwgt.data(),  &_wgtflag,
    //        //            &_numflag, &_ndims, _xyz.data(), &_ncon, &_nparts, _tpwgts.data(),
    //        //            _ubvec, _options, &_edgecut, _part.data(), &_mpi_comm
    //        //        );

    //    }


    // check if parmetis can be used
    bool nb_void_sst = 0;
    for( unsigned i = 0; i + 1 < node_off.size(); ++i  )
    if ( nb_void_sst ) {
        if ( nb_void_sst != node_off.size() - 2 ) {
            // -> centralize data
            TODO;
        }
        // -> call metis
        TODO;
        return;
    }

    std::vector<idx_t>  _vtxdist     ( node_off.begin(), node_off.end() );
    std::vector<idx_t>  _xadj        ( edge_indices.begin(), edge_indices.end() );
    std::vector<idx_t>  _adjncy      ( edge_values.begin(), edge_values.end() );
    std::vector<idx_t>  _adjwgt      ( edge_costs.begin(), edge_costs.end() );
    idx_t              *_vwgt        = NULL;
    //idx_t              *_adjwgt      = NULL;
    idx_t               _wgtflag     = 1;      // => weight on edges
    idx_t               _numflag     = 0;      // => C style numbering
    idx_t               _ndims       = dim;    //
    std::vector<real_t> _xyz         ( xyz.begin(), xyz.end() );
    idx_t               _ncon        = 1;      //
    idx_t               _nparts      = size();
    std::vector<real_t> _tpwgts      ( size(), 1.0 / size() );
    real_t              _ubvec[]     = { 1.05 };
    idx_t               _options[10] = { 0 };
    idx_t               _edgecut     = 0;
    std::vector<idx_t>  _part        ( edge_indices.size() - 1, 17 );
    MPI_Comm            _mpi_comm    = MPI_COMM_WORLD;
    real_t              _itr         = 1.0;

    if ( _adjwgt.empty() ) _adjwgt.push_back( 0 );
    if ( _adjncy.empty() ) _adjncy.push_back( 0 );
    if ( _part  .empty() ) _part  .push_back( 0 );
    if ( _xyz   .empty() ) _xyz   .push_back( 0 );

    if ( full_redistribution ) {
        //        ParMETIS_V3_PartKway(
        //            _vtxdist.data(), _xadj.data(), _adjncy.data(), _vwgt, _adjwgt.data(),  &_wgtflag,
        //            &_numflag, &_ncon, &_nparts, _tpwgts.data(),
        //            _ubvec, _options, &_edgecut, _part.data(), &_mpi_comm
        //        );
        ParMETIS_V3_PartGeomKway(
            _vtxdist.data(), _xadj.data(), _adjncy.data(), _vwgt, _adjwgt.data(),  &_wgtflag,
            &_numflag, &_ndims, _xyz.data(), &_ncon, &_nparts, _tpwgts.data(),
            _ubvec, _options, &_edgecut, _part.data(), &_mpi_comm
        );
    } else {
        ParMETIS_V3_AdaptiveRepart(
            _vtxdist.data(), _xadj.data(), _adjncy.data(), _vwgt, 0, _adjwgt.data(),
            &_wgtflag, &_numflag, &_ncon, &_nparts, _tpwgts.data(), _ubvec,
            &_itr, _options, &_edgecut, _part.data(), &_mpi_comm
        );
    }

    partition = { _part.begin(), _part.end() };
}

void MpiInst::partition( std::vector<int> &partition, const std::vector<std::size_t> &node_off, const std::vector<double> &xyz, int dim ) {
    std::vector<idx_t> _node_off( node_off.begin(), node_off.end() );
    std::vector<real_t> _xyz( xyz.begin(), xyz.end() );
    std::vector<idx_t> _partition( xyz.size() / dim, 17 );
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    idx_t _dim = dim;

    ParMETIS_V3_PartGeom( _node_off.data(), &_dim, _xyz.data(), _partition.data(), &mpi_comm );

    partition.resize( _partition.size() );
    for( std::size_t i = 0; i < partition.size(); ++i )
        partition[ i ] = _partition[ i ];
}

template<class T>
void MpiInst::_cross_sends( std::vector<std::vector<T> > &dst, const std::vector<std::vector<T> > &src ) {
    MPI_Datatype data_type = mpi_datatype( S<T>() );
    dst.resize( size() );

    // send asynchronously all the `src` content
    std::vector<MPI_Request> wr_requests( size() );
    for( int mpi_rank = 0; mpi_rank < size(); ++mpi_rank )
        if ( mpi_rank != rank() )
            MPI_Isend( src[ mpi_rank ].data(), src[ mpi_rank ].size(), data_type, mpi_rank, 0, MPI_COMM_WORLD, &wr_requests[ mpi_rank ] );

    // wait for the size, and read the content
    MPI_Status status;
    for( int mpi_rank = 0; mpi_rank < size(); ++mpi_rank ) {
        if ( mpi_rank != rank() ) {
            int count = probe_size( mpi_rank );
            dst[ mpi_rank ].resize( count / sizeof( T ) );
            MPI_Recv( dst[ mpi_rank ].data(), dst[ mpi_rank ].size(), data_type, mpi_rank, 0, MPI_COMM_WORLD, &status );
        }
    }

    //    std::vector<int> done( size(), false );
    //    for( int nb_remaining_reads = size() - 1; nb_remaining_reads; ) {
    //        for( int mpi_rank = 0; mpi_rank < size(); ++mpi_rank ) {
    //            if ( mpi_rank != rank() ) {
    //                if ( done[ mpi_rank ] )
    //                    continue;

    //                int flag;
    //                MPI_Iprobe( mpi_rank, 0, MPI_COMM_WORLD, &flag, &status );
    //                if ( flag ) {
    //                    int count;
    //                    --nb_remaining_reads;
    //                    done[ mpi_rank ] = true;
    //                    MPI_Get_count( &status, MPI_CHAR, &count );
    //                    dst[ mpi_rank ].resize( count / sizeof( std::size_t ) );
    //                    MPI_Recv( dst[ mpi_rank ].data(), dst[ mpi_rank ].size(), data_type, mpi_rank, 0, MPI_COMM_WORLD, &status );
    //                }
    //            }
    //        }
    //    }

    //  copy local data
    dst[ mpi->rank() ] = src[ mpi->rank() ];

    // wait for requests to be finished
    for( int mpi_rank = 0; mpi_rank < size(); ++mpi_rank )
        if ( mpi_rank != rank() )
            MPI_Wait( &wr_requests[ mpi_rank ], &status );
}

void MpiInst::cross_sends( std::vector<std::vector<std::uint8_t >> &dst, const std::vector<std::vector<std::uint8_t >> &src ) { _cross_sends( dst, src ); }
void MpiInst::cross_sends( std::vector<std::vector<std::uint32_t>> &dst, const std::vector<std::vector<std::uint32_t>> &src ) { _cross_sends( dst, src ); }
void MpiInst::cross_sends( std::vector<std::vector<std::uint64_t>> &dst, const std::vector<std::vector<std::uint64_t>> &src ) { _cross_sends( dst, src ); }

void MpiInst::barrier() {
    MPI_Barrier( MPI_COMM_WORLD );
}

