// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef LIBMESH_PARALLEL_IMPLEMENTATION_H
#define LIBMESH_PARALLEL_IMPLEMENTATION_H

// Local includes
#include "parallel.h"
#include "libmesh_logging.h"

// C++ includes
#include <iterator> // iterator_traits

namespace libMesh {
namespace Parallel {

// First declare StandardType specializations so we can use them in anonymous
// helper functions later

#ifdef LIBMESH_HAVE_MPI

#define LIBMESH_STANDARD_TYPE(cxxtype,mpitype)                          \
  template<>                                                            \
  class StandardType<cxxtype> : public DataType                         \
  {                                                                     \
  public:                                                               \
    explicit                                                            \
      StandardType(const cxxtype * = libmesh_nullptr) : DataType(mpitype) {} \
  }

#else

#define LIBMESH_STANDARD_TYPE(cxxtype,mpitype)                          \
  template<>                                                            \
  class StandardType<cxxtype> : public DataType                         \
  {                                                                     \
  public:                                                               \
    explicit                                                            \
      StandardType(const cxxtype * = libmesh_nullptr) : DataType() {}   \
  }

#endif

#define LIBMESH_INT_TYPE(cxxtype,mpitype)                               \
  LIBMESH_STANDARD_TYPE(cxxtype,mpitype);                               \
                                                                        \
  template<>                                                            \
  struct Attributes<cxxtype>                                            \
  {                                                                     \
    static const bool has_min_max = true;                               \
    static void set_lowest(cxxtype & x) { x = std::numeric_limits<cxxtype>::min(); } \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::max(); } \
  }

#define LIBMESH_FLOAT_TYPE(cxxtype,mpitype)                             \
  LIBMESH_STANDARD_TYPE(cxxtype,mpitype);                               \
                                                                        \
  template<>                                                            \
  struct Attributes<cxxtype>                                            \
  {                                                                     \
    static const bool has_min_max = true;                               \
    static void set_lowest(cxxtype & x) { x = -std::numeric_limits<cxxtype>::infinity(); } \
    static void set_highest(cxxtype & x) { x = std::numeric_limits<cxxtype>::infinity(); } \
  }

#define LIBMESH_CONTAINER_TYPE(cxxtype)                                 \
  template<typename T>                                                  \
  struct Attributes<cxxtype<T> >                                        \
  {                                                                     \
    static const bool has_min_max = Attributes<T>::has_min_max;         \
    static void set_lowest(cxxtype<T> & x) {                            \
      for (typename cxxtype<T>::iterator i = x.begin(); i != x.end(); ++i) \
        Attributes<T>::set_lowest(*i); }                                \
    static void set_highest(cxxtype<T> & x) {                           \
      for (typename cxxtype<T>::iterator i = x.begin(); i != x.end(); ++i) \
        Attributes<T>::set_highest(*i); }                               \
  }


LIBMESH_INT_TYPE(char,MPI_CHAR);
#if MPI_VERSION > 1
LIBMESH_INT_TYPE(signed char,MPI_SIGNED_CHAR);
#endif
LIBMESH_INT_TYPE(unsigned char,MPI_UNSIGNED_CHAR);
LIBMESH_INT_TYPE(short int,MPI_SHORT);
LIBMESH_INT_TYPE(unsigned short int,MPI_UNSIGNED_SHORT);
LIBMESH_INT_TYPE(int,MPI_INT);
LIBMESH_INT_TYPE(unsigned int,MPI_UNSIGNED);
LIBMESH_INT_TYPE(long,MPI_LONG);
LIBMESH_INT_TYPE(long long,MPI_LONG_LONG_INT);
LIBMESH_INT_TYPE(unsigned long,MPI_UNSIGNED_LONG);
#if MPI_VERSION > 1 || !defined(LIBMESH_HAVE_MPI)
LIBMESH_INT_TYPE(unsigned long long,MPI_UNSIGNED_LONG_LONG);
#else
// MPI 1.0 did not have an unsigned long long type, so we have to use
// MPI_UNSIGNED_LONG in this case.  If "unsigned long" and "unsigned
// long long" are different sizes on your system, we detect this and
// throw an error in dbg mode rather than communicating values
// incorrectly.
template<>
class StandardType<unsigned long long> : public DataType
{
public:
  explicit
  StandardType(const cxxtype * = libmesh_nullptr) :
    DataType(MPI_UNSIGNED_LONG)
  {
    libmesh_assert_equal_to(sizeof(unsigned long long),
                            sizeof(unsigned long));
  }
};

template<>
struct Attributes<unsigned long long>
{
  static const bool has_min_max = true;
  static void set_lowest(unsigned long long & x) { x = std::numeric_limits<unsigned long long>::min(); }
  static void set_highest(unsigned long long & x) { x = std::numeric_limits<unsigned long long>::max(); }
};
#endif

LIBMESH_FLOAT_TYPE(float,MPI_FLOAT);
LIBMESH_FLOAT_TYPE(double,MPI_DOUBLE);
LIBMESH_FLOAT_TYPE(long double,MPI_LONG_DOUBLE);
LIBMESH_CONTAINER_TYPE(std::set);
LIBMESH_CONTAINER_TYPE(std::vector);

template<typename T1, typename T2>
class StandardType<std::pair<T1, T2> > : public DataType
{
public:
  explicit
  StandardType(const std::pair<T1, T2> * example = libmesh_nullptr) {
    // We need an example for MPI_Address to use
    static const std::pair<T1, T2> p;
    if (!example)
      example = &p;

    // _static_type never gets freed, but it only gets committed once
    // per T, so it's not a *huge* memory leak...
    static data_type _static_type;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
#ifdef LIBMESH_HAVE_MPI

        // Get the sub-data-types, and make sure they live long enough
        // to construct the derived type
        StandardType<T1> d1(&example->first);
        StandardType<T2> d2(&example->second);

#if MPI_VERSION == 1

        // Use MPI_LB and MPI_UB here to workaround potential bugs from
        // nested MPI_LB and MPI_UB in the specifications of d1 and/or d2:
        // https://github.com/libMesh/libmesh/issues/631
        MPI_Datatype types[] = { MPI_LB, (data_type)d1, (data_type)d2, MPI_UB };
        int blocklengths[] = {1,1,1,1};
        MPI_Aint displs[4];

        libmesh_call_mpi
          (MPI_Address (const_cast<std::pair<T1,T2> *>(example),
                        &displs[0]));
        libmesh_call_mpi
          (MPI_Address (const_cast<T1*>(&example->first),
                        &displs[1]));
        libmesh_call_mpi
          (MPI_Address (const_cast<T2*>(&example->second),
                        &displs[2]));
        libmesh_call_mpi
          (MPI_Address (const_cast<std::pair<T1,T2> *>(example+1),
                        &displs[3]));

        displs[1] -= displs[0];
        displs[2] -= displs[0];
        displs[3] -= displs[0];
        displs[0] = 0;

        libmesh_call_mpi
          (MPI_Type_struct (4, blocklengths, displs, types,
                            &_static_type));
#else
        MPI_Datatype types[] = { (data_type)d1, (data_type)d2 };
        int blocklengths[] = {1,1};
        MPI_Aint displs[2], start;

        libmesh_call_mpi
          (MPI_Get_address (const_cast<std::pair<T1,T2> *>(example),
                            &start));
        libmesh_call_mpi
          (MPI_Get_address (const_cast<T1*>(&example->first),
                            &displs[0]));
        libmesh_call_mpi
          (MPI_Get_address (const_cast<T2*>(&example->second),
                            &displs[1]));
        displs[0] -= start;
        displs[1] -= start;

        // create a prototype structure
        MPI_Datatype tmptype;
        libmesh_call_mpi
          (MPI_Type_create_struct (2, blocklengths, displs, types,
                                   &tmptype));
        libmesh_call_mpi
          (MPI_Type_commit (&tmptype));

        // resize the structure type to account for padding, if any
        libmesh_call_mpi
          (MPI_Type_create_resized (tmptype, 0,
                                    sizeof(std::pair<T1,T2>),
                                    &_static_type));
#endif

        libmesh_call_mpi
          (MPI_Type_commit (&_static_type));
#endif // LIBMESH_HAVE_MPI

        _is_initialized = true;
      }

    _datatype = _static_type;
  }

  // Make sure not to free our singleton
  ~StandardType() {}
};

template<typename T>
class StandardType<std::complex<T> > : public DataType
{
public:
  explicit
  StandardType(const std::complex<T> * /*example*/ = libmesh_nullptr) :
    DataType(StandardType<T>(libmesh_nullptr), 2) {}

  ~StandardType() { this->free(); }
};

} // namespace Parallel

} // namespace libMesh


// Anonymous namespace for helper functions
namespace {

// Internal helper function to create vector<something_useable> from
// vector<bool> for compatibility with MPI bitwise operations
template <typename T>
inline void pack_vector_bool(const std::vector<bool> & vec_in,
                             std::vector<T> & vec_out)
{
  unsigned int data_bits = 8*sizeof(T);
  std::size_t in_size = vec_in.size();
  std::size_t out_size = in_size/data_bits + ((in_size%data_bits)?1:0);
  vec_out.clear();
  vec_out.resize(out_size);
  for (std::size_t i=0; i != in_size; ++i)
    {
      std::size_t index = i/data_bits;
      std::size_t offset = i%data_bits;
      vec_out[index] += (vec_in[i]?1:0) << offset;
    }
}

// Internal helper function to create vector<bool> from
// vector<something usable> for compatibility with MPI byte
// operations
template <typename T>
inline void unpack_vector_bool(const std::vector<T> & vec_in,
                               std::vector<bool> & vec_out)
{
  unsigned int data_bits = 8*sizeof(T);
  // We need the output vector to already be properly sized
  std::size_t out_size = vec_out.size();
  libmesh_assert_equal_to
    (out_size/data_bits + (out_size%data_bits?1:0), vec_in.size());

  for (std::size_t i=0; i != out_size; ++i)
    {
      std::size_t index = i/data_bits;
      std::size_t offset = i%data_bits;
      vec_out[i] = vec_in[index] << (data_bits-1-offset) >> (data_bits-1);
    }
}


#ifdef LIBMESH_HAVE_MPI
// We use a helper function here to avoid ambiguity when calling
// send_receive of (vector<vector<T>>,vector<vector<T>>)
template <typename T1, typename T2>
inline void send_receive_vec_of_vec(const unsigned int dest_processor_id,
                                    const std::vector<std::vector<T1> > & send,
                                    const unsigned int source_processor_id,
                                    std::vector<std::vector<T2> > & recv,
                                    const libMesh::Parallel::MessageTag & send_tag,
                                    const libMesh::Parallel::MessageTag & recv_tag,
                                    const libMesh::Parallel::Communicator & comm)
{
  LOG_SCOPE("send_receive()", "Parallel");

  if (dest_processor_id   == comm.rank() &&
      source_processor_id == comm.rank())
    {
      recv = send;
      return;
    }

  // temporary buffers - these will be sized in bytes
  // and manipulated with MPI_Pack and friends
  std::vector<char> sendbuf, recvbuf;

  // figure out how many bytes we need to pack all the data
  int packedsize=0, sendsize=0;

  // The outer buffer size
  libmesh_call_mpi
    (MPI_Pack_size (1,
                    libMesh::Parallel::StandardType<unsigned int>(),
                    comm.get(),
                    &packedsize));

  sendsize += packedsize;

  for (std::size_t i=0; i<send.size(); i++)
    {
      // The size of the ith inner buffer
      libmesh_call_mpi
        (MPI_Pack_size (1,
                        libMesh::Parallel::StandardType<unsigned int>(),
                        comm.get(),
                        &packedsize));

      sendsize += packedsize;

      // The data for each inner buffer
      libmesh_call_mpi
        (MPI_Pack_size (libMesh::cast_int<int>(send[i].size()),
                        libMesh::Parallel::StandardType<T1>
                        (send[i].empty() ? libmesh_nullptr : &send[i][0]),
                        comm.get(),
                        &packedsize));

      sendsize += packedsize;
    }

  libmesh_assert (sendsize /* should at least be 1! */);
  sendbuf.resize (sendsize);

  // Pack the send buffer
  int pos=0;

  // ... the size of the outer buffer
  sendsize = libMesh::cast_int<int>(send.size());

  libmesh_call_mpi
    (MPI_Pack (&sendsize, 1,
               libMesh::Parallel::StandardType<unsigned int>(),
               &sendbuf[0], libMesh::cast_int<int>(sendbuf.size()),
               &pos, comm.get()));

  for (std::size_t i=0; i<send.size(); i++)
    {
      // ... the size of the ith inner buffer
      sendsize = libMesh::cast_int<int>(send[i].size());

      libmesh_call_mpi
        (MPI_Pack (&sendsize, 1, libMesh::Parallel::StandardType<unsigned int>(),
                   &sendbuf[0], libMesh::cast_int<int>(sendbuf.size()), &pos,
                   comm.get()));

      // ... the contents of the ith inner buffer
      if (!send[i].empty())
        libmesh_call_mpi
          (MPI_Pack (const_cast<T1*>(&send[i][0]),
                     libMesh::cast_int<int>(send[i].size()),
                     libMesh::Parallel::StandardType<T1>(&send[i][0]),
                     &sendbuf[0],
                     libMesh::cast_int<int>(sendbuf.size()), &pos,
                     comm.get()));
    }

  libmesh_assert_equal_to (static_cast<unsigned int>(pos), sendbuf.size());

  libMesh::Parallel::Request request;

  comm.send (dest_processor_id, sendbuf, MPI_PACKED, request, send_tag);

  comm.receive (source_processor_id, recvbuf, MPI_PACKED, recv_tag);

  // Unpack the received buffer
  libmesh_assert (!recvbuf.empty());
  pos=0;
  libmesh_call_mpi
    (MPI_Unpack (&recvbuf[0], libMesh::cast_int<int>(recvbuf.size()), &pos,
                 &sendsize, 1, libMesh::Parallel::StandardType<unsigned int>(),
                 comm.get()));

  // ... size the outer buffer
  recv.resize (sendsize);

  for (std::size_t i=0; i<recv.size(); i++)
    {
      libmesh_call_mpi
        (MPI_Unpack (&recvbuf[0],
                     libMesh::cast_int<int>(recvbuf.size()), &pos,
                     &sendsize, 1,
                     libMesh::Parallel::StandardType<unsigned int>(),
                     comm.get()));

      // ... size the inner buffer
      recv[i].resize (sendsize);

      // ... unpack the inner buffer if it is not empty
      if (!recv[i].empty())
        libmesh_call_mpi
          (MPI_Unpack (&recvbuf[0],
                       libMesh::cast_int<int>(recvbuf.size()), &pos,
                       &recv[i][0],
                       libMesh::cast_int<int>(recv[i].size()),
                       libMesh::Parallel::StandardType<T2>(&recv[i][0]),
                       comm.get()));
    }

  request.wait();
}

#endif // LIBMESH_HAVE_MPI

} // Anonymous namespace



namespace libMesh
{

namespace Parallel
{

/*
 * A reference to the default libMesh communicator.  This is now
 * deprecated - instead of libMesh::Parallel::Communicator_World use
 * libMesh::CommWorld
 */
#ifdef LIBMESH_DISABLE_COMMWORLD
extern FakeCommunicator & Communicator_World;
#else
extern Communicator & Communicator_World;
#endif


/**
 * Helper function for range packing
 */
template <typename Context, typename Iter>
inline std::size_t packed_range_size (const Context * context,
                                      Iter range_begin,
                                      const Iter range_end)
{
  typedef typename std::iterator_traits<Iter>::value_type T;

  std::size_t buffer_size = 0;
  for (Iter range_count = range_begin;
       range_count != range_end;
       ++range_count)
    {
      buffer_size += Parallel::Packing<T>::packable_size(*range_count, context);
    }
  return buffer_size;
}


/**
 * Helper function for range packing
 */
template <typename Context, typename buffertype, typename Iter>
inline Iter pack_range (const Context * context,
                        Iter range_begin,
                        const Iter range_end,
                        std::vector<buffertype> & buffer,
                        // When we serialize into buffers, we need to use large buffers to optimize MPI
                        // bandwidth, but not so large as to risk allocation failures.  max_buffer_size
                        // is measured in number of buffer type entries; number of bytes may be 4 or 8
                        // times larger depending on configuration.
                        std::size_t approx_buffer_size)
{
  typedef typename std::iterator_traits<Iter>::value_type T;

  // Count the total size of and preallocate buffer for efficiency.
  // Prepare to stop early if the buffer would be too large.
  std::size_t buffer_size = 0;
  Iter range_stop = range_begin;
  for (; range_stop != range_end && buffer_size < approx_buffer_size;
       ++range_stop)
    {
      std::size_t next_buffer_size =
        Parallel::Packing<T>::packable_size(*range_stop, context);
      buffer_size += next_buffer_size;
    }
  buffer.reserve(buffer.size() + buffer_size);

  // Pack the objects into the buffer
  for (; range_begin != range_stop; ++range_begin)
    {
#ifndef NDEBUG
      std::size_t old_size = buffer.size();
#endif

      Parallel::Packing<T>::pack
        (*range_begin, back_inserter(buffer), context);

#ifndef NDEBUG
      unsigned int my_packable_size =
        Parallel::Packing<T>::packable_size(*range_begin, context);
      unsigned int my_packed_size =
        Parallel::Packing<T>::packed_size (buffer.begin() + old_size);
      libmesh_assert_equal_to (my_packable_size, my_packed_size);
      libmesh_assert_equal_to (buffer.size(), old_size + my_packable_size);
#endif
    }

  return range_stop;
}



/**
 * Helper function for range unpacking
 */
template <typename Context, typename buffertype,
          typename OutputIter, typename T>
inline void unpack_range (const std::vector<buffertype> & buffer,
                          Context * context,
                          OutputIter out_iter,
                          const T * /* output_type */)
{
  // Loop through the buffer and unpack each object, returning the
  // object pointer via the output iterator
  typename std::vector<buffertype>::const_iterator
    next_object_start = buffer.begin();

  while (next_object_start < buffer.end())
    {
      *out_iter++ = Parallel::Packing<T>::unpack(next_object_start, context);
      next_object_start +=
        Parallel::Packing<T>::packed_size(next_object_start);
    }

  // We should have used up the exact amount of data in the buffer
  libmesh_assert (next_object_start == buffer.end());
}


inline Communicator::Communicator () :
#ifdef LIBMESH_HAVE_MPI
  _communicator(MPI_COMM_SELF),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _I_duped_it(false) {}

inline Communicator::Communicator (const communicator & comm) :
#ifdef LIBMESH_HAVE_MPI
  _communicator(MPI_COMM_SELF),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _I_duped_it(false)
{
  this->assign(comm);
}

inline Communicator::~Communicator ()
{
  this->clear();
}

#ifdef LIBMESH_HAVE_MPI
inline void Communicator::split(int color, int key, Communicator & target) const
{
  target.clear();
  MPI_Comm newcomm;
  libmesh_call_mpi
    (MPI_Comm_split(this->get(), color, key, &newcomm));

  target.assign(newcomm);
  target._I_duped_it = true;
  target.send_mode(this->send_mode());
}
#else
inline void Communicator::split(int, int, Communicator & target) const
{
  target.assign(this->get());
}
#endif

inline void Communicator::duplicate(const Communicator & comm)
{
  this->duplicate(comm._communicator);
  this->send_mode(comm.send_mode());
}

#ifdef LIBMESH_HAVE_MPI
inline void Communicator::duplicate(const communicator & comm)
{
  if (_communicator != MPI_COMM_NULL)
    {
      libmesh_call_mpi
        (MPI_Comm_dup(comm, &_communicator));

      _I_duped_it = true;
    }
  this->assign(_communicator);
}
#else
inline void Communicator::duplicate(const communicator &) { }
#endif

inline void Communicator::clear() {
#ifdef LIBMESH_HAVE_MPI
  if (_I_duped_it)
    {
      libmesh_assert (_communicator != MPI_COMM_NULL);
      libmesh_call_mpi
        (MPI_Comm_free(&_communicator));

      _communicator = MPI_COMM_NULL;
    }
  _I_duped_it = false;
#endif
}

inline Communicator & Communicator::operator= (const communicator & comm)
{
  this->clear();
  this->assign(comm);
  return *this;
}

// Disallowed copy constructor
inline Communicator::Communicator (const Communicator &) :
#ifdef LIBMESH_HAVE_MPI
  _communicator(MPI_COMM_NULL),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _I_duped_it(false)
{
  libmesh_not_implemented();
}

inline void Communicator::assign(const communicator & comm)
{
  _communicator = comm;
#ifdef LIBMESH_HAVE_MPI
  if (_communicator != MPI_COMM_NULL)
    {
      int i;
      libmesh_call_mpi
        (MPI_Comm_size(_communicator, &i));

      libmesh_assert_greater_equal (i, 0);
      _size = static_cast<unsigned int>(i);

      libmesh_call_mpi
        (MPI_Comm_rank(_communicator, &i));

      libmesh_assert_greater_equal (i, 0);
      _rank = static_cast<unsigned int>(i);
    }
  else
    {
      _rank = 0;
      _size = 1;
    }
#endif
  _send_mode = DEFAULT;
}



inline Status::Status () :
  _status(),
  _datatype()
{}

inline Status::Status (const data_type & type) :
  _status(),
  _datatype(type)
{}

inline Status::Status (const status & stat) :
  _status(stat),
  _datatype()
{}

inline Status::Status (const status & stat,
                       const data_type & type) :
  _status(stat),
  _datatype(type)
{}

inline Status::Status (const Status & stat) :
  _status(stat._status),
  _datatype(stat._datatype)
{}

inline Status::Status (const Status    & stat,
                       const data_type & type) :
  _status(stat._status),
  _datatype(type)
{}

inline int Status::source () const
{
#ifdef LIBMESH_HAVE_MPI
  return _status.MPI_SOURCE;
#else
  return 0;
#endif
}

inline int Status::tag () const
{
#ifdef LIBMESH_HAVE_MPI
  return _status.MPI_TAG;
#else
  libmesh_not_implemented();
  return 0;
#endif
}

#ifdef LIBMESH_HAVE_MPI
inline unsigned int Status::size (const data_type & type) const
{
  int msg_size;
  libmesh_call_mpi
    (MPI_Get_count (const_cast<MPI_Status*>(&_status), type,
                    &msg_size));

  libmesh_assert_greater_equal (msg_size, 0);
  return msg_size;
}
#else
inline unsigned int Status::size (const data_type &) const
{
  libmesh_not_implemented();
  return 0;
}
#endif

inline unsigned int Status::size () const
{ return this->size (this->datatype()); }



inline Request::Request () :
#ifdef LIBMESH_HAVE_MPI
  _request(MPI_REQUEST_NULL),
#else
  _request(),
#endif
  post_wait_work(libmesh_nullptr)
{}

inline Request::Request (const request & r) :
  _request(r),
  post_wait_work(libmesh_nullptr)
{}

inline Request::Request (const Request & other) :
  _request(other._request),
  post_wait_work(other.post_wait_work)
{
  if (other._prior_request.get())
    _prior_request = UniquePtr<Request>
      (new Request(*other._prior_request.get()));

  // operator= should behave like a shared pointer
  if (post_wait_work)
    post_wait_work->second++;
}

inline void Request::cleanup()
{
  if (post_wait_work)
    {
      // Decrement the use count
      post_wait_work->second--;

      if (!post_wait_work->second)
        {
#ifdef DEBUG
          // If we're done using this request, then we'd better have
          // done the work we waited for
          for (std::vector<PostWaitWork *>::iterator i =
                 post_wait_work->first.begin();
               i != post_wait_work->first.end(); ++i)
            libmesh_assert(!(*i));
#endif
          delete post_wait_work;
          post_wait_work = libmesh_nullptr;
        }
    }
}

inline Request & Request::operator = (const Request & other)
{
  this->cleanup();
  _request = other._request;
  post_wait_work = other.post_wait_work;

  if (other._prior_request.get())
    _prior_request = UniquePtr<Request>
      (new Request(*other._prior_request.get()));

  // operator= should behave like a shared pointer
  if (post_wait_work)
    post_wait_work->second++;

  return *this;
}

inline Request & Request::operator = (const request & r)
{
  this->cleanup();
  _request = r;
  post_wait_work = libmesh_nullptr;
  return *this;
}

inline Request::~Request () {
  this->cleanup();
}

inline Status Request::wait ()
{
  LOG_SCOPE("wait()", "Parallel::Request");

  if (_prior_request.get())
    _prior_request->wait();

  Status stat;
#ifdef LIBMESH_HAVE_MPI
  libmesh_call_mpi
    (MPI_Wait (&_request, stat.get()));
#endif
  if (post_wait_work)
    for (std::vector<PostWaitWork *>::iterator i =
           post_wait_work->first.begin();
         i != post_wait_work->first.end(); ++i)
      {
        // The user should never try to give us NULL work or try
        // to wait() twice.
        libmesh_assert (*i);
        (*i)->run();
        delete (*i);
        *i = libmesh_nullptr;
      }

  return stat;
}

inline bool Request::test ()
{
#ifdef LIBMESH_HAVE_MPI
  int val=0;

  // MPI_STATUS_IGNORE is from MPI-2; using it with some versions of
  // MPICH may cause a crash:
  // https://bugzilla.mcs.anl.gov/globus/show_bug.cgi?id=1798
#if MPI_VERSION > 1
  libmesh_call_mpi
    (MPI_Test (&_request, &val, MPI_STATUS_IGNORE));
#else
  MPI_Status stat;
  libmesh_call_mpi
    (MPI_Test (&_request, &val, &stat));
#endif

  if (val)
    {
      libmesh_assert          (_request == MPI_REQUEST_NULL);
      libmesh_assert_equal_to (val, 1);
    }

  return val;
#else
  return true;
#endif
}

#ifdef LIBMESH_HAVE_MPI
inline bool Request::test (status & stat)
{
  int val=0;

  libmesh_call_mpi
    (MPI_Test (&_request, &val, &stat));

  return val;
}
#else
inline bool Request::test (status &)
{
  return true;
}
#endif

inline void Request::add_prior_request(const Request & req)
{
  // We're making a chain of prior requests, not a tree
  libmesh_assert(!req._prior_request.get());

  Request * new_prior_req = new Request(req);

  // new_prior_req takes ownership of our existing _prior_request
  new_prior_req->_prior_request.reset(this->_prior_request.release());

  // Our _prior_request now manages the new resource we just set up
  this->_prior_request.reset(new_prior_req);
}

inline void Request::add_post_wait_work(PostWaitWork * work)
{
  if (!post_wait_work)
    post_wait_work = new
      std::pair<std::vector <PostWaitWork * >, unsigned int>
      (std::vector <PostWaitWork * >(), 1);
  post_wait_work->first.push_back(work);
}



/**
 * Pause execution until all processors reach a certain point.
 */
#ifdef LIBMESH_HAVE_MPI
inline void Communicator::barrier () const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("barrier()", "Parallel");
      libmesh_call_mpi(MPI_Barrier (this->get()));
    }
}
#else
inline void Communicator::barrier () const {}
#endif


// legacy e.g. Paralell::send() methods, requires
// Communicator_World
#ifndef LIBMESH_DISABLE_COMMWORLD
inline void barrier (const Communicator & comm = Communicator_World)
{
  comm.barrier();
}

template <typename T>
inline bool verify(const T & r,
                   const Communicator & comm = Communicator_World)
{ return comm.verify(r); }

template <typename T>
inline void min(T & r,
                const Communicator & comm = Communicator_World)
{ comm.min(r); }

template <typename T, typename U>
inline void minloc(T & r,
                   U & min_id,
                   const Communicator & comm = Communicator_World)
{ comm.minloc(r, min_id); }

template <typename T>
inline void max(T & r,
                const Communicator & comm = Communicator_World)
{ comm.max(r); }

template <typename T, typename U>
inline void maxloc(T & r,
                   U & max_id,
                   const Communicator & comm = Communicator_World)
{ comm.maxloc(r, max_id); }

template <typename T>
inline void sum(T & r,
                const Communicator & comm = Communicator_World)
{ comm.sum(r); }

template <typename T>
inline void set_union(T & data, const unsigned int root_id,
                      const Communicator & comm = Communicator_World)
{ comm.set_union(data, root_id); }

template <typename T>
inline void set_union(T & data,
                      const Communicator & comm = Communicator_World)
{ comm.set_union(data); }

inline status probe (const unsigned int src_processor_id,
                     const MessageTag & tag=any_tag,
                     const Communicator & comm = Communicator_World)
{ return comm.probe(src_processor_id, tag); }

template <typename T>
inline void send (const unsigned int dest_processor_id,
                  const T & data,
                  const MessageTag & tag=no_tag,
                  const Communicator & comm = Communicator_World)
{ comm.send(dest_processor_id, data, tag); }

template <typename T>
inline void send (const unsigned int dest_processor_id,
                  const T & data,
                  Request & req,
                  const MessageTag & tag=no_tag,
                  const Communicator & comm = Communicator_World)
{ comm.send(dest_processor_id, data, req, tag); }

template <typename T>
inline void send (const unsigned int dest_processor_id,
                  const T & data,
                  const DataType & type,
                  const MessageTag & tag=no_tag,
                  const Communicator & comm = Communicator_World)
{ comm.send(dest_processor_id, data, type, tag); }

template <typename T>
inline void send (const unsigned int dest_processor_id,
                  const T & data,
                  const DataType & type,
                  Request & req,
                  const MessageTag & tag=no_tag,
                  const Communicator & comm = Communicator_World)
{ comm.send(dest_processor_id, data, type, req, tag); }


template <typename Context, typename Iter>
inline void send_packed_range (const unsigned int dest_processor_id,
                               const Context * context,
                               Iter range_begin,
                               const Iter range_end,
                               const MessageTag & tag=no_tag,
                               const Communicator & comm = Communicator_World)
{ comm.send_packed_range(dest_processor_id, context, range_begin, range_end, tag); }


template <typename Context, typename Iter>
inline void send_packed_range (const unsigned int dest_processor_id,
                               const Context * context,
                               Iter range_begin,
                               const Iter range_end,
                               Request & req,
                               const MessageTag & tag=no_tag,
                               const Communicator & comm = Communicator_World)
{ comm.send_packed_range(dest_processor_id, context, range_begin, range_end, req, tag); }


template <typename T>
inline void nonblocking_send (const unsigned int dest_processor_id,
                              T & buf,
                              const DataType & type,
                              Request & r,
                              const MessageTag & tag=no_tag,
                              const Communicator & comm = Communicator_World)
{ comm.send (dest_processor_id, buf, type, r, tag); }

template <typename T>
inline void nonblocking_send (const unsigned int dest_processor_id,
                              T & buf,
                              Request & r,
                              const MessageTag & tag=no_tag,
                              const Communicator & comm = Communicator_World)
{ comm.send (dest_processor_id, buf, r, tag); }

template <typename T>
inline Status receive (const unsigned int src_processor_id,
                       T & buf,
                       const MessageTag & tag=any_tag,
                       const Communicator & comm = Communicator_World)
{ return comm.receive (src_processor_id, buf, tag); }

template <typename T>
inline void receive (const unsigned int src_processor_id,
                     T & buf,
                     Request & req,
                     const MessageTag & tag=any_tag,
                     const Communicator & comm = Communicator_World)
{ comm.receive (src_processor_id, buf, req, tag); }

template <typename T>
inline Status receive (const unsigned int src_processor_id,
                       T & buf,
                       const DataType & type,
                       const MessageTag & tag=any_tag,
                       const Communicator & comm = Communicator_World)
{ return comm.receive (src_processor_id, buf, type, tag); }

template <typename T>
inline void receive (const unsigned int src_processor_id,
                     T & buf,
                     const DataType & type,
                     Request & req,
                     const MessageTag & tag=any_tag,
                     const Communicator & comm = Communicator_World)
{ comm.receive (src_processor_id, buf, type, req, tag); }

template <typename Context, typename OutputIter, typename T>
inline void receive_packed_range (const unsigned int src_processor_id,
                                  Context * context,
                                  OutputIter out_iter,
                                  const T * output_type,
                                  const MessageTag & tag=any_tag,
                                  const Communicator & comm = Communicator_World)
{
  comm.receive_packed_range (src_processor_id, context, out_iter,
                             output_type, tag);
}

// template <typename Context, typename OutputIter>
// inline void receive_packed_range (const unsigned int src_processor_id,
//                                   Context * context,
//                                   OutputIter out_iter,
//                                   Request & req,
//                                   const MessageTag & tag=any_tag,
//                                   const Communicator & comm = Communicator_World)
// { comm.receive_packed_range (src_processor_id, context, out_iter, req, tag); }

template <typename T>
inline void nonblocking_receive (const unsigned int src_processor_id,
                                 T & buf,
                                 const DataType & type,
                                 Request & r,
                                 const MessageTag & tag=any_tag,
                                 const Communicator & comm = Communicator_World)
{ comm.receive (src_processor_id, buf, type, r, tag); }

template <typename T>
inline void nonblocking_receive (const unsigned int src_processor_id,
                                 T & buf,
                                 Request & r,
                                 const MessageTag & tag=any_tag,
                                 const Communicator & comm = Communicator_World)
{ comm.receive (src_processor_id, buf, r, tag); }

template <typename T1, typename T2>
inline void send_receive(const unsigned int dest_processor_id,
                         T1 & send,
                         const unsigned int source_processor_id,
                         T2 & recv,
                         const MessageTag & send_tag = no_tag,
                         const MessageTag & recv_tag = any_tag,
                         const Communicator & comm = Communicator_World)
{ comm.send_receive(dest_processor_id, send, source_processor_id, recv,
                    send_tag, recv_tag); }

template <typename Context1, typename RangeIter,
          typename Context2, typename OutputIter, typename T>
inline void send_receive_packed_range(const unsigned int dest_processor_id,
                                      const Context1 * context1,
                                      RangeIter send_begin,
                                      const RangeIter send_end,
                                      const unsigned int source_processor_id,
                                      Context2 * context2,
                                      OutputIter out_iter,
                                      const T * output_type,
                                      const MessageTag & send_tag = no_tag,
                                      const MessageTag & recv_tag = any_tag,
                                      const Communicator & comm = Communicator_World)
{
  comm.send_receive_packed_range(dest_processor_id, context1,
                                 send_begin, send_end,
                                 source_processor_id, context2,
                                 out_iter, output_type,
                                 send_tag, recv_tag);
}

template <typename T1, typename T2>
inline void send_receive(const unsigned int dest_processor_id,
                         T1 & send,
                         const DataType & type1,
                         const unsigned int source_processor_id,
                         T2 & recv,
                         const DataType & type2,
                         const MessageTag & send_tag = no_tag,
                         const MessageTag & recv_tag = any_tag,
                         const Communicator & comm = Communicator_World)
{ comm.send_receive(dest_processor_id, send, type1, source_processor_id,
                    recv, type2, send_tag, recv_tag); }

template <typename T>
inline void gather(const unsigned int root_id,
                   T send,
                   std::vector<T> & recv,
                   const Communicator & comm = Communicator_World)
{ comm.gather(root_id, send, recv); }

template <typename T>
inline void gather(const unsigned int root_id,
                   std::vector<T> & r,
                   const Communicator & comm = Communicator_World)
{ comm.gather(root_id, r); }

template <typename T>
inline void allgather(T send,
                      std::vector<T> & recv,
                      const Communicator & comm = Communicator_World)
{ comm.allgather(send, recv); }

template <typename T>
inline void allgather(std::vector<T> & r,
                      const bool identical_buffer_sizes = false,
                      const Communicator & comm = Communicator_World)
{ comm.allgather(r, identical_buffer_sizes); }

template <typename Context, typename Iter, typename OutputIter>
inline void gather_packed_range (const unsigned int root_id,
                                 Context * context,
                                 Iter range_begin,
                                 const Iter range_end,
                                 OutputIter out_iter,
                                 const Communicator & comm = Communicator_World)
{ comm.gather_packed_range(root_id, context, range_begin, range_end, out_iter); }

template <typename Context, typename Iter, typename OutputIter>
inline void allgather_packed_range (Context * context,
                                    Iter range_begin,
                                    const Iter range_end,
                                    OutputIter out_iter,
                                    const Communicator & comm = Communicator_World)
{ comm.allgather_packed_range(context, range_begin, range_end, out_iter); }

template <typename T>
inline void alltoall(std::vector<T> & r,
                     const Communicator & comm = Communicator_World)
{ comm.alltoall(r); }

template <typename T>
inline void broadcast(T & data,
                      const unsigned int root_id=0,
                      const Communicator & comm = Communicator_World)
{ comm.broadcast(data, root_id); }

template <typename Context, typename OutputContext, typename Iter, typename OutputIter>
inline void broadcast_packed_range (const Context * context1,
                                    Iter range_begin,
                                    const Iter range_end,
                                    OutputContext * context2,
                                    OutputIter out_iter,
                                    const unsigned int root_id = 0,
                                    const Communicator & comm = Communicator_World)
{ comm.broadcast_packed_range(context1, range_begin, range_end, context2, out_iter, root_id); }

#endif // #ifndef LIBMESH_DISABLE_COMMWORLD

//-----------------------------------------------------------------------
// Parallel members

inline
MessageTag::~MessageTag()
{
  if (_comm)
    _comm->dereference_unique_tag(_tagvalue);
}


inline
MessageTag::MessageTag(const MessageTag & other)
  : _tagvalue(other._tagvalue), _comm(other._comm)
{
  if (_comm)
    _comm->reference_unique_tag(_tagvalue);
}


inline
MessageTag Communicator::get_unique_tag(int tagvalue) const
{
  if (used_tag_values.count(tagvalue))
    {
      // Get the largest value in the used values, and pick one
      // larger
      tagvalue = used_tag_values.rbegin()->first+1;
      libmesh_assert(!used_tag_values.count(tagvalue));
    }
  used_tag_values[tagvalue] = 1;

  // #ifndef NDEBUG
  //   // Make sure everyone called get_unique_tag and make sure
  //   // everyone got the same value
  //   int maxval = tagvalue;
  //   this->max(maxval);
  //   libmesh_assert_equal_to (tagvalue, maxval);
  // #endif

  return MessageTag(tagvalue, this);
}


inline
void Communicator::reference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  libmesh_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]++;
}


inline
void Communicator::dereference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  libmesh_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]--;
  // If we don't have any more outstanding references, we
  // don't even need to keep this tag in our "used" set.
  if (!used_tag_values[tagvalue])
    used_tag_values.erase(tagvalue);
}


#ifdef LIBMESH_HAVE_MPI
template<>
inline data_type dataplusint_type<short int>() { return MPI_SHORT_INT; }

template<>
inline data_type dataplusint_type<int>() { return MPI_2INT; }

template<>
inline data_type dataplusint_type<long>() { return MPI_LONG_INT; }

template<>
inline data_type dataplusint_type<float>() { return MPI_FLOAT_INT; }

template<>
inline data_type dataplusint_type<double>() { return MPI_DOUBLE_INT; }

template<>
inline data_type dataplusint_type<long double>() { return MPI_LONG_DOUBLE_INT; }

template <typename T>
inline bool Communicator::verify(const T & r) const
{
  if (this->size() > 1 && Attributes<T>::has_min_max == true)
    {
      T tempmin = r, tempmax = r;
      this->min(tempmin);
      this->max(tempmax);
      bool verified = (r == tempmin) &&
        (r == tempmax);
      this->min(verified);
      return verified;
    }

#ifdef LIBMESH_HAVE_CXX11
  static_assert(Attributes<T>::has_min_max,
                "Tried to verify an unverifiable type");
#endif

  return true;
}



template <>
inline bool Communicator::verify(const bool & r) const
{
  const unsigned char rnew = r;
  return this->verify(rnew);
}



template <typename T>
inline bool Communicator::semiverify(const T * r) const
{
  if (this->size() > 1 && Attributes<T>::has_min_max == true)
    {
      T tempmin, tempmax;
      if (r)
        tempmin = tempmax = *r;
      else
        {
          Attributes<T>::set_highest(tempmin);
          Attributes<T>::set_lowest(tempmax);
        }
      this->min(tempmin);
      this->max(tempmax);
      bool invalid = r && ((*r != tempmin) ||
                           (*r != tempmax));
      this->max(invalid);
      return !invalid;
    }

#ifdef LIBMESH_HAVE_CXX11
  static_assert(Attributes<T>::has_min_max,
                "Tried to semiverify an unverifiable type");
#endif

  return true;
}



template <>
inline bool Communicator::semiverify(const bool * r) const
{
  if (r)
    {
      const unsigned char rnew = *r;
      return this->semiverify(&rnew);
    }

  const unsigned char * rptr = libmesh_nullptr;
  return this->semiverify(rptr);
}



template <typename T>
inline bool Communicator::semiverify(const std::vector<T> * r) const
{
  if (this->size() > 1 && Attributes<T>::has_min_max == true)
    {
      std::size_t rsize = r ? r->size() : 0;
      std::size_t * psize = r ? &rsize : libmesh_nullptr;

      if (!this->semiverify(psize))
        return false;

      this->max(rsize);

      std::vector<T> tempmin, tempmax;
      if (r)
        {
          tempmin = tempmax = *r;
        }
      else
        {
          tempmin.resize(rsize);
          tempmax.resize(rsize);
          Attributes<std::vector<T> >::set_highest(tempmin);
          Attributes<std::vector<T> >::set_lowest(tempmax);
        }
      this->min(tempmin);
      this->max(tempmax);
      bool invalid = r && ((*r != tempmin) ||
                           (*r != tempmax));
      this->max(invalid);
      return !invalid;
    }

#ifdef LIBMESH_HAVE_CXX11
  static_assert(Attributes<T>::has_min_max,
                "Tried to semiverify a vector of an unverifiable type");
#endif

  return true;
}



inline bool Communicator::verify(const std::string & r) const
{
  if (this->size() > 1)
    {
      // Cannot use <char> since MPI_MIN is not
      // strictly defined for chars!
      std::vector<short int> temp; temp.reserve(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        temp.push_back(r[i]);
      return this->verify(temp);
    }
  return true;
}



inline bool Communicator::semiverify(const std::string * r) const
{
  if (this->size() > 1)
    {
      std::size_t rsize = r ? r->size() : 0;
      std::size_t * psize = r ? &rsize : libmesh_nullptr;

      if (!this->semiverify(psize))
        return false;

      this->max(rsize);

      // Cannot use <char> since MPI_MIN is not
      // strictly defined for chars!
      std::vector<short int> temp (rsize);
      if (r)
        {
          temp.reserve(rsize);
          for (std::size_t i=0; i != rsize; ++i)
            temp.push_back((*r)[i]);
        }

      std::vector<short int> * ptemp = r ? &temp: libmesh_nullptr;

      return this->semiverify(ptemp);
    }
  return true;
}



template <typename T>
inline void Communicator::min(T & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("min(scalar)", "Parallel");

      T temp = r;
      libmesh_call_mpi
        (MPI_Allreduce (&temp, &r, 1, StandardType<T>(&temp), MPI_MIN,
                        this->get()));
    }
}


inline void Communicator::min(bool & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("min(bool)", "Parallel");

      unsigned int tempsend = r;
      unsigned int temp;

      libmesh_call_mpi
        (MPI_Allreduce (&tempsend, &temp, 1,
                        StandardType<unsigned int>(), MPI_MIN,
                        this->get()));
      r = temp;
    }
}


template <typename T>
inline void Communicator::min(std::vector<T> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("min(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<T> temp(r);
      libmesh_call_mpi
        (MPI_Allreduce (&temp[0], &r[0], cast_int<int>(r.size()),
                        StandardType<T>(&temp[0]), MPI_MIN,
                        this->get()));
    }
}


inline void Communicator::min(std::vector<bool> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("min(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<unsigned int> ruint;
      pack_vector_bool(r, ruint);
      std::vector<unsigned int> temp(ruint.size());
      libmesh_call_mpi
        (MPI_Allreduce (&ruint[0], &temp[0],
                        cast_int<int>(ruint.size()),
                        StandardType<unsigned int>(), MPI_BAND,
                        this->get()));
      unpack_vector_bool(temp, r);
    }
}


template <typename T>
inline void Communicator::minloc(T & r,
                                 unsigned int & min_id) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("minloc(scalar)", "Parallel");

      DataPlusInt<T> data_in;
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<T> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1, dataplusint_type<T>(),
                        MPI_MINLOC, this->get()));
      r = data_out.val;
      min_id = data_out.rank;
    }
  else
    min_id = this->rank();
}


inline void Communicator::minloc(bool & r,
                                 unsigned int & min_id) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("minloc(bool)", "Parallel");

      DataPlusInt<int> data_in;
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type<int>(),
                        MPI_MINLOC, this->get()));
      r = data_out.val;
      min_id = data_out.rank;
    }
  else
    min_id = this->rank();
}


template <typename T>
inline void Communicator::minloc(std::vector<T> & r,
                                 std::vector<unsigned int> & min_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("minloc(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<T> > data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<T> > data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (&data_in[0], &data_out[0],
                        cast_int<int>(r.size()),
                        dataplusint_type<T>(), MPI_MINLOC,
                        this->get()));
      for (std::size_t i=0; i != r.size(); ++i)
        {
          r[i]      = data_out[i].val;
          min_id[i] = data_out[i].rank;
        }
    }
  else if (!r.empty())
    {
      for (std::size_t i=0; i != r.size(); ++i)
        min_id[i] = this->rank();
    }
}


inline void Communicator::minloc(std::vector<bool> & r,
                                 std::vector<unsigned int> & min_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("minloc(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<int> > data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<int> > data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (&data_in[0], &data_out[0],
                        cast_int<int>(r.size()),
                        StandardType<int>(), MPI_MINLOC,
                        this->get()));
      for (std::size_t i=0; i != r.size(); ++i)
        {
          r[i]      = data_out[i].val;
          min_id[i] = data_out[i].rank;
        }
    }
  else if (!r.empty())
    {
      for (std::size_t i=0; i != r.size(); ++i)
        min_id[i] = this->rank();
    }
}


template <typename T>
inline void Communicator::max(T & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("max(scalar)", "Parallel");

      T temp;
      libmesh_call_mpi
        (MPI_Allreduce (&r, &temp, 1, StandardType<T>(&r), MPI_MAX,
                        this->get()));
      r = temp;
    }
}


inline void Communicator::max(bool & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("max(bool)", "Parallel");

      unsigned int tempsend = r;
      unsigned int temp;
      libmesh_call_mpi
        (MPI_Allreduce (&tempsend, &temp, 1,
                        StandardType<unsigned int>(), MPI_MAX,
                        this->get()));
      r = temp;
    }
}


template <typename T>
inline void Communicator::max(std::vector<T> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("max(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<T> temp(r);
      libmesh_call_mpi
        (MPI_Allreduce (&temp[0], &r[0], cast_int<int>(r.size()),
                        StandardType<T>(&temp[0]), MPI_MAX,
                        this->get()));
    }
}


inline void Communicator::max(std::vector<bool> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("max(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<unsigned int> ruint;
      pack_vector_bool(r, ruint);
      std::vector<unsigned int> temp(ruint.size());
      libmesh_call_mpi
        (MPI_Allreduce (&ruint[0], &temp[0],
                        cast_int<int>(ruint.size()),
                        StandardType<unsigned int>(), MPI_BOR,
                        this->get()));
      unpack_vector_bool(temp, r);
    }
}


template <typename T>
inline void Communicator::maxloc(T & r,
                                 unsigned int & max_id) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("maxloc(scalar)", "Parallel");

      DataPlusInt<T> data_in;
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<T> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type<T>(),
                        MPI_MAXLOC, this->get()));
      r = data_out.val;
      max_id = data_out.rank;
    }
  else
    max_id = this->rank();
}


inline void Communicator::maxloc(bool & r,
                                 unsigned int & max_id) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("maxloc(bool)", "Parallel");

      DataPlusInt<int> data_in;
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type<int>(),
                        MPI_MAXLOC, this->get()));
      r = data_out.val;
      max_id = data_out.rank;
    }
  else
    max_id = this->rank();
}


template <typename T>
inline void Communicator::maxloc(std::vector<T> & r,
                                 std::vector<unsigned int> & max_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("maxloc(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<T> > data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<T> > data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (&data_in[0], &data_out[0],
                        cast_int<int>(r.size()),
                        dataplusint_type<T>(), MPI_MAXLOC,
                        this->get()));
      for (std::size_t i=0; i != r.size(); ++i)
        {
          r[i]      = data_out[i].val;
          max_id[i] = data_out[i].rank;
        }
    }
  else if (!r.empty())
    {
      for (std::size_t i=0; i != r.size(); ++i)
        max_id[i] = this->rank();
    }
}


inline void Communicator::maxloc(std::vector<bool> & r,
                                 std::vector<unsigned int> & max_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("maxloc(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<int> > data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<int> > data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (&data_in[0], &data_out[0],
                        cast_int<int>(r.size()),
                        StandardType<int>(), MPI_MAXLOC,
                        this->get()));
      for (std::size_t i=0; i != r.size(); ++i)
        {
          r[i]      = data_out[i].val;
          max_id[i] = data_out[i].rank;
        }
    }
  else if (!r.empty())
    {
      for (std::size_t i=0; i != r.size(); ++i)
        max_id[i] = this->rank();
    }
}


template <typename T>
inline void Communicator::sum(T & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("sum()", "Parallel");

      T temp = r;
      libmesh_call_mpi
        (MPI_Allreduce (&temp, &r, 1, StandardType<T>(&temp), MPI_SUM,
                        this->get()));
    }
}


template <typename T>
inline void Communicator::sum(std::vector<T> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<T> temp(r);
      libmesh_call_mpi
        (MPI_Allreduce (&temp[0], &r[0], cast_int<int>(r.size()),
                        StandardType<T>(&temp[0]), MPI_SUM,
                        this->get()));
    }
}


// We still do function overloading for complex sums - in a perfect
// world we'd have a StandardSumOp to go along with StandardType...
template <typename T>
inline void Communicator::sum(std::complex<T> & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("sum()", "Parallel");

      std::complex<T> temp(r);
      libmesh_call_mpi
        (MPI_Allreduce (&temp, &r, 2, StandardType<T>(), MPI_SUM,
                        this->get()));
    }
}


template <typename T>
inline void Communicator::sum(std::vector<std::complex<T> > & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<std::complex<T> > temp(r);
      libmesh_call_mpi
        (MPI_Allreduce (&temp[0], &r[0], cast_int<int>(r.size() * 2),
                        StandardType<T>(libmesh_nullptr), MPI_SUM, this->get()));
    }
}


template <typename T>
inline void Communicator::set_union(std::set<T> & data,
                                    const unsigned int root_id) const
{
  std::vector<T> vecdata(data.begin(), data.end());
  this->gather(root_id, vecdata);
  if (this->rank() == root_id)
    data.insert(vecdata.begin(), vecdata.end());
}



template <typename T>
inline void Communicator::set_union(std::set<T> & data) const
{
  std::vector<T> vecdata(data.begin(), data.end());
  this->allgather(vecdata, false);
  data.insert(vecdata.begin(), vecdata.end());
}



template <typename T1, typename T2>
inline void Communicator::set_union(std::map<T1,T2> & data,
                                    const unsigned int root_id) const
{
  std::vector<std::pair<T1,T2> > vecdata(data.begin(), data.end());
  this->gather(root_id, vecdata);
  if (this->rank() == root_id)
    data.insert(vecdata.begin(), vecdata.end());
}



template <typename T1, typename T2>
inline void Communicator::set_union(std::map<T1,T2> & data) const
{
  std::vector<std::pair<T1,T2> > vecdata(data.begin(), data.end());
  this->allgather(vecdata, false);
  data.insert(vecdata.begin(), vecdata.end());
}



inline status Communicator::probe (const unsigned int src_processor_id,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("probe()", "Parallel");

  status stat;

  libmesh_call_mpi
    (MPI_Probe (src_processor_id, tag.value(), this->get(), &stat));

  return stat;
}

template<typename T>
inline Status Communicator::packed_range_probe (const unsigned int src_processor_id,
                                                const MessageTag & tag,
                                                bool & flag) const
{
  START_LOG("packed_range_probe()", "Parallel");

  libmesh_experimental();

  Status stat((StandardType<typename Packing<T>::buffer_type>()));

  int int_flag;

  libmesh_call_mpi(MPI_Iprobe(src_processor_id,
                              tag.value(),
                              this->get(),
                              &int_flag,
                              stat.get()));

  flag = int_flag;

  STOP_LOG("packed_range_probe()", "Parallel");

  return stat;
}


template<typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::basic_string<T> & buf,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = buf.empty() ? libmesh_nullptr : const_cast<T *>(buf.data());

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Ssend : MPI_Send) (dataptr,
                             cast_int<int>(buf.size()),
                             StandardType<T>(dataptr),
                             dest_processor_id,
                             tag.value(),
                             this->get()));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::basic_string<T> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = buf.empty() ? libmesh_nullptr : const_cast<T *>(buf.data());

  std::cerr<<"Sending: "<<buf.size()<<std::endl;

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (dataptr,
                               cast_int<int>(buf.size()),
                               StandardType<T>(dataptr),
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const T & buf,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = const_cast<T*> (&buf);

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Ssend : MPI_Send) (dataptr,
                             1,
                             StandardType<T>(dataptr),
                             dest_processor_id,
                             tag.value(),
                             this->get()));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const T & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = const_cast<T*>(&buf);

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (dataptr,
                               1,
                               StandardType<T>(dataptr),
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T> & buf,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), req, tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T> & buf,
                                const DataType & type,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  std::vector<T> vecbuf(buf.begin(), buf.end());
  this->send(dest_processor_id, vecbuf, type, tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T> & buf,
                                const DataType & type,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  // Allocate temporary buffer on the heap so it lives until after
  // the non-blocking send completes
  std::vector<T> * vecbuf =
    new std::vector<T>(buf.begin(), buf.end());

  // Make the Request::wait() handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T> >(vecbuf));

  this->send(dest_processor_id, *vecbuf, type, req, tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T> & buf,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? libmesh_nullptr : &buf.front()), tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? libmesh_nullptr : &buf.front()), req, tag);
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T> & buf,
                                const DataType & type,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Ssend : MPI_Send) (buf.empty() ? libmesh_nullptr : const_cast<T*>(&buf[0]),
                             cast_int<int>(buf.size()),
                             type,
                             dest_processor_id,
                             tag.value(),
                             this->get()));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T> & buf,
                                const DataType & type,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (buf.empty() ? libmesh_nullptr : const_cast<T*>(&buf[0]),
                               cast_int<int>(buf.size()),
                               type,
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));
}


template <typename Context, typename Iter>
inline void Communicator::send_packed_range (const unsigned int dest_processor_id,
                                             const Context * context,
                                             Iter range_begin,
                                             const Iter range_end,
                                             const MessageTag & tag) const
{
  // We will serialize variable size objects from *range_begin to
  // *range_end as a sequence of plain data (e.g. ints) in this buffer
  typedef typename std::iterator_traits<Iter>::value_type T;

  std::size_t total_buffer_size =
    Parallel::packed_range_size (context, range_begin, range_end);

  this->send(dest_processor_id, total_buffer_size, tag);

#ifdef DEBUG
  std::size_t used_buffer_size = 0;
#endif

  while (range_begin != range_end)
    {
      libmesh_assert_greater (std::distance(range_begin, range_end), 0);

      std::vector<typename Parallel::Packing<T>::buffer_type> buffer;

      const Iter next_range_begin = Parallel::pack_range
        (context, range_begin, range_end, buffer);

      libmesh_assert_greater (std::distance(range_begin, next_range_begin), 0);

      range_begin = next_range_begin;

#ifdef DEBUG
      used_buffer_size += buffer.size();
#endif

      // Blocking send of the buffer
      this->send(dest_processor_id, buffer, tag);
    }

#ifdef DEBUG
  libmesh_assert_equal_to(used_buffer_size, total_buffer_size);
#endif
}


template <typename Context, typename Iter>
inline void Communicator::send_packed_range (const unsigned int dest_processor_id,
                                             const Context * context,
                                             Iter range_begin,
                                             const Iter range_end,
                                             Request & req,
                                             const MessageTag & tag) const
{
  // Allocate a buffer on the heap so we don't have to free it until
  // after the Request::wait()
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  std::size_t total_buffer_size =
    Parallel::packed_range_size (context, range_begin, range_end);

  // That local variable will be gone soon; we need a send buffer that
  // will stick around.  I heard you like buffering so I put a buffer
  // for your buffer size so you can buffer the size of your buffer.
  std::size_t * total_buffer_size_buffer = new std::size_t;
  *total_buffer_size_buffer = total_buffer_size;

  // Delete the buffer size's buffer when we're done
  Request intermediate_req = request();
  intermediate_req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::size_t>(total_buffer_size_buffer));
  this->send(dest_processor_id, *total_buffer_size_buffer, intermediate_req, tag);

  // And don't finish up the full request until we're done with its
  // dependencies
  req.add_prior_request(intermediate_req);

#ifdef DEBUG
  std::size_t used_buffer_size = 0;
#endif

  while (range_begin != range_end)
    {
      libmesh_assert_greater (std::distance(range_begin, range_end), 0);

      std::vector<buffer_t> * buffer = new std::vector<buffer_t>();

      const Iter next_range_begin =
        Parallel::pack_range(context, range_begin, range_end,
                             *buffer);

      libmesh_assert_greater (std::distance(range_begin, next_range_begin), 0);

      range_begin = next_range_begin;

#ifdef DEBUG
      used_buffer_size += buffer->size();
#endif

      Request next_intermediate_req;

      Request * my_req = (range_begin == range_end) ? &req : &next_intermediate_req;

      // Make the Request::wait() handle deleting the buffer
      my_req->add_post_wait_work
        (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t> >
         (buffer));

      // Non-blocking send of the buffer
      this->send(dest_processor_id, *buffer, *my_req, tag);

      if (range_begin != range_end)
        req.add_prior_request(*my_req);
    }
}







template <typename Context, typename Iter>
inline void Communicator::nonblocking_send_packed_range (const unsigned int dest_processor_id,
                                                         const Context * context,
                                                         Iter range_begin,
                                                         const Iter range_end,
                                                         Request & req,
                                                         const MessageTag & tag) const
{
  libmesh_experimental();

  // Allocate a buffer on the heap so we don't have to free it until
  // after the Request::wait()
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  if (range_begin != range_end)
    {
      std::vector<buffer_t> * buffer = new std::vector<buffer_t>();

      range_begin =
        Parallel::pack_range(context,
                             range_begin,
                             range_end,
                             *buffer,
                             // MPI-2 can only use integers for size
                             std::numeric_limits<int>::max());

      if (range_begin != range_end)
        libmesh_error_msg("Non-blocking packed range sends cannot exceed " << std::numeric_limits<int>::max() << "in size");

      // Make the Request::wait() handle deleting the buffer
      req.add_post_wait_work
        (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t> >
         (buffer));

      // Non-blocking send of the buffer
      this->send(dest_processor_id, *buffer, req, tag);
    }
}


template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::basic_string<T> & buf,
                                     const MessageTag & tag) const
{
  std::vector<T> tempbuf;  // Officially C++ won't let us get a
                           // modifiable array from a string

  Status stat = this->receive(src_processor_id, tempbuf, tag);
  buf.assign(tempbuf.begin(), tempbuf.end());
  return stat;
}



template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::basic_string<T> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  // Officially C++ won't let us get a modifiable array from a
  // string, and we can't even put one on the stack for the
  // non-blocking case.
  std::vector<T> * tempbuf = new std::vector<T>();

  // We can clear the string, but the Request::wait() will need to
  // handle copying our temporary buffer to it
  buf.clear();

  req.add_post_wait_work
    (new Parallel::PostWaitCopyBuffer<std::vector<T>,
     std::back_insert_iterator<std::basic_string<T> > >
     (tempbuf, std::back_inserter(buf)));

  // Make the Request::wait() then handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T> >(tempbuf));

  this->receive(src_processor_id, tempbuf, req, tag);
}



template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     T & buf,
                                     const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  // Get the status of the message, explicitly provide the
  // datatype so we can later query the size
  Status stat(this->probe(src_processor_id, tag), StandardType<T>(&buf));

  libmesh_call_mpi
    (MPI_Recv (&buf, 1, StandardType<T>(&buf), src_processor_id,
               tag.value(), this->get(), stat.get()));

  return stat;
}



template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   T & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  libmesh_call_mpi
    (MPI_Irecv (&buf, 1, StandardType<T>(&buf), src_processor_id,
                tag.value(), this->get(), req.get()));
}



template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::set<T> & buf,
                                     const MessageTag & tag) const
{
  return this->receive
    (src_processor_id, buf,
     StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), tag);
}



/*
 * No non-blocking receives of std::set until we figure out how to
 * resize the temporary buffer
 */
#if 0
template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::set<T> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  this->receive (src_processor_id, buf,
                 StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), req, tag);
}
#endif // 0



template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::set<T> & buf,
                                     const DataType & type,
                                     const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  std::vector<T> vecbuf;
  Status stat = this->receive(src_processor_id, vecbuf, type, tag);
  buf.clear();
  buf.insert(vecbuf.begin(), vecbuf.end());

  return stat;
}



/*
 * No non-blocking receives of std::set until we figure out how to
 * resize the temporary buffer
 */
#if 0
template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::set<T> & buf,
                                   const DataType & type,
                                   Request & req,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  // Allocate temporary buffer on the heap so it lives until after
  // the non-blocking send completes
  std::vector<T> * vecbuf = new std::vector<T>();

  // We can clear the set, but the Request::wait() will need to
  // handle copying our temporary buffer to it
  buf.clear();

  req.add_post_wait_work
    (new Parallel::PostWaitCopyBuffer<std::vector<T>,
     std::insert_iterator<std::set<T> > >
     (*vecbuf, std::inserter(buf,buf.end())));

  // Make the Request::wait() then handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T> >(vecbuf));

  this->receive(src_processor_id, *vecbuf, type, req, tag);
}
#endif // 0



template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<T> & buf,
                                     const MessageTag & tag) const
{
  return this->receive
    (src_processor_id, buf,
     StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), tag);
}



template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<T> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  this->receive (src_processor_id, buf,
                 StandardType<T>(buf.empty() ? libmesh_nullptr : &(*buf.begin())), req, tag);
}



template <typename T>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<T> & buf,
                                     const DataType & type,
                                     const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  // Get the status of the message, explicitly provide the
  // datatype so we can later query the size
  Status stat(this->probe(src_processor_id, tag), type);

  buf.resize(stat.size());

  // Use stat.source() and stat.tag() in the receive - if
  // src_processor_id is or tag is "any" then we want to be sure we
  // try to receive the same message we just probed.
  libmesh_call_mpi
    (MPI_Recv (buf.empty() ? libmesh_nullptr : &buf[0],
               cast_int<int>(buf.size()), type, stat.source(),
               stat.tag(), this->get(), stat.get()));

  libmesh_assert_equal_to (stat.size(), buf.size());

  return stat;
}



template <typename T>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<T> & buf,
                                   const DataType & type,
                                   Request & req,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  libmesh_call_mpi
    (MPI_Irecv (buf.empty() ? libmesh_nullptr : &buf[0],
                cast_int<int>(buf.size()), type, src_processor_id,
                tag.value(), this->get(), req.get()));
}


template <typename Context, typename OutputIter, typename T>
inline void Communicator::receive_packed_range (const unsigned int src_processor_id,
                                                Context * context,
                                                OutputIter out_iter,
                                                const T * output_type,
                                                const MessageTag & tag) const
{
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  // Receive serialized variable size objects as sequences of buffer_t
  std::size_t total_buffer_size = 0;
  Status stat = this->receive(src_processor_id, total_buffer_size, tag);

  // Use stat.source() and stat.tag() in subsequent receives - if
  // src_processor_id is or tag is "any" then we want to be sure we
  // try to receive messages all corresponding to the same send.

  std::size_t received_buffer_size = 0;
  while (received_buffer_size < total_buffer_size)
    {
      std::vector<buffer_t> buffer;
      this->receive(stat.source(), buffer, MessageTag(stat.tag()));
      received_buffer_size += buffer.size();
      Parallel::unpack_range
        (buffer, context, out_iter, output_type);
    }
}



// template <typename Context, typename OutputIter>
// inline void Communicator::receive_packed_range (const unsigned int src_processor_id,
//                                                 Context * context,
//                                                 OutputIter out_iter,
//                                                 Request & req,
//                                                 const MessageTag & tag) const
// {
//   typedef typename std::iterator_traits<OutputIter>::value_type T;
//   typedef typename Parallel::Packing<T>::buffer_type buffer_t;
//
//   // Receive serialized variable size objects as a sequence of
//   // buffer_t.
//   // Allocate a buffer on the heap so we don't have to free it until
//   // after the Request::wait()
//   std::vector<buffer_t> * buffer = new std::vector<buffer_t>();
//   this->receive(src_processor_id, *buffer, req, tag);
//
//   // Make the Request::wait() handle unpacking the buffer
//   req.add_post_wait_work
//     (new Parallel::PostWaitUnpackBuffer<std::vector<buffer_t>, Context, OutputIter>
//      (buffer, context, out_iter));
//
//   // Make the Request::wait() then handle deleting the buffer
//   req.add_post_wait_work
//     (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t> >(buffer));
// }

template <typename Context, typename OutputIter, typename T>
inline void Communicator::nonblocking_receive_packed_range (const unsigned int src_processor_id,
                                                            Context * context,
                                                            OutputIter out,
                                                            const T * /* output_type */,
                                                            Request & req,
                                                            Status & stat,
                                                            const MessageTag & tag) const
{
  libmesh_experimental();

  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  // Receive serialized variable size objects as a sequence of
  // buffer_t.
  // Allocate a buffer on the heap so we don't have to free it until
  // after the Request::wait()
  std::vector<buffer_t> * buffer = new std::vector<buffer_t>(stat.size());
  this->receive(src_processor_id, *buffer, req, tag);

  // Make the Request::wait() handle unpacking the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitUnpackBuffer<std::vector<buffer_t>, Context, OutputIter, T>(*buffer, context, out));

  // Make the Request::wait() then handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t> >(buffer));
}



template <typename T1, typename T2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T1> & sendvec,
                                       const DataType & type1,
                                       const unsigned int source_processor_id,
                                       std::vector<T2> & recv,
                                       const DataType & type2,
                                       const MessageTag & send_tag,
                                       const MessageTag & recv_tag) const
{
  LOG_SCOPE("send_receive()", "Parallel");

  if (dest_processor_id   == this->rank() &&
      source_processor_id == this->rank())
    {
      recv = sendvec;
      return;
    }

  Parallel::Request req;

  this->send (dest_processor_id, sendvec, type1, req, send_tag);

  this->receive (source_processor_id, recv, type2, recv_tag);

  req.wait();
}



template <typename T1, typename T2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const T1 & sendvec,
                                       const unsigned int source_processor_id,
                                       T2 & recv,
                                       const MessageTag & send_tag,
                                       const MessageTag & recv_tag) const
{
  LOG_SCOPE("send_receive()", "Parallel");

  if (dest_processor_id   == this->rank() &&
      source_processor_id == this->rank())
    {
      recv = sendvec;
      return;
    }

  // MPI_STATUS_IGNORE is from MPI-2; using it with some versions of
  // MPICH may cause a crash:
  // https://bugzilla.mcs.anl.gov/globus/show_bug.cgi?id=1798
#if MPI_VERSION > 1
  libmesh_call_mpi
    (MPI_Sendrecv(const_cast<T1*>(&sendvec), 1, StandardType<T1>(&sendvec),
                  dest_processor_id, send_tag.value(), &recv, 1,
                  StandardType<T2>(&recv), source_processor_id,
                  recv_tag.value(), this->get(), MPI_STATUS_IGNORE));
#else
  MPI_Status stat;
  libmesh_call_mpi
    (MPI_Sendrecv(const_cast<T1*>(&sendvec), 1, StandardType<T1>(&sendvec),
                  dest_processor_id, send_tag.value(), &recv, 1,
                  StandardType<T2>(&recv), source_processor_id,
                  recv_tag.value(), this->get(), &stat));
#endif
}



// This is both a declaration and definition for a new overloaded
// function template, so we have to re-specify the default
// arguments.
//
// We specialize on the T1==T2 case so that we can handle
// send_receive-to-self with a plain copy rather than going through
// MPI.
template <typename T>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<T> & recv,
                                       const MessageTag & send_tag,
                                       const MessageTag & recv_tag) const
{
  if (dest_processor_id   == this->rank() &&
      source_processor_id == this->rank())
    {
      LOG_SCOPE("send_receive()", "Parallel");
      recv = sendvec;
      return;
    }

  const T* example = sendvec.empty() ?
    (recv.empty() ? libmesh_nullptr : &recv[0]) : &sendvec[0];

  // Call the user-defined type version with automatic
  // type conversion based on template argument:
  this->send_receive (dest_processor_id, sendvec,
                      StandardType<T>(example),
                      source_processor_id, recv,
                      StandardType<T>(example),
                      send_tag, recv_tag);
}


// This is both a declaration and definition for a new overloaded
// function template, so we have to re-specify the default arguments
template <typename T1, typename T2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T1> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<T2> & recv,
                                       const MessageTag & send_tag,
                                       const MessageTag & recv_tag) const
{
  // Call the user-defined type version with automatic
  // type conversion based on template argument:
  this->send_receive (dest_processor_id, sendvec,
                      StandardType<T1>(sendvec.empty() ? libmesh_nullptr : &sendvec[0]),
                      source_processor_id, recv,
                      StandardType<T2>(recv.empty() ? libmesh_nullptr : &recv[0]),
                      send_tag, recv_tag);
}




template <typename T1, typename T2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<std::vector<T1> > & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<std::vector<T2> > & recv,
                                       const MessageTag & /* send_tag */,
                                       const MessageTag & /* recv_tag */) const
{
  // FIXME - why aren't we honoring send_tag and recv_tag here?
  send_receive_vec_of_vec
    (dest_processor_id, sendvec, source_processor_id, recv,
     no_tag, any_tag, *this);
}



// This is both a declaration and definition for a new overloaded
// function template, so we have to re-specify the default arguments
template <typename T>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<std::vector<T> > & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<std::vector<T> > & recv,
                                       const MessageTag & /* send_tag */,
                                       const MessageTag & /* recv_tag */) const
{
  // FIXME - why aren't we honoring send_tag and recv_tag here?
  send_receive_vec_of_vec
    (dest_processor_id, sendvec, source_processor_id, recv,
     no_tag, any_tag, *this);
}




template <typename Context1, typename RangeIter, typename Context2,
          typename OutputIter, typename T>
inline void
Communicator::send_receive_packed_range (const unsigned int dest_processor_id,
                                         const Context1 * context1,
                                         RangeIter send_begin,
                                         const RangeIter send_end,
                                         const unsigned int source_processor_id,
                                         Context2 * context2,
                                         OutputIter out_iter,
                                         const T * output_type,
                                         const MessageTag & send_tag,
                                         const MessageTag & recv_tag) const
{
  LOG_SCOPE("send_receive()", "Parallel");

  Parallel::Request req;

  this->send_packed_range (dest_processor_id, context1, send_begin, send_end,
                           req, send_tag);

  this->receive_packed_range (source_processor_id, context2, out_iter,
                              output_type, recv_tag);

  req.wait();
}



template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 const T & sendval,
                                 std::vector<T> & recv) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->rank() == root_id)
    recv.resize(this->size());

  if (this->size() > 1)
    {
      LOG_SCOPE("gather()", "Parallel");

      StandardType<T> send_type(&sendval);

      libmesh_call_mpi
        (MPI_Gather(const_cast<T*>(&sendval), 1, send_type,
                    recv.empty() ? libmesh_nullptr : &recv[0], 1, send_type,
                    root_id, this->get()));
    }
  else
    recv[0] = sendval;
}



template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 const std::basic_string<T> & sendval,
                                 std::vector<std::basic_string<T> > & recv,
                                 const bool identical_buffer_sizes) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->rank() == root_id)
    recv.resize(this->size());

  if (this->size() > 1)
    {
      LOG_SCOPE ("gather()","Parallel");

      std::vector<int>
        sendlengths  (this->size(), 0),
        displacements(this->size(), 0);

      const int mysize = static_cast<int>(sendval.size());

      if (identical_buffer_sizes)
        sendlengths.assign(this->size(), mysize);
      else
        // first comm step to determine buffer sizes from all processors
        this->gather(root_id, mysize, sendlengths);

      // Find the total size of the final array and
      // set up the displacement offsets for each processor
      unsigned int globalsize = 0;
      for (unsigned int i=0; i < this->size(); ++i)
        {
          displacements[i] = globalsize;
          globalsize += sendlengths[i];
        }

      // monolithic receive buffer
      std::string r;
      if (this->rank() == root_id)
        r.resize(globalsize, 0);

      // and get the data from the remote processors.
      libmesh_call_mpi
        (MPI_Gatherv (const_cast<T*>(&sendval[0]),
                      mysize, StandardType<T>(),
                      this->rank() == root_id ? &r[0] : libmesh_nullptr,
                      &sendlengths[0], &displacements[0],
                      StandardType<T>(), root_id, this->get()));

      // slice receive buffer up
      if (this->rank() == root_id)
        for (unsigned int i=0; i != this->size(); ++i)
          recv[i] = r.substr(displacements[i], sendlengths[i]);
    }
  else
    recv[0] = sendval;
}



template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 std::vector<T> & r) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  std::vector<int>
    sendlengths  (this->size(), 0),
    displacements(this->size(), 0);

  const int mysize = static_cast<int>(r.size());
  this->allgather(mysize, sendlengths);

  LOG_SCOPE("gather()", "Parallel");

  // Find the total size of the final array and
  // set up the displacement offsets for each processor.
  unsigned int globalsize = 0;
  for (unsigned int i=0; i != this->size(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

  // Check for quick return
  if (globalsize == 0)
    return;

  // copy the input buffer
  std::vector<T> r_src(r);

  // now resize it to hold the global data
  // on the receiving processor
  if (root_id == this->rank())
    r.resize(globalsize);

  // and get the data from the remote processors
  libmesh_call_mpi
    (MPI_Gatherv (r_src.empty() ? libmesh_nullptr : &r_src[0], mysize,
                  StandardType<T>(), r.empty() ? libmesh_nullptr : &r[0],
                  &sendlengths[0], &displacements[0],
                  StandardType<T>(), root_id, this->get()));
}


template <typename T>
inline void Communicator::allgather(const T & sendval,
                                    std::vector<T> & recv) const
{
  LOG_SCOPE ("allgather()","Parallel");

  libmesh_assert(this->size());
  recv.resize(this->size());

  unsigned int comm_size = this->size();
  if (comm_size > 1)
    {
      StandardType<T> send_type(&sendval);

      libmesh_call_mpi
        (MPI_Allgather (const_cast<T*>(&sendval), 1, send_type, &recv[0], 1,
                        send_type, this->get()));
    }
  else if (comm_size > 0)
    recv[0] = sendval;
}



template <typename T>
inline void Communicator::allgather(const std::basic_string<T> & sendval,
                                    std::vector<std::basic_string<T> > & recv,
                                    const bool identical_buffer_sizes) const
{
  LOG_SCOPE ("allgather()","Parallel");

  libmesh_assert(this->size());
  recv.assign(this->size(), "");

  // serial case
  if (this->size() < 2)
    {
      recv.resize(1);
      recv[0] = sendval;
      return;
    }

  std::vector<int>
    sendlengths  (this->size(), 0),
    displacements(this->size(), 0);

  const int mysize = static_cast<int>(sendval.size());

  if (identical_buffer_sizes)
    sendlengths.assign(this->size(), mysize);
  else
    // first comm step to determine buffer sizes from all processors
    this->allgather(mysize, sendlengths);

  // Find the total size of the final array and
  // set up the displacement offsets for each processor
  unsigned int globalsize = 0;
  for (unsigned int i=0; i != this->size(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

  // Check for quick return
  if (globalsize == 0)
    return;

  // monolithic receive buffer
  std::string r(globalsize, 0);

  // and get the data from the remote processors.
  libmesh_call_mpi
    (MPI_Allgatherv (const_cast<T*>(mysize ? &sendval[0] : libmesh_nullptr),
                     mysize, StandardType<T>(),
                     &r[0], &sendlengths[0], &displacements[0],
                     StandardType<T>(), this->get()));

  // slice receive buffer up
  for (unsigned int i=0; i != this->size(); ++i)
    recv[i] = r.substr(displacements[i], sendlengths[i]);
}



template <typename T>
inline void Communicator::allgather(std::vector<T> & r,
                                    const bool identical_buffer_sizes) const
{
  if (this->size() < 2)
    return;

  LOG_SCOPE("allgather()", "Parallel");

  if (identical_buffer_sizes)
    {
      if (r.empty())
        return;

      libmesh_assert(this->verify(r.size()));

      std::vector<T> r_src(r.size()*this->size());
      r_src.swap(r);
      StandardType<T> send_type(&r_src[0]);

      libmesh_call_mpi
        (MPI_Allgather (&r_src[0], cast_int<int>(r_src.size()),
                        send_type, &r[0], cast_int<int>(r_src.size()),
                        send_type, this->get()));
      // libmesh_assert(this->verify(r));
      return;
    }

  std::vector<int>
    sendlengths  (this->size(), 0),
    displacements(this->size(), 0);

  const int mysize = static_cast<int>(r.size());
  this->allgather(mysize, sendlengths);

  // Find the total size of the final array and
  // set up the displacement offsets for each processor.
  unsigned int globalsize = 0;
  for (unsigned int i=0; i != this->size(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += sendlengths[i];
    }

  // Check for quick return
  if (globalsize == 0)
    return;

  // copy the input buffer
  std::vector<T> r_src(globalsize);
  r_src.swap(r);

  StandardType<T> send_type(&r[0]);

  // and get the data from the remote processors.
  // Pass NULL if our vector is empty.
  libmesh_call_mpi
    (MPI_Allgatherv (r_src.empty() ? libmesh_nullptr : &r_src[0], mysize,
                     send_type, &r[0], &sendlengths[0],
                     &displacements[0], send_type, this->get()));
}



template <typename T>
void Communicator::scatter(const std::vector<T> & data,
                           T & recv,
                           const unsigned int root_id) const
{
  libmesh_assert_less (root_id, this->size());

  // Do not allow the root_id to scatter a NULL vector.
  // That would leave recv in an indeterminate state.
  libmesh_assert (this->rank() != root_id || this->size() == data.size());

  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      recv = data[0];
      return;
    }

  LOG_SCOPE("scatter()", "Parallel");

  T * data_ptr = const_cast<T*>(data.empty() ? libmesh_nullptr : &data[0]);

  libmesh_call_mpi
    (MPI_Scatter (data_ptr, 1, StandardType<T>(data_ptr),
                  &recv, 1, StandardType<T>(&recv), root_id, this->get()));
}



template <typename T>
void Communicator::scatter(const std::vector<T> & data,
                           std::vector<T> & recv,
                           const unsigned int root_id) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      recv.assign(data.begin(), data.end());
      return;
    }

  LOG_SCOPE("scatter()", "Parallel");

  int recv_buffer_size;
  if (this->rank() == root_id)
    {
      libmesh_assert(data.size() % this->size() == 0);
      recv_buffer_size = data.size() / this->size();
    }

  this->broadcast(recv_buffer_size);
  recv.resize(recv_buffer_size);

  T * data_ptr = const_cast<T*>(data.empty() ? libmesh_nullptr : &data[0]);
  T * recv_ptr = recv.empty() ? libmesh_nullptr : &recv[0];

  libmesh_call_mpi
    (MPI_Scatter (data_ptr, recv_buffer_size, StandardType<T>(data_ptr),
                  recv_ptr, recv_buffer_size, StandardType<T>(recv_ptr), root_id, this->get()));
}



template <typename T>
void Communicator::scatter(const std::vector<T> & data,
                           const std::vector<int> counts,
                           std::vector<T> & recv,
                           const unsigned int root_id) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      libmesh_assert (counts.size() == this->size());
      recv.assign(data.begin(), data.begin() + counts[0]);
      return;
    }

  std::vector<int> displacements(this->size(), 0);
  if (root_id == this->rank())
    {
      libmesh_assert(counts.size() == this->size());

      // Create a displacements vector from the incoming counts vector
      unsigned int globalsize = 0;
      for (unsigned int i=0; i < this->size(); ++i)
        {
          displacements[i] = globalsize;
          globalsize += counts[i];
        }

      libmesh_assert(data.size() == globalsize);
    }

  LOG_SCOPE("scatter()", "Parallel");

  // Scatter the buffer sizes to size remote buffers
  int recv_buffer_size;
  this->scatter(counts, recv_buffer_size, root_id);
  recv.resize(recv_buffer_size);

  T * data_ptr = const_cast<T*>(data.empty() ? libmesh_nullptr : &data[0]);
  int * count_ptr = const_cast<int*>(counts.empty() ? libmesh_nullptr : &counts[0]);
  T * recv_ptr = recv.empty() ? libmesh_nullptr : &recv[0];

  // Scatter the non-uniform chunks
  libmesh_call_mpi
    (MPI_Scatterv (data_ptr, count_ptr, &displacements[0], StandardType<T>(data_ptr),
                   recv_ptr, recv_buffer_size, StandardType<T>(recv_ptr), root_id, this->get()));
}



template <typename T>
void Communicator::scatter(const std::vector<std::vector<T> > & data,
                           std::vector<T> & recv,
                           const unsigned int root_id,
                           const bool identical_buffer_sizes) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      libmesh_assert (data.size() == this->size());
      recv.assign(data[0].begin(), data[0].end());
      return;
    }

  std::vector<T> stacked_data;
  std::vector<int> counts;

  if (root_id == this->rank())
    {
      libmesh_assert (data.size() == this->size());

      if (!identical_buffer_sizes)
        counts.resize(this->size());

      for (std::size_t i=0; i < data.size(); ++i)
        {
          if (!identical_buffer_sizes)
            counts[i] = data[i].size();
#ifndef NDEBUG
          else
            // Check that buffer sizes are indeed equal
            libmesh_assert(!i || data[i-1].size() == data[i].size());
#endif
          std::copy(data[i].begin(), data[i].end(), std::back_inserter(stacked_data));
        }
    }

  if (identical_buffer_sizes)
    this->scatter(stacked_data, recv, root_id);
  else
    this->scatter(stacked_data, counts, recv, root_id);
}



template <typename T>
inline void Communicator::alltoall(std::vector<T> & buf) const
{
  if (this->size() < 2 || buf.empty())
    return;

  LOG_SCOPE("alltoall()", "Parallel");

  // the per-processor size.  this is the same for all
  // processors using MPI_Alltoall, could be variable
  // using MPI_Alltoallv
  const int size_per_proc =
    cast_int<int>(buf.size()/this->size());

  libmesh_assert_equal_to (buf.size()%this->size(), 0);

  libmesh_assert(this->verify(size_per_proc));

  std::vector<T> tmp(buf);

  StandardType<T> send_type(&tmp[0]);

  libmesh_call_mpi
    (MPI_Alltoall (&tmp[0], size_per_proc, send_type, &buf[0],
                   size_per_proc, send_type, this->get()));
}



template <typename T>
inline void Communicator::broadcast (T & data, const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  // Spread data to remote processors.
  libmesh_call_mpi
    (MPI_Bcast (&data, 1, StandardType<T>(&data), root_id,
                this->get()));
}


template <>
inline void Communicator::broadcast (bool & data, const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  // We don't want to depend on MPI-2 or C++ MPI, so we don't have
  // MPI::BOOL available
  char char_data = data;

  // Spread data to remote processors.
  libmesh_call_mpi
    (MPI_Bcast (&char_data, 1, StandardType<char>(&char_data),
                root_id, this->get()));

  data = char_data;
}


template <typename T>
inline void Communicator::broadcast (std::basic_string<T> & data,
                                     const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  std::size_t data_size = data.size();
  this->broadcast(data_size, root_id);

  std::vector<T> data_c(data_size);
#ifndef NDEBUG
  std::string orig(data);
#endif

  if (this->rank() == root_id)
    for (std::size_t i=0; i<data.size(); i++)
      data_c[i] = data[i];

  this->broadcast (data_c, root_id);

  data.assign(data_c.begin(), data_c.end());

#ifndef NDEBUG
  if (this->rank() == root_id)
    libmesh_assert_equal_to (data, orig);
#endif
}



template <typename T>
inline void Communicator::broadcast (std::vector<T> & data,
                                     const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  // and get the data from the remote processors.
  // Pass NULL if our vector is empty.
  T * data_ptr = data.empty() ? libmesh_nullptr : &data[0];

  libmesh_call_mpi
    (MPI_Bcast (data_ptr, cast_int<int>(data.size()),
                StandardType<T>(data_ptr), root_id, this->get()));
}


template <typename T>
inline void Communicator::broadcast (std::vector<std::basic_string<T> > & data,
                                     const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  std::size_t bufsize=0;
  if (root_id == this->rank())
    {
      for (std::size_t i=0; i<data.size(); ++i)
        bufsize += data[i].size() + 1;  // Add one for the string length word
    }
  this->broadcast(bufsize, root_id);

  // Here we use unsigned int to store up to 32-bit characters
  std::vector<unsigned int> temp; temp.reserve(bufsize);
  // Pack the strings
  if (root_id == this->rank())
    {
      for (std::size_t i=0; i<data.size(); ++i)
        {
          temp.push_back(cast_int<unsigned int>(data[i].size()));
          for (std::size_t j=0; j != data[i].size(); ++j)
            /**
             * The strings will be packed in one long array with the size of each
             * string preceeding the actual characters
             */
            temp.push_back(data[i][j]);
        }
    }
  else
    temp.resize(bufsize);

  // broad cast the packed strings
  this->broadcast(temp, root_id);

  // Unpack the strings
  if (root_id != this->rank())
    {
      data.clear();
      std::vector<unsigned int>::const_iterator iter = temp.begin();
      while (iter != temp.end())
        {
          std::size_t curr_len = *iter++;
          data.push_back(std::string(iter, iter+curr_len));
          iter += curr_len;
        }
    }
}




template <typename T>
inline void Communicator::broadcast (std::set<T> & data,
                                     const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  std::vector<T> vecdata;
  if (this->rank() == root_id)
    vecdata.assign(data.begin(), data.end());

  std::size_t vecsize = vecdata.size();
  this->broadcast(vecsize, root_id);
  if (this->rank() != root_id)
    vecdata.resize(vecsize);

  this->broadcast(vecdata, root_id);
  if (this->rank() != root_id)
    {
      data.clear();
      data.insert(vecdata.begin(), vecdata.end());
    }
}



template <typename T1, typename T2>
inline void Communicator::broadcast(std::map<T1, T2> & data,
                                    const unsigned int root_id) const
{
  if (this->size() == 1)
    {
      libmesh_assert (!this->rank());
      libmesh_assert (!root_id);
      return;
    }

  libmesh_assert_less (root_id, this->size());

  LOG_SCOPE("broadcast()", "Parallel");

  std::size_t data_size=data.size();
  this->broadcast(data_size, root_id);

  std::vector<T1> pair_first; pair_first.reserve(data_size);
  std::vector<T2> pair_second; pair_first.reserve(data_size);

  if (root_id == this->rank())
    {
      for (typename std::map<T1, T2>::const_iterator it = data.begin();
           it != data.end(); ++it)
        {
          pair_first.push_back(it->first);
          pair_second.push_back(it->second);
        }
    }
  else
    {
      pair_first.resize(data_size);
      pair_second.resize(data_size);
    }

  this->broadcast(pair_first, root_id);
  this->broadcast(pair_second, root_id);

  libmesh_assert(pair_first.size() == pair_first.size());

  if (this->rank() != root_id)
    {
      data.clear();
      for (std::size_t i=0; i<pair_first.size(); ++i)
        data[pair_first[i]] = pair_second[i];
    }
}



template <typename Context, typename OutputContext,
          typename Iter, typename OutputIter>
inline void Communicator::broadcast_packed_range(const Context * context1,
                                                 Iter range_begin,
                                                 const Iter range_end,
                                                 OutputContext * context2,
                                                 OutputIter out_iter,
                                                 const unsigned int root_id) const
{
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  do
    {
      // We will serialize variable size objects from *range_begin to
      // *range_end as a sequence of ints in this buffer
      std::vector<buffer_t> buffer;

      if (this->rank() == root_id)
        range_begin = Parallel::pack_range
          (context1, range_begin, range_end, buffer);

      // this->broadcast(vector) requires the receiving vectors to
      // already be the appropriate size
      std::size_t buffer_size = buffer.size();
      this->broadcast (buffer_size, root_id);

      // We continue until there's nothing left to broadcast
      if (!buffer_size)
        break;

      buffer.resize(buffer_size);

      // Broadcast the packed data
      this->broadcast (buffer, root_id);

      if (this->rank() != root_id)
        Parallel::unpack_range
          (buffer, context2, out_iter, (T*)libmesh_nullptr);
    } while (true);  // break above when we reach buffer_size==0
}


#else // LIBMESH_HAVE_MPI

template <typename T>
inline bool Communicator::verify(const T &) const { return true; }

template <typename T>
inline bool Communicator::semiverify(const T *) const { return true; }

template <typename T>
inline void Communicator::min(T &) const {}

template <typename T>
inline void Communicator::minloc(T &, unsigned int & min_id) const { min_id = 0; }

template <typename T>
inline void Communicator::minloc(std::vector<T> & r, std::vector<unsigned int> & min_id) const
{ for (std::size_t i=0; i!= r.size(); ++i) min_id[i] = 0; }

template <typename T>
inline void Communicator::max(T &) const {}

template <typename T>
inline void Communicator::maxloc(T &, unsigned int & max_id) const { max_id = 0; }

template <typename T>
inline void Communicator::maxloc(std::vector<T> & r, std::vector<unsigned int> & max_id) const
{ for (std::size_t i=0; i!= r.size(); ++i) max_id[i] = 0; }

template <typename T>
inline void Communicator::sum(T &) const {}

template <typename T>
inline void Communicator::set_union(T &) const {}

template <typename T>
inline void Communicator::set_union(T &, const unsigned int root_id) const
{ libmesh_assert_equal_to(root_id, 0); }

/**
 * We do not currently support probes on one processor without MPI.
 */
inline status Communicator::probe (const unsigned int,
                                   const MessageTag &) const
{ libmesh_not_implemented(); status s; return s; }

/**
 * We do not currently support sends on one processor without MPI.
 */
template <typename T>
inline void Communicator::send (const unsigned int,
                                const T &,
                                const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename T>
inline void Communicator::send (const unsigned int,
                                const T &,
                                Request &,
                                const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename T>
inline void Communicator::send (const unsigned int,
                                const T &,
                                const DataType &,
                                const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename T>
inline void Communicator::send (const unsigned int,
                                const T &,
                                const DataType &,
                                Request &,
                                const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename Context, typename Iter>
inline void Communicator::send_packed_range(const unsigned int,
                                            const Context *,
                                            Iter,
                                            const Iter,
                                            const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename Context, typename Iter>
inline void Communicator::send_packed_range (const unsigned int,
                                             const Context *,
                                             Iter,
                                             const Iter,
                                             Request &,
                                             const MessageTag &) const
{ libmesh_not_implemented(); }

/**
 * We do not currently support receives on one processor without MPI.
 */
template <typename T>
inline Status Communicator::receive (const unsigned int,
                                     T &,
                                     const MessageTag &) const
{ libmesh_not_implemented(); return Status(); }

template <typename T>
inline void Communicator::receive(const unsigned int,
                                  T &,
                                  Request &,
                                  const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename T>
inline Status Communicator::receive(const unsigned int,
                                    T &,
                                    const DataType &,
                                    const MessageTag &) const
{ libmesh_not_implemented(); return Status(); }

template <typename T>
inline void Communicator::receive(const unsigned int,
                                  T &,
                                  const DataType &,
                                  Request &,
                                  const MessageTag &) const
{ libmesh_not_implemented(); }

template <typename Context, typename OutputIter, typename T>
inline void
Communicator::receive_packed_range(const unsigned int,
                                   Context *,
                                   OutputIter,
                                   const T *,
                                   const MessageTag &) const
{ libmesh_not_implemented(); }

// template <typename Context, typename OutputIter>
// inline void Communicator::receive_packed_range(const unsigned int, Context *, OutputIter, Request &, const MessageTag &) const
// { libmesh_not_implemented(); }

/**
 * Send-receive data from one processor.
 */
template <typename T1, typename T2>
inline void Communicator::send_receive (const unsigned int send_tgt,
                                        const T1 & send_val,
                                        const unsigned int recv_source,
                                        T2 & recv_val,
                                        const MessageTag &,
                                        const MessageTag &) const
{
  libmesh_assert_equal_to (send_tgt, 0);
  libmesh_assert_equal_to (recv_source, 0);
  recv_val = send_val;
}

/**
 * Send-receive range-of-pointers from one processor.
 *
 * If you call this without MPI you might be making a mistake, but
 * we'll support it.
 */
template <typename Context1, typename RangeIter,
          typename Context2, typename OutputIter, typename T>
inline void
Communicator::send_receive_packed_range (const unsigned int dest_processor_id,
                                         const Context1 * context1,
                                         RangeIter send_begin,
                                         const RangeIter send_end,
                                         const unsigned int source_processor_id,
                                         Context2 * context2,
                                         OutputIter out_iter,
                                         const T * output_type,
                                         const MessageTag &,
                                         const MessageTag &) const
{
  // This makes no sense on one processor unless we're deliberately
  // sending to ourself.
  libmesh_assert_equal_to(dest_processor_id, 0);
  libmesh_assert_equal_to(source_processor_id, 0);

  // On one processor, we just need to pack the range and then unpack
  // it again.
  typedef typename std::iterator_traits<RangeIter>::value_type T1;
  typedef typename Parallel::Packing<T1>::buffer_type buffer_t;

  while (send_begin != send_end)
    {
      libmesh_assert_greater (std::distance(send_begin, send_end), 0);

      // We will serialize variable size objects from *range_begin to
      // *range_end as a sequence of ints in this buffer
      std::vector<buffer_t> buffer;

      const RangeIter next_send_begin = Parallel::pack_range
        (context1, send_begin, send_end, buffer);

      libmesh_assert_greater (std::distance(send_begin, next_send_begin), 0);

      send_begin = next_send_begin;

      Parallel::unpack_range
        (buffer, context2, out_iter, output_type);
    }
}

/**
 * Gather-to-root on one processor.
 */
template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 const T & send_val,
                                 std::vector<T> & recv_val) const
{
  libmesh_assert_equal_to (root_id, 0);
  recv_val.resize(1);
  recv_val[0] = send_val;
}

template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 const std::basic_string<T> & sendval,
                                 std::vector<std::basic_string<T> > & recv,
                                 const bool identical_buffer_sizes) const
{
  libmesh_assert_equal_to (root_id, 0);
  recv.resize(1);
  recv[0] = sendval;
}

template <typename T>
inline void Communicator::gather(const unsigned int root_id,
                                 std::vector<T> &) const
{ libmesh_assert_equal_to(root_id, 0); }

template <typename T>
inline void Communicator::allgather(const T & send_val,
                                    std::vector<T> & recv_val) const
{
  recv_val.resize(1);
  recv_val[0] = send_val;
}

template <typename T>
inline void Communicator::allgather(std::vector<T> &,
                                    const bool) const {}

template <typename T>
inline void Communicator::scatter(const std::vector<T> & data,
                                  T & recv,
                                  const unsigned int root_id) const
{
  libmesh_assert_equal_to (root_id, 0);
  recv = data[0];
}


template <typename T>
inline void Communicator::scatter(const std::vector<T> & data,
                                  std::vector<T> & recv,
                                  const unsigned int root_id) const
{
  libmesh_assert_equal_to (root_id, 0);
  recv.assign(data.begin(), data.end());
}


template <typename T>
inline void Communicator::scatter(const std::vector<T> & data,
                                  const std::vector<int> counts,
                                  std::vector<T> & recv,
                                  const unsigned int root_id) const
{
  libmesh_assert_equal_to (root_id, 0);
  libmesh_assert_equal_to (counts.size(), 1);
  recv.assign(data.begin(), data.begin() + counts[0]);
}


template <typename T>
inline void Communicator::scatter(const std::vector<std::vector<T> > & data,
                                  std::vector<T> & recv,
                                  const unsigned int root_id,
                                  const bool identical_buffer_sizes) const
{
  libmesh_assert_equal_to (root_id, 0);
  libmesh_assert_equal_to (data.size(), 1);
  recv.assign(data[0].begin(), data[0].end());
}



template <typename T>
inline void Communicator::alltoall(std::vector<T> &) const {}

template <typename T>
inline void Communicator::broadcast (T &,
                                     const unsigned int root_id) const
{ libmesh_assert_equal_to(root_id, 0); }

#endif // LIBMESH_HAVE_MPI

// Some of our methods are implemented indirectly via other
// MPI-encapsulated methods and the implementation works with or
// without MPI.

template <typename Context, typename Iter, typename OutputIter>
inline void Communicator::gather_packed_range(const unsigned int root_id,
                                              Context * context,
                                              Iter range_begin,
                                              const Iter range_end,
                                              OutputIter out_iter) const
{
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  bool nonempty_range = (range_begin != range_end);
  this->max(nonempty_range);

  while (nonempty_range)
    {
      // We will serialize variable size objects from *range_begin to
      // *range_end as a sequence of ints in this buffer
      std::vector<buffer_t> buffer;

      range_begin = Parallel::pack_range
        (context, range_begin, range_end, buffer);

      this->gather(root_id, buffer);

      Parallel::unpack_range
        (buffer, context, out_iter, (T*)(libmesh_nullptr));

      nonempty_range = (range_begin != range_end);
      this->max(nonempty_range);
    }
}


template <typename Context, typename Iter, typename OutputIter>
inline void Communicator::allgather_packed_range(Context * context,
                                                 Iter range_begin,
                                                 const Iter range_end,
                                                 OutputIter out_iter) const
{
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  bool nonempty_range = (range_begin != range_end);
  this->max(nonempty_range);

  while (nonempty_range)
    {
      // We will serialize variable size objects from *range_begin to
      // *range_end as a sequence of ints in this buffer
      std::vector<buffer_t> buffer;

      range_begin = Parallel::pack_range
        (context, range_begin, range_end, buffer);

      this->allgather(buffer, false);

      libmesh_assert(buffer.size());

      Parallel::unpack_range
        (buffer, context, out_iter, (T*)libmesh_nullptr);

      nonempty_range = (range_begin != range_end);
      this->max(nonempty_range);
    }
}


} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_PARALLEL_IMPLEMENTATION_H
