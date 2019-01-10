// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <complex>
#include <cstddef>
#include <cstring> // memcpy
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <type_traits>

namespace libMesh {
namespace Parallel {

#ifdef LIBMESH_HAVE_MPI

/**
 * Templated function to return the appropriate MPI datatype
 * for use with built-in C types when combined with an int
 */
template <typename T>
inline data_type dataplusint_type();

#endif // LIBMESH_HAVE_MPI

/**
 * Types combined with an int
 */
template <typename T>
class DataPlusInt
{
public:
  T val;
  int rank;
};

} // namespace Parallel

} // namespace libMesh


// Anonymous namespace for helper functions
namespace {

// Internal helper function to create vector<something_usable> from
// vector<bool> for compatibility with MPI bitwise operations
template <typename T, typename A1, typename A2>
inline void pack_vector_bool(const std::vector<bool,A1> & vec_in,
                             std::vector<T,A2> & vec_out)
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
template <typename T, typename A1, typename A2>
inline void unpack_vector_bool(const std::vector<T,A1> & vec_in,
                               std::vector<bool,A2> & vec_out)
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
template <typename T1, typename T2, typename A1, typename A2, typename A3, typename A4>
inline void send_receive_vec_of_vec(const unsigned int dest_processor_id,
                                    const std::vector<std::vector<T1,A1>,A2> & send,
                                    const unsigned int source_processor_id,
                                    std::vector<std::vector<T2,A3>,A4> & recv,
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

  libMesh::Parallel::Request request;
  comm.send (dest_processor_id, send, request, send_tag);
  comm.receive (source_processor_id, recv, recv_tag);
  request.wait();
}

#endif // LIBMESH_HAVE_MPI

} // Anonymous namespace



namespace libMesh
{

namespace Parallel
{

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



inline status Communicator::probe (const unsigned int src_processor_id,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("probe()", "Parallel");

  status stat;

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  libmesh_call_mpi
    (MPI_Probe (src_processor_id, tag.value(), this->get(), &stat));

  return stat;
}

template<typename T>
inline Status Communicator::packed_range_probe (const unsigned int src_processor_id,
                                                const MessageTag & tag,
                                                bool & flag) const
{
  LOG_SCOPE("packed_range_probe()", "Parallel");

  libmesh_experimental();

  Status stat((StandardType<typename Packing<T>::buffer_type>()));

  int int_flag;

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  libmesh_call_mpi(MPI_Iprobe(src_processor_id,
                              tag.value(),
                              this->get(),
                              &int_flag,
                              stat.get()));

  flag = int_flag;

  return stat;
}


template<typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::basic_string<T> & buf,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = buf.empty() ? nullptr : const_cast<T *>(buf.data());

  libmesh_assert_less(dest_processor_id, this->size());

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

  T * dataptr = buf.empty() ? nullptr : const_cast<T *>(buf.data());

  libmesh_assert_less(dest_processor_id, this->size());

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (dataptr,
                               cast_int<int>(buf.size()),
                               StandardType<T>(dataptr),
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));

  // The MessageTag should stay registered for the Request lifetime
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceTag(tag));
}



template <typename T>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const T & buf,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  T * dataptr = const_cast<T*> (&buf);

  libmesh_assert_less(dest_processor_id, this->size());

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

  libmesh_assert_less(dest_processor_id, this->size());

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (dataptr,
                               1,
                               StandardType<T>(dataptr),
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));

  // The MessageTag should stay registered for the Request lifetime
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceTag(tag));
}



template <typename T, typename C, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T,C,A> & buf,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), tag);
}



template <typename T, typename C, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T,C,A> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), req, tag);
}



template <typename T, typename C, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T,C,A> & buf,
                                const DataType & type,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  std::vector<T> vecbuf(buf.begin(), buf.end());
  this->send(dest_processor_id, vecbuf, type, tag);
}



template <typename T, typename C, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::set<T,C,A> & buf,
                                const DataType & type,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  // Allocate temporary buffer on the heap so it lives until after
  // the non-blocking send completes
  std::vector<T> * vecbuf =
    new std::vector<T,A>(buf.begin(), buf.end());

  // Make the Request::wait() handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T,A>>(vecbuf));

  this->send(dest_processor_id, *vecbuf, type, req, tag);
}



template <typename T, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T,A> & buf,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? nullptr : &buf.front()), tag);
}



template <typename T, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T,A> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>(buf.empty() ? nullptr : &buf.front()), req, tag);
}



template <typename T, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T,A> & buf,
                                const DataType & type,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Ssend : MPI_Send) (buf.empty() ? nullptr : const_cast<T*>(buf.data()),
                             cast_int<int>(buf.size()),
                             type,
                             dest_processor_id,
                             tag.value(),
                             this->get()));
}



template <typename T, typename A>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<T,A> & buf,
                                const DataType & type,
                                Request & req,
                                const MessageTag & tag) const
{
  LOG_SCOPE("send()", "Parallel");

  libmesh_assert_less(dest_processor_id, this->size());

  libmesh_call_mpi
    (((this->send_mode() == SYNCHRONOUS) ?
      MPI_Issend : MPI_Isend) (buf.empty() ? nullptr : const_cast<T*>(buf.data()),
                               cast_int<int>(buf.size()),
                               type,
                               dest_processor_id,
                               tag.value(),
                               this->get(),
                               req.get()));

  // The MessageTag should stay registered for the Request lifetime
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceTag(tag));
}



template <typename T, typename A1, typename A2>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<std::vector<T,A1>,A2> & buf,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>((buf.empty() || buf.front().empty()) ?
                             nullptr : &(buf.front().front())), tag);
}



template <typename T, typename A1, typename A2>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<std::vector<T,A1>,A2> & buf,
                                Request & req,
                                const MessageTag & tag) const
{
  this->send(dest_processor_id, buf,
             StandardType<T>((buf.empty() || buf.front().empty()) ?
                             nullptr : &(buf.front().front())), req, tag);
}



template <typename T, typename A1, typename A2>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<std::vector<T,A1>,A2> & buf,
                                const DataType & type,
                                const MessageTag & tag) const
{
  // We'll avoid redundant code (at the cost of using heap rather
  // than stack buffer allocation) by reusing the non-blocking send
  Request req;
  this->send(dest_processor_id, buf, type, req, tag);
  req.wait();
}



template <typename T, typename A1, typename A2>
inline void Communicator::send (const unsigned int dest_processor_id,
                                const std::vector<std::vector<T,A1>,A2> & send_vecs,
                                const DataType & type,
                                Request & req,
                                const MessageTag & tag) const
{
  // temporary buffer - this will be sized in bytes
  // and manipulated with MPI_Pack
  std::vector<char> * sendbuf = new std::vector<char>();

  // figure out how many bytes we need to pack all the data
  int packedsize=0;

  // The outer buffer size
  libmesh_call_mpi
    (MPI_Pack_size (1,
                    libMesh::Parallel::StandardType<unsigned int>(),
                    this->get(),
                    &packedsize));

  int sendsize = packedsize;

  const std::size_t n_vecs = send_vecs.size();

  for (std::size_t i = 0; i != n_vecs; ++i)
    {
      // The size of the ith inner buffer
      libmesh_call_mpi
        (MPI_Pack_size (1,
                        libMesh::Parallel::StandardType<unsigned int>(),
                        this->get(),
                        &packedsize));

      sendsize += packedsize;

      // The data for each inner buffer
      libmesh_call_mpi
        (MPI_Pack_size (libMesh::cast_int<int>(send_vecs[i].size()), type,
                        this->get(), &packedsize));

      sendsize += packedsize;
    }

  libmesh_assert (sendsize /* should at least be 1! */);
  sendbuf->resize (sendsize);

  // Pack the send buffer
  int pos=0;

  // ... the size of the outer buffer
  const int mpi_n_vecs = libMesh::cast_int<int>(n_vecs);

  libmesh_call_mpi
    (MPI_Pack (&mpi_n_vecs, 1,
               libMesh::Parallel::StandardType<unsigned int>(),
               sendbuf->data(), sendsize, &pos, this->get()));

  for (std::size_t i = 0; i != n_vecs; ++i)
    {
      // ... the size of the ith inner buffer
      const int subvec_size = libMesh::cast_int<int>(send_vecs[i].size());

      libmesh_call_mpi
        (MPI_Pack (&subvec_size, 1, libMesh::Parallel::StandardType<unsigned int>(),
                   sendbuf->data(), sendsize, &pos, this->get()));

      // ... the contents of the ith inner buffer
      if (!send_vecs[i].empty())
        libmesh_call_mpi
          (MPI_Pack (const_cast<T*>(send_vecs[i].data()),
                     libMesh::cast_int<int>(subvec_size), type,
                     sendbuf->data(), sendsize, &pos, this->get()));
    }

  libmesh_assert_equal_to (pos, sendsize);

  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<char>> (sendbuf));

  this->send (dest_processor_id, *sendbuf, MPI_PACKED, req, tag);
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
        (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t>>
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
        (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t>>
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
     std::back_insert_iterator<std::basic_string<T>>>
     (tempbuf, std::back_inserter(buf)));

  // Make the Request::wait() then handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T>>(tempbuf));

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

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

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

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  libmesh_call_mpi
    (MPI_Irecv (&buf, 1, StandardType<T>(&buf), src_processor_id,
                tag.value(), this->get(), req.get()));

  // The MessageTag should stay registered for the Request lifetime
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceTag(tag));
}



template <typename T, typename C, typename A>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::set<T,C,A> & buf,
                                     const MessageTag & tag) const
{
  return this->receive
    (src_processor_id, buf,
     StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), tag);
}



/*
 * No non-blocking receives of std::set until we figure out how to
 * resize the temporary buffer
 */
#if 0
template <typename T, typename C, typename A>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::set<T,C,A> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  this->receive (src_processor_id, buf,
                 StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), req, tag);
}
#endif // 0



template <typename T, typename C, typename A>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::set<T,C,A> & buf,
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
template <typename T, typename C, typename A>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::set<T,C,A> & buf,
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
     std::insert_iterator<std::set<T,C,A>>>
     (*vecbuf, std::inserter(buf,buf.end())));

  // Make the Request::wait() then handle deleting the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitDeleteBuffer<std::vector<T>>(vecbuf));

  this->receive(src_processor_id, *vecbuf, type, req, tag);
}
#endif // 0



template <typename T, typename A>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<T,A> & buf,
                                     const MessageTag & tag) const
{
  return this->receive
    (src_processor_id, buf,
     StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), tag);
}



template <typename T, typename A>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<T,A> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  this->receive (src_processor_id, buf,
                 StandardType<T>(buf.empty() ? nullptr : &(*buf.begin())), req, tag);
}



template <typename T, typename A>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<T,A> & buf,
                                     const DataType & type,
                                     const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  // Get the status of the message, explicitly provide the
  // datatype so we can later query the size
  Status stat(this->probe(src_processor_id, tag), type);

  buf.resize(stat.size());

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  // Use stat.source() and stat.tag() in the receive - if
  // src_processor_id is or tag is "any" then we want to be sure we
  // try to receive the same message we just probed.
  libmesh_call_mpi
    (MPI_Recv (buf.empty() ? nullptr : buf.data(),
               cast_int<int>(buf.size()), type, stat.source(),
               stat.tag(), this->get(), stat.get()));

  libmesh_assert_equal_to (stat.size(), buf.size());

  return stat;
}



template <typename T, typename A>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<T,A> & buf,
                                   const DataType & type,
                                   Request & req,
                                   const MessageTag & tag) const
{
  LOG_SCOPE("receive()", "Parallel");

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  libmesh_call_mpi
    (MPI_Irecv (buf.empty() ? nullptr : buf.data(),
                cast_int<int>(buf.size()), type, src_processor_id,
                tag.value(), this->get(), req.get()));

  // The MessageTag should stay registered for the Request lifetime
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceTag(tag));
}



template <typename T, typename A1, typename A2>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<std::vector<T,A1>,A2> & buf,
                                     const MessageTag & tag) const
{
  return this->receive
    (src_processor_id, buf,
     StandardType<T>((buf.empty() || buf.front().empty()) ?
                     nullptr : &(buf.front().front())), tag);
}



template <typename T, typename A1, typename A2>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<std::vector<T,A1>,A2> & buf,
                                   Request & req,
                                   const MessageTag & tag) const
{
  this->receive (src_processor_id, buf,
                 StandardType<T>((buf.empty() || buf.front().empty()) ?
                                 nullptr : &(buf.front().front())), req, tag);
}



template <typename T, typename A1, typename A2>
inline Status Communicator::receive (const unsigned int src_processor_id,
                                     std::vector<std::vector<T,A1>,A2> & recv,
                                     const DataType & type,
                                     const MessageTag & tag) const
{
  // temporary buffer - this will be sized in bytes
  // and manipulated with MPI_Unpack
  std::vector<char> recvbuf;

  Status stat = this->receive (src_processor_id, recvbuf, MPI_PACKED, tag);

  // We should at least have one header datum, for outer vector size
  libmesh_assert (!recvbuf.empty());

  // Unpack the received buffer
  int bufsize = libMesh::cast_int<int>(recvbuf.size());
  int recvsize, pos=0;
  libmesh_call_mpi
    (MPI_Unpack (recvbuf.data(), bufsize, &pos,
                 &recvsize, 1, libMesh::Parallel::StandardType<unsigned int>(),
                 this->get()));

  // ... size the outer buffer
  recv.resize (recvsize);

  const std::size_t n_vecs = recvsize;
  for (std::size_t i = 0; i != n_vecs; ++i)
    {
      int subvec_size;

      libmesh_call_mpi
        (MPI_Unpack (recvbuf.data(), bufsize, &pos,
                     &subvec_size, 1,
                     libMesh::Parallel::StandardType<unsigned int>(),
                     this->get()));

      // ... size the inner buffer
      recv[i].resize (subvec_size);

      // ... unpack the inner buffer if it is not empty
      if (!recv[i].empty())
        libmesh_call_mpi
          (MPI_Unpack (recvbuf.data(), bufsize, &pos, recv[i].data(),
                       subvec_size, type, this->get()));
    }

  return stat;
}



// FIXME - non-blocking receive of vector-of-vectors is currently unimplemented
/*
template <typename T, typename A1, typename A2>
inline void Communicator::receive (const unsigned int src_processor_id,
                                   std::vector<std::vector<T,A1>,A2> & buf,
                                   const DataType & type,
                                   Request & req,
                                   const MessageTag & tag) const
{
}
*/


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
//     (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t>>(buffer));
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
    (new Parallel::PostWaitDeleteBuffer<std::vector<buffer_t>>(buffer));
}



template <typename T1, typename T2, typename A1, typename A2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T1,A1> & sendvec,
                                       const DataType & type1,
                                       const unsigned int source_processor_id,
                                       std::vector<T2,A2> & recv,
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

  libmesh_assert_less(dest_processor_id, this->size());
  libmesh_assert(source_processor_id < this->size() ||
                 source_processor_id == any_source);

  // MPI_STATUS_IGNORE is from MPI-2; using it with some versions of
  // MPICH may cause a crash:
  // https://bugzilla.mcs.anl.gov/globus/show_bug.cgi?id=1798
  libmesh_call_mpi
    (MPI_Sendrecv(const_cast<T1*>(&sendvec), 1, StandardType<T1>(&sendvec),
                  dest_processor_id, send_tag.value(), &recv, 1,
                  StandardType<T2>(&recv), source_processor_id,
                  recv_tag.value(), this->get(), MPI_STATUS_IGNORE));
}



// This is both a declaration and definition for a new overloaded
// function template, so we have to re-specify the default
// arguments.
//
// We specialize on the T1==T2 case so that we can handle
// send_receive-to-self with a plain copy rather than going through
// MPI.
template <typename T, typename A>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T,A> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<T,A> & recv,
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
    (recv.empty() ? nullptr : recv.data()) : sendvec.data();

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
template <typename T1, typename T2, typename A1, typename A2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<T1,A1> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<T2,A2> & recv,
                                       const MessageTag & send_tag,
                                       const MessageTag & recv_tag) const
{
  // Call the user-defined type version with automatic
  // type conversion based on template argument:
  this->send_receive (dest_processor_id, sendvec,
                      StandardType<T1>(sendvec.empty() ? nullptr : sendvec.data()),
                      source_processor_id, recv,
                      StandardType<T2>(recv.empty() ? nullptr : recv.data()),
                      send_tag, recv_tag);
}




template <typename T1, typename T2, typename A1, typename A2, typename A3, typename A4>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<std::vector<T1,A1>,A2> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<std::vector<T2,A3>,A4> & recv,
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
template <typename T, typename A1, typename A2>
inline void Communicator::send_receive(const unsigned int dest_processor_id,
                                       const std::vector<std::vector<T,A1>,A2> & sendvec,
                                       const unsigned int source_processor_id,
                                       std::vector<std::vector<T,A1>,A2> & recv,
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



template <typename Context, typename Iter>
inline void Communicator::nonblocking_send_packed_range (const unsigned int dest_processor_id,
                                                         const Context * context,
                                                         Iter range_begin,
                                                         const Iter range_end,
                                                         Request & req,
                                                         std::shared_ptr<std::vector<typename Parallel::Packing<typename std::iterator_traits<Iter>::value_type>::buffer_type>> & buffer,
                                                         const MessageTag & tag) const
{
  libmesh_experimental();

  // Allocate a buffer on the heap so we don't have to free it until
  // after the Request::wait()
  typedef typename std::iterator_traits<Iter>::value_type T;
  typedef typename Parallel::Packing<T>::buffer_type buffer_t;

  if (range_begin != range_end)
    {
      if (buffer == nullptr)
        buffer = std::make_shared<std::vector<buffer_t>>();
      else
        buffer->clear();

      range_begin =
        Parallel::pack_range(context,
                             range_begin,
                             range_end,
                             *buffer,
                             // MPI-2 can only use integers for size
                             std::numeric_limits<int>::max());

      if (range_begin != range_end)
        libmesh_error_msg("Non-blocking packed range sends cannot exceed " << std::numeric_limits<int>::max() << "in size");

      // Make it dereference the shared pointer (possibly freeing the buffer)
      req.add_post_wait_work
        (new Parallel::PostWaitDereferenceSharedPtr<std::vector<buffer_t>>(buffer));

      // Non-blocking send of the buffer
      this->send(dest_processor_id, *buffer, req, tag);
    }
}



template <typename T, typename A>
inline void Communicator::allgather(const std::basic_string<T> & sendval,
                                    std::vector<std::basic_string<T>,A> & recv,
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
    (MPI_Allgatherv (const_cast<T*>(mysize ? sendval.data() : nullptr),
                     mysize, StandardType<T>(),
                     &r[0], sendlengths.data(), displacements.data(),
                     StandardType<T>(), this->get()));

  // slice receive buffer up
  for (unsigned int i=0; i != this->size(); ++i)
    recv[i] = r.substr(displacements[i], sendlengths[i]);
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

  libmesh_assert_less(root_id, this->size());

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



template <typename T, typename A>
inline void Communicator::broadcast (std::vector<T,A> & data,
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
  // Pass nullptr if our vector is empty.
  T * data_ptr = data.empty() ? nullptr : data.data();

  libmesh_assert_less(root_id, this->size());

  libmesh_call_mpi
    (MPI_Bcast (data_ptr, cast_int<int>(data.size()),
                StandardType<T>(data_ptr), root_id, this->get()));
}


template <typename T, typename A>
inline void Communicator::broadcast (std::vector<std::basic_string<T>,A> & data,
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
             * string preceding the actual characters
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
      typename std::vector<unsigned int>::const_iterator iter = temp.begin();
      while (iter != temp.end())
        {
          std::size_t curr_len = *iter++;
          data.push_back(std::string(iter, iter+curr_len));
          iter += curr_len;
        }
    }
}




template <typename T, typename C, typename A>
inline void Communicator::broadcast (std::set<T,C,A> & data,
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



template <typename T1, typename T2, typename C, typename A>
inline void Communicator::broadcast(std::map<T1,T2,C,A> & data,
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
      for (const auto & pr : data)
        {
          pair_first.push_back(pr.first);
          pair_second.push_back(pr.second);
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
          (buffer, context2, out_iter, (T*)nullptr);
    } while (true);  // break above when we reach buffer_size==0
}



template <typename Context, typename OutputIter, typename T>
inline void Communicator::nonblocking_receive_packed_range (const unsigned int src_processor_id,
                                                            Context * context,
                                                            OutputIter out,
                                                            const T * /*output_type*/,
                                                            Request & req,
                                                            Status & stat,
                                                            std::shared_ptr<std::vector<typename Parallel::Packing<T>::buffer_type>> & buffer,
                                                            const MessageTag & tag) const
{
  libmesh_experimental();

  // If they didn't pass in a buffer - let's make one
  if (buffer == nullptr)
    buffer = std::make_shared<std::vector<typename Parallel::Packing<T>::buffer_type>>();
  else
    buffer->clear();

  // Receive serialized variable size objects as a sequence of
  // buffer_t.
  // Allocate a buffer on the heap so we don't have to free it until
  // after the Request::wait()
  buffer->resize(stat.size());
  this->receive(src_processor_id, *buffer, req, tag);

  // Make the Request::wait() handle unpacking the buffer
  req.add_post_wait_work
    (new Parallel::PostWaitUnpackBuffer<std::vector<typename Parallel::Packing<T>::buffer_type>, Context, OutputIter, T>(*buffer, context, out));

  // Make it dereference the shared pointer (possibly freeing the buffer)
  req.add_post_wait_work
    (new Parallel::PostWaitDereferenceSharedPtr<std::vector<typename Parallel::Packing<T>::buffer_type>>(buffer));
}


#else // LIBMESH_HAVE_MPI

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
Communicator::send_receive_packed_range
  (const unsigned int libmesh_dbg_var(dest_processor_id),
   const Context1 * context1,
   RangeIter send_begin,
   const RangeIter send_end,
   const unsigned int libmesh_dbg_var(source_processor_id),
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

#endif // LIBMESH_HAVE_MPI

// Some of our methods are implemented indirectly via other
// MPI-encapsulated methods and the implementation works with or
// without MPI.
//
// Other methods have a "this->size() == 1" shortcut which still
// applies in the without-MPI case, and the libmesh_call_mpi macro
// becomes a no-op so wrapped MPI methods still compile without MPI

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

  static_assert(Attributes<T>::has_min_max,
                "Tried to verify an unverifiable type");

  return true;
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

  static_assert(Attributes<T>::has_min_max,
                "Tried to semiverify an unverifiable type");

  return true;
}



inline bool Communicator::verify(const bool & r) const
{
  const unsigned char rnew = r;
  return this->verify(rnew);
}



inline bool Communicator::semiverify(const bool * r) const
{
  if (r)
    {
      const unsigned char rnew = *r;
      return this->semiverify(&rnew);
    }

  const unsigned char * rptr = nullptr;
  return this->semiverify(rptr);
}



template <typename T, typename A>
inline bool Communicator::semiverify(const std::vector<T,A> * r) const
{
  if (this->size() > 1 && Attributes<T>::has_min_max == true)
    {
      std::size_t rsize = r ? r->size() : 0;
      std::size_t * psize = r ? &rsize : nullptr;

      if (!this->semiverify(psize))
        return false;

      this->max(rsize);

      std::vector<T,A> tempmin, tempmax;
      if (r)
        {
          tempmin = tempmax = *r;
        }
      else
        {
          tempmin.resize(rsize);
          tempmax.resize(rsize);
          Attributes<std::vector<T,A>>::set_highest(tempmin);
          Attributes<std::vector<T,A>>::set_lowest(tempmax);
        }
      this->min(tempmin);
      this->max(tempmax);
      bool invalid = r && ((*r != tempmin) ||
                           (*r != tempmax));
      this->max(invalid);
      return !invalid;
    }

  static_assert(Attributes<T>::has_min_max,
                "Tried to semiverify a vector of an unverifiable type");

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
      std::size_t * psize = r ? &rsize : nullptr;

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

      std::vector<short int> * ptemp = r ? &temp: nullptr;

      return this->semiverify(ptemp);
    }
  return true;
}



template <typename T>
inline void Communicator::min(T & libmesh_mpi_var(r)) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("min(scalar)", "Parallel");

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &r, 1,
                        StandardType<T>(&r), OpFunction<T>::min(),
                        this->get()));
    }
}



inline void Communicator::min(bool & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("min(bool)", "Parallel");

      unsigned int temp = r;
      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &temp, 1,
                        StandardType<unsigned int>(),
                        OpFunction<unsigned int>::min(),
                        this->get()));
      r = temp;
    }
}


template <typename T, typename A>
inline void Communicator::min(std::vector<T,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("min(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, r.data(),
                        cast_int<int>(r.size()),
                        StandardType<T>(r.data()),
                        OpFunction<T>::min(),
                        this->get()));
    }
}


template <typename A>
inline void Communicator::min(std::vector<bool,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("min(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<unsigned int> ruint;
      pack_vector_bool(r, ruint);
      std::vector<unsigned int> temp(ruint.size());
      libmesh_call_mpi
        (MPI_Allreduce (ruint.data(), temp.data(),
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
      libmesh_ignore(data_in); // unused ifndef LIBMESH_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &data_in, 1, dataplusint_type<T>(),
                        OpFunction<T>::min_location(), this->get()));
      r = data_in.val;
      min_id = data_in.rank;
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
      libmesh_ignore(data_in); // unused ifndef LIBMESH_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type<int>(),
                        OpFunction<int>::min_location(), this->get()));
      r = data_out.val;
      min_id = data_out.rank;
    }
  else
    min_id = this->rank();
}


template <typename T, typename A1, typename A2>
inline void Communicator::minloc(std::vector<T,A1> & r,
                                 std::vector<unsigned int,A2> & min_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("minloc(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<T>> data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<T>> data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (data_in.data(), data_out.data(),
                        cast_int<int>(r.size()),
                        dataplusint_type<T>(),
                        OpFunction<T>::min_location(), this->get()));
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


template <typename A1, typename A2>
inline void Communicator::minloc(std::vector<bool,A1> & r,
                                 std::vector<unsigned int,A2> & min_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("minloc(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<int>> data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<int>> data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (data_in.data(), data_out.data(),
                        cast_int<int>(r.size()),
                        StandardType<int>(),
                        OpFunction<int>::min_location(), this->get()));
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
inline void Communicator::max(T & libmesh_mpi_var(r)) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("max(scalar)", "Parallel");

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &r, 1,
                        StandardType<T>(&r),
                        OpFunction<T>::max(),
                        this->get()));
    }
}


inline void Communicator::max(bool & r) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("max(bool)", "Parallel");

      unsigned int temp = r;
      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &temp, 1,
                        StandardType<unsigned int>(),
                        OpFunction<unsigned int>::max(),
                        this->get()));
      r = temp;
    }
}


template <typename T, typename A>
inline void Communicator::max(std::vector<T,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("max(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, r.data(),
                        cast_int<int>(r.size()),
                        StandardType<T>(r.data()),
                        OpFunction<T>::max(),
                        this->get()));
    }
}


template <typename A>
inline void Communicator::max(std::vector<bool,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("max(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<unsigned int> ruint;
      pack_vector_bool(r, ruint);
      std::vector<unsigned int> temp(ruint.size());
      libmesh_call_mpi
        (MPI_Allreduce (ruint.data(), temp.data(),
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
      libmesh_ignore(data_in); // unused ifndef LIBMESH_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &data_in, 1,
                        dataplusint_type<T>(),
                        OpFunction<T>::max_location(),
                        this->get()));
      r = data_in.val;
      max_id = data_in.rank;
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
      libmesh_ignore(data_in); // unused ifndef LIBMESH_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out;
      libmesh_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type<int>(),
                        OpFunction<int>::max_location(),
                        this->get()));
      r = data_out.val;
      max_id = data_out.rank;
    }
  else
    max_id = this->rank();
}


template <typename T, typename A1, typename A2>
inline void Communicator::maxloc(std::vector<T,A1> & r,
                                 std::vector<unsigned int,A2> & max_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("maxloc(vector)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<T>> data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<T>> data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (data_in.data(), data_out.data(),
                        cast_int<int>(r.size()),
                        dataplusint_type<T>(),
                        OpFunction<T>::max_location(),
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


template <typename A1, typename A2>
inline void Communicator::maxloc(std::vector<bool,A1> & r,
                                 std::vector<unsigned int,A2> & max_id) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("maxloc(vector<bool>)", "Parallel");

      libmesh_assert(this->verify(r.size()));

      std::vector<DataPlusInt<int>> data_in(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        {
          data_in[i].val  = r[i];
          data_in[i].rank = this->rank();
        }
      std::vector<DataPlusInt<int>> data_out(r.size());
      libmesh_call_mpi
        (MPI_Allreduce (data_in.data(), data_out.data(),
                        cast_int<int>(r.size()),
                        StandardType<int>(),
                        OpFunction<int>::max_location(),
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
inline void Communicator::sum(T & libmesh_mpi_var(r)) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &r, 1,
                        StandardType<T>(&r),
                        OpFunction<T>::sum(),
                        this->get()));
    }
}


template <typename T, typename A>
inline void Communicator::sum(std::vector<T,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_assert(this->verify(r.size()));

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, r.data(),
                        cast_int<int>(r.size()),
                        StandardType<T>(r.data()),
                        OpFunction<T>::sum(),
                        this->get()));
    }
}


// We still do function overloading for complex sums - in a perfect
// world we'd have a StandardSumOp to go along with StandardType...
template <typename T>
inline void Communicator::sum(std::complex<T> & libmesh_mpi_var(r)) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &r, 2,
                        StandardType<T>(),
                        OpFunction<T>::sum(),
                        this->get()));
    }
}


template <typename T, typename A>
inline void Communicator::sum(std::vector<std::complex<T>,A> & r) const
{
  if (this->size() > 1 && !r.empty())
    {
      LOG_SCOPE("sum()", "Parallel");

      libmesh_assert(this->verify(r.size()));

      libmesh_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, r.data(),
                        cast_int<int>(r.size() * 2),
                        StandardType<T>(nullptr),
                        OpFunction<T>::sum(), this->get()));
    }
}



template <typename T, typename C, typename A>
inline void Communicator::set_union(std::set<T,C,A> & data,
                                    const unsigned int root_id) const
{
  if (this->size() > 1)
    {
      std::vector<T> vecdata(data.begin(), data.end());
      this->gather(root_id, vecdata);
      if (this->rank() == root_id)
        data.insert(vecdata.begin(), vecdata.end());
    }
}



template <typename T, typename C, typename A>
inline void Communicator::set_union(std::set<T,C,A> & data) const
{
  if (this->size() > 1)
    {
      std::vector<T> vecdata(data.begin(), data.end());
      this->allgather(vecdata, false);
      data.insert(vecdata.begin(), vecdata.end());
    }
}



template <typename T1, typename T2, typename C, typename A>
inline void Communicator::set_union(std::map<T1,T2,C,A> & data,
                                    const unsigned int root_id) const
{
  if (this->size() > 1)
    {
      std::vector<std::pair<T1,T2>> vecdata(data.begin(), data.end());
      this->gather(root_id, vecdata);
      if (this->rank() == root_id)
        data.insert(vecdata.begin(), vecdata.end());
    }
}



template <typename T1, typename T2, typename C, typename A>
inline void Communicator::set_union(std::map<T1,T2,C,A> & data) const
{
  if (this->size() > 1)
    {
      std::vector<std::pair<T1,T2>> vecdata(data.begin(), data.end());
      this->allgather(vecdata, false);
      data.insert(vecdata.begin(), vecdata.end());
    }
}



template <typename T, typename A>
inline void Communicator::gather(const unsigned int root_id,
                                 const T & sendval,
                                 std::vector<T,A> & recv) const
{
  libmesh_assert_less (root_id, this->size());

  if (this->rank() == root_id)
    recv.resize(this->size());

  if (this->size() > 1)
    {
      LOG_SCOPE("gather()", "Parallel");

      StandardType<T> send_type(&sendval);

      libmesh_assert_less(root_id, this->size());

      libmesh_call_mpi
        (MPI_Gather(const_cast<T*>(&sendval), 1, send_type,
                    recv.empty() ? nullptr : recv.data(), 1, send_type,
                    root_id, this->get()));
    }
  else
    recv[0] = sendval;
}



template <typename T, typename A>
inline void Communicator::gather(const unsigned int root_id,
                                 std::vector<T,A> & r) const
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
  std::vector<T,A> r_src(r);

  // now resize it to hold the global data
  // on the receiving processor
  if (root_id == this->rank())
    r.resize(globalsize);

  libmesh_assert_less(root_id, this->size());

  // and get the data from the remote processors
  libmesh_call_mpi
    (MPI_Gatherv (r_src.empty() ? nullptr : r_src.data(), mysize,
                  StandardType<T>(), r.empty() ? nullptr : r.data(),
                  sendlengths.data(), displacements.data(),
                  StandardType<T>(), root_id, this->get()));
}



template <typename T, typename A>
inline void Communicator::gather(const unsigned int root_id,
                                 const std::basic_string<T> & sendval,
                                 std::vector<std::basic_string<T>,A> & recv,
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

      libmesh_assert_less(root_id, this->size());

      // and get the data from the remote processors.
      libmesh_call_mpi
        (MPI_Gatherv (const_cast<T*>(sendval.data()),
                      mysize, StandardType<T>(),
                      this->rank() == root_id ? &r[0] : nullptr,
                      sendlengths.data(), displacements.data(),
                      StandardType<T>(), root_id, this->get()));

      // slice receive buffer up
      if (this->rank() == root_id)
        for (unsigned int i=0; i != this->size(); ++i)
          recv[i] = r.substr(displacements[i], sendlengths[i]);
    }
  else
    recv[0] = sendval;
}



template <typename T, typename A>
inline void Communicator::allgather(const T & sendval,
                                    std::vector<T,A> & recv) const
{
  LOG_SCOPE ("allgather()","Parallel");

  libmesh_assert(this->size());
  recv.resize(this->size());

  unsigned int comm_size = this->size();
  if (comm_size > 1)
    {
      StandardType<T> send_type(&sendval);

      libmesh_call_mpi
        (MPI_Allgather (const_cast<T*>(&sendval), 1, send_type, recv.data(), 1,
                        send_type, this->get()));
    }
  else if (comm_size > 0)
    recv[0] = sendval;
}



template <typename T, typename A>
inline void Communicator::allgather(std::vector<T,A> & r,
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

      std::vector<T,A> r_src(r.size()*this->size());
      r_src.swap(r);
      StandardType<T> send_type(r_src.data());

      libmesh_call_mpi
        (MPI_Allgather (r_src.data(), cast_int<int>(r_src.size()),
                        send_type, r.data(), cast_int<int>(r_src.size()),
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
  std::vector<T,A> r_src(globalsize);
  r_src.swap(r);

  StandardType<T> send_type(r.data());

  // and get the data from the remote processors.
  // Pass nullptr if our vector is empty.
  libmesh_call_mpi
    (MPI_Allgatherv (r_src.empty() ? nullptr : r_src.data(), mysize,
                     send_type, r.data(), sendlengths.data(),
                     displacements.data(), send_type, this->get()));
}



template <typename T, typename A>
inline void Communicator::allgather(std::vector<std::basic_string<T>,A> & r,
                                    const bool identical_buffer_sizes) const
{
  if (this->size() < 2)
    return;

  LOG_SCOPE("allgather()", "Parallel");

  if (identical_buffer_sizes)
    {
      libmesh_assert(this->verify(r.size()));

      // identical_buffer_sizes doesn't buy us much since we have to
      // communicate the lengths of strings within each buffer anyway
      if (r.empty())
        return;
    }

  // Concatenate the input buffer into a send buffer, and keep track
  // of input string lengths
  std::vector<int> mystrlengths (r.size());
  std::vector<T> concat_src;

  int myconcatsize = 0;
  for (unsigned int i=0; i != r.size(); ++i)
    {
      int stringlen = cast_int<int>(r[i].size());
      mystrlengths[i] = stringlen;
      myconcatsize += stringlen;
    }
  concat_src.reserve(myconcatsize);
  for (unsigned int i=0; i != r.size(); ++i)
    concat_src.insert
      (concat_src.end(), r[i].begin(), r[i].end());

  // Get the string lengths from all other processors
  std::vector<int> strlengths = mystrlengths;
  this->allgather(strlengths, identical_buffer_sizes);

  // We now know how many strings we'll be receiving
  r.resize(strlengths.size());

  // Get the concatenated data sizes from all other processors
  std::vector<int> concat_sizes;
  this->allgather(myconcatsize, concat_sizes);

  // Find the total size of the final concatenated array and
  // set up the displacement offsets for each processor.
  std::vector<int> displacements(this->size(), 0);
  unsigned int globalsize = 0;
  for (unsigned int i=0; i != this->size(); ++i)
    {
      displacements[i] = globalsize;
      globalsize += concat_sizes[i];
    }

  // Check for quick return
  if (globalsize == 0)
    return;

  // Get the concatenated data from the remote processors.
  // Pass nullptr if our vector is empty.
  std::vector<T> concat(globalsize);

  // We may have concat_src.empty(), but we know concat has at least
  // one element we can use as an example for StandardType
  StandardType<T> send_type(concat.data());

  libmesh_call_mpi
    (MPI_Allgatherv (concat_src.empty() ?
                     nullptr : concat_src.data(), myconcatsize,
                     send_type, concat.data(), concat_sizes.data(),
                     displacements.data(), send_type, this->get()));

  // Finally, split concatenated data into strings
  const T * begin = concat.data();
  for (unsigned int i=0; i != r.size(); ++i)
    {
      const T * end = begin + strlengths[i];
      r[i].assign(begin, end);
      begin = end;
    }
}



template <typename T, typename A>
void Communicator::scatter(const std::vector<T,A> & data,
                           T & recv,
                           const unsigned int root_id) const
{
  libmesh_ignore(root_id); // Only needed for MPI and/or dbg/devel
  libmesh_assert_less (root_id, this->size());

  // Do not allow the root_id to scatter a nullptr vector.
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

  T * data_ptr = const_cast<T*>(data.empty() ? nullptr : data.data());
  libmesh_ignore(data_ptr); // unused ifndef LIBMESH_HAVE_MPI

  libmesh_assert_less(root_id, this->size());

  libmesh_call_mpi
    (MPI_Scatter (data_ptr, 1, StandardType<T>(data_ptr),
                  &recv, 1, StandardType<T>(&recv), root_id, this->get()));
}



template <typename T, typename A>
void Communicator::scatter(const std::vector<T,A> & data,
                           std::vector<T,A> & recv,
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

  int recv_buffer_size = 0;
  if (this->rank() == root_id)
    {
      libmesh_assert(data.size() % this->size() == 0);
      recv_buffer_size = cast_int<int>(data.size() / this->size());
    }

  this->broadcast(recv_buffer_size);
  recv.resize(recv_buffer_size);

  T * data_ptr = const_cast<T*>(data.empty() ? nullptr : data.data());
  T * recv_ptr = recv.empty() ? nullptr : recv.data();
  libmesh_ignore(data_ptr, recv_ptr); // unused ifndef LIBMESH_HAVE_MPI

  libmesh_assert_less(root_id, this->size());

  libmesh_call_mpi
    (MPI_Scatter (data_ptr, recv_buffer_size, StandardType<T>(data_ptr),
                  recv_ptr, recv_buffer_size, StandardType<T>(recv_ptr), root_id, this->get()));
}



template <typename T, typename A1, typename A2>
void Communicator::scatter(const std::vector<T,A1> & data,
                           const std::vector<int,A2> counts,
                           std::vector<T,A1> & recv,
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

  std::vector<int,A2> displacements(this->size(), 0);
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
  int recv_buffer_size = 0;
  this->scatter(counts, recv_buffer_size, root_id);
  recv.resize(recv_buffer_size);

  T * data_ptr = const_cast<T*>(data.empty() ? nullptr : data.data());
  int * count_ptr = const_cast<int*>(counts.empty() ? nullptr : counts.data());
  T * recv_ptr = recv.empty() ? nullptr : recv.data();
  libmesh_ignore(data_ptr, count_ptr, recv_ptr); // unused ifndef LIBMESH_HAVE_MPI

  libmesh_assert_less(root_id, this->size());

  // Scatter the non-uniform chunks
  libmesh_call_mpi
    (MPI_Scatterv (data_ptr, count_ptr, displacements.data(), StandardType<T>(data_ptr),
                   recv_ptr, recv_buffer_size, StandardType<T>(recv_ptr), root_id, this->get()));
}



template <typename T, typename A1, typename A2>
void Communicator::scatter(const std::vector<std::vector<T,A1>,A2> & data,
                           std::vector<T,A1> & recv,
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

  std::vector<T,A1> stacked_data;
  std::vector<int> counts;

  if (root_id == this->rank())
    {
      libmesh_assert (data.size() == this->size());

      if (!identical_buffer_sizes)
        counts.resize(this->size());

      for (std::size_t i=0; i < data.size(); ++i)
        {
          if (!identical_buffer_sizes)
            counts[i] = cast_int<int>(data[i].size());
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



template <typename T, typename A>
inline void Communicator::alltoall(std::vector<T,A> & buf) const
{
  if (this->size() < 2 || buf.empty())
    return;

  LOG_SCOPE("alltoall()", "Parallel");

  // the per-processor size.  this is the same for all
  // processors using MPI_Alltoall, could be variable
  // using MPI_Alltoallv
  const int size_per_proc =
    cast_int<int>(buf.size()/this->size());
  libmesh_ignore(size_per_proc);

  libmesh_assert_equal_to (buf.size()%this->size(), 0);

  libmesh_assert(this->verify(size_per_proc));

  StandardType<T> send_type(buf.data());

  libmesh_call_mpi
    (MPI_Alltoall (MPI_IN_PLACE, size_per_proc, send_type, buf.data(),
                   size_per_proc, send_type, this->get()));
}



template <typename T>
inline void Communicator::broadcast (T & libmesh_mpi_var(data), const unsigned int root_id) const
{
  libmesh_ignore(root_id); // Only needed for MPI and/or dbg/devel
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
        (buffer, context, out_iter, (T*)(nullptr));

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
        (buffer, context, out_iter, (T*)nullptr);

      nonempty_range = (range_begin != range_end);
      this->max(nonempty_range);
    }
}



template <typename T, typename A>
inline bool Communicator::possibly_receive (unsigned int & src_processor_id,
                                            std::vector<T,A> & buf,
                                            const DataType & type,
                                            Request & req,
                                            const MessageTag & tag) const
{
  LOG_SCOPE("possibly_receive()", "Parallel");

  Status stat(type);

  int int_flag = 0;

  libmesh_assert(src_processor_id < this->size() ||
                 src_processor_id == any_source);

  libmesh_call_mpi(MPI_Iprobe(src_processor_id,
                              tag.value(),
                              this->get(),
                              &int_flag,
                              stat.get()));

  if (int_flag)
  {
    buf.resize(stat.size());

    src_processor_id = stat.source();

    libmesh_call_mpi
      (MPI_Irecv (buf.data(),
                  cast_int<int>(buf.size()),
                  type,
                  src_processor_id,
                  tag.value(),
                  this->get(),
                  req.get()));

    // The MessageTag should stay registered for the Request lifetime
    req.add_post_wait_work
      (new Parallel::PostWaitDereferenceTag(tag));
  }

  return int_flag;
}


} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_PARALLEL_IMPLEMENTATION_H
