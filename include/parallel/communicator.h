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


#ifndef LIBMESH_COMMUNICATOR_H
#define LIBMESH_COMMUNICATOR_H

// Parallel includes
#include "libmesh/data_type.h"
#include "libmesh/message_tag.h"
#include "libmesh/request.h"
#include "libmesh/status.h"
#include "libmesh/post_wait_dereference_shared_ptr.h"

// libMesh Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <map>
#include <string>
#include <vector>

// C++ includes needed for parallel_communicator_specializations
//
// These could be forward declarations if only that wasn't illegal
#include <numeric> // complex
#include <set>

namespace libMesh
{

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 *
 * For MPI 1.1 compatibility, temporary buffers are used
 * instead of MPI 2's MPI_IN_PLACE
 */
namespace Parallel
{

template <typename T>
class Packing;

#ifdef LIBMESH_HAVE_MPI

//-------------------------------------------------------------------
/**
 * Communicator object for talking with subsets of processors
 */
typedef MPI_Comm communicator;

/**
 * Info object used by some MPI-3 methods
 */
typedef MPI_Info info;

/**
 * Processor id meaning "Accept from any source"
 */
const unsigned int any_source =
  static_cast<unsigned int>(MPI_ANY_SOURCE);

#else

// These shouldn't actually be needed, but must be
// unique types for function overloading to work
// properly.
typedef int communicator; // Must match petsc-nompi definition

typedef int info;

const unsigned int any_source=0;

#endif // LIBMESH_HAVE_MPI

//-------------------------------------------------------------------
/**
 * Encapsulates the MPI_Comm object.  Allows the size of the group
 * and this process's position in the group to be determined.
 *
 * Methods of this object are the preferred way to perform
 * distributed-memory parallel operations.
 */
class Communicator
{
  // Basic operations:
public:

  /**
   * Default Constructor.
   */
  Communicator ();

  /*
   * Constructor from MPI_Comm
   */
  explicit Communicator (const communicator & comm);

  /*
   * Don't use copy construction or assignment, just copy by reference
   * or pointer - it's too hard to keep a common used_tag_values if
   * each communicator is shared by more than one Communicator
   */
  Communicator (const Communicator &) = delete;
  Communicator & operator= (const Communicator &) = delete;

  /*
   * Move constructor and assignment operator
   */
  Communicator (Communicator &&) = default;
  Communicator & operator= (Communicator &&) = default;

  /*
   * NON-VIRTUAL destructor
   */
  ~Communicator ();

  /*
   * Create a new communicator between some subset of \p this,
   * based on specified "color"
   */
  void split(int color, int key, Communicator & target) const;

  /*
   * Create a new communicator between some subset of \p this,
   * based on specified split "type" (e.g. MPI_COMM_TYPE_SHARED)
   */
  void split_by_type(int split_type, int key, info i, Communicator & target) const;

  /*
   * Create a new duplicate of \p this communicator
   */
  void duplicate(const Communicator & comm);

  /*
   * Create a new duplicate of an MPI communicator
   */
  void duplicate(const communicator & comm);

  communicator & get() { return _communicator; }

  const communicator & get() const { return _communicator; }

  /**
   * Get a tag that is unique to this Communicator.  A requested tag
   * value may be provided.  If no request is made then an automatic
   * unique tag value will be generated; such usage of
   * get_unique_tag() must be done on every processor in a consistent
   * order.
   *
   * \note If people are also using magic numbers or copying
   * raw communicators around then we can't guarantee the tag is
   * unique to this MPI_Comm.
   *
   * \note Leaving \p tagvalue unspecified is recommended in most
   * cases.  Manually selecting tag values is dangerous, as tag values may be
   * freed and reselected earlier than expected in asynchronous
   * communication algorithms.
   */
  MessageTag get_unique_tag(int tagvalue = MessageTag::invalid_tag) const;

  /**
   * Reference an already-acquired tag, so that we know it will
   * be dereferenced multiple times before we can re-release it.
   */
  void reference_unique_tag(int tagvalue) const;

  /**
   * Dereference an already-acquired tag, and see if we can
   * re-release it.
   */
  void dereference_unique_tag(int tagvalue) const;

  /**
   * Free and reset this communicator
   */
  void clear();

  Communicator & operator= (const communicator & comm);

  processor_id_type rank() const { return _rank; }

  processor_id_type size() const { return _size; }

  /**
   * Whether to use default or synchronous sends?
   */
  enum SendMode { DEFAULT=0, SYNCHRONOUS };

private:

  /**
   * Utility function for setting our member variables from an MPI
   * communicator
   */
  void assign(const communicator & comm);

  communicator  _communicator;
  processor_id_type _rank, _size;
  SendMode _send_mode;

  // mutable used_tag_values and tag_queue - not thread-safe, but then
  // Parallel:: isn't thread-safe in general.
  mutable std::map<int, unsigned int> used_tag_values;
  mutable int _next_tag;

  int _max_tag;

  // Keep track of duplicate/split operations so we know when to free
  bool _I_duped_it;

  // Communication operations:
public:

  /**
   * Explicitly sets the \p SendMode type used for send operations.
   */
  void send_mode (const SendMode sm) { _send_mode = sm; }

  /**
   * Gets the user-requested SendMode.
   */
  SendMode send_mode() const { return _send_mode; }

  /**
   * Pause execution until all processors reach a certain point.
   */
  void barrier () const;

  /**
   * Start a barrier that doesn't block
   */
  void nonblocking_barrier (Request & req) const;

  /**
   * Verify that a local variable has the same value on all processors.
   * Containers must have the same value in every entry.
   */
  template <typename T>
  bool verify(const T & r) const;

  /**
   * Verify that a local pointer points to the same value on all
   * processors where it is not nullptr.
   * Containers must have the same value in every entry.
   */
  template <typename T>
  bool semiverify(const T * r) const;

  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors.  Containers are replaced element-wise.
   */
  template <typename T>
  void min(T & r) const;

  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors, returning the minimum rank of a processor
   * which originally held the minimum value.
   */
  template <typename T>
  void minloc(T & r,
              unsigned int & min_id) const;

  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors.  Set each \p min_id entry to
   * the minimum rank where a corresponding minimum was found.
   */
  template <typename T, typename A1, typename A2>
  void minloc(std::vector<T,A1> & r,
              std::vector<unsigned int,A2> & min_id) const;

  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors.  Containers are replaced element-wise.
   */
  template <typename T>
  void max(T & r) const;

  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors, returning the minimum rank of a processor
   * which originally held the maximum value.
   */
  template <typename T>
  void maxloc(T & r,
              unsigned int & max_id) const;

  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors.  Set each \p min_id entry to
   * the minimum rank where a corresponding maximum was found.
   */
  template <typename T, typename A1, typename A2>
  void maxloc(std::vector<T,A1> & r,
              std::vector<unsigned int,A2> & max_id) const;

  /**
   * Take a local variable and replace it with the sum of it's values
   * on all processors.  Containers are replaced element-wise.
   */
  template <typename T>
  void sum(T & r) const;

  /**
   * Take a container of local variables on each processor, and
   * collect their union over all processors, replacing the set on
   * processor 0.
   */
  template <typename T>
  void set_union(T & data, const unsigned int root_id) const;

  /**
   * Take a container of local variables on each processor, and
   * replace it with their union over all processors.
   */
  template <typename T>
  void set_union(T & data) const;

  /**
   * Blocking message probe.  Allows information about a message to be
   * examined before the message is actually received.
   */
  status probe (const unsigned int src_processor_id,
                const MessageTag & tag=any_tag) const;

  /**
   * Non-Blocking message probe for a packed range message.
   * Allows information about a message to be
   * examined before the message is actually received.
   *
   * Template type must match the object type that will be in
   * the packed range
   *
   * \param src_processor_id The processor the message is expected from or Parallel::any_source
   * \param tag The message tag or Parallel::any_tag
   * \param flag Output.  True if a message exists.  False otherwise.
   */
  template <typename T>
  Status packed_range_probe (const unsigned int src_processor_id,
                             const MessageTag & tag,
                             bool & flag) const;

  /**
   * Blocking-send to one processor with data-defined type.
   */
  template <typename T>
  void send (const unsigned int dest_processor_id,
             const T & buf,
             const MessageTag & tag=no_tag) const;

  /**
   * Nonblocking-send to one processor with data-defined type.
   */
  template <typename T>
  void send (const unsigned int dest_processor_id,
             const T & buf,
             Request & req,
             const MessageTag & tag=no_tag) const;

  /**
   * Blocking-send to one processor with user-defined type.
   *
   * If \p T is a container, container-of-containers, etc., then
   * \p type should be the DataType of the underlying fixed-size
   * entries in the container(s).
   */
  template <typename T>
  void send (const unsigned int dest_processor_id,
             const T & buf,
             const DataType & type,
             const MessageTag & tag=no_tag) const;

  /**
   * Nonblocking-send to one processor with user-defined type.
   *
   * If \p T is a container, container-of-containers, etc., then
   * \p type should be the DataType of the underlying fixed-size
   * entries in the container(s).
   */
  template <typename T>
  void send (const unsigned int dest_processor_id,
             const T & buf,
             const DataType & type,
             Request & req,
             const MessageTag & tag=no_tag) const;

  /**
   * Blocking-receive from one processor with data-defined type.
   */
  template <typename T>
  Status receive (const unsigned int dest_processor_id,
                  T & buf,
                  const MessageTag & tag=any_tag) const;

  /**
   * Nonblocking-receive from one processor with data-defined type.
   */
  template <typename T>
  void receive (const unsigned int dest_processor_id,
                T & buf,
                Request & req,
                const MessageTag & tag=any_tag) const;

  /**
   * Blocking-receive from one processor with user-defined type.
   *
   * If \p T is a container, container-of-containers, etc., then
   * \p type should be the DataType of the underlying fixed-size
   * entries in the container(s).
   */
  template <typename T>
  Status receive (const unsigned int dest_processor_id,
                  T & buf,
                  const DataType & type,
                  const MessageTag & tag=any_tag) const;

  /**
   * Nonblocking-receive from one processor with user-defined type.
   *
   * If \p T is a container, container-of-containers, etc., then
   * \p type should be the DataType of the underlying fixed-size
   * entries in the container(s).
   */
  template <typename T>
  void receive (const unsigned int dest_processor_id,
                T & buf,
                const DataType & type,
                Request & req,
                const MessageTag & tag=any_tag) const;

  /**
   * Nonblocking-receive from one processor with user-defined type.
   *
   * Checks to see if a message can be received from the
   * src_processor_id .  If so, it starts a non-blocking
   * receive using the passed in request and returns true
   *
   * Otherwise - if there is no message to receive it returns false
   *
   * Note: The buf does NOT need to properly sized before this call
   * this will resize the buffer automatically
   *
   * If \p T is a container, container-of-containers, etc., then
   * \p type should be the DataType of the underlying fixed-size
   * entries in the container(s).
   *
   * @param src_processor_id The pid to receive from or "any".
   * will be set to the actual src being receieved from
   * @param buf THe buffer to receive into
   * @param type The intrinsic datatype to receive
   * @param req The request to use
   * @param tag The tag to use
   */
  template <typename T, typename A>
  bool possibly_receive (unsigned int & src_processor_id,
                         std::vector<T,A> & buf,
                         const DataType & type,
                         Request & req,
                         const MessageTag & tag) const;

  /**
   * Blocking-send range-of-pointers to one processor.  This
   * function does not send the raw pointers, but rather constructs
   * new objects at the other end whose contents match the objects
   * pointed to by the sender.
   *
   * void Parallel::pack(const T *, vector<int> & data, const Context *)
   * is used to serialize type T onto the end of a data vector.
   *
   * unsigned int Parallel::packable_size(const T *, const Context *) is
   * used to allow data vectors to reserve memory, and for additional
   * error checking
   */
  template <typename Context, typename Iter>
  void send_packed_range (const unsigned int dest_processor_id,
                          const Context * context,
                          Iter range_begin,
                          const Iter range_end,
                          const MessageTag & tag=no_tag) const;

  /**
   * Nonblocking-send range-of-pointers to one processor.  This
   * function does not send the raw pointers, but rather constructs
   * new objects at the other end whose contents match the objects
   * pointed to by the sender.
   *
   * void Parallel::pack(const T *, vector<int> & data, const Context *)
   * is used to serialize type T onto the end of a data vector.
   *
   * unsigned int Parallel::packable_size(const T *, const Context *) is
   * used to allow data vectors to reserve memory, and for additional
   * error checking
   */
  template <typename Context, typename Iter>
  void send_packed_range (const unsigned int dest_processor_id,
                          const Context * context,
                          Iter range_begin,
                          const Iter range_end,
                          Request & req,
                          const MessageTag & tag=no_tag) const;

  /**
   * Similar to the above Nonblocking send_packed_range with a few important differences:
   *
   * 1. The total size of the packed buffer MUST be less than std::numeric_limits<int>::max()
   * 2. Only _one_ message is generated
   * 3. On the receiving end the message should be tested for using Communicator::packed_range_probe()
   * 4. The message must be received by Communicator::nonblocking_receive_packed_range()
   */
  template <typename Context, typename Iter>
  void nonblocking_send_packed_range (const unsigned int dest_processor_id,
                                      const Context * context,
                                      Iter range_begin,
                                      const Iter range_end,
                                      Request & req,
                                      const MessageTag & tag=no_tag) const;


  /**
   * Similar to the above Nonblocking send_packed_range with a few important differences:
   *
   * 1. The total size of the packed buffer MUST be less than std::numeric_limits<int>::max()
   * 2. Only _one_ message is generated
   * 3. On the receiving end the message should be tested for using Communicator::packed_range_probe()
   * 4. The message must be received by Communicator::nonblocking_receive_packed_range()
   */
  template <typename Context, typename Iter>
  void nonblocking_send_packed_range (const unsigned int dest_processor_id,
                                      const Context * context,
                                      Iter range_begin,
                                      const Iter range_end,
                                      Request & req,
                                      std::shared_ptr<std::vector<typename Parallel::Packing<typename std::iterator_traits<Iter>::value_type>::buffer_type>> & buffer,
                                      const MessageTag & tag=no_tag) const;

  /**
   * Blocking-receive range-of-pointers from one processor.  This
   * function does not receive raw pointers, but rather constructs new
   * objects whose contents match the objects pointed to by the
   * sender.
   *
   * The objects will be of type
   * T = iterator_traits<OutputIter>::value_type.
   *
   * Using std::back_inserter as the output iterator allows receive to
   * fill any container type.  Using libMesh::null_output_iterator
   * allows the receive to be dealt with solely by Parallel::unpack(),
   * for objects whose unpack() is written so as to not leak memory
   * when used in this fashion.
   *
   * A future version of this method should be created to preallocate
   * memory when receiving vectors...
   *
   * void Parallel::unpack(vector<int>::iterator in, T ** out, Context *)
   * is used to unserialize type T, typically into a new
   * heap-allocated object whose pointer is returned as *out.
   *
   * unsigned int Parallel::packed_size(const T *,
   *                                    vector<int>::const_iterator)
   * is used to advance to the beginning of the next object's data.
   */
  template <typename Context, typename OutputIter, typename T>
  void receive_packed_range (const unsigned int dest_processor_id,
                             Context * context,
                             OutputIter out,
                             const T * output_type, // used only to infer T
                             const MessageTag & tag=any_tag) const;

  /**
   * Non-Blocking-receive range-of-pointers from one processor.
   *
   * This is meant to receive messages from nonblocking_send_packed_range
   *
   * Similar in design to the above receive_packed_range.  However,
   * this version requires a Request and a Status.
   *
   * The Status must be a positively tested Status for a message of this
   * type (i.e. a message _does_ exist).  It should most likely be generated by
   * Communicator::packed_range_probe.
   */
  template <typename Context, typename OutputIter, typename T>
  void nonblocking_receive_packed_range (const unsigned int src_processor_id,
                                         Context * context,
                                         OutputIter out,
                                         const T * output_type,
                                         Request & req,
                                         Status & stat,
                                         const MessageTag & tag=any_tag) const;

  /**
   * Non-Blocking-receive range-of-pointers from one processor.
   *
   * This is meant to receive messages from nonblocking_send_packed_range
   *
   * Similar in design to the above receive_packed_range.  However,
   * this version requires a Request and a Status.
   *
   * The Status must be a positively tested Status for a message of this
   * type (i.e. a message _does_ exist).  It should most likely be generated by
   * Communicator::packed_range_probe.
   */
  template <typename Context, typename OutputIter, typename T>
  void nonblocking_receive_packed_range (const unsigned int src_processor_id,
                                         Context * context,
                                         OutputIter out,
                                         const T * output_type,
                                         Request & req,
                                         Status & stat,
                                         std::shared_ptr<std::vector<typename Parallel::Packing<T>::buffer_type>> & buffer,
                                         const MessageTag & tag=any_tag
                                         ) const;

  /**
   * Send data \p send to one processor while simultaneously receiving
   * other data \p recv from a (potentially different) processor.
   */
  template <typename T1, typename T2>
  void send_receive(const unsigned int dest_processor_id,
                    const T1 & send,
                    const unsigned int source_processor_id,
                    T2 & recv,
                    const MessageTag & send_tag = no_tag,
                    const MessageTag & recv_tag = any_tag) const;

  /**
   * Send a range-of-pointers to one processor while simultaneously receiving
   * another range from a (potentially different) processor.  This
   * function does not send or receive raw pointers, but rather constructs
   * new objects at each receiver whose contents match the objects
   * pointed to by the sender.
   *
   * The objects being sent will be of type
   * T1 = iterator_traits<RangeIter>::value_type, and the objects
   * being received will be of type
   * T2 = iterator_traits<OutputIter>::value_type
   *
   * void Parallel::pack(const T1*, vector<int> & data, const Context1*)
   * is used to serialize type T1 onto the end of a data vector.
   *
   * Using std::back_inserter as the output iterator allows
   * send_receive to fill any container type.  Using
   * libMesh::null_output_iterator allows the receive to be dealt with
   * solely by Parallel::unpack(), for objects whose unpack() is
   * written so as to not leak memory when used in this fashion.
   *
   * A future version of this method should be created to preallocate
   * memory when receiving vectors...
   *
   * void Parallel::unpack(vector<int>::iterator in, T2** out, Context *)
   * is used to unserialize type T2, typically into a new
   * heap-allocated object whose pointer is returned as *out.
   *
   * unsigned int Parallel::packable_size(const T1*, const Context1*)
   * is used to allow data vectors to reserve memory, and for
   * additional error checking.
   *
   * unsigned int Parallel::packed_size(const T2*,
   *                                    vector<int>::const_iterator)
   * is used to advance to the beginning of the next object's data.
   */
  template <typename Context1, typename RangeIter, typename Context2,
            typename OutputIter, typename T>
  void send_receive_packed_range(const unsigned int dest_processor_id,
                                 const Context1 * context1,
                                 RangeIter send_begin,
                                 const RangeIter send_end,
                                 const unsigned int source_processor_id,
                                 Context2 * context2,
                                 OutputIter out,
                                 const T * output_type, // used only to infer T
                                 const MessageTag & send_tag = no_tag,
                                 const MessageTag & recv_tag = any_tag) const;

  /**
   * Send data \p send to one processor while simultaneously receiving
   * other data \p recv from a (potentially different) processor, using
   * a user-specified MPI Dataype.
   */
  template <typename T1, typename T2>
  void send_receive(const unsigned int dest_processor_id,
                    const T1 & send,
                    const DataType & type1,
                    const unsigned int source_processor_id,
                    T2 & recv,
                    const DataType & type2,
                    const MessageTag & send_tag = no_tag,
                    const MessageTag & recv_tag = any_tag) const;

  /**
   * Take a vector of length comm.size(), and on processor root_id fill in
   * recv[processor_id] = the value of send on processor processor_id
   */
  template <typename T, typename A>
  inline void gather(const unsigned int root_id,
                     const T & send,
                     std::vector<T,A> & recv) const;

  /**
   * The gather overload for string types has an optional
   * identical_buffer_sizes optimization for when all strings are the
   * same length.
   */
  template <typename T, typename A>
  inline void gather(const unsigned int root_id,
                     const std::basic_string<T> & send,
                     std::vector<std::basic_string<T>,A> & recv,
                     const bool identical_buffer_sizes=false) const;

  /**
   * Take a vector of local variables and expand it on processor root_id
   * to include values from all processors
   *
   * This handles the case where the lengths of the vectors may vary.
   * Specifically, this function transforms this:
   * \verbatim
   * Processor 0: [ ... N_0 ]
   * Processor 1: [ ....... N_1 ]
   * ...
   * Processor M: [ .. N_M]
   * \endverbatim
   *
   * into this:
   *
   * \verbatim
   * [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
   * \endverbatim
   *
   * on processor root_id. This function is collective and therefore
   * must be called by all processors in the Communicator.
   */
  template <typename T, typename A>
  inline void gather(const unsigned int root_id,
                     std::vector<T,A> & r) const;

  /**
   * Take a vector of length \p this->size(), and fill in
   * \p recv[processor_id] = the value of \p send on that processor
   */
  template <typename T, typename A>
  inline void allgather(const T & send,
                        std::vector<T,A> & recv) const;

  /**
   * The allgather overload for string types has an optional
   * identical_buffer_sizes optimization for when all strings are the
   * same length.
   */
  template <typename T, typename A>
  inline void allgather(const std::basic_string<T> & send,
                        std::vector<std::basic_string<T>,A> & recv,
                        const bool identical_buffer_sizes=false) const;

  /**
   * Take a vector of local variables and expand it to include
   * values from all processors. By default, each processor is
   * allowed to have its own unique input buffer length. If
   * it is known that all processors have the same input sizes
   * additional communication can be avoided.
   *
   * Specifically, this function transforms this:
   * \verbatim
   * Processor 0: [ ... N_0 ]
   * Processor 1: [ ....... N_1 ]
   * ...
   * Processor M: [ .. N_M]
   * \endverbatim
   *
   * into this:
   *
   * \verbatim
   * [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
   * \endverbatim
   *
   * on each processor. This function is collective and therefore
   * must be called by all processors in the Communicator.
   */
  template <typename T, typename A>
  inline void allgather(std::vector<T,A> & r,
                        const bool identical_buffer_sizes = false) const;

  /**
   * AllGather overload for vectors of string types
   */
  template <typename T, typename A>
  inline void allgather(std::vector<std::basic_string<T>,A> & r,
                        const bool identical_buffer_sizes = false) const;

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and scatter the ith item to the ith
   * processor in the communicator. The result is saved into recv.
   */
  template <typename T, typename A>
  inline void scatter(const std::vector<T,A> & data,
                      T & recv,
                      const unsigned int root_id=0) const;

  /**
   * Take a vector of local variables and scatter the ith equal-sized chunk
   * to the ith processor in the communicator. The data size must be a
   * multiple of the communicator size. The result is saved into recv buffer.
   * The recv buffer does not have to be sized prior to this operation.
   */
  template <typename T, typename A>
  inline void scatter(const std::vector<T,A> & data,
                      std::vector<T,A> & recv,
                      const unsigned int root_id=0) const;

  /**
   * Take a vector of local variables and scatter the ith variable-sized chunk
   * to the ith processor in the communicator. The counts vector should contain
   * the number of items for each processor. The result is saved into recv buffer.
   * The recv buffer does not have to be sized prior to this operation.
   */
  template <typename T, typename A1, typename A2>
  inline void scatter(const std::vector<T,A1> & data,
                      const std::vector<int,A2> counts,
                      std::vector<T,A1> & recv,
                      const unsigned int root_id=0) const;

  /**
   * Take a vector of vectors and scatter the ith inner vector
   * to the ith processor in the communicator. The result is saved into recv buffer.
   * The recv buffer does not have to be sized prior to this operation.
   */
  template <typename T, typename A1, typename A2>
  inline void scatter(const std::vector<std::vector<T,A1>,A2> & data,
                      std::vector<T,A1> & recv,
                      const unsigned int root_id=0,
                      const bool identical_buffer_sizes=false) const;

  //-------------------------------------------------------------------
  /**
   * Take a range of local variables, combine it with ranges from all
   * processors, and write the output to the output iterator on rank root.
   */
  template <typename Context, typename Iter, typename OutputIter>
  inline void gather_packed_range (const unsigned int root_id,
                                   Context * context,
                                   Iter range_begin,
                                   const Iter range_end,
                                   OutputIter out) const;

  /**
   * Take a range of local variables, combine it with ranges from all
   * processors, and write the output to the output iterator.
   */
  template <typename Context, typename Iter, typename OutputIter>
  inline void allgather_packed_range (Context * context,
                                      Iter range_begin,
                                      const Iter range_end,
                                      OutputIter out) const;

  /**
   * Effectively transposes the input vector across all processors.
   * The jth entry on processor i is replaced with the ith entry
   * from processor j.
   */
  template <typename T, typename A>
  inline void alltoall(std::vector<T,A> & r) const;

  /**
   * Take a local value and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor initiating the broadcast.
   * If \p data is a vector, the user is responsible for resizing it
   * on all processors, except in the case when \p data is a vector
   * of strings.
   */
  template <typename T>
  inline void broadcast(T & data, const unsigned int root_id=0) const;

  /**
   * Blocking-broadcast range-of-pointers to one processor.  This
   * function does not send the raw pointers, but rather constructs
   * new objects at the other end whose contents match the objects
   * pointed to by the sender.
   *
   * void Parallel::pack(const T *, vector<int> & data, const Context *)
   * is used to serialize type T onto the end of a data vector.
   *
   * unsigned int Parallel::packable_size(const T *, const Context *) is
   * used to allow data vectors to reserve memory, and for additional
   * error checking
   *
   * unsigned int Parallel::packed_size(const T *,
   *                                    vector<int>::const_iterator)
   * is used to advance to the beginning of the next object's data.
   */
  template <typename Context, typename OutputContext, typename Iter, typename OutputIter>
  inline void broadcast_packed_range (const Context * context1,
                                      Iter range_begin,
                                      const Iter range_end,
                                      OutputContext * context2,
                                      OutputIter out,
                                      const unsigned int root_id = 0) const;

  /**
   * C++ doesn't let us partially specialize functions (we're really
   * just doing operator overloading on them), so more-specialized
   * versions of the functions above need to be redeclared.  We'll
   * do so in a separate file so that users don't have to look at
   * the redundancy.
   */
#include "libmesh/parallel_communicator_specializations"

}; // class Communicator


} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_COMMUNICATOR_H
