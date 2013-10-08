// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_PARALLEL_H
#define LIBMESH_PARALLEL_H

// Local includes
#include "libmesh/libmesh_common.h" // libmesh_assert, libmesh_cast_int
#include "libmesh/libmesh_logging.h"

// C++ includes
#include <cstddef>
#include <climits>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace libMesh
{


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only
#ifndef NDEBUG
  #define parallel_only() do { \
    libmesh_deprecated(); \
    libmesh_assert(CommWorld.verify(std::string(__FILE__).size())); \
    libmesh_assert(CommWorld.verify(std::string(__FILE__))); \
    libmesh_assert(CommWorld.verify(__LINE__)); } while (0)
#else
  #define parallel_only()  ((void) 0)
#endif

#undef libmesh_parallel_only
#ifndef NDEBUG
  #define libmesh_parallel_only(comm_obj) do { \
    libmesh_assert(comm_obj.verify(std::string(__FILE__).size())); \
    libmesh_assert(comm_obj.verify(std::string(__FILE__))); \
    libmesh_assert(comm_obj.verify(__LINE__)); } while (0)
#else
  #define libmesh_parallel_only(comm_obj)  ((void) 0)
#endif

// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only_on
#ifndef NDEBUG
  #define parallel_only_on(comm_arg) do { \
    libmesh_deprecated(); \
    libmesh_assert(CommWorld.verify(std::string(__FILE__).size(), comm_arg)); \
    libmesh_assert(CommWorld.verify(std::string(__FILE__), comm_arg)); \
    libmesh_assert(CommWorld.verify(__LINE__), comm_arg); } while (0)
#else
  #define parallel_only_on(comm_arg)  ((void) 0)
#endif

#undef libmesh_parallel_only_on
#ifndef NDEBUG
  #define libmesh_parallel_only_on(comm_obj,comm_arg) do {	\
    libmesh_assert(comm_obj.verify(std::string(__FILE__).size(), comm_arg)); \
    libmesh_assert(comm_obj.verify(std::string(__FILE__), comm_arg)); \
    libmesh_assert(comm_obj.verify(__LINE__), comm_arg); } while (0)
#else
  #define libmesh_parallel_only_on(comm_obj,comm_arg)  ((void) 0)
#endif

/**
 * The Parallel namespace is for wrapper functions
 * for common general parallel synchronization tasks.
 *
 * For MPI 1.1 compatibility, temporary buffers are used
 * instead of MPI 2's MPI_IN_PLACE
 */
namespace Parallel
{
  //-------------------------------------------------------------------
  /**
   * Forward declarations of classes we will define later.
   */
  class Communicator;
  class DataType;
  class Request;
  class Status;

#ifdef LIBMESH_HAVE_MPI
  //-------------------------------------------------------------------
  /**
   * Data types for communication
   */
  typedef MPI_Datatype data_type;

  /**
   * Request object for non-blocking I/O
   */
  typedef MPI_Request request;

  /**
   * Status object for querying messages
   */
  typedef MPI_Status status;

  /**
   * Communicator object for talking with subsets of processors
   */
  typedef MPI_Comm communicator;

  /**
   * Templated function to return the appropriate MPI datatype
   * for use with built-in C types when combined with an int
   */
  template <typename T>
  inline data_type dataplusint_type();

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

  /**
   * Accept from any source
   */
  const unsigned int any_source =
    static_cast<unsigned int>(MPI_ANY_SOURCE);


#else

  // These shouldn't actually be needed, but must be
  // unique types for function overloading to work
  // properly.
  struct data_type    { /* unsigned int t; */ };
  struct request      { /* unsigned int r; */ };
  struct status       { /* unsigned int s; */ };
  typedef int communicator; // Must match petsc-nompi definition

  const unsigned int any_source=0;
#endif // LIBMESH_HAVE_MPI



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI tag integers.
   */
  class MessageTag
  {
  public:

    /**
     * Invalid tag, to allow for default construction.
     */
    static const int invalid_tag = INT_MIN;

    /**
     * Explicit constructor, to discourage using "magic numbers"
     * as tags.  Communicator::get_unique_tag is recommended instead.
     */
    explicit MessageTag(int tagvalue = invalid_tag)
      : _tagvalue(tagvalue), _comm(NULL) {}

    /**
     * Copy constructor.  Helps Communicator do reference counting on
     * unique tags
     */
    MessageTag(const MessageTag& other);

    /**
     * Destructor.  Helps Communicator do reference counting on unique
     * tags
     */
    ~MessageTag();

    int value() const {
      return _tagvalue;
    }

  private:
    int _tagvalue;
    const Communicator *_comm;

    // Constructor for reference-counted unique tags
    MessageTag(int tagvalue, const Communicator *comm)
      : _tagvalue(tagvalue), _comm(comm) {}

    // Let Communicator handle the reference counting
    friend class Communicator;
  };


  //-------------------------------------------------------------------
  /**
   * Default message tag ids
   */
#ifdef LIBMESH_HAVE_MPI
  const MessageTag any_tag = MessageTag(MPI_ANY_TAG);
#else
  const MessageTag any_tag = MessageTag(-1);
#endif

  const MessageTag no_tag = MessageTag(0);


  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Datatype.
   */
  class DataType
  {
  public:
    DataType () : _datatype() {}

    DataType (const DataType &other) :
      _datatype(other._datatype)
    {}

    DataType (const data_type &type) :
      _datatype(type)
    {}

#ifdef LIBMESH_HAVE_MPI
    DataType (const DataType &other, unsigned int count)
    {
      MPI_Type_contiguous(count, other._datatype, &_datatype);
      this->commit();
    }
#else
    DataType (const DataType &, unsigned int)
    {
    }
#endif

    DataType & operator = (const DataType &other)
    { _datatype = other._datatype; return *this; }

    DataType & operator = (const data_type &type)
    { _datatype = type; return *this; }

    operator const data_type & () const
    { return _datatype; }

    operator data_type & ()
    { return _datatype; }

//     operator data_type const * () const
//     { return &_datatype; }

//     operator data_type * ()
//     { return &_datatype; }

    void commit ()
    {
#ifdef LIBMESH_HAVE_MPI
      MPI_Type_commit (&_datatype);
#endif
    }

    void free ()
    {
#ifdef LIBMESH_HAVE_MPI
      MPI_Type_free (&_datatype);
#endif
    }

  protected:

    data_type _datatype;
  };


  //-------------------------------------------------------------------
  /**
   * Templated class to provide the appropriate MPI datatype
   * for use with built-in C types or simple C++ constructions.
   *
   * More complicated data types may need to provide a pointer-to-T so
   * that we can use MPI_Address without constructing a new T.
   */
  template <typename T>
  class StandardType : public DataType
  {
  /*
   * The unspecialized class is useless, so we make its constructor
   * private to catch mistakes at compile-time rather than link-time.
   * Specializations should have a public constructor of the same
   * form.
   */
  private:
    StandardType(const T* example = NULL);
  };

  /*
   * The unspecialized class gives default, lowest-common-denominator
   * attributes, for values which can't be used with Parallel min/max.
   * Specialized classes can set this to true, and should define
   * the lowest and highest values possible for the type.
   */
  template<typename T>
  struct Attributes
  {
    static const bool has_min_max = false;
    static void set_lowest(T&) {}
    static void set_highest(T&) {}
  };



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Status struct.  Allows the source and size
   * of the message to be determined.
   */
  class Status
  {
  public:
    Status ();

    explicit Status (const data_type &type);

    explicit Status (const status &status);

    Status (const status    &status,
	    const data_type &type);

    Status (const Status &status);

    Status (const Status    &status,
	    const data_type &type);

    status * get() { return &_status; }

    status const * get() const { return &_status; }

    int source () const;

    int tag () const;

    data_type& datatype () { return _datatype; }

    const data_type& datatype () const { return _datatype; }

    unsigned int size (const data_type &type) const;

    unsigned int size () const;

  private:

    status    _status;
    data_type _datatype;
  };


  //-------------------------------------------------------------------
  /**
   * A class that can be subclassed to allow other code to
   * perform work after a MPI_Wait succeeds
   */
  struct PostWaitWork {
    virtual ~PostWaitWork() {}

    virtual void run() {}
  };


  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Request
   */
  class Request
  {
  public:
    Request ();

    Request (const request &r);

    Request (const Request &other);

    void cleanup();

    Request & operator = (const Request &other);

    Request & operator = (const request &r);

    ~Request ();

    request* get() { return &_request; }

    const request* get() const { return &_request; }

    Status wait ();

    bool test ();

    bool test (status &status);

    void add_post_wait_work(PostWaitWork* work);

  private:
    request _request;

    // post_wait_work->first is a vector of work to do after a wait
    // finishes; post_wait_work->second is a reference count so that
    // Request objects will behave roughly like a shared_ptr and be
    // usable in STL containers
    std::pair<std::vector <PostWaitWork* >, unsigned int>* post_wait_work;
  };

  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline Status wait (Request &r) { return r.wait(); }

  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline void wait (std::vector<Request> &r)
  { for (unsigned int i=0; i<r.size(); i++) r[i].wait(); }


  /**
   * Define the data type to be used for data arrays whne encoding
   * a potentially-variable-size object of type T.
   */
  template <typename T>
  struct BufferType;

  /**
   * Encoding non-const T* always uses the same buffer type as
   * encoding const T*, so other classes only have to specialize the
   * latter.
   */
  template <typename T>
  struct BufferType<T*> {
    typedef typename BufferType<const T*>::type type;
  };

  /**
   * Encode a potentially-variable-size object at the end of a data
   * array.
   *
   * Parallel::pack() has no default implementation, and must be
   * specialized for each class which is to be communicated via packed
   * ranges.
   */
  template <typename T, typename buffertype, typename Context>
  void pack(const T* object,
	    typename std::vector<buffertype>& data,
	    const Context* context);

  /**
   * Output the number of integers required to encode a
   * potentially-variable-size object into a data array.
   *
   * Parallel::packable_size() has no default implementation, and must
   * be specialized for each class which is to be communicated via
   * packed ranges.
   */
  template <typename T, typename Context>
  unsigned int packable_size(const T*, const Context*);

  /**
   * Output the number of integers that were used to encode the next
   * variable-size object in the data array.
   *
   * Parallel::packed_size() has no default implementation, and must
   * be specialized for each class which is to be communicated via
   * packed ranges.
   *
   * The output of this method should be based *only* on the data
   * array; the T* argument is solely for function specialization.
   */
  template <typename T, typename BufferIter>
  unsigned int packed_size(const T*,
			   BufferIter);

  /**
   * Decode a potentially-variable-size object from a subsequence of a
   * data array.
   *
   * Parallel::unpack() has no default implementation, and must be
   * specialized for each class which is to be communicated via packed
   * ranges.
   */
  template <typename T, typename BufferIter, typename Context>
  void unpack(BufferIter in, T** out, Context* ctx);

  /**
   * Decode a range of potentially-variable-size objects from a data
   * array.
   */
  template <typename Context, typename buffertype, typename OutputIter>
  inline void unpack_range (const typename std::vector<buffertype>& buffer,
		            Context *context,
		            OutputIter out);

  /**
   * Encode a range of potentially-variable-size objects to a data
   * array.
   */
  template <typename Context, typename buffertype, typename Iter>
  inline void pack_range (const Context *context,
		          Iter range_begin,
		          const Iter range_end,
                          typename std::vector<buffertype>& buffer);

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
    explicit Communicator (const communicator &comm);

    /*
     * NON-VIRTUAL destructor
     */
    ~Communicator ();

    /*
     * Create a new communicator between some subset of \p this
     */
    void split(int color, int key, Communicator &target);

    /*
     * Create a new duplicate of \p this communicator
     */
    void duplicate(const Communicator &comm);

    /*
     * Create a new duplicate of an MPI communicator
     */
    void duplicate(const communicator &comm);

    communicator& get() { return _communicator; }

    const communicator& get() const { return _communicator; }

    /**
     * Get a tag that is unique to this Communicator.  Note that if
     * people are also using magic numbers or copying communicators
     * around then we can't guarantee the tag is unique to this
     * MPI_Comm.
     */
    MessageTag get_unique_tag(int tagvalue) const;

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

    Communicator& operator= (const communicator &comm);

    unsigned int rank() const { return _rank; }

    unsigned int size() const { return _size; }

    /**
     * Whether to use default or synchronous sends?
     */
    enum SendMode { DEFAULT=0, SYNCHRONOUS };

  private:

    // Don't use the copy constructor, just copy by reference or
    // pointer - it's too hard to keep a common used_tag_values if
    // each communicator is shared by more than one Communicator
    explicit Communicator (const Communicator &);

    /**
     * Utility function for setting our member variables from an MPI
     * communicator
     */
    void assign(const communicator &comm);

    communicator  _communicator;
    unsigned int  _rank, _size;
    SendMode _send_mode;

    // mutable used_tag_values - not thread-safe, but then Parallel::
    // isn't thread-safe in general.
    mutable std::map<int, unsigned int> used_tag_values;
    bool          _I_duped_it;

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
     * Verify that a local variable has the same value on all processors.
     * Containers must have the same value in every entry.
     */
    template <typename T>
    bool verify(const T &r) const;

    /**
     * Verify that a local pointer points to the same value on all
     * processors where it is not NULL.
     * Containers must have the same value in every entry.
     */
    template <typename T>
    bool semiverify(const T *r) const;

    /**
     * Take a local variable and replace it with the minimum of it's values
     * on all processors.  Containers are replaced element-wise.
     */
    template <typename T>
    void min(T &r) const;

    /**
     * Take a local variable and replace it with the minimum of it's values
     * on all processors, returning the minimum rank of a processor
     * which originally held the minimum value.
     */
    template <typename T>
    void minloc(T &r,
                unsigned int &min_id) const;

    /**
     * Take a vector of local variables and replace each entry with the minimum
     * of it's values on all processors.  Set each \p min_id entry to
     * the minimum rank where a corresponding minimum was found.
     */
    template <typename T>
    void minloc(std::vector<T> &r,
                std::vector<unsigned int> &min_id) const;

    /**
     * Take a local variable and replace it with the maximum of it's values
     * on all processors.  Containers are replaced element-wise.
     */
    template <typename T>
    void max(T &r) const;

    /**
     * Take a local variable and replace it with the maximum of it's values
     * on all processors, returning the minimum rank of a processor
     * which originally held the maximum value.
     */
    template <typename T>
    void maxloc(T &r,
                unsigned int &max_id) const;

    /**
     * Take a vector of local variables and replace each entry with the maximum
     * of it's values on all processors.  Set each \p min_id entry to
     * the minimum rank where a corresponding maximum was found.
     */
    template <typename T>
    void maxloc(std::vector<T> &r,
                std::vector<unsigned int> &max_id) const;

    /**
     * Take a local variable and replace it with the sum of it's values
     * on all processors.  Containers are replaced element-wise.
     */
    template <typename T>
    void sum(T &r) const;

    /**
     * Take a container of local variables on each processor, and
     * collect their union over all processors, replacing the set on
     * processor 0.
     */
    template <typename T>
    void set_union(T &data, const unsigned int root_id) const;

    /**
     * Take a container of local variables on each processor, and
     * replace it with their union over all processors.
     */
    template <typename T>
    void set_union(T &data) const;

    /**
     * Blocking message probe.  Allows information about a message to be
     * examined before the message is actually received.
     */
    status probe (const unsigned int src_processor_id,
                  const MessageTag &tag=any_tag) const;

    /**
     * Blocking-send to one processor with data-defined type.
     */
    template <typename T>
    void send (const unsigned int dest_processor_id,
               T &buf,
               const MessageTag &tag=no_tag) const;

    /**
     * Nonblocking-send to one processor with data-defined type.
     */
    template <typename T>
    void send (const unsigned int dest_processor_id,
               T &buf,
               Request &req,
               const MessageTag &tag=no_tag) const;

    /**
     * Blocking-send to one processor with user-defined type.
     */
    template <typename T>
    void send (const unsigned int dest_processor_id,
               T &buf,
               const DataType &type,
               const MessageTag &tag=no_tag) const;

    /**
     * Nonblocking-send to one processor with user-defined type.
     */
    template <typename T>
    void send (const unsigned int dest_processor_id,
               T &buf,
               const DataType &type,
               Request &req,
               const MessageTag &tag=no_tag) const;

    /**
     * Blocking-receive from one processor with data-defined type.
     */
    template <typename T>
    Status receive (const unsigned int dest_processor_id,
                    T &buf,
                    const MessageTag &tag=any_tag) const;

    /**
     * Nonblocking-receive from one processor with data-defined type.
     */
    template <typename T>
    void receive (const unsigned int dest_processor_id,
                  T &buf,
                  Request &req,
                  const MessageTag &tag=any_tag) const;

    /**
     * Blocking-receive from one processor with user-defined type.
     */
    template <typename T>
    Status receive (const unsigned int dest_processor_id,
                    T &buf,
                    const DataType &type,
                    const MessageTag &tag=any_tag) const;

    /**
     * Nonblocking-receive from one processor with user-defined type.
     */
    template <typename T>
    void receive (const unsigned int dest_processor_id,
                  T &buf,
                  const DataType &type,
                  Request &req,
                  const MessageTag &tag=any_tag) const;

    /**
     * Blocking-send range-of-pointers to one processor.  This
     * function does not send the raw pointers, but rather constructs
     * new objects at the other end whose contents match the objects
     * pointed to by the sender.
     *
     * void Parallel::pack(const T*, vector<int>& data, const Context*)
     * is used to serialize type T onto the end of a data vector.
     *
     * unsigned int Parallel::packable_size(const T*, const Context*) is
     * used to allow data vectors to reserve memory, and for additional
     * error checking
     */
    template <typename Context, typename Iter>
    void send_packed_range (const unsigned int dest_processor_id,
                            const Context *context,
                            Iter range_begin,
                            const Iter range_end,
                            const MessageTag &tag=no_tag) const;

    /**
     * Nonblocking-send range-of-pointers to one processor.  This
     * function does not send the raw pointers, but rather constructs
     * new objects at the other end whose contents match the objects
     * pointed to by the sender.
     *
     * void Parallel::pack(const T*, vector<int>& data, const Context*)
     * is used to serialize type T onto the end of a data vector.
     *
     * unsigned int Parallel::packable_size(const T*, const Context*) is
     * used to allow data vectors to reserve memory, and for additional
     * error checking
     */
    template <typename Context, typename Iter>
    void send_packed_range (const unsigned int dest_processor_id,
                            const Context *context,
                            Iter range_begin,
                            const Iter range_end,
                            Request &req,
                            const MessageTag &tag=no_tag) const;

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
     * void Parallel::unpack(vector<int>::iterator in, T** out, Context*)
     * is used to unserialize type T, typically into a new
     * heap-allocated object whose pointer is returned as *out.
     *
     * unsigned int Parallel::packed_size(const T*,
     *                                    vector<int>::const_iterator)
     * is used to advance to the beginning of the next object's data.
     */
    template <typename Context, typename OutputIter>
    void receive_packed_range (const unsigned int dest_processor_id,
                               Context *context,
                               OutputIter out,
                               const MessageTag &tag=any_tag) const;

    /**
     * Nonblocking-receive range-of-pointers from one processor.  This
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
     * void Parallel::unpack(vector<int>::iterator in, T** out, Context*)
     * is used to unserialize type T, typically into a new
     * heap-allocated object whose pointer is returned as *out.
     *
     * unsigned int Parallel::packed_size(const T*,
     *                                    vector<int>::const_iterator)
     * is used to advance to the beginning of the next object's data.
     */
    template <typename Context, typename OutputIter>
    void receive_packed_range (const unsigned int dest_processor_id,
                               Context *context,
                               OutputIter out,
                               Request &req,
                               const MessageTag &tag=any_tag) const;

    /**
     * Send data \p send to one processor while simultaneously receiving
     * other data \p recv from a (potentially different) processor.
     */
    template <typename T1, typename T2>
    void send_receive(const unsigned int dest_processor_id,
                      T1 &send,
                      const unsigned int source_processor_id,
                      T2 &recv,
                      const MessageTag &send_tag = no_tag,
                      const MessageTag &recv_tag = any_tag) const;

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
     * void Parallel::pack(const T1*, vector<int>& data, const Context1*)
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
     * void Parallel::unpack(vector<int>::iterator in, T2** out, Context*)
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
    template <typename Context1, typename RangeIter, typename Context2, typename OutputIter>
    void send_receive_packed_range(const unsigned int dest_processor_id,
                                   const Context1* context1,
                                   RangeIter send_begin,
                                   const RangeIter send_end,
                                   const unsigned int source_processor_id,
                                   Context2* context2,
                                   OutputIter out,
                                   const MessageTag &send_tag = no_tag,
                                   const MessageTag &recv_tag = any_tag) const;

    /**
     * Send data \p send to one processor while simultaneously receiving
     * other data \p recv from a (potentially different) processor, using
     * a user-specified MPI Dataype.
     */
    template <typename T1, typename T2>
    void send_receive(const unsigned int dest_processor_id,
                      T1 &send,
                      const DataType &type1,
                      const unsigned int source_processor_id,
                      T2 &recv,
                      const DataType &type2,
                      const MessageTag &send_tag = no_tag,
                      const MessageTag &recv_tag = any_tag) const;

    /**
     * Take a vector of length comm.size(), and on processor root_id fill in
     * recv[processor_id] = the value of send on processor processor_id
     */
    template <typename T>
    inline void gather(const unsigned int root_id,
                       T send,
                       std::vector<T> &recv) const;

    /**
     * Take a vector of local variables and expand it on processor root_id
     * to include values from all processors
     *
     *
     * This handles the
     * case where the lengths of the vectors may vary.
     * Specifically, this function transforms this:
     \verbatim
      Processor 0: [ ... N_0 ]
      Processor 1: [ ....... N_1 ]
        ...
      Processor M: [ .. N_M]
     \endverbatim
     *
     * into this:
     *
     \verbatim
     [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
     \endverbatim
     *
     * on processor root_id. This function is collective and therefore
     * must be called by all processors in the Communicator.
     */
    template <typename T>
    inline void gather(const unsigned int root_id,
                       std::vector<T> &r) const;

    /**
     * Take a vector of length \p this->size(), and fill in
     * \p recv[processor_id] = the value of \p send on that processor
     */
    template <typename T>
    inline void allgather(T send,
                          std::vector<T> &recv) const;


    /**
     * Take a vector of local variables and expand it to include
     * values from all processors. By default, each processor is
     * allowed to have its own unique input buffer length. If
     * it is known that all processors have the same input sizes
     * additional communication can be avoided.
     *
     * Specifically, this function transforms this:
     \verbatim
      Processor 0: [ ... N_0 ]
      Processor 1: [ ....... N_1 ]
        ...
      Processor M: [ .. N_M]
     \endverbatim
     *
     * into this:
     *
     \verbatim
     [ [ ... N_0 ] [ ....... N_1 ] ... [ .. N_M] ]
     \endverbatim
     *
     * on each processor. This function is collective and therefore
     * must be called by all processors in the Communicator.
     */
    template <typename T>
    inline void allgather(std::vector<T> &r,
                          const bool identical_buffer_sizes = false) const;

    //-------------------------------------------------------------------
    /**
     * Take a range of local variables, combine it with ranges from all
     * processors, and write the output to the output iterator.
     */

    template <typename Context, typename Iter, typename OutputIter>
    inline void allgather_packed_range (Context *context,
                                        Iter range_begin,
                                        const Iter range_end,
                                        OutputIter out) const;

    /**
     * Effectively transposes the input vector across all processors.
     * The jth entry on processor i is replaced with the ith entry
     * from processor j.
     */
    template <typename T>
    inline void alltoall(std::vector<T> &r) const;

    /**
     * Take a local value and broadcast it to all processors.
     * Optionally takes the \p root_id processor, which specifies
     * the processor initiating the broadcast.
     * If \p data is a vector, the user is responsible for resizing it
     * on all processors, except in the case when \p data is a vector
     * of strings.
     */
    template <typename T>
    inline void broadcast(T &data, const unsigned int root_id=0) const;

    /**
     * Blocking-broadcast range-of-pointers to one processor.  This
     * function does not send the raw pointers, but rather constructs
     * new objects at the other end whose contents match the objects
     * pointed to by the sender.
     *
     * void Parallel::pack(const T*, vector<int>& data, const Context*)
     * is used to serialize type T onto the end of a data vector.
     *
     * unsigned int Parallel::packable_size(const T*, const Context*) is
     * used to allow data vectors to reserve memory, and for additional
     * error checking
     *
     * unsigned int Parallel::packed_size(const T*,
     *                                    vector<int>::const_iterator)
     * is used to advance to the beginning of the next object's data.
     */
    template <typename Context, typename OutputContext, typename Iter, typename OutputIter>
    inline void broadcast_packed_range (const Context *context1,
                                        Iter range_begin,
                                        const Iter range_end,
                                        OutputContext *context2,
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

  // FakeCommunicator for debugging inappropriate CommWorld uses
  class FakeCommunicator
  {
    operator Communicator& () {
      libmesh_error();
      static Communicator temp;
      return temp;
    }
  };

  // PostWaitWork specialization for copying from temporary to
  // output containers
  template <typename Container, typename OutputIter>
  struct PostWaitCopyBuffer : public PostWaitWork {
    PostWaitCopyBuffer(const Container& buffer, const OutputIter out)
      : _buf(buffer), _out(out) {}

    virtual void run() { std::copy(_buf.begin(), _buf.end(), _out); }

    private:
    const Container& _buf;
    OutputIter _out;
  };

  // PostWaitWork specialization for unpacking received buffers.
  template <typename Container, typename Context, typename OutputIter>
  struct PostWaitUnpackBuffer : public PostWaitWork {
    PostWaitUnpackBuffer(const Container& buffer, Context *context, OutputIter out) :
      _buf(buffer), _context(context), _out(out) {}

    virtual void run() { Parallel::unpack_range(_buf, _context, _out); }

    private:
    const Container& _buf;
    Context *_context;
    OutputIter _out;
  };


  // PostWaitWork specialization for freeing no-longer-needed buffers.
  template <typename Container>
  struct PostWaitDeleteBuffer : public PostWaitWork {
    PostWaitDeleteBuffer(Container* buffer) : _buf(buffer) {}

    virtual void run() { delete _buf; }

    private:
    Container* _buf;
  };

} // namespace Parallel

  /**
   * The default libMesh communicator.
   *
   * If this communicator is disabled, we also disable it as a default
   * argument to functions which accept a default communicator
   * argument.  This should expose implicit uses of the default
   * communicator as compile-time rather than run-time errors.
   *
   * The macro LIBMESH_CAN_DEFAULT_TO_COMMWORLD effects this
   * functionality; it is empty (and so leaves arguments with no
   * default value) if the default is disabled, and it sets something
   * equal to the default otherwise.
   */
#ifdef LIBMESH_DISABLE_COMMWORLD
  extern Parallel::FakeCommunicator CommWorld;
  #define LIBMESH_CAN_DEFAULT_TO_COMMWORLD
#else
  extern Parallel::Communicator CommWorld;
  #define LIBMESH_CAN_DEFAULT_TO_COMMWORLD = libMesh::CommWorld
#endif

} // namespace libMesh

// Define all the implementations separately; users might want to look
// through this file for APIs, and it's long enough already.

#include "libmesh/parallel_implementation.h"

#endif // LIBMESH_PARALLEL_H
