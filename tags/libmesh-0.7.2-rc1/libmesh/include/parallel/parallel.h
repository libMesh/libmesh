// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#ifndef __parallel_h__
#define __parallel_h__

// System includes
#include <set>
#include <string>
#include <vector>
#include <map>

// Local includes
#include "libmesh_common.h" // for libmesh_assert
#include "libmesh_logging.h"

namespace libMesh
{


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only
#ifndef NDEBUG
  #define parallel_only() do { \
    libmesh_assert(Parallel::verify(std::string(__FILE__).size())); \
    libmesh_assert(Parallel::verify(std::string(__FILE__))); \
    libmesh_assert(Parallel::verify(__LINE__)); } while (0)
#else
  #define parallel_only()  ((void) 0)
#endif

// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only_on
#ifndef NDEBUG
  #define parallel_only_on(comm_arg) do { \
    libmesh_assert(Parallel::verify(std::string(__FILE__).size(), comm_arg)); \
    libmesh_assert(Parallel::verify(std::string(__FILE__), comm_arg)); \
    libmesh_assert(Parallel::verify(__LINE__), comm_arg); } while (0)
#else
  #define parallel_only_on(comm_arg)  ((void) 0)
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
  const int any_source=MPI_ANY_SOURCE;

  
#else

  // These shouldn't actually be needed, but must be 
  // unique types for function overloading to work 
  // properly.
  struct data_type    { /* unsigned int t; */ };
  struct request      { /* unsigned int r; */ };
  struct status       { /* unsigned int s; */ };
  struct communicator { /* unsigned int s; */ };

  const int any_source=0;
#endif // LIBMESH_HAVE_MPI



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI tag integers.
   */
  class MessageTag
  {
  public:
    /**
     * Explicit constructor, to discourage using "magic numbers"
     * as tags.  Communicator::get_unique_tag is recommended instead.
     */
    explicit MessageTag(int tagvalue)
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
   * attributes.  More specialized classes can override them.
   */
  template<typename T>
  struct Attributes
  {
    const static bool has_min_max = false;
  };

  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Comm object.  Allows the size of the group
   * and this process's position in the group to be determined.
   */
  class Communicator
  {
  public:
    Communicator () :
#ifdef LIBMESH_HAVE_MPI
      _communicator(MPI_COMM_NULL),
      _rank(0),
      _size(0),
#else
      _rank(0), 
      _size(1), 
#endif
      used_tag_values(),
      _I_duped_it(false) {}

    /*
     * Constructor from MPI_Comm
     */
    explicit Communicator (const communicator &comm) : 
#ifdef LIBMESH_HAVE_MPI
      _communicator(MPI_COMM_NULL),
      _rank(0),
      _size(0),
#else
      _rank(0), 
      _size(1), 
#endif
      used_tag_values(),
      _I_duped_it(false)
    {
      this->assign(comm);
    }

    ~Communicator () {
      this->clear();
    }

    /*
     * Create a new communicator between some subset of \p this
     */
#ifdef LIBMESH_HAVE_MPI
    void split(int color, int key, Communicator &target) {
      MPI_Comm_split(this->get(), color, key, &target.get());
    }
#else
    void split(int, int, Communicator &target) {
      target.assign(this->get());
    }
#endif

    void duplicate(const Communicator &comm) {
      this->duplicate(comm._communicator);
    }

#ifdef LIBMESH_HAVE_MPI
    void duplicate(const communicator &comm) {
      if (_communicator != MPI_COMM_NULL) 
        {
          MPI_Comm_dup(comm, &_communicator);
          _I_duped_it = true;
        }
      this->assign(_communicator);
    }
#else
    void duplicate(const communicator &) {
    }
#endif

    communicator& get() {
      return _communicator;
    }

    const communicator& get() const {
      return _communicator;
    }

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

    void clear() {
#ifdef LIBMESH_HAVE_MPI
      if (_I_duped_it)
        {
          libmesh_assert(_communicator != MPI_COMM_NULL);
          MPI_Comm_free(&_communicator);
          _communicator = MPI_COMM_NULL;
        }
      _I_duped_it = false;
#endif
    }

    Communicator& operator= (const communicator &comm) {
      this->clear();
      this->assign(comm);
      return *this;
    }

    unsigned int rank() const {
      return _rank;
    }

    unsigned int size() const {
      return _size;
    }

  private:

    // Don't use the copy constructor, just copy by reference or
    // pointer - it's too hard to keep a common used_tag_values if
    // each communicator is shared by more than one Communicator
    explicit Communicator (const Communicator &) :
#ifdef LIBMESH_HAVE_MPI
      _communicator(MPI_COMM_NULL),
      _rank(0),
      _size(0),
#else
      _rank(0), 
      _size(1), 
#endif
      used_tag_values(),
      _I_duped_it(false)
    {
      libmesh_error();
    }

    void assign(const communicator &comm) {
      _communicator = comm;
#ifdef LIBMESH_HAVE_MPI
      if (_communicator != MPI_COMM_NULL)
        {
          int i;
          MPI_Comm_size(_communicator, &i);
          libmesh_assert (i >= 0);
          _size = static_cast<unsigned int>(i);

          MPI_Comm_rank(_communicator, &i);
          libmesh_assert (i >= 0);
          _rank = static_cast<unsigned int>(i);
        }
      else
        {
          _rank = 0;
          _size = 0;
        }
#endif
    }

    communicator  _communicator;
    unsigned int  _rank, _size;

    // mutable used_tag_values - not thread-safe, but then Parallel::
    // isn't thread-safe in general.
    mutable std::map<int, unsigned int> used_tag_values;
    bool          _I_duped_it;
  };



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Status struct.  Allows the source and size 
   * of the message to be determined.
   */
  class Status
  {
  public:
    Status () :
      _status(),
      _datatype()           
    {}
    
    Status (const data_type &type) :
      _status(),
      _datatype(type)           
    {}

    Status (const status &status) :
      _status(status),
      _datatype()           
    {}

    Status (const status    &status,
	    const data_type &type) :
      _status(status),
      _datatype(type)           
    {}

    Status (const Status &status) :
      _status(status._status),
      _datatype(status._datatype)
    {}

    Status (const Status    &status,
	    const data_type &type) :
      _status(status._status),
      _datatype(type)
    {}
    
    operator status * () 
    { return &_status; }

    operator status const * () const
    { return &_status; }
    
//     operator status & ()
//     { return _status; }

//     operator const status & () const
//     { return _status; }

    int source () const
    {
#ifdef LIBMESH_HAVE_MPI 
      return _status.MPI_SOURCE; 
#else
      return 0;
#endif
    }    

    int tag () const
    {
#ifdef LIBMESH_HAVE_MPI 
      return _status.MPI_TAG; 
#else
      libmesh_error();
      return 0;
#endif
    }

    data_type& datatype () 
    { return _datatype; }

    const data_type& datatype () const
    { return _datatype; }
  
#ifdef LIBMESH_HAVE_MPI
    unsigned int size (const data_type &type) const
    {
      int msg_size;
      MPI_Get_count (const_cast<MPI_Status*>(&_status), type, &msg_size);
      libmesh_assert (msg_size >= 0);
      return msg_size;
    }
#else
    unsigned int size (const data_type &) const
    {
      libmesh_error();
      return 0;
    }
#endif

    unsigned int size () const
    { return this->size (this->datatype()); }

  private:

    status    _status;
    data_type _datatype;
  };



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Request
   */
  class Request
  {
  public:
    Request () :
#ifdef LIBMESH_HAVE_MPI
      _request(MPI_REQUEST_NULL)
#else
      _request()
#endif
    {}

    Request (const request &r) :
      _request(r)
    {}

    Request & operator = (const Request &other) 
    { _request = other._request; return *this; }

    Request & operator = (const request &r)
    { _request = r; return *this; }

    ~Request () {}

    operator const request & () const
    { return _request; }

    operator request & () 
    { return _request; }

    request* get()
    { return &_request; }

    void free (void)
    {
#ifdef LIBMESH_HAVE_MPI
      if (_request != MPI_REQUEST_NULL)
	MPI_Request_free (&_request);
#endif
    }

    const request* get() const
    { return &_request; }

    status wait ()
    {
      status status;
#ifdef LIBMESH_HAVE_MPI
      MPI_Wait (&_request, &status);
#endif
      return status;
    }

    bool test ()
    {
#ifdef LIBMESH_HAVE_MPI
      int val=0;

      MPI_Test (&_request,
		&val,
		MPI_STATUS_IGNORE);
      if (val)
	{
	  libmesh_assert (_request == MPI_REQUEST_NULL);
	  libmesh_assert (val == 1);
	}

      return val;
#else
      return true;
#endif
    }

#ifdef LIBMESH_HAVE_MPI
    bool test (status &status)
    {
      int val=0;

      MPI_Test (&_request,
		&val,
		&status);

      return val;
    }
#else
    bool test (status &)
    {
      return true;
    }
#endif


  private:

    request _request;    
  };


  
  //-------------------------------------------------------------------
  /* 
   * The default libMesh communicator
   */
  extern Communicator Communicator_World;



  //-------------------------------------------------------------------
  /**
   * Pause execution until all processors reach a certain point.
   */
#ifdef LIBMESH_HAVE_MPI
  inline void barrier (const Communicator &comm = Communicator_World)
  {
    if (libMesh::n_processors() > 1)
      {
	START_LOG("barrier()", "Parallel");
    
	MPI_Barrier (comm.get());

	STOP_LOG("barrier()", "Parallel");
      }
  }    
#else
  inline void barrier (const Communicator & = Communicator_World) {}
#endif
  
  //-------------------------------------------------------------------
  /**
   * Verify that a local variable has the same value on all processors
   */
  template <typename T>
  inline bool verify(const T &r,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors
   */
  template <typename T>
  inline void min(T &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors
   */
  template <typename T>
  inline void min(std::vector<T> &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors
   */
  template <typename T>
  inline void minloc(T &r,
                     unsigned int &min_id,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors
   */
  template <typename T>
  inline void minloc(std::vector<T> &r,
                     std::vector<unsigned int> &min_id,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors
   */
  template <typename T>
  inline void max(T &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors
   */
  template <typename T>
  inline void max(std::vector<T> &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors, returning the minimum rank of a processor
   * which originally held the maximum value.
   */
  template <typename T>
  inline void maxloc(T &r,
                     unsigned int &max_id,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors
   */
  template <typename T>
  inline void maxloc(std::vector<T> &r,
                     std::vector<unsigned int> &max_id,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the sum of it's values
   * on all processors
   */
  template <typename T>
  inline void sum(T &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the sum of
   * it's values on all processors
   */
  template <typename T>
  inline void sum(std::vector<T> &r,
                  const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and collect their
   * union over all processors, replacing the set on processor 0.
   */
  template <typename T>
  inline void set_union(std::set<T> &data, const unsigned int root_id,
                        const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and replace it
   * with their union over all processors.
   */
  template <typename T>
  inline void set_union(std::set<T> &data,
                        const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and collect their
   * union over all processors, replacing the set on processor 0.
   */
  template <typename T1, typename T2>
  inline void set_union(std::map<T1,T2> &data, const unsigned int root_id,
                        const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and replace it
   * with their union over all processors.
   */
  template <typename T1, typename T2>
  inline void set_union(std::map<T1, T2> &data,
                        const Communicator &comm = Communicator_World);


  //-------------------------------------------------------------------
  /**
   * Blocking message probe.  Allows information about a message to be 
   * examined before the message is actually received.
   */
  inline status probe (const int src_processor_id,
		       const MessageTag &tag=any_tag,
                       const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    Request &req,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  StandardType<T>(buf.empty() ? NULL : &buf[0]),
	  tag,
          comm);
  }

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    Request &req,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  StandardType<T>(buf.empty() ? NULL : &buf[0]),
	  req,
	  tag,
          comm);
  }


  //-------------------------------------------------------------------
  /**
   * Blocking-send set to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send set to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    Request &req,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send set to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  StandardType<T>(buf.empty() ? NULL : &(*buf.begin())),
	  tag,
          comm);
  }


  //-------------------------------------------------------------------
  /**
   * Nonblocking-send set to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    Request &req,
		    const MessageTag &tag=no_tag,
                    const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  StandardType<T>(buf.empty() ? NULL : &(*buf.begin())),
	  req,
	  tag,
          comm);
  }


  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void nonblocking_send (const unsigned int dest_processor_id,
		                std::vector<T> &buf,
		                const DataType &type,
		                Request &r,
		                const MessageTag &tag=no_tag,
                                const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  type,
	  r,
	  tag,
          comm);
  }

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor.
   */
  template <typename T>
  inline void nonblocking_send (const unsigned int dest_processor_id,
		                std::vector<T> &buf,
		                Request &r,
		                const MessageTag &tag=no_tag,
                                const Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  StandardType<T>(buf.empty() ? NULL : &buf[0]),
	  r,
	  tag,
          comm);
  }


  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor with user-defined type.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<T> &buf,
		         const DataType &type,
		         const MessageTag &tag=any_tag,
                         const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor with user-defined type.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<T> &buf,
		       const DataType &type,
		       Request &req,
		       const MessageTag &tag=any_tag,
                       const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<T> &buf,
		         const MessageTag &tag=any_tag,
                         const Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
	            StandardType<T>(buf.empty() ? NULL : &buf[0]),
		    tag,
                    comm);
  }


  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<T> &buf,
		       Request &req,
		       const MessageTag &tag=any_tag,
                       const Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     StandardType<T>(buf.empty() ? NULL : &buf[0]),
	     req,
	     tag,
             comm);
  }


  //-------------------------------------------------------------------
  /**
   * Blocking-receive set from one processor with user-defined type.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::set<T> &buf,
		         const DataType &type,
		         const MessageTag &tag=any_tag,
                         const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive set from one processor with user-defined type.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<T> &buf,
		       const DataType &type,
		       Request &req,
		       const MessageTag &tag=any_tag,
                       const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive set from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::set<T> &buf,
		         const MessageTag &tag=any_tag,
                         const Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
	            StandardType<T>(buf.empty() ? NULL : &(*buf.begin())),
		    tag,
                    comm);
  }


  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive set from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<T> &buf,
		       Request &req,
		       const MessageTag &tag=any_tag,
                       const Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     StandardType<T>(buf.empty() ? NULL : &(*buf.begin())),
	     req,
	     tag,
             comm);
  }


  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor with user-defined type
   */
  template <typename T>
  inline void nonblocking_receive (const int src_processor_id,
		                   std::vector<T> &buf,
				   const DataType &type,
		                   Request &r,
		                   const MessageTag &tag=any_tag,
                                   const Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     type,
	     r,
	     tag,
             comm);
  }

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor.
   */
  template <typename T>
  inline void nonblocking_receive (const int src_processor_id,
		                   std::vector<T> &buf,
		                   Request &r,
		                   const MessageTag &tag=any_tag,
                                   const Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     StandardType<T>(buf.empty() ? NULL : &buf[0]),
	     r,
	     tag,
             comm);
  }

  
  //-------------------------------------------------------------------
  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline status wait (request &r);
  
  //-------------------------------------------------------------------
  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline void wait (std::vector<request> &r);
  
  //-------------------------------------------------------------------
  /**
   * Wait for a non-blocking send or receive to finish
   */
  inline void wait (std::vector<Request> &r)
  { for (unsigned int i=0; i<r.size(); i++) r[i].wait(); }
  
  //-------------------------------------------------------------------
  /**
   * Send vector \p send to one processor while simultaneously receiving
   * another vector \p recv from a (potentially different) processor.
   */
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
                           T1 &send,
			   const unsigned int source_processor_id,
                           T2 &recv,
                           const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Send vector \p send to one processor while simultaneously receiving
   * another vector \p recv from a (potentially different) processor using
   * a user-specified MPI Dataype.
   */
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
                           T1 &send,
			   const DataType &type1,
			   const unsigned int source_processor_id,
                           T2 &recv,
			   const DataType &type2,
                           const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length comm.size(), and on processor root_id fill in
   * recv[processor_id] = the value of send on processor processor_id
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it on processor root_id 
   * to include values from all processors
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::vector<T> &r,
                     const Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length comm.size(), and fill in 
   * \p recv[processor_id] = the value of \p send on that processor
   */
  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        const Communicator &comm = Communicator_World);


  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it to include 
   * values from all processors. By default, each processor is 
   * allowed to have its own unique input buffer length. If 
   * it is known that all processors have the same input sizes
   * additional communication can be avoided.
   */
  template <typename T>
  inline void allgather(std::vector<T> &r,
			const bool identical_buffer_sizes = false,
                        const Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Effectively transposes the input vector across all processors.  
   * The jth entry on processor i is replaced with the ith entry 
   * from processor j.
   */
  template <typename T>
  inline void alltoall(std::vector<T> &r,
                       const Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Take a local value and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(T &data, const unsigned int root_id=0,
                        const Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Take a local vector and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.  The user is responsible
   * for appropriately sizing the input buffer on all processors.
   */
  template <typename T>
  inline void broadcast(std::vector<T> &data, const unsigned int root_id=0,
                        const Communicator &comm = Communicator_World);


  //-------------------------------------------------------------------
  /**
   * Take a local set and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(std::set<T> &data, const unsigned int root_id=0,
                        const Communicator &comm = Communicator_World);


  //-----------------------------------------------------------------------
  // Parallel members

  inline
  MessageTag::~MessageTag()
  {
    if (_comm)
      _comm->dereference_unique_tag(_tagvalue);
  }


  inline
  MessageTag::MessageTag(const MessageTag &other)
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

#ifndef NDEBUG
    // Make sure everyone called get_unique_tag and make sure
    // everyone got the same value
    int maxval = tagvalue;
    Parallel::max(maxval);
    libmesh_assert(tagvalue == maxval);
#endif

    return MessageTag(tagvalue, this);
  }


  inline
  void Communicator::reference_unique_tag(int tagvalue) const
  {
    // This has better be an already-acquired tag.
    libmesh_assert(used_tag_values.count(tagvalue));

    used_tag_values[tagvalue]++;
  }


  inline
  void Communicator::dereference_unique_tag(int tagvalue) const
  {
    // This has better be an already-acquired tag.
    libmesh_assert(used_tag_values.count(tagvalue));

    used_tag_values[tagvalue]--;
    // If we don't have any more outstanding references, we
    // don't even need to keep this tag in our "used" set.
    if (!used_tag_values[tagvalue])
      used_tag_values.erase(tagvalue);
  }


  // Internal helper function to create vector<something_useable> from
  // vector<bool> for compatibility with MPI bitwise operations
  template <typename T>
  inline void pack_vector_bool(const std::vector<bool> &in,
			       std::vector<T> &out)
  {
    unsigned int data_bits = 8*sizeof(T);
    unsigned int in_size = in.size();
    unsigned int out_size = in_size/data_bits + (in_size%data_bits?1:0);
    out.clear();
    out.resize(out_size);
    for (unsigned int i=0; i != in_size; ++i)
      {
        unsigned int index = i/data_bits;
        unsigned int offset = i%data_bits;
        out[index] += (in[i]?1:0) << offset;
      }
  }

  // Internal helper function to create vector<bool> from
  // vector<something usable> for compatibility with MPI byte
  // operations
  template <typename T>
  inline void unpack_vector_bool(const std::vector<T> &in,
				 std::vector<bool> &out)
  {
    unsigned int data_bits = 8*sizeof(T);
    // We need the output vector to already be properly sized
    unsigned int out_size = out.size();
    libmesh_assert(out_size/data_bits + (out_size%data_bits?1:0) == in.size());

    for (unsigned int i=0; i != out_size; ++i)
      {
        unsigned int index = i/data_bits;
        unsigned int offset = i%data_bits;
        out[i] = in[index] << (data_bits-1-offset) >> (data_bits-1);
      }
  }

#ifdef LIBMESH_HAVE_MPI

#define STANDARD_TYPE(cxxtype,mpitype) \
  template<> \
  class StandardType<cxxtype> : public DataType \
  { \
  public: \
    StandardType(const cxxtype* example = NULL) : DataType(mpitype) \
      { /* make compiler think 'example' is used */ libmesh_ignore(example); } \
  }; \
 \
  template<> \
  struct Attributes<cxxtype> \
  { \
    const static bool has_min_max = true; \
  }

#else

#define STANDARD_TYPE(cxxtype,mpitype) \
  template<> \
  class StandardType<cxxtype> : public DataType \
  { \
  public: \
    StandardType(const cxxtype* example = NULL) : DataType() \
      { /* make compiler think 'example' is used */ libmesh_ignore(example); } \
  }; \
 \
  template<> \
  struct Attributes<cxxtype> \
  { \
    const static bool has_min_max = true; \
  }

#endif

  STANDARD_TYPE(char,MPI_CHAR);
  STANDARD_TYPE(unsigned char,MPI_UNSIGNED_CHAR);
  STANDARD_TYPE(short int,MPI_SHORT);
  STANDARD_TYPE(unsigned short int,MPI_UNSIGNED_SHORT);
  STANDARD_TYPE(int,MPI_INT);
  STANDARD_TYPE(unsigned int,MPI_UNSIGNED);
  STANDARD_TYPE(long,MPI_LONG);
  STANDARD_TYPE(unsigned long,MPI_UNSIGNED_LONG);
  STANDARD_TYPE(float,MPI_FLOAT);
  STANDARD_TYPE(double,MPI_DOUBLE);
  STANDARD_TYPE(long double,MPI_LONG_DOUBLE);

  // We'd love to do a singleton pattern on derived data types, rather
  // than commit, free, commit, free, ad infinitum... but it's a
  // little tricky when we T1 and T2 are undefined.
  template<typename T1, typename T2>
  class StandardType<std::pair<T1, T2> > : public DataType
  {
  public:
    StandardType(const std::pair<T1, T2> *example = NULL) {
      // We need an example for MPI_Address to use
      libmesh_assert(example);

#ifdef LIBMESH_HAVE_MPI
      // Get the sub-data-types, and make sure they live long enough
      // to construct the derived type
      StandardType<T1> d1(&example->first);
      StandardType<T2> d2(&example->second);
      MPI_Datatype types[] = { (data_type)d1, (data_type)d2 };
      int blocklengths[] = {1,1};

      MPI_Aint displs[2];
      MPI_Address (const_cast<T1*>(&example->first), &displs[0]);
      MPI_Address (const_cast<T2*>(&example->second), &displs[1]);
      displs[1] -= displs[0];
      displs[0] = 0;

#if MPI_VERSION > 1
      MPI_Type_create_struct (2, blocklengths, displs, types, &_datatype);
#else
      MPI_Type_struct (2, blocklengths, displs, types, &_datatype);
#endif // #if MPI_VERSION > 1
      MPI_Type_commit (&_datatype);
#endif // LIBMESH_HAVE_MPI
    }

    inline ~StandardType() { this->free(); }
  };

  template<typename T>
  class StandardType<std::complex<T> > : public DataType
  {
  public:
    inline StandardType(const std::complex<T> *example = NULL) :
      DataType(StandardType<T>(example ? &(example->real()) : NULL), 2) {}

    inline ~StandardType() { this->free(); }
  };

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
  inline bool verify(const T &r,
                     const Communicator &comm)
  {
    if (comm.size() > 1 && Attributes<T>::has_min_max == true)
      {
	T tempmin = r, tempmax = r;
	Parallel::min(tempmin, comm);
	Parallel::max(tempmax, comm);
	bool verified = (r == tempmin) &&
	                (r == tempmax);
	Parallel::min(verified, comm);
	return verified;
      }
    return true;
  }


  template <>
  inline bool verify(const std::string & r,
                     const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	// Cannot use <char> since MPI_MIN is not
	// strictly defined for chars!
	std::vector<short int> temp; temp.reserve(r.size());
	for (unsigned int i=0; i != r.size(); ++i)
	  temp.push_back(r[i]);
	return Parallel::verify(temp, comm);
      }
    return true;
  }


  template <typename T>
  inline void min(T &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("min(scalar)", "Parallel");
    
	T temp = r;
	MPI_Allreduce (&temp,
		       &r,
		       1,
		       StandardType<T>(&temp),
		       MPI_MIN,
		       comm.get());

	STOP_LOG("min(scalar)", "Parallel");
      }
  }


  template <>
  inline void min(bool &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("min(bool)", "Parallel");
    
	unsigned int tempsend = r;
	unsigned int temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       StandardType<unsigned int>(),
		       MPI_MIN,
		       comm.get());
	r = temp;

	STOP_LOG("min(bool)", "Parallel");
      }
  }


  template <typename T>
  inline void min(std::vector<T> &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("min(vector)", "Parallel");

        libmesh_assert(verify(r.size(), comm));
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       StandardType<T>(&temp[0]),
		       MPI_MIN,
		       comm.get());

	STOP_LOG("min(vector)", "Parallel");
      }
  }


  template <>
  inline void min(std::vector<bool> &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("min(vector<bool>)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
        std::vector<unsigned int> ruint;
        pack_vector_bool(r, ruint);
	std::vector<unsigned int> temp(ruint.size());
	MPI_Allreduce (&ruint[0],
		       &temp[0],
		       ruint.size(),
		       StandardType<unsigned int>(),
		       MPI_BAND,
		       comm.get());
        unpack_vector_bool(temp, r);

	STOP_LOG("min(vector<bool>)", "Parallel");
      }
  }


  template <typename T>
  inline void minloc(T &r,
                     unsigned int &min_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("minloc(scalar)", "Parallel");
    
	DataPlusInt<T> in;
        in.val = r;
        in.rank = comm.rank();
	DataPlusInt<T> out;
	MPI_Allreduce (&in,
		       &out,
		       1,
		       dataplusint_type<T>(),
		       MPI_MINLOC,
		       comm.get());
	r = out.val;
        min_id = out.rank;

	STOP_LOG("minloc(scalar)", "Parallel");
      }
    else
      min_id = comm.rank();
  }


  template <>
  inline void minloc(bool &r,
                     unsigned int &min_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("minloc(bool)", "Parallel");
    
	DataPlusInt<int> in;
        in.val = r;
        in.rank = comm.rank();
	DataPlusInt<int> out;
	MPI_Allreduce (&in,
		       &out,
		       1,
		       dataplusint_type<int>(),
		       MPI_MINLOC,
		       comm.get());
	r = out.val;
        min_id = out.rank;

	STOP_LOG("minloc(bool)", "Parallel");
      }
    else
      min_id = comm.rank();
  }


  template <typename T>
  inline void minloc(std::vector<T> &r,
                     std::vector<unsigned int> &min_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("minloc(vector)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<DataPlusInt<T> > in(r.size());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            in[i].val  = r[i];
            in[i].rank = comm.rank();
          }
	std::vector<DataPlusInt<T> > out(r.size());
	MPI_Allreduce (&in[0],
		       &out[0],
		       r.size(),
		       dataplusint_type<T>(),
		       MPI_MINLOC,
		       comm.get());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            r[i]      = out[i].val;
            min_id[i] = out[i].rank;
          }

	STOP_LOG("minloc(vector)", "Parallel");
      }
    else if (!r.empty())
      {
        for (unsigned int i=0; i != r.size(); ++i)
          min_id[i] = comm.rank();
      }
  }


  template <>
  inline void minloc(std::vector<bool> &r,
                     std::vector<unsigned int> &min_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("minloc(vector<bool>)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<DataPlusInt<int> > in(r.size());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            in[i].val  = r[i];
            in[i].rank = comm.rank();
          }
	std::vector<DataPlusInt<int> > out(r.size());
	MPI_Allreduce (&in[0],
		       &out[0],
		       r.size(),
		       StandardType<int>(),
		       MPI_MINLOC,
		       comm.get());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            r[i]      = out[i].val;
            min_id[i] = out[i].rank;
          }

	STOP_LOG("minloc(vector<bool>)", "Parallel");
      }
    else if (!r.empty())
      {
        for (unsigned int i=0; i != r.size(); ++i)
          min_id[i] = comm.rank();
      }
  }


  template <typename T>
  inline void max(T &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("max(scalar)", "Parallel");
    
	T temp;
	MPI_Allreduce (&r,
		       &temp,
		       1,
		       StandardType<T>(&r),
		       MPI_MAX,
		       comm.get());
	r = temp;

	STOP_LOG("max(scalar)", "Parallel");
      }
  }


  template <>
  inline void max(bool &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("max(bool)", "Parallel");
    
	unsigned int tempsend = r;
	unsigned int temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       StandardType<unsigned int>(),
		       MPI_MAX,
		       comm.get());
	r = temp;

	STOP_LOG("max(bool)", "Parallel");
      }
  }


  template <typename T>
  inline void max(std::vector<T> &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("max(vector)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       StandardType<T>(&temp[0]),
		       MPI_MAX,
		       comm.get());

	STOP_LOG("max(vector)", "Parallel");
      }
  }


  template <>
  inline void max(std::vector<bool> &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("max(vector<bool>)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
        std::vector<unsigned int> ruint;
        pack_vector_bool(r, ruint);
	std::vector<unsigned int> temp(ruint.size());
	MPI_Allreduce (&ruint[0],
		       &temp[0],
		       ruint.size(),
		       StandardType<unsigned int>(),
		       MPI_BOR,
		       comm.get());
        unpack_vector_bool(temp, r);

	STOP_LOG("max(vector<bool>)", "Parallel");
      }
  }


  template <typename T>
  inline void maxloc(T &r,
                     unsigned int &max_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("maxloc(scalar)", "Parallel");
    
	DataPlusInt<T> in;
        in.val = r;
        in.rank = comm.rank();
	DataPlusInt<T> out;
	MPI_Allreduce (&in,
		       &out,
		       1,
		       dataplusint_type<T>(),
		       MPI_MAXLOC,
		       comm.get());
	r = out.val;
        max_id = out.rank;

	STOP_LOG("maxloc(scalar)", "Parallel");
      }
    else
      max_id = comm.rank();
  }


  template <>
  inline void maxloc(bool &r,
                     unsigned int &max_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("maxloc(bool)", "Parallel");
    
	DataPlusInt<int> in;
        in.val = r;
        in.rank = comm.rank();
	DataPlusInt<int> out;
	MPI_Allreduce (&in,
		       &out,
		       1,
		       dataplusint_type<int>(),
		       MPI_MAXLOC,
		       comm.get());
	r = out.val;
        max_id = out.rank;

	STOP_LOG("maxloc(bool)", "Parallel");
      }
    else
      max_id = comm.rank();
  }


  template <typename T>
  inline void maxloc(std::vector<T> &r,
                     std::vector<unsigned int> &max_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("maxloc(vector)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<DataPlusInt<T> > in(r.size());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            in[i].val  = r[i];
            in[i].rank = comm.rank();
          }
	std::vector<DataPlusInt<T> > out(r.size());
	MPI_Allreduce (&in[0],
		       &out[0],
		       r.size(),
		       dataplusint_type<T>(),
		       MPI_MAXLOC,
		       comm.get());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            r[i]      = out[i].val;
            max_id[i] = out[i].rank;
          }

	STOP_LOG("maxloc(vector)", "Parallel");
      }
    else if (!r.empty())
      {
        for (unsigned int i=0; i != r.size(); ++i)
          max_id[i] = comm.rank();
      }
  }


  template <>
  inline void maxloc(std::vector<bool> &r,
                     std::vector<unsigned int> &max_id,
                     const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("maxloc(vector<bool>)", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<DataPlusInt<int> > in(r.size());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            in[i].val  = r[i];
            in[i].rank = comm.rank();
          }
	std::vector<DataPlusInt<int> > out(r.size());
	MPI_Allreduce (&in[0],
		       &out[0],
		       r.size(),
		       StandardType<int>(),
		       MPI_MAXLOC,
		       comm.get());
        for (unsigned int i=0; i != r.size(); ++i)
          {
            r[i]      = out[i].val;
            max_id[i] = out[i].rank;
          }

	STOP_LOG("maxloc(vector<bool>)", "Parallel");
      }
    else if (!r.empty())
      {
        for (unsigned int i=0; i != r.size(); ++i)
          max_id[i] = comm.rank();
      }
  }


  template <typename T>
  inline void sum(T &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("sum()", "Parallel");
    
	T temp = r;
	MPI_Allreduce (&temp,
		       &r,
		       1,
		       StandardType<T>(&temp),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  template <typename T>
  inline void sum(std::vector<T> &r,
                  const Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("sum()", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       StandardType<T>(&temp[0]),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  // We still do function overloading for complex sums - in a perfect
  // world we'd have a StandardSumOp to go along with StandardType...
  template <typename T>
  inline void sum(std::complex<T> &r,
                  const Communicator &comm = Communicator_World)
  {
    if (comm.size() > 1)
      {
	START_LOG("sum()", "Parallel");
    
	std::complex<T> temp(r);
	MPI_Allreduce (&temp,
		       &r,
		       2,
		       StandardType<T>(&(r.real())),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  template <typename T>
  inline void sum(std::vector<std::complex<T> > &r,
                  const Communicator &comm = Communicator_World)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("sum()", "Parallel");
    
        libmesh_assert(verify(r.size(), comm));
    
	std::vector<std::complex<T> > temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size() * 2,
		       StandardType<T>(&(r[0].real())),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  template <typename T>
  inline void set_union(std::set<T> &data, const unsigned int root_id,
                        const Communicator &comm)
  {
    std::vector<T> vecdata(data.begin(), data.end());
    Parallel::gather(root_id, vecdata, comm);
    if (comm.rank() == root_id)
      data.insert(vecdata.begin(), vecdata.end());
  }



  template <typename T>
  inline void set_union(std::set<T> &data,
                        const Communicator &comm)
  {
    std::vector<T> vecdata(data.begin(), data.end());
    Parallel::allgather(vecdata, false, comm);
    data.insert(vecdata.begin(), vecdata.end());
  }



  template <typename T1, typename T2>
  inline void set_union(std::map<T1,T2> &data, const unsigned int root_id,
                        const Communicator &comm)
  {
    std::vector<std::pair<T1,T2> > vecdata(data.begin(), data.end());
    Parallel::gather(root_id, vecdata, comm);
    if (comm.rank() == root_id)
      data.insert(vecdata.begin(), vecdata.end());
  }



  template <typename T1, typename T2>
  inline void set_union(std::map<T1,T2> &data,
                        const Communicator &comm)
  {
    std::vector<std::pair<T1,T2> > vecdata(data.begin(), data.end());
    Parallel::allgather(vecdata, false, comm);
    data.insert(vecdata.begin(), vecdata.end());
  }



  inline status probe (const int src_processor_id,
		       const MessageTag &tag,
                       const Communicator &comm)
  {
    START_LOG("probe()", "Parallel");

    status status;

    MPI_Probe (src_processor_id, 
	       tag.value(), 
	       comm.get(), 
	       &status);

    STOP_LOG("probe()", "Parallel");

    return status;
  }  



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    const MessageTag &tag,
                    const Communicator &comm)
  {    
    START_LOG("send()", "Parallel");

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Send (buf.empty() ? NULL : &buf[0],
		buf.size(),
		type,
		dest_processor_id,
		tag.value(),
		comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);    
    
    STOP_LOG("send()", "Parallel");
  }


  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    Request &req,
		    const MessageTag &tag,
                    const Communicator &comm)
  {    
    START_LOG("send()", "Parallel");

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Isend (buf.empty() ? NULL : &buf[0],
		 buf.size(),
		 type,
		 dest_processor_id,
		 tag.value(),
		 comm.get(),
		 req.get());    
    libmesh_assert (ierr == MPI_SUCCESS);    

    STOP_LOG("send()", "Parallel");
  }


  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<T> &buf,
		         const DataType &type,
		         const MessageTag &tag,
                         const Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

    // Get the status of the message, explicitly provide the
    // datatype so we can later query the size
    Status status(Parallel::probe(src_processor_id, tag, comm), type);

    buf.resize(status.size());
    
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Recv (buf.empty() ? NULL : &buf[0],
		buf.size(),
		type,
		src_processor_id,
		tag.value(),
		comm.get(),
		status);
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("receive()", "Parallel");

    return status;
  }



  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<T> &buf,
		       const DataType &type,
		       Request &req,
		       const MessageTag &tag,
                       const Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Irecv (buf.empty() ? NULL : &buf[0],
		 buf.size(),
		 type,
		 src_processor_id,
		 tag.value(),
		 comm.get(),
		 req.get());    
    libmesh_assert (ierr == MPI_SUCCESS);    

    STOP_LOG("receive()", "Parallel");
  }



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    const MessageTag &tag,
                    const Communicator &comm)
  {    
    START_LOG("send()", "Parallel");

    std::vector<T> vecbuf(buf.begin(), buf.end());
    // Use Parallel::send() so we get its specialization(s)
    Parallel::send(dest_processor_id, vecbuf, type, tag, comm);

    STOP_LOG("send()", "Parallel");
  }



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    Request &req,
		    const MessageTag &tag,
                    const Communicator &comm)
  {    
    START_LOG("send()", "Parallel");

    std::vector<T> vecbuf(buf.begin(), buf.end());
    // Use Parallel::send() so we get its specialization(s)
    Parallel::send(dest_processor_id, vecbuf, type, req, tag, comm);

    STOP_LOG("send()", "Parallel");
  }



  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::set<T> &buf,
		         const DataType &type,
		         const MessageTag &tag,
                         const Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

    std::vector<T> vecbuf;
    // Use Parallel::receive() so we get its specialization(s)
    Status status = Parallel::receive(src_processor_id, vecbuf, type, tag, comm);
    buf.clear();
    buf.insert(vecbuf.begin(), vecbuf.end());

    STOP_LOG("receive()", "Parallel");

    return status;
  }



  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<T> &buf,
		       const DataType &type,
		       Request &req,
		       const MessageTag &tag,
                       const Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

    std::vector<T> vecbuf;
    // Use Parallel::receive() so we get its specialization(s)
    Parallel::receive(src_processor_id, vecbuf, type, req, tag, comm);
    buf.clear();
    buf.insert(vecbuf.begin(), vecbuf.end());

    STOP_LOG("receive()", "Parallel");
  }



  inline status wait (request &r)
  {
    START_LOG("wait()", "Parallel");

    status status;

    MPI_Wait (&r, &status);

    STOP_LOG("wait()", "Parallel");

    return status;
  }

  

  inline void wait (std::vector<request> &r)
  {
    START_LOG("wait()", "Parallel");
    
    MPI_Waitall (r.size(), r.empty() ? NULL : &r[0], MPI_STATUSES_IGNORE);

    STOP_LOG("wait()", "Parallel");
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<T1> &send,
			   const DataType &type1,
			   const unsigned int source_processor_id,
			   std::vector<T2> &recv,
			   const DataType &type2,
                           const Communicator &comm = Communicator_World)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }

    Parallel::Request request;
    
    Parallel::nonblocking_send (dest_processor_id,
				send,
				type1,
				request,
				MessageTag(321),
                                comm);
    
    Parallel::receive (source_processor_id,
		       recv,
		       type2,
		       MessageTag(321),
                       comm);
    
    Parallel::wait (request);
    
    STOP_LOG("send_receive()", "Parallel");
  }
  

  
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
			   T1 &send,
			   const unsigned int source_processor_id,
			   T2 &recv,
                           const Communicator &comm)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }

    MPI_Sendrecv(&send, 1, StandardType<T1>(&send),
		 dest_processor_id, 0,
		 &recv, 1, StandardType<T2>(&recv),
		 source_processor_id, 0,
		 comm.get(),
		 MPI_STATUS_IGNORE);

    STOP_LOG("send_receive()", "Parallel");
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<T1> &send,
			   const unsigned int source_processor_id,
			   std::vector<T2> &recv,
                           const Communicator &comm = Communicator_World)
  {
    // Call the user-defined type version with automatic 
    // type conversion based on template argument:    
    send_receive (dest_processor_id,
		  send,
		  StandardType<T1>(send.empty() ? NULL : &send[0]),
		  source_processor_id,
		  recv,
		  StandardType<T2>(recv.empty() ? NULL : &recv[0]),
                  comm);
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T1, typename T2>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<std::vector<T1> > &send,
			   const unsigned int source_processor_id,
			   std::vector<std::vector<T2> > &recv,
                           const Communicator &comm = Communicator_World)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }
    
    // temporary buffers - these will be sized in bytes
    // and manipulated with MPI_Pack and friends
    std::vector<char> sendbuf, recvbuf;

    // figure out how many bytes we need to pack all the data
    int packedsize=0, sendsize=0;

    // The outer buffer size
    MPI_Pack_size (1,           
		   StandardType<unsigned int>(),
		   comm.get(),
		   &packedsize);
    sendsize += packedsize;

    for (unsigned int i=0; i<send.size(); i++)
      {
	// The size of the ith inner buffer
	MPI_Pack_size (1,           
		       StandardType<unsigned int>(),
		       comm.get(),
		       &packedsize);
	sendsize += packedsize;

	// The data for each inner buffer
	MPI_Pack_size (send[i].size(),
		       StandardType<T1>(send[i].empty() ? NULL : &send[i][0]),
		       comm.get(),
		       &packedsize);
	sendsize += packedsize;
      }
    
    libmesh_assert (sendsize /* should at least be 1! */);
    sendbuf.resize (sendsize);

    // Pack the send buffer
    int pos=0;

    // ... the size of the outer buffer
    sendsize = send.size();
    MPI_Pack (&sendsize, 1, StandardType<unsigned int>(),
	      &sendbuf[0], sendbuf.size(), &pos,
	      comm.get());
    
    for (unsigned int i=0; i<send.size(); i++)
      {
	// ... the size of the ith inner buffer
	sendsize = send[i].size();
	MPI_Pack (&sendsize, 1, StandardType<unsigned int>(),
		  &sendbuf[0], sendbuf.size(), &pos,
		  comm.get());
	
	// ... the contents of the ith inner buffer
	if (!send[i].empty())
	  MPI_Pack (&send[i][0], send[i].size(),
		    StandardType<T1>(&send[i][0]),
		    &sendbuf[0], sendbuf.size(), &pos,
		    comm.get());
      }

    libmesh_assert (static_cast<unsigned int>(pos) == sendbuf.size());

    Parallel::Request request;

    Parallel::MessageTag tag = comm.get_unique_tag(123);

    Parallel::nonblocking_send (dest_processor_id,
				sendbuf,
				MPI_PACKED,
				request,
				tag,
                                comm);

    Parallel::receive (source_processor_id,
		       recvbuf,
		       MPI_PACKED,
		       tag,
                       comm);

    // Unpack the received buffer
    libmesh_assert (!recvbuf.empty());
    pos=0;
    MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		&sendsize, 1, StandardType<unsigned int>(),
		comm.get());
    
    // ... size the outer buffer
    recv.resize (sendsize);

    for (unsigned int i=0; i<recv.size(); i++)
      {
	MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		    &sendsize, 1, StandardType<unsigned int>(),
		    comm.get());
    
	// ... size the inner buffer
	recv[i].resize (sendsize);

	// ... unpack the inner buffer if it is not empty
	if (!recv[i].empty())
	  MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		      &recv[i][0], recv[i].size(),
                      StandardType<T2>(&recv[i][0]),
		      comm.get());
      }

    Parallel::wait (request);

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     const Communicator &comm)
  {
    libmesh_assert(root_id < comm.size());

    if (comm.rank() == root_id)
      recv.resize(comm.size());
    
    if (comm.size() > 1)
      {
        START_LOG("gather()", "Parallel");

        StandardType<T> send_type(&send);

	MPI_Gather(&send,
		   1,
		   send_type,
		   recv.empty() ? NULL : &recv[0],
		   1,
		   send_type,
		   root_id,
		   comm.get());

        STOP_LOG("gather()", "Parallel");
      }
    else
      recv[0] = send;
  }



  /**
   * This function provides a convenient method
   * for combining vectors from each processor into one
   * contiguous chunk on one processor.  This handles the 
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
		     std::vector<T> &r,
                     const Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

    libmesh_assert(root_id < comm.size());
    
    std::vector<int>
      sendlengths  (comm.size(), 0),
      displacements(comm.size(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths, comm);
    
    START_LOG("gather()", "Parallel");

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0; 
    for (unsigned int i=0; i != comm.size(); ++i)
      {
	displacements[i] = globalsize;
	globalsize += sendlengths[i];
      }

    // Check for quick return
    if (globalsize == 0)
      {
	STOP_LOG("gather()", "Parallel");
	return;
      }

    // copy the input buffer
    std::vector<T> r_src(r);

    // now resize it to hold the global data
    // on the receiving processor
    if (root_id == comm.rank())
      r.resize(globalsize);

    // and get the data from the remote processors
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, StandardType<T>(),
		   r.empty() ? NULL : &r[0], &sendlengths[0],
		   &displacements[0], StandardType<T>(),
		   root_id,
		   comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        const Communicator &comm)
  {
    START_LOG ("allgather()","Parallel");
    
    libmesh_assert(comm.size());
    recv.resize(comm.size());
    
    unsigned int comm_size = comm.size();
    if (comm_size > 1)
      {
        StandardType<T> send_type(&send);

	MPI_Allgather (&send,
		       1,
		       send_type,
		       &recv[0],
		       1, 
		       send_type,
		       comm.get());
      }
    else if (comm_size > 0)
      recv[0] = send;

    STOP_LOG ("allgather()","Parallel");
  }



  /**
   * This function provides a convenient method
   * for combining vectors from each processor into one
   * contiguous chunk.  This handles the case where the
   * lengths of the vectors may vary.  Specifically, this
   * function transforms this:
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
			const bool identical_buffer_sizes,
                        const Communicator &comm)
  {
    if (comm.size() < 2)
      return;

    START_LOG("allgather()", "Parallel");

    if (identical_buffer_sizes)
      {
        if (r.empty())
          return;

        libmesh_assert(verify(r.size(), comm));
    
	std::vector<T> r_src(r.size()*comm.size());
	r_src.swap(r);
        StandardType<T> send_type(&r_src[0]);

	MPI_Allgather (&r_src[0],
		       r_src.size(),
		       send_type,
		       &r[0],
		       r_src.size(), 
		       send_type,
		       comm.get());
        libmesh_assert(verify(r));
	STOP_LOG("allgather()", "Parallel");
	return;
      }

    std::vector<int>
      sendlengths  (comm.size(), 0),
      displacements(comm.size(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths, comm);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0; 
    for (unsigned int i=0; i != comm.size(); ++i)
      {
	displacements[i] = globalsize;
	globalsize += sendlengths[i];
      }

    // Check for quick return
    if (globalsize == 0)
      {
	STOP_LOG("allgather()", "Parallel");
	return;
      }

    // copy the input buffer
    std::vector<T> r_src(globalsize);
    r_src.swap(r);

    StandardType<T> send_type(&r[0]);

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, send_type,
		      &r[0], &sendlengths[0],
		      &displacements[0], send_type, comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("allgather()", "Parallel");
  }


  /**
   * Replaces the input buffer with the result of MPI_Alltoall.
   * The vector size must be of te form N*n_procs, where N is 
   * the number of elements to be sent/received from each 
   * processor.
   */
  template <typename T>
  inline void alltoall(std::vector<T> &buf,
                       const Communicator &comm)
  {
    if (comm.size() < 2 || buf.empty())
      return;

    START_LOG("alltoall()", "Parallel");

    // the per-processor size.  this is the same for all
    // processors using MPI_Alltoall, could be variable
    // using MPI_Alltoallv
    const unsigned int size_per_proc = 
      buf.size()/comm.size();

    libmesh_assert (buf.size()%comm.size() == 0);

    libmesh_assert(verify(size_per_proc, comm));
    
    std::vector<T> tmp(buf);

    StandardType<T> send_type(&tmp[0]);
    
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Alltoall (&tmp[0],
		    size_per_proc,
		    send_type,
		    &buf[0],
		    size_per_proc,
		    send_type,
		    comm.get());
    libmesh_assert (ierr == MPI_SUCCESS);
    
    STOP_LOG("alltoall()", "Parallel");
  }



  template <typename T>
  inline void broadcast (T &data, const unsigned int root_id,
                         const Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

    libmesh_assert(root_id < comm.size());
    
    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Bcast (&data, 1, StandardType<T>(&data), root_id, comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <>
  inline void broadcast (std::string &data, const unsigned int root_id,
                         const Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

    libmesh_assert(root_id < comm.size());
    
    START_LOG("broadcast()", "Parallel");

    unsigned int data_size = data.size();
    Parallel::broadcast(data_size, root_id, comm);
    
    std::vector<char> data_c(data_size);
    std::string orig(data);

    if (comm.rank() == root_id)
      for(unsigned int i=0; i<data.size(); i++)
	data_c[i] = data[i];

    Parallel::broadcast (data_c, root_id, comm);
    
    data.clear(); data.reserve(data_c.size());
    for(unsigned int i=0; i<data_c.size(); i++)
      data.push_back(data_c[i]);
    
    if (comm.rank() == root_id)
      libmesh_assert(data == orig);

    STOP_LOG("broadcast()", "Parallel");
  }



  template <typename T>
  inline void broadcast (std::vector<T> &data, const unsigned int root_id,
                         const Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

    libmesh_assert(root_id < comm.size());
    
    START_LOG("broadcast()", "Parallel");

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    T *data_ptr = data.empty() ? NULL : &data[0];

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Bcast (data_ptr, data.size(), StandardType<T>(data_ptr),
		 root_id, comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T>
  inline void broadcast (std::set<T> &data, const unsigned int root_id,
                         const Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

    libmesh_assert(root_id < comm.size());
    
    START_LOG("broadcast()", "Parallel");

    std::vector<T> vecdata;
    if (comm.rank() == root_id)
      vecdata.assign(data.begin(), data.end());

    unsigned int vecsize = vecdata.size();
    Parallel::broadcast(vecsize, root_id, comm);
    if (comm.rank() != root_id)
      vecdata.resize(vecsize);

    Parallel::broadcast(vecdata, root_id, comm);
    if (comm.rank() != root_id)
    {
      data.clear();
      data.insert(vecdata.begin(), vecdata.end());
    }

    STOP_LOG("broadcast()", "Parallel");
  }

#else // LIBMESH_HAVE_MPI

  template <typename T>
  inline bool verify(const T &, const Communicator&) { return true; }

  template <typename T>
  inline void min(T &, const Communicator&) {}

  template <typename T>
  inline void min(std::vector<T> &, const Communicator&) {}

  template <typename T>
  inline void minloc(T &, unsigned int &min_id, const Communicator&) { min_id = 0; }

  template <typename T>
  inline void minloc(std::vector<T> &r, std::vector<unsigned int> &min_id, const Communicator&)
    { for (unsigned int i=0; i!= r.size(); ++i) min_id[i] = 0; }

  template <typename T>
  inline void max(T &, const Communicator&) {}

  template <typename T>
  inline void max(std::vector<T> &, const Communicator&) {}

  template <typename T>
  inline void maxloc(T &, unsigned int &max_id, const Communicator&) { max_id = 0; }

  template <typename T>
  inline void maxloc(std::vector<T> &r, std::vector<unsigned int> &max_id, const Communicator&)
    { for (unsigned int i=0; i!= r.size(); ++i) max_id[i] = 0; }

  template <typename T>
  inline void sum(T &, const Communicator&) {}

  template <typename T>
  inline void sum(std::vector<T> &, const Communicator&) {}

  template <typename T>
  inline void set_union(std::set<T> &, const unsigned int, const Communicator&) {}

  template <typename T>
  inline void set_union(std::set<T> &, const Communicator&) {}

  template <typename T1, typename T2>
  inline void set_union(std::map<T1,T2> &, const unsigned int, const Communicator &) {}

  template <typename T1, typename T2>
  inline void set_union(std::map<T1,T2> &, const Communicator &) {}

 //-------------------------------------------------------------------
  /**
   * Blocking message probe.  Allows information about a message to be 
   * examined before the message is actually received.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  inline status probe (const int,
		       const MessageTag&,
                       const Communicator&)
  { libmesh_error(); status status; return status; }

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor with user-defined type.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  template <typename T>
  inline void send (const unsigned int,
		    std::vector<T> &,
		    const DataType &,
		    const MessageTag &,
                    const Communicator&)
  { libmesh_error(); }

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor with user-defined type.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  template <typename T>
  inline void send (const unsigned int,
		    std::vector<T> &,
		    const DataType &,
		    Request &,
		    const MessageTag &,
                    const Communicator&)
  { libmesh_error(); }

  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor with user-defined type.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  template <typename T>
  inline Status receive (const int,
		         std::vector<T> &,
		         const DataType &,
		         const MessageTag &,
                         const Communicator&)
  { libmesh_error(); return Status(); }

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor with user-defined type.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  template <typename T>
  inline void receive (const int,
		       std::vector<T> &,
		       const DataType &,
		       Request &,
		       const MessageTag &,
                       const Communicator&)
  { libmesh_error(); }


//   // on one processor a blocking probe can only be used to 
//   // test a nonblocking send, which we don't really support
//   inline status probe (const int,
// 		       const int)
//   { libmesh_error(); status status; return status; }


//   // Blocking sends don't make sense on one processor
//   template <typename T>
//   inline void send (const unsigned int,
// 		    std::vector<T> &,
// 		    const unsigned int) { libmesh_error(); }

//   template <typename T>
//   inline void nonblocking_send (const unsigned int,
// 		                std::vector<T> &,
// 		                Request &,
// 		                const int) {}

//   // Blocking receives don't make sense on one processor
//   template <typename T>
//   inline Status receive (const int,
// 		         std::vector<T> &,
// 		         const int) { libmesh_error(); return Status(); }

//   template <typename T>
//   inline void nonblocking_receive (const int,
// 		                   std::vector<T> &,
// 		                   Request &,
// 		                   const int) {}
  
  inline status wait (request &) { status status; return status; }
  
  inline void wait (std::vector<request> &) {}
  
  template <typename T1, typename T2>
  inline void send_receive (const unsigned int send_tgt,
			    T1 &send,
			    const unsigned int recv_source,
			    T2 &recv,
                            const Communicator&)
  {
    libmesh_assert (send_tgt == recv_source);
    recv = send;
  }

  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     const Communicator&)
  {
    libmesh_assert (!root_id);
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void gather(const unsigned int, std::vector<T> &,
                     const Communicator&) {}

  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        const Communicator&)
  {
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void allgather(std::vector<T> &, const bool,
                        const Communicator&) {}

  template <typename T>
  inline void alltoall(std::vector<T> &, const Communicator&) {}

  template <typename T>
    inline void broadcast (T &, const unsigned int, const Communicator&) {}

  template <typename T>
    inline void broadcast (std::vector<T> &, const unsigned int, const Communicator&) {}

  template <typename T>
    inline void broadcast(std::set<T> &, const unsigned int, const Communicator &) {}

#endif // LIBMESH_HAVE_MPI


}

} // namespace libMesh

#endif // #define __parallel_h__
