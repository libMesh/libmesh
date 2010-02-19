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

// Local includes
#include "libmesh_common.h" // for Real
#include "libmesh_logging.h"


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only
#ifndef NDEBUG
  #define parallel_only() do { libmesh_assert(Parallel::verify(std::string(__FILE__))); libmesh_assert(Parallel::verify(__LINE__)); } while (0)
#else
  #define parallel_only()
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
   * for use with built-in C types
   */
  template <typename T>
  inline data_type datatype();

  /**
   * Default message tag id
   */
  const int any_tag=MPI_ANY_TAG;

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

  template <typename T>
  inline data_type datatype() { return data_type(); }

  const int any_tag=-1;
  const int any_source=0;
#endif // LIBMESH_HAVE_MPI



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Datatype.
   */
  class DataType
  {
  public:
    DataType () {}

    DataType (const DataType &other) :
      _datatype(other._datatype)
    {}

    DataType (const data_type &type) :
      _datatype(type)
    {}

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
    
  private:

    data_type _datatype;    
  };



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Comm object.  Allows the size of the group
   * and this process's position in the group to be determined.
   */
  class Communicator
  {
  public:
#ifdef LIBMESH_HAVE_MPI
    Communicator () : _communicator(MPI_COMM_NULL), _rank(0), _size(0) {}
#else
    Communicator () : _rank(0), _size(1) {}
#endif

    Communicator (const communicator &comm) {
      this->assign(comm);
    }

    ~Communicator () {
      this->clear();
    }

    communicator& get() {
      return _communicator;
    }

    void clear() {
#ifdef LIBMESH_HAVE_MPI
      if (_communicator != MPI_COMM_NULL)
        {
          MPI_Comm_free(&_communicator);
          _communicator = MPI_COMM_NULL;
        }
#endif
    }

    void operator= (const communicator &comm) {
      this->clear();
      this->assign(comm);
    }

    unsigned int rank() const {
      return _rank;
    }

    unsigned int size() const {
      return _size;
    }

  private:

    void assign(const communicator &comm) {
#ifdef LIBMESH_HAVE_MPI
      if (comm != MPI_COMM_NULL)
        {
          MPI_Comm_dup(comm, &_communicator);

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
          _communicator = MPI_COMM_NULL;
          _rank = 0;
          _size = 0;
        }
#endif
    }

    communicator _communicator;
    unsigned int _rank, _size;
  };



  //-------------------------------------------------------------------
  /**
   * Encapsulates the MPI_Status struct.  Allows the source and size 
   * of the message to be determined.
   */
  class Status
  {
  public:
    Status () {}
    
    Status (const data_type &type) :
    _datatype(type)           
    {}

    Status (const status &status) :
      _status(status)
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
    Request () 
    {
#ifdef LIBMESH_HAVE_MPI
      _request = MPI_REQUEST_NULL;
#endif
    }

    Request (const request &r) :
      _request(r)
    {}

    Request & operator = (const Request &other) 
    { _request = other._request; return *this; }

    Request & operator = (const request &r)
    { _request = r; return *this; }

    ~Request ()
    {
#ifdef LIBMESH_HAVE_MPI
      // explicitly free this request if not 
      // done so already, otherwise this would
      // be a memory leak!
      if (_request != MPI_REQUEST_NULL)
	MPI_Request_free (&_request);
#endif
    }

    operator const request & () const
    { return _request; }

    operator request & () 
    { return _request; }

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
  inline void barrier (Communicator &comm = Communicator_World)
  {
#ifdef LIBMESH_HAVE_MPI
    MPI_Barrier (comm.get());
#endif
  }    
  
  //-------------------------------------------------------------------
  /**
   * Verify that a local variable has the same value on all processors
   */
  template <typename T>
  inline bool verify(const T &r,
                     Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors
   */
  template <typename T>
  inline void min(T &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors
   */
  template <typename T>
  inline void min(std::vector<T> &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors
   */
  template <typename T>
  inline void max(T &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors
   */
  template <typename T>
  inline void max(std::vector<T> &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the sum of it's values
   * on all processors
   */
  template <typename T>
  inline void sum(T &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the sum of
   * it's values on all processors
   */
  template <typename T>
  inline void sum(std::vector<T> &r,
                  Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and collect their
   * union over all processors, replacing the set on processor 0.
   */
  template <typename T>
  inline void set_union(std::set<T> &data, const unsigned int root_id,
                        Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take set of local variables on each processor, and replace it
   * with their union over all processors.
   */
  template <typename T>
  inline void set_union(std::set<T> &data,
                        Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking message probe.  Allows information about a message to be 
   * examined before the message is actually received.
   */
  inline status probe (const int src_processor_id,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    const int tag=0,
                    Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send vector to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send vector to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
	  tag,
          comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<std::complex<T> > &buf,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
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
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
	  req,
	  tag,
          comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<std::complex<T> > &buf,
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
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
		    const int tag=0,
                    Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-send set to one processor with user-defined type.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-send set to one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
	  tag,
          comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<std::complex<T> > &buf,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
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
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
	  req,
	  tag,
          comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<std::complex<T> > &buf,
		    request &req,
		    const int tag=0,
                    Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
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
		                request &r,
		                const int tag=0,
                                Communicator &comm = Communicator_World)
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
		                request &r,
		                const int tag=0,
                                Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
	  r,
	  tag,
          comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void nonblocking_send (const unsigned int dest_processor_id,
		                std::vector<std::complex<T> > &buf,
		                request &r,
		                const int tag=0,
                                Communicator &comm = Communicator_World)
  {
    send (dest_processor_id,
	  buf,
	  datatype<T>(),
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
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive vector from one processor with user-defined type.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<T> &buf,
		       const DataType &type,
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive vector from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<T> &buf,
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
		    datatype<T>(),
		    tag,
                    comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<std::complex<T> > &buf,
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
		    datatype<T>(),
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
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
	     req,
	     tag,
             comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<std::complex<T> > &buf,
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
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
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Nonblocking-receive set from one processor with user-defined type.
   */
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<T> &buf,
		       const DataType &type,
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Blocking-receive set from one processor where the communication type 
   * is inferred from the template argument.
   */
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::set<T> &buf,
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
		    datatype<T>(),
		    tag,
                    comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::set<std::complex<T> > &buf,
		         const int tag=any_tag,
                         Communicator &comm = Communicator_World)
  {
    return receive (src_processor_id,
		    buf,
		    datatype<T>(),
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
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
	     req,
	     tag);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<std::complex<T> > &buf,
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
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
		                   request &r,
		                   const int tag=any_tag,
                                   Communicator &comm = Communicator_World)
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
		                   request &r,
		                   const int tag=any_tag,
                                   Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
	     r,
	     tag,
             comm);
  }

  // Function overloading for std::complex<>
  template <typename T>
  inline void nonblocking_receive (const int src_processor_id,
		                   std::vector<std::complex<T> > &buf,
		                   request &r,
		                   const int tag=any_tag,
                                   Communicator &comm = Communicator_World)
  {
    receive (src_processor_id,
	     buf,
	     datatype<T>(),
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
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
			   const unsigned int source_processor_id,
                           T &recv,
                           Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Send vector \p send to one processor while simultaneously receiving
   * another vector \p recv from a (potentially different) processor using
   * a user-specified MPI Dataype.
   */
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
			   const unsigned int source_processor_id,
                           T &recv,
			   const DataType &type,
                           Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length comm.size(), and on processor root_id fill in
   * recv[processor_id] = the value of send on processor processor_id
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it on processor root_id 
   * to include values from all processors
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::vector<T> &r,
                     Communicator &comm = Communicator_World);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length comm.size(), and fill in 
   * \p recv[processor_id] = the value of \p send on that processor
   */
  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        Communicator &comm = Communicator_World);


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
                        Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Effectively transposes the input vector across all processors.  
   * The jth entry on processor i is replaced with the ith entry 
   * from processor j.
   */
  template <typename T>
  inline void alltoall(std::vector<T> &r,
                       Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Take a local value and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(T &data, const unsigned int root_id=0,
                        Communicator &comm = Communicator_World);



  //-------------------------------------------------------------------
  /**
   * Take a local vector and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.  The user is responsible
   * for appropriately sizing the input buffer on all processors.
   */
  template <typename T>
  inline void broadcast(std::vector<T> &data, const unsigned int root_id=0,
                        Communicator &comm = Communicator_World);


  //-------------------------------------------------------------------
  /**
   * Take a local set and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(std::set<T> &data, const unsigned int root_id=0,
                        Communicator &comm = Communicator_World);


  //-----------------------------------------------------------------------
  // Parallel members

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
 template<>
 inline data_type datatype<char>() { return MPI_CHAR; }

 template<>
 inline data_type datatype<unsigned char>() { return MPI_UNSIGNED_CHAR; }

  template<>
  inline data_type datatype<short int>() { return MPI_SHORT; }

  template<>
  inline data_type datatype<unsigned short int>() { return MPI_UNSIGNED_SHORT; }

  template<>
  inline data_type datatype<int>() { return MPI_INT; }

  template<>
  inline data_type datatype<unsigned int>() { return MPI_UNSIGNED; }

  template<>
  inline data_type datatype<long>() { return MPI_LONG; }

  template<>
  inline data_type datatype<unsigned long>() { return MPI_UNSIGNED_LONG; }

  template<>
  inline data_type datatype<float>() { return MPI_FLOAT; }

  template<>
  inline data_type datatype<double>() { return MPI_DOUBLE; }

  template<>
  inline data_type datatype<long double>() { return MPI_LONG_DOUBLE; }

  template <typename T>
  inline bool verify(const T &r,
                     Communicator &comm)
  {
    if (comm.size() > 1)
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
                     Communicator &comm)
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
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("min()", "Parallel");
    
	T temp = r;
	MPI_Allreduce (&temp,
		       &r,
		       1,
		       datatype<T>(),
		       MPI_MIN,
		       comm.get());

	STOP_LOG("min()", "Parallel");
      }
  }


  template <>
  inline void min(bool &r,
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("min()", "Parallel");
    
	unsigned int tempsend = r;
	unsigned int temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       datatype<unsigned int>(),
		       MPI_MIN,
		       comm.get());
	r = temp;

	STOP_LOG("min()", "Parallel");
      }
  }


  template <typename T>
  inline void min(std::vector<T> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("min()", "Parallel");
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       datatype<T>(),
		       MPI_MIN,
		       comm.get());

	STOP_LOG("min()", "Parallel");
      }
  }


  template <>
  inline void min(std::vector<bool> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("min()", "Parallel");
    
        std::vector<unsigned int> ruint;
        pack_vector_bool(r, ruint);
	std::vector<unsigned int> temp(ruint.size());
	MPI_Allreduce (&ruint[0],
		       &temp[0],
		       ruint.size(),
		       datatype<unsigned int>(),
		       MPI_BAND,
		       comm.get());
        unpack_vector_bool(temp, r);

	STOP_LOG("min()", "Parallel");
      }
  }


  template <typename T>
  inline void max(T &r,
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("max()", "Parallel");
    
	T temp;
	MPI_Allreduce (&r,
		       &temp,
		       1,
		       datatype<T>(),
		       MPI_MAX,
		       comm.get());
	r = temp;

	STOP_LOG("max()", "Parallel");
      }
  }


  template <>
  inline void max(bool &r,
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("max()", "Parallel");
    
	unsigned int tempsend = r;
	unsigned int temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       datatype<unsigned int>(),
		       MPI_MAX,
		       comm.get());
	r = temp;

	STOP_LOG("max()", "Parallel");
      }
  }


  template <typename T>
  inline void max(std::vector<T> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("max()", "Parallel");
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       datatype<T>(),
		       MPI_MAX,
		       comm.get());

	STOP_LOG("max()", "Parallel");
      }
  }


  template <>
  inline void max(std::vector<bool> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("max()", "Parallel");
    
        std::vector<unsigned int> ruint;
        pack_vector_bool(r, ruint);
	std::vector<unsigned int> temp(ruint.size());
	MPI_Allreduce (&ruint[0],
		       &temp[0],
		       ruint.size(),
		       datatype<unsigned int>(),
		       MPI_BOR,
		       comm.get());
        unpack_vector_bool(temp, r);

	STOP_LOG("max()", "Parallel");
      }
  }


  template <typename T>
  inline void sum(T &r,
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("sum()", "Parallel");
    
	T temp = r;
	MPI_Allreduce (&temp,
		       &r,
		       1,
		       datatype<T>(),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  template <typename T>
  inline void sum(std::vector<T> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("sum()", "Parallel");
    
	std::vector<T> temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size(),
		       datatype<T>(),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }



  template <typename T>
  inline void sum(std::complex<T> &r,
                  Communicator &comm)
  {
    if (comm.size() > 1)
      {
	START_LOG("sum()", "Parallel");
    
	std::complex<T> temp(r);
	MPI_Allreduce (&temp,
		       &r,
		       2,
		       datatype<T>(),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }


  template <typename T>
  inline void sum(std::vector<std::complex<T> > &r,
                  Communicator &comm)
  {
    if (comm.size() > 1 && !r.empty())
      {
	START_LOG("sum()", "Parallel");
    
	std::vector<std::complex<T> > temp(r);
	MPI_Allreduce (&temp[0],
		       &r[0],
		       r.size() * 2,
		       datatype<T>(),
		       MPI_SUM,
		       comm.get());

	STOP_LOG("sum()", "Parallel");
      }
  }



  template <typename T>
  inline void set_union(std::set<T> &data, const unsigned int root_id,
                        Communicator &comm)
  {
    std::vector<T> vecdata(data.begin(), data.end());
    Parallel::gather(root_id, vecdata);
    if (comm.rank() == root_id)
      data.insert(vecdata.begin(), vecdata.end());
  }



  template <typename T>
  inline void set_union(std::set<T> &data,
                        Communicator &comm)
  {
    Parallel::set_union(data, 0, comm);
    broadcast(data);
  }



  inline status probe (const int src_processor_id,
		       const int tag,
                       Communicator &comm)
  {
    START_LOG("probe()", "Parallel");

    status status;

    MPI_Probe (src_processor_id, 
	       tag, 
	       comm.get(), 
	       &status);

    STOP_LOG("probe()", "Parallel");

    return status;
  }  



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    const int tag,
                    Communicator &comm)
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
		tag,
		comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);    
    
    STOP_LOG("send()", "Parallel");
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default arguments
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<std::complex<T> > &buf,
		    const DataType &type,
		    const int tag=0,
                    Communicator &comm=Communicator_World)
  {    
    START_LOG("send()", "Parallel");

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Send (buf.empty() ? NULL : &buf[0],
		buf.size() * 2,
		type,
		dest_processor_id,
		tag,
		comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);    
    
    STOP_LOG("send()", "Parallel");
  }



  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<T> &buf,
		    const DataType &type,
		    request &req,
		    const int tag,
                    Communicator &comm)
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
		 tag,
		 comm.get(),
		 &req);    
    libmesh_assert (ierr == MPI_SUCCESS);    

    STOP_LOG("send()", "Parallel");
  }


  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default arguments
  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::vector<std::complex<T> > &buf,
		    const DataType &type,
		    request &req,
		    const int tag=0,
                    Communicator &comm=Communicator_World)
  {
    START_LOG("send()", "Parallel");

#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =	  
#endif
      MPI_ISend (buf.empty() ? NULL : &buf[0],
		 buf.size() * 2,
		 type,
		 dest_processor_id,
		 tag,
		 comm.get(),
                 &req);    
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("send()", "Parallel");
  }



  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<T> &buf,
		         const DataType &type,
		         const int tag,
                         Communicator &comm)
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
		tag,
		comm.get(),
		status);
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("receive()", "Parallel");

    return status;
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default arguments
  template <typename T>
  inline Status receive (const int src_processor_id,
		         std::vector<std::complex<T> > &buf,
		         const DataType &type,
		         const int tag=any_tag,
                         Communicator &comm=Communicator_World)
  {
    START_LOG("receive()", "Parallel");

    // Get the status of the message, explicitly provide the
    // datatype so we can later query the size
    Status status(Parallel::probe(src_processor_id, tag, comm), datatype<T>());

    libmesh_assert(!(status.size()%2));
    buf.resize(status.size()/2);
    
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =	  
#endif
      MPI_Recv (buf.empty() ? NULL : &buf[0],
		buf.size() * 2,
		type,
		src_processor_id,
		tag,
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
		       request &req,
		       const int tag,
                       Communicator &comm)
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
		 tag,
		 comm.get(),
		 &req);    
    libmesh_assert (ierr == MPI_SUCCESS);    

    STOP_LOG("receive()", "Parallel");
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default arguments
  template <typename T>
  inline void receive (const int src_processor_id,
		       std::vector<std::complex<T> > &buf,
		       const DataType &type,
		       request &req,
		       const int tag=any_tag,
                       Communicator &comm=Communicator_World)
  {
    START_LOG("receive()", "Parallel");

    const int ierr =	  
      MPI_Irecv (buf.empty() ? NULL : &buf[0],
		buf.size() * 2,
		type,
		src_processor_id,
		tag,
		comm.get(),
		&req);

    libmesh_assert (ierr == MPI_SUCCESS);
    
    STOP_LOG("receive()", "Parallel");

    return;
  }

  

  template <typename T>
  inline void send (const unsigned int dest_processor_id,
		    std::set<T> &buf,
		    const DataType &type,
		    const int tag,
                    Communicator &comm)
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
		    request &req,
		    const int tag,
                    Communicator &comm)
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
		         const int tag,
                         Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

    std::vector<T> vecbuf;
    // Use Parallel::receive() so we get its specialization(s)
    Status status = Parallel::receive(src_processor_id, vecbuf, type, tag, comm);
    buf.assign(vecbuf.begin(), vecbuf.end());

    STOP_LOG("receive()", "Parallel");

    return status;
  }



  template <typename T>
  inline void receive (const int src_processor_id,
		       std::set<T> &buf,
		       const DataType &type,
		       request &req,
		       const int tag,
                       Communicator &comm)
  {
    START_LOG("receive()", "Parallel");

    std::vector<T> vecbuf;
    // Use Parallel::receive() so we get its specialization(s)
    Parallel::receive(src_processor_id, vecbuf, type, req, tag, comm);
    buf.assign(vecbuf.begin(), vecbuf.end());

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
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<T> &send,
			   const unsigned int source_processor_id,
			   std::vector<T> &recv,
			   const DataType &type,
                           Communicator &comm = Communicator_World)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }

    Parallel::request request;
    
    Parallel::nonblocking_send (dest_processor_id,
				send,
				type,
				request,
				/* tag = */ 321,
                                comm);
    
    Parallel::receive (source_processor_id,
		       recv,
		       type,
		       /* tag = */ 321,
                       comm);
    
    Parallel::wait (request);
    
    STOP_LOG("send_receive()", "Parallel");
  }
  

  
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
			   T &send,
			   const unsigned int source_processor_id,
			   T &recv,
                           Communicator &comm)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }

    MPI_Sendrecv(&send, 1, datatype<T>(),
		 dest_processor_id, 0,
		 &recv, 1, datatype<T>(),
		 source_processor_id, 0,
		 comm.get(),
		 MPI_STATUS_IGNORE);

    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::complex<T> &send,
			   const unsigned int source_processor_id,
			   std::complex<T> &recv,
                           Communicator &comm)
  {
    START_LOG("send_receive()", "Parallel");

    if (dest_processor_id   == comm.rank() &&
	source_processor_id == comm.rank())
      {
	recv = send;
	STOP_LOG("send_receive()", "Parallel");
	return;
      }

    MPI_Sendrecv(&send, 2, datatype<T>(),
		 dest_processor_id, 0,
		 &recv, 2, datatype<T>(),
		 source_processor_id, 0,
		 comm.get(),
		 MPI_STATUS_IGNORE);

    STOP_LOG("send_receive()", "Parallel");
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<T> &send,
			   const unsigned int source_processor_id,
			   std::vector<T> &recv,
                           Communicator &comm = Communicator_World)
  {
    // Call the user-defined type version with automatic 
    // type conversion based on template argument:    
    send_receive (dest_processor_id,
		  send,
		  source_processor_id,
		  recv,
		  datatype<T>(),
                  comm);
  }



  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
			   std::vector<std::vector<T> > &send,
			   const unsigned int source_processor_id,
			   std::vector<std::vector<T> > &recv,
                           Communicator &comm = Communicator_World)
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
		   datatype<unsigned int>(),
		   comm.get(),
		   &packedsize);
    sendsize += packedsize;

    for (unsigned int i=0; i<send.size(); i++)
      {
	// The size of the ith inner buffer
	MPI_Pack_size (1,           
		       datatype<unsigned int>(),
		       comm.get(),
		       &packedsize);
	sendsize += packedsize;

	// The data for each inner buffer
	MPI_Pack_size (send[i].size(),
		       datatype<T>(),
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
    MPI_Pack (&sendsize, 1, datatype<unsigned int>(),
	      &sendbuf[0], sendbuf.size(), &pos,
	      comm.get());
    
    for (unsigned int i=0; i<send.size(); i++)
      {
	// ... the size of the ith inner buffer
	sendsize = send[i].size();
	MPI_Pack (&sendsize, 1, datatype<unsigned int>(),
		  &sendbuf[0], sendbuf.size(), &pos,
		  comm.get());
	
	// ... the contents of the ith inner buffer
	if (!send[i].empty())
	  MPI_Pack (&send[i][0], send[i].size(), datatype<T>(),
		    &sendbuf[0], sendbuf.size(), &pos,
		    comm.get());
      }

    libmesh_assert (static_cast<unsigned int>(pos) == sendbuf.size());

    Parallel::request request;

    Parallel::nonblocking_send (dest_processor_id,
				sendbuf,
				MPI_PACKED,
				request,
				/* tag = */ 123,
                                comm);

    Parallel::receive (source_processor_id,
		       recvbuf,
		       MPI_PACKED,
		       /* tag = */ 123,
                       comm);

    // Unpack the received buffer
    libmesh_assert (!recvbuf.empty());
    pos=0;
    MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		&sendsize, 1, datatype<unsigned int>(),
		comm.get());
    
    // ... size the outer buffer
    recv.resize (sendsize);

    for (unsigned int i=0; i<recv.size(); i++)
      {
	MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		    &sendsize, 1, datatype<unsigned int>(),
		    comm.get());
    
	// ... size the inner buffer
	recv[i].resize (sendsize);

	// ... unpack the inner buffer if it is not empty
	if (!recv[i].empty())
	  MPI_Unpack (&recvbuf[0], recvbuf.size(), &pos,
		      &recv[i][0], recv[i].size(), datatype<T>(),
		      comm.get());
      }

    Parallel::wait (request);

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     Communicator &comm)
  {
    libmesh_assert(root_id < comm.size());

    if (comm.rank() == root_id)
      recv.resize(comm.size());
    
    if (comm.size() > 1)
      {
        START_LOG("gather()", "Parallel");

	MPI_Gather(&send,
		   1,
		   datatype<T>(),
		   recv.empty() ? NULL : &recv[0],
		   1,
		   datatype<T>(),
		   root_id,
		   comm.get());

        STOP_LOG("gather()", "Parallel");
      }
    else
      recv[0] = send;
  }



  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::complex<T> send,
		     std::vector<std::complex<T> > &recv,
                     Communicator &comm)
  {
    libmesh_assert(root_id < comm.size());

    if (comm.rank() == root_id)
      recv.resize(comm.size());
    
    if (comm.size() > 1)
      {
        START_LOG("gather()", "Parallel");

	MPI_Gather(&send,
		   2,
		   datatype<T>(),
		   recv.empty() ? NULL : &recv[0],
		   2,
		   datatype<T>(),
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
                     Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank()==root_id);
	return;
      }
    
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
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		   r.empty() ? NULL :  &r[0], &sendlengths[0],
		   &displacements[0], datatype<T>(),
		   root_id,
		   comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::vector<std::complex<T> > &r,
                     Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank()==root_id);
	return;
      }
    
    std::vector<int>
      sendlengths  (comm.size(), 0),
      displacements(comm.size(), 0);

    const int mysize = r.size() * 2;
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

    // Make temporary buffers for the input/output data
    std::vector<std::complex<T> > r_src(r);

    // now resize r to hold the global data
    // on the receiving processor
    if (root_id == comm.rank())
      r.resize(globalsize/2);

    // and get the data from the remote processors
    const int ierr =
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		   r.empty() ? NULL : &r[0], &sendlengths[0],
		   &displacements[0], datatype<T>(),
		   root_id, comm.get());
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        Communicator &comm)
  {
    START_LOG ("allgather()","Parallel");
    
    recv.resize(comm.size());
    
    if (comm.size() > 1)
      {
	MPI_Allgather (&send,
		       1,
		       datatype<T>(),
		       &recv[0],
		       1, 
		       datatype<T>(),
		       comm.get());
      }
    else
      recv[0] = send;

    STOP_LOG ("allgather()","Parallel");
  }



  template <typename T>
  inline void allgather(std::complex<T> send,
			std::vector<std::complex<T> > &recv,
                        Communicator &comm)
  {
    START_LOG ("allgather()","Parallel");

    recv.resize(comm.size());

    if (comm.size() > 1)
      {
	MPI_Allgather (&send,
		       2,
		       datatype<T>(),
		       &recv[0],
		       2, 
		       datatype<T>(),
		       comm.get());
      }
    else
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
                        Communicator &comm)
  {
    if (comm.size() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    if (identical_buffer_sizes)
      {
	std::vector<T> r_src(r.size()*comm.size());
	r_src.swap(r);
	MPI_Allgather (r_src.empty() ? NULL : &r_src[0],
		       r_src.size(),
		       datatype<T>(),
		       r.empty() ? NULL : &r[0],
		       r_src.size(), 
		       datatype<T>(),
		       comm.get());
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

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		      r.empty()     ? NULL : &r[0],     &sendlengths[0],
		      &displacements[0], datatype<T>(), comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("allgather()", "Parallel");
  }


  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default arguments
  template <typename T>
  inline void allgather(std::vector<std::complex<T> > &r,
			const bool identical_buffer_sizes = false,
                        Communicator &comm = Communicator_World)
  {
    if (comm.size() == 1)
      return;

    START_LOG("allgather()", "Parallel");

    if (identical_buffer_sizes)
      {
	std::vector<std::complex<T> > r_src(r.size()*comm.size());
	r_src.swap(r);
	MPI_Allgather (r_src.empty() ? NULL : &r_src[0],
		       2*r_src.size(),
		       datatype<T>(),
		       r.empty() ? NULL : &r[0],
		       2*r_src.size(), 
		       datatype<T>(),
		       comm.get());
	STOP_LOG("allgather()", "Parallel");
	return;
      }

    std::vector<int>
      sendlengths  (comm.size(), 0),
      displacements(comm.size(), 0);

    const int mysize = r.size() * 2;
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
    std::vector<std::complex<T> > r_src(globalsize);
    r_src.swap(r);

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		      r.empty()     ? NULL : &r[0],     &sendlengths[0],
		      &displacements[0], datatype<T>(),
		      comm.get());
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
                       Communicator &comm)
  {
    if (comm.size() == 1)
      return;

    START_LOG("alltoall()", "Parallel");

    // the per-processor size.  this is the same for all
    // processors using MPI_Alltoall, could be variable
    // using MPI_Alltoallv
    const unsigned int size_per_proc = 
      buf.size()/comm.size();

    libmesh_assert (buf.size()%comm.size() == 0);

    std::vector<T> tmp(buf);
    
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Alltoall (tmp.empty() ? NULL : &tmp[0],
		    size_per_proc,
		    datatype<T>(),
		    buf.empty() ? NULL : &buf[0],
		    size_per_proc,
		    datatype<T>(),
		    comm.get());
    libmesh_assert (ierr == MPI_SUCCESS);
    
    STOP_LOG("alltoall()", "Parallel");
  }



  template <typename T>
  inline void broadcast (T &data, const unsigned int root_id,
                         Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank() == root_id);
	return;
      }
    
    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Bcast (&data, 1, datatype<T>(), root_id, comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T>
  inline void broadcast (std::complex<T> &data, const unsigned int root_id=0,
                         Communicator &comm = Communicator_World)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank() == root_id);
	return;
      }
    
    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
    const int ierr =
      MPI_Bcast (&data, 2, datatype<T>(), root_id, comm.get());
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }



  template <>
  inline void broadcast (std::string &data, const unsigned int root_id,
                         Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (!comm.rank());
	libmesh_assert (!root_id);
	return;
      }

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
                         Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank() == root_id);
	return;
      }

    START_LOG("broadcast()", "Parallel");

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
#ifndef NDEBUG
    // Only catch the return value when asserts are active.
    const int ierr =
#endif
      MPI_Bcast (data.empty() ? NULL : &data[0], data.size(), datatype<T>(),
		 root_id, comm.get());

    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  // This is both a declaration and definition for a new overloaded
  // function template, so we have to re-specify the default argument
  template <typename T>
  inline void broadcast (std::vector<std::complex<T> > &data,
			 const unsigned int root_id=0,
                         Communicator &comm = Communicator_World)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank() == root_id);
	return;
      }

    START_LOG("broadcast()", "Parallel");

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Bcast (data.empty() ? NULL : &data[0], data.size() * 2, datatype<T>(),
		 root_id, comm.get());
    libmesh_assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


  template <typename T>
  inline void broadcast (std::set<T> &data, const unsigned int root_id,
                         Communicator &comm)
  {
    if (comm.size() == 1)
      {
	libmesh_assert (comm.rank() == root_id);
	return;
      }

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
      data.assign(vecdata.begin(), vecdata.end());

    STOP_LOG("broadcast()", "Parallel");
  }

#else // LIBMESH_HAVE_MPI

  template <typename T>
  inline bool verify(const T &, Communicator&) { return true; }

  template <typename T>
  inline void min(T &, Communicator&) {}

  template <typename T>
  inline void min(std::vector<T> &, Communicator&) {}

  template <typename T>
  inline void max(T &, Communicator&) {}

  template <typename T>
  inline void max(std::vector<T> &, Communicator&) {}

  template <typename T>
  inline void sum(T &, Communicator&) {}

  template <typename T>
  inline void sum(std::vector<T> &, Communicator&) {}

  template <typename T>
  inline void set_union(std::set<T> &, const unsigned int, Communicator&) {}

  template <typename T>
  inline void set_union(std::set<T> &, Communicator&) {}

 //-------------------------------------------------------------------
  /**
   * Blocking message probe.  Allows information about a message to be 
   * examined before the message is actually received.
   *
   * we do not currently support this operation on one processor without MPI.
   */
  inline status probe (const int,
		       const int,
                       Communicator&)
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
		    const int,
                    Communicator&)
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
		    request &,
		    const int,
                    Communicator&)
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
		         const int,
                         Communicator&)
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
		       request &,
		       const int,
                       Communicator&)
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
// 		                request &,
// 		                const int) {}

//   // Blocking receives don't make sense on one processor
//   template <typename T>
//   inline Status receive (const int,
// 		         std::vector<T> &,
// 		         const int) { libmesh_error(); return Status(); }

//   template <typename T>
//   inline void nonblocking_receive (const int,
// 		                   std::vector<T> &,
// 		                   request &,
// 		                   const int) {}
  
  inline status wait (request &) { status status; return status; }
  
  inline void wait (std::vector<request> &) {}
  
  template <typename T>
  inline void send_receive (const unsigned int send_tgt,
			    T &send,
			    const unsigned int recv_source,
			    T &recv,
                            Communicator&)
  {
    libmesh_assert (send_tgt == recv_source);
    recv = send;
  }

  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv,
                     Communicator&)
  {
    libmesh_assert (!root_id);
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void gather(const unsigned int, std::vector<T> &,
                     Communicator&) {}

  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv,
                        Communicator&)
  {
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void allgather(std::vector<T> &, const bool,
                        Communicator&) {}

  template <typename T>
  inline void alltoall(std::vector<T> &, Communicator&) {}

  template <typename T>
    inline void broadcast (T &, const unsigned int, Communicator&) {}

  template <typename T>
    inline void broadcast (std::vector<T> &, const unsigned int, Communicator&) {}

#endif // LIBMESH_HAVE_MPI


}

#endif // #define __parallel_h__
