// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <string>
#include <vector>

// Local includes
#include "libmesh_common.h" // for Real
#include "libmesh_logging.h"


// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#undef parallel_only
#ifdef DEBUG
  #define parallel_only() { assert(Parallel::verify(std::string(__FILE__))); assert(Parallel::verify(__LINE__)); }
#else
  #define parallel_only() { }
#endif

// ------------------------------------------------------------
// The Parallel namespace is for wrapper functions
// for common general parallel synchronization tasks
//
// For MPI 1.1 compatibility, temporary buffers are used
// instead of MPI 2's MPI_IN_PLACE

namespace Parallel
{
#ifdef HAVE_MPI
  //-------------------------------------------------------------------
  /**
   * Templated function to return the appropriate MPI datatype
   * for use with built-in C types
   */
  template <typename T>
  inline MPI_Datatype datatype();
#endif // HAVE_MPI

  //-------------------------------------------------------------------
  /**
   * Verify that a local variable has the same value on all processors
   */
  template <typename T>
  inline bool verify(const T &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the minimum of it's values
   * on all processors
   */
  template <typename T>
  inline void min(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the minimum
   * of it's values on all processors
   */
  template <typename T>
  inline void min(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the maximum of it's values
   * on all processors
   */
  template <typename T>
  inline void max(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the maximum
   * of it's values on all processors
   */
  template <typename T>
  inline void max(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Take a local variable and replace it with the sum of it's values
   * on all processors
   */
  template <typename T>
  inline void sum(T &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and replace each entry with the sum of
   * it's values on all processors
   */
  template <typename T>
  inline void sum(std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Send vector send to one processor while simultaneously receiving
   * another vector recv from a (potentially different) processor.
   */
  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
			   const unsigned int source_processor_id,
                           T &recv);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length n_processors, and on processor root_id fill in
   * recv[processor_id] = the value of send on processor processor_id
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv);

  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it on processor root_id 
   * to include values from all processors
   */
  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::vector<T> &r);

  //-------------------------------------------------------------------
  /**
   * Take a vector of length n_processors, and fill in recv[processor_id] = the
   * value of send on that processor
   */
  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv);


  //-------------------------------------------------------------------
  /**
   * Take a vector of local variables and expand it to include 
   * values from all processors
   */
  template <typename T>
  inline void allgather(std::vector<T> &r);



  //-------------------------------------------------------------------
  /**
   * Take a local value and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.
   */
  template <typename T>
  inline void broadcast(T &data, const unsigned int root_id=0);



  //-------------------------------------------------------------------
  /**
   * Take a local vector and broadcast it to all processors.
   * Optionally takes the \p root_id processor, which specifies
   * the processor intiating the broadcast.  The user is responsible
   * for appropriately sizing the input buffer on all processors.
   */
  template <typename T>
  inline void broadcast(std::vector<T> &data, const unsigned int root_id=0);



  //-----------------------------------------------------------------------
  // Parallel members

#ifdef HAVE_MPI
  template<>
  inline MPI_Datatype datatype<char>() { return MPI_CHAR; }

  template<>
  inline MPI_Datatype datatype<unsigned char>() { return MPI_UNSIGNED_CHAR; }

  template<>
  inline MPI_Datatype datatype<short int>() { return MPI_SHORT; }

  template<>
  inline MPI_Datatype datatype<unsigned short int>() { return MPI_UNSIGNED_SHORT; }

  template<>
  inline MPI_Datatype datatype<int>() { return MPI_INT; }

  template<>
  inline MPI_Datatype datatype<unsigned int>() { return MPI_UNSIGNED; }

  template<>
  inline MPI_Datatype datatype<long>() { return MPI_LONG; }

  template<>
  inline MPI_Datatype datatype<unsigned long>() { return MPI_UNSIGNED_LONG; }

  template<>
  inline MPI_Datatype datatype<float>() { return MPI_FLOAT; }

  template<>
  inline MPI_Datatype datatype<double>() { return MPI_DOUBLE; }

  template<>
  inline MPI_Datatype datatype<long double>() { return MPI_LONG_DOUBLE; }


  template <typename T>
  inline bool verify(const T &r)
  {
    if (libMesh::n_processors() > 1)
      {
	T tempmin = r, tempmax = r;
	Parallel::min(tempmin);
	Parallel::max(tempmax);
	bool verified = (r == tempmin) &&
	                (r == tempmax);
	Parallel::min(verified);
      }
    return true;
  }


  template <>
  inline bool verify(const std::string & r)
  {
    if (libMesh::n_processors() > 1)
      {
	std::vector<char> temp(r.size());
	for (unsigned int i=0; i != r.size(); ++i)
	  temp.push_back(r[i]);
	return Parallel::verify(temp);
      }
    return true;
  }


  template <typename T>
  inline void min(T &r)
  {
    if (libMesh::n_processors() > 1)
      {
	T temp;
	MPI_Allreduce (&r,
		       &temp,
		       1,
		       datatype<T>(),
		       MPI_MIN,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <>
  inline void min(bool &r)
  {
    if (libMesh::n_processors() > 1)
      {
	unsigned char tempsend = r;
	unsigned char temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       datatype<unsigned char>(),
		       MPI_MIN,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void min(std::vector<T> &r)
  {
    if (libMesh::n_processors() > 1)
      {
	std::vector<T> temp(r.size());
	MPI_Allreduce (&r[0],
		       &temp[0],
		       r.size(),
		       datatype<T>(),
		       MPI_MIN,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void max(T &r)
  {
    if (libMesh::n_processors() > 1)
      {
	T temp;
	MPI_Allreduce (&r,
		       &temp,
		       1,
		       datatype<T>(),
		       MPI_MAX,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <>
  inline void max(bool &r)
  {
    if (libMesh::n_processors() > 1)
      {
	unsigned char tempsend = r;
	unsigned char temp;
	MPI_Allreduce (&tempsend,
		       &temp,
		       1,
		       datatype<unsigned char>(),
		       MPI_MAX,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void max(std::vector<T> &r)
  {
    if (libMesh::n_processors() > 1)
      {
	std::vector<T> temp(r.size());
	MPI_Allreduce (&r[0],
		       &temp[0],
		       r.size(),
		       datatype<T>(),
		       MPI_MAX,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void sum(T &r)
  {
    if (libMesh::n_processors() > 1)
      {
	T temp;
	MPI_Allreduce (&r,
		       &temp,
		       1,
		       datatype<T>(),
		       MPI_SUM,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void sum(std::vector<T> &r)
  {
    if (libMesh::n_processors() > 1 && !r.empty())
      {
	std::vector<T> temp(r.size());
	MPI_Allreduce (&r[0],
		       &temp[0],
		       r.size(),
		       datatype<T>(),
		       MPI_SUM,
		       libMesh::COMM_WORLD);
	r = temp;
      }
  }


  template <typename T>
  inline void sum(std::complex<T> &r)
  {
    if (libMesh::n_processors() > 1)
      {
	T tempinput[2], tempoutput[2];
	tempinput[0] = r.real();
	tempinput[1] = r.imag();
	MPI_Allreduce (&tempinput,
		       &tempoutput,
		       2,
		       datatype<T>(),
		       MPI_SUM,
		       libMesh::COMM_WORLD);
	r.real() = tempoutput[0];
	r.imag() = tempoutput[1];
      }
  }


  template <typename T>
  inline void sum(std::vector<std::complex<T> > &r)
  {
    if (libMesh::n_processors() > 1 && !r.empty())
      {
	std::vector<T> temprealinput(r.size()),
	  tempimaginput(r.size()),
	  temprealoutput(r.size()),
	  tempimagoutput(r.size());
	for (unsigned int i=0; i != r.size(); ++i)
	  {
	    temprealinput[i] = r[i].real();
	    tempimaginput[i] = r[i].imag();
	  }
	MPI_Allreduce (&temprealinput[0],
		       &temprealoutput[0],
		       r.size(),
		       datatype<T>(),
		       MPI_SUM,
		       libMesh::COMM_WORLD);
	MPI_Allreduce (&tempimaginput[0],
		       &tempimagoutput[0],
		       r.size(),
		       datatype<T>(),
		       MPI_SUM,
		       libMesh::COMM_WORLD);
	for (unsigned int i=0; i != r.size(); ++i)
	  {
	    r[i].real() = temprealoutput[i];
	    r[i].imag() = tempimagoutput[i];
	  }
      }
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           T &send,
			   const unsigned int source_processor_id,
                           T &recv)
  {
    START_LOG("send_receive()", "Parallel");

    MPI_Status status;
    MPI_Sendrecv(&send, 1, datatype<T>(),
		 dest_processor_id, 0,
		 &recv, 1, datatype<T>(),
		 source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);

    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<T> &send,
			   const unsigned int source_processor_id,
                           std::vector<T> &recv)
  {
    START_LOG("send_receive()", "Parallel");

    // Trade buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
		 dest_processor_id, 0,
		 &recvsize, 1, datatype<unsigned int>(),
		 source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);

    recv.resize(recvsize);

    MPI_Sendrecv(sendsize ? &send[0] : NULL, sendsize, datatype<T>(),
		 dest_processor_id, 0,
		 recvsize ? &recv[0] : NULL, recvsize, datatype<T>(),
		 source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);
    
    STOP_LOG("send_receive()", "Parallel");
  }



  template <typename T>
  inline void send_receive(const unsigned int dest_processor_id,
                           std::vector<std::vector<T> > &send,
			   const unsigned int source_processor_id,
                           std::vector<std::vector<T> > &recv)
  {
    START_LOG("send_receive()", "Parallel");

    // Trade outer buffer sizes first
    unsigned int sendsize = send.size(), recvsize;
    MPI_Status status;
    MPI_Sendrecv(&sendsize, 1, datatype<unsigned int>(),
		 dest_processor_id, 0,
		 &recvsize, 1, datatype<unsigned int>(),
		 source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);

    recv.resize(recvsize);

    // Trade inner buffer sizes next
    std::vector<unsigned int> sendsizes(sendsize), recvsizes(recvsize);
    unsigned int sendsizesum = 0, recvsizesum = 0;
    for (unsigned int i = 0; i != sendsize; ++i)
      {
        sendsizes[i] = send[i].size();
	sendsizesum += sendsizes[i];
      }

    MPI_Sendrecv(sendsize ? &sendsizes[0] : NULL, sendsize, 
		 datatype<unsigned int>(), dest_processor_id, 0,
		 recvsize ? &recvsizes[0] : NULL, recvsize, 
		 datatype<unsigned int>(), source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);

    for (unsigned int i = 0; i != recvsize; ++i)
      {
	recvsizesum += recvsizes[i];
      }

    // Build temporary buffers third
    // We can't do multiple Sendrecv calls instead because send.size() may
    // differ on different processors
    std::vector<T> senddata(sendsizesum), recvdata(recvsizesum);

    // Fill the temporary send buffer
    typename std::vector<T>::iterator out = senddata.begin();
    for (unsigned int i = 0; i != sendsize; ++i)
      {
	out = std::copy(send[i].begin(), send[i].end(), out);
      }
    assert(out == senddata.end());
    
    MPI_Sendrecv(sendsizesum ? &senddata[0] : NULL, sendsizesum,
		 datatype<T>(), dest_processor_id, 0,
		 recvsizesum ? &recvdata[0] : NULL, recvsizesum,
		 datatype<T>(), source_processor_id, 0,
		 libMesh::COMM_WORLD,
		 &status);

    // Empty the temporary recv buffer
    typename std::vector<T>::iterator in = recvdata.begin();
    for (unsigned int i = 0; i != recvsize; ++i)
      {
	std::copy(in, in + recvsizes[i], recv[i].begin());
	in += recvsizes[i];
      }
    assert(in == recvdata.end());

    STOP_LOG("send_receive()", "Parallel");
  }


  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv)
  {
    assert(root_id < libMesh::n_processors());

    if (libMesh::processor_id() == root_id)
      recv.resize(libMesh::n_processors());
    
    if (libMesh::n_processors() > 1)
      {
	MPI_Gather(&send,
		   1,
		   datatype<T>(),
		   &recv[0],
		   1,
		   datatype<T>(),
		   root_id,
		   libMesh::COMM_WORLD);
      }
    else
      recv[0] = send;
  }



  template <typename T>
  inline void gather(const unsigned int root_id,
		     std::complex<T> send,
		     std::vector<std::complex<T> > &recv)
  {
    assert(root_id < libMesh::n_processors());

    std::vector<T> temprealoutput, tempimagoutput;

    if (libMesh::processor_id() == root_id)
      {
        temprealoutput.resize(libMesh::n_processors());
        tempimagoutput.resize(libMesh::n_processors());
        recv.resize(libMesh::n_processors());
      }
    
    if (libMesh::n_processors() > 1)
      {
	T realinput = send.real(),
	  imaginput = send.imag();

	MPI_Gather(&realinput,
		   1,
		   datatype<T>(),
		   &temprealoutput[0],
		   1,
		   datatype<T>(),
		   root_id,
		   libMesh::COMM_WORLD);

	MPI_Gather(&imaginput,
		   1,
		   datatype<T>(),
		   &tempimagoutput[0],
		   1,
		   datatype<T>(),
		   root_id,
		   libMesh::COMM_WORLD);

	for (unsigned int i=0; i != recv.size(); ++i)
	  {
	    recv[i].real() = temprealoutput[i];
	    recv[i].imag() = tempimagoutput[i];
	  }
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
   * must be called by all processors.
   */
  template <typename T>
  void gather(const unsigned int root_id,
	      std::vector<T> &r)
  {
    if (libMesh::n_processors() == 1)
      return;
    
    START_LOG("gather()", "Parallel");
    
    std::vector<int>
      sendlengths  (libMesh::n_processors(), 0),
      displacements(libMesh::n_processors(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0; 
    for (unsigned int i=0; i != libMesh::n_processors(); ++i)
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
    if (root_id == libMesh::processor_id())
      r.resize(globalsize);

    // and get the data from the remote processors
    const int ierr =
      MPI_Gatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		   &r[0], &sendlengths[0], &displacements[0], datatype<T>(),
		   root_id,
		   libMesh::COMM_WORLD);

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("gather()", "Parallel");
  }


  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv)
  {
    //assert(recv.size() == libMesh::n_processors());
    recv.resize(libMesh::n_processors());
    
    if (libMesh::n_processors() > 1)
      {
	MPI_Allgather (&send,
		       1,
		       datatype<T>(),
		       &recv[0],
		       1, 
		       datatype<T>(),
		       libMesh::COMM_WORLD);
      }
    else
      recv[0] = send;
  }



  template <typename T>
  inline void allgather(std::complex<T> send,
			std::vector<std::complex<T> > &recv)
  {
    assert(recv.size() == libMesh::n_processors());

    if (libMesh::n_processors() > 1)
      {
	std::vector<T> temprealoutput(recv.size()),
	  tempimagoutput(recv.size());
	T realinput = send.real(),
	  imaginput = send.imag();

	MPI_Allgather (&realinput,
		       1,
		       datatype<T>(),
		       &temprealoutput[0],
		       1, 
		       libMesh::COMM_WORLD);

	MPI_Allgather (&imaginput,
		       1,
		       datatype<T>(),
		       &tempimagoutput[0],
		       1, 
		       libMesh::COMM_WORLD);

	for (unsigned int i=0; i != recv.size(); ++i)
	  {
	    recv[i].real() = temprealoutput[i];
	    recv[i].imag() = tempimagoutput[i];
	  }
      }
    else
      recv[0] = send;
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
   * must be called by all processors.
   */
  template <typename T>
  void allgather(std::vector<T> &r)
  {
    if (libMesh::n_processors() == 1)
      return;
    
    START_LOG("allgather()", "Parallel");

    std::vector<int>
      sendlengths  (libMesh::n_processors(), 0),
      displacements(libMesh::n_processors(), 0);

    const int mysize = r.size();
    Parallel::allgather(mysize, sendlengths);

    // Find the total size of the final array and
    // set up the displacement offsets for each processor.
    unsigned int globalsize = 0; 
    for (unsigned int i=0; i != libMesh::n_processors(); ++i)
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
    std::vector<T> r_src(r);

    // now resize it to hold the global data
    r.resize(globalsize);

    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Allgatherv (r_src.empty() ? NULL : &r_src[0], mysize, datatype<T>(),
		      r.empty()     ? NULL : &r[0],     &sendlengths[0], &displacements[0], datatype<T>(),
		      libMesh::COMM_WORLD);

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("allgather()", "Parallel");
  }



  template <typename T>
  void broadcast (T &data, const unsigned int root_id)
  {
    if (libMesh::n_processors() == 1)
      return;
    
    START_LOG("broadcast()", "Parallel");

    // Spread data to remote processors.
    const int ierr =
      MPI_Bcast (&data, 1, datatype<T>(), root_id, libMesh::COMM_WORLD);

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }



  template <typename T>
  void broadcast (std::vector<T> &data, const unsigned int root_id)
  {
    if (libMesh::n_processors() == 1)
      return;
    
    START_LOG("broadcast()", "Parallel");


    // and get the data from the remote processors.
    // Pass NULL if our vector is empty.
    const int ierr =
      MPI_Bcast (data.empty() ? NULL : &data[0], data.size(), datatype<T>(), root_id, libMesh::COMM_WORLD);

    assert (ierr == MPI_SUCCESS);

    STOP_LOG("broadcast()", "Parallel");
  }


#else // HAVE_MPI

  template <typename T>
  inline bool verify(const T &) { return true; }

  template <typename T>
  inline void min(T &) {}

  template <typename T>
  inline void min(std::vector<T> &) {}

  template <typename T>
  inline void max(T &) {}

  template <typename T>
  inline void max(std::vector<T> &) {}

  template <typename T>
  inline void sum(T &) {}

  template <typename T>
  inline void sum(std::vector<T> &) {}

  template <typename T>
  inline void send_receive(const unsigned int,
                           T &send,
			   const unsigned int,
                           T &recv)
  {
    recv = send;
  }

  template <typename T>
  inline void gather(const unsigned int root_id,
		     T send,
		     std::vector<T> &recv)
  {
    assert (!root_id);
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void gather(const unsigned int, std::vector<T> &) {}

  template <typename T>
  inline void allgather(T send,
			std::vector<T> &recv)
  {
    recv.resize(1);
    recv[0] = send;
  }

  template <typename T>
  inline void allgather(std::vector<T> &) {}

  template <typename T>
    inline void broadcast (std::vector<T> &, const unsigned int =0) {}

#endif // HAVE_MPI


}

#endif // #define __parallel_h__
