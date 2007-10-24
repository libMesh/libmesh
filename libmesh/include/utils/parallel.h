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
#include <vector>

// Local includes
#include "libmesh_common.h" // for Real


// ------------------------------------------------------------
// The Parallel namespace is for wrapper functions
// for common general parallel synchronization tasks

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
  void vector_union(std::vector<T> &r);



  //-----------------------------------------------------------------------
  // Parallel members

#ifdef HAVE_MPI
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
    if (libMesh::n_processors() > 1)
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
    if (libMesh::n_processors() > 1)
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
  void vector_union(std::vector<T> &r)
  {
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

    // copy the input buffer
    std::vector<T> r_src(r);

    // now resize it to hold the global data
    r.resize(globalsize);

    // and get the data from the remote processors
    const int ierr =
      MPI_Allgatherv (&r_src[0], mysize, datatype<T>(),
		      &r[0], &sendlengths[0], &displacements[0], datatype<T>(),
		      libMesh::COMM_WORLD);

    assert (ierr == MPI_SUCCESS);

  }



#else // HAVE_MPI

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
  inline void vector_union(std::vector<T> &r) {}

#endif // HAVE_MPI


}

#endif // #define __parallel_h__
