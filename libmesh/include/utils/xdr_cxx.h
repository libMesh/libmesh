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



#ifndef __xdr_cxx_h__
#define __xdr_cxx_h__

// Local includes
#include "libmesh_common.h"
#include "libmesh.h"
#include "enum_xdr_mode.h"
#include "auto_ptr.h"

// C++ includes
#ifdef LIBMESH_HAVE_XDR
#  include <rpc/rpc.h>
#endif

#include <iosfwd>
#include <vector>
#include <string>
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
# include <complex>
#endif

#define xdr_MAX_STRING_LENGTH 256

#ifndef SINGLE_PRECISION
#define xdr_REAL xdr_double
#else
#define xdr_REAL xdr_float
#endif



//--------------------------------------------------------------
// Xdr class definition

/**
 * This class implements a C++ interface to the XDR 
 * (eXternal Data Representation) format.  XDR is useful for
 * creating platform-independent binary files.  This class was
 * created to handle equation system output as a replacement for
 * XdrIO since that is somewhat limited.
 */

class Xdr
{
  
public:

  /**
   * Constructor.  Takes the filename and the mode.
   * Valid modes are ENCODE, DECODE, READ, and WRITE.
   */
  Xdr (const std::string& name="", const libMeshEnums::XdrMODE m=UNKNOWN);

  /**
   * Destructor.  Closes the file if it is open.
   */
  ~Xdr ();

  /**
   * Opens the file.
   */ 
  void open (const std::string& name);

  /**
   * Closes the file if it is open.
   */ 
  void close();

  /**
   * Returns true if the Xdr file is open, false
   * if it is closed.
   */ 
  bool is_open() const;

  /**
   * Returns true if the file is opened in a reading
   * state, false otherwise.
   */
  bool reading() const { return ((mode == DECODE) || (mode == READ)); }

  /**
   * Returns true if the file is opened in a writing
   * state, false otherwise.
   */
  bool writing() const { return ((mode == ENCODE) || (mode == WRITE)); }

  /**
   * Returns the mode used to access the file.  Valid modes
   * are ENCODE, DECODE, READ, or WRITE.
   */
  XdrMODE access_mode () const { return mode; }

  // Data access methods

  /**
   * Inputs or outputs a single integer.
   */
  void data(int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (int& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (int& a) { libmesh_assert (reading()); data(a); return *this; }
  
  /**
   * Inputs or outputs a single unsigned integer.
   */
  void data(unsigned int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned int& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned int& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single short integer.
   */
  void data(short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (short int& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (short int& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single unsigned short integer.
   */
  void data(unsigned short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned short int& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned short int& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single float.
   */
  void data(float& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (float& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (float& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single double.
   */
  void data(double& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (double& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (double& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single long double, but in double precision.
   */
  void data(long double& a, const char* comment="") 
    { double ad = a;
      data(ad, comment); }

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (long double& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (long double& a) { libmesh_assert (reading()); data(a); return *this; }


#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a single complex<double>.
   */
  void data(std::complex<double>& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::complex<double>& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::complex<double>& a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single complex<long double>, but in double
   * precision.
   */
  void data(std::complex<long double>& a, const char* comment="")
    { std::complex<double> ad (a);
      data(ad, comment); }

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::complex<long double>& a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::complex<long double>& a) { libmesh_assert (reading()); data(a); return *this; }

#endif


  /**
   * Inputs or outputs a vector of integers.
   */
  void data(std::vector<int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<int>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<int>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of unsigned integers.
   */
  void data(std::vector<unsigned int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<unsigned int>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<unsigned int>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of short integers.
   */
  void data(std::vector<short int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<short int>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<short int>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of unsigned short integers.
   */
  void data(std::vector<unsigned short int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<unsigned short int>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<unsigned short int>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of floats.
   */
  void data(std::vector<float>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<float>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<float>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of doubles.
   */
  void data(std::vector<double>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<double>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<double>& v) { libmesh_assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of long doubles, but in double
   * precision.
   */
  void data(std::vector<long double>& v, const char* comment="")
    { std::vector<double> vd(v.size());
      for (unsigned int i = 0; i != v.size(); ++i)
	vd[i] = static_cast<double>(v[i]);
      data(vd, comment);
      v.resize(vd.size());
      for (unsigned int i = 0; i != vd.size(); ++i)
	v[i] = static_cast<long double>(vd[i]);
    }

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<long double>& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<long double>& v) { libmesh_assert (reading()); data(v); return *this; }


#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a vector of complex<double>.
   */
  void data(std::vector< std::complex<double> >& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector< std::complex<double> >& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector< std::complex<double> >& v) { libmesh_assert (reading()); data(v); return *this; }

#endif


  /**
   * Inputs or outputs a single string.
   */
  void data(std::string& s, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::string& v) { libmesh_assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::string& v) { libmesh_assert (reading()); data(v); return *this; }


  /**
   * Inputs or outputs a raw data stream.
   */
  template <typename T>
  void data_stream (T *val, const unsigned int len, const unsigned int line_break=libMesh::invalid_uint);

  /**
   * Writes or reads (ignores) a comment line.
   */   
  void comment (std::string &);
  

private:

  /**
   * The mode used for accessing the file.
   */ 
  const XdrMODE mode;

  /**
   * The file name
   */
  std::string file_name;

#ifdef LIBMESH_HAVE_XDR
  
  /**
   * Pointer to the standard @p xdr
   * struct.  See the standard
   * header file rpc/rpc.h
   * for more information.
   */
  XDR* xdrs;

  /**
   * File pointer.
   */
  FILE* fp;
  
#endif

  /**
   * The input file stream.
   */
  AutoPtr<std::istream> in;

  /**
   * The output file stream.
   */
  AutoPtr<std::ostream> out;

  /**
   * A buffer to put comment strings into.
   */
  const int comm_len;
  char comm[xdr_MAX_STRING_LENGTH];  

  /**
   * Are we reading/writing bzipped or gzipped files?
   */
  bool gzipped_file;
  bool bzipped_file;
};



#endif
