// $Id: xdr_cxx.h,v 1.3 2003-01-20 17:06:13 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_common.h"

// C++ includes
#ifdef HAVE_RPC_RPC_H
#  include <rpc/rpc.h>
#endif

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#ifdef USE_COMPLEX_NUMBERS
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
   * Defines the method used to access the file.
   */  
  enum XdrMODE {UNKNOWN = -1, ENCODE=0, DECODE, WRITE, READ};

  /**
   * Constructor.  Takes the filename and the mode.
   * Valid modes are ENCODE, DECODE, READ, and WRITE.
   */
  Xdr (const std::string& name="", const XdrMODE m=UNKNOWN);

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
  bool reading() const { return ((mode == DECODE) || (mode == READ)); };

  /**
   * Returns true if the file is opened in a writing
   * state, false otherwise.
   */
  bool writing() const { return ((mode == ENCODE) || (mode == WRITE)); };

  /**
   * Returns the mode used to access the file.  Valid modes
   * are ENCODE, DECODE, READ, or WRITE.
   */
  XdrMODE access_mode() const { return mode; };

  // Data access methods

  /**
   * Inputs or outputs a single integer.
   */
  void data(int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (int& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (int& a) { assert (reading()); data(a); return *this; };
  
  /**
   * Inputs or outputs a single unsigned integer.
   */
  void data(unsigned int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned int& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned int& a) { assert (reading()); data(a); return *this; };

  /**
   * Inputs or outputs a single short integer.
   */
  void data(short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (short int& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (short int& a) { assert (reading()); data(a); return *this; };

  /**
   * Inputs or outputs a single unsigned short integer.
   */
  void data(unsigned short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned short int& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned short int& a) { assert (reading()); data(a); return *this; };

  /**
   * Inputs or outputs a single float.
   */
  void data(float& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (float& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (float& a) { assert (reading()); data(a); return *this; };

  /**
   * Inputs or outputs a single double.
   */
  void data(double& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (double& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (double& a) { assert (reading()); data(a); return *this; };


#ifdef USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a single complex<float>.
   */
  void data(std::complex<float>& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::complex<float>& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::complex<float>& a) { assert (reading()); data(a); return *this; };

  /**
   * Inputs or outputs a single complex<double>.
   */
  void data(std::complex<double>& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::complex<double>& a) { assert (writing()); data(a); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::complex<double>& a) { assert (reading()); data(a); return *this; };


#endif

  /**
   * Inputs or outputs a vector of integers.
   */
  void data(std::vector<int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<int>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<int>& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of unsigned integers.
   */
  void data(std::vector<unsigned int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<unsigned int>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<unsigned int>& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of short integers.
   */
  void data(std::vector<short int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<short int>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<short int>& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of unsigned short integers.
   */
  void data(std::vector<unsigned short int>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<unsigned short int>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<unsigned short int>& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of floats.
   */
  void data(std::vector<float>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<float>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<float>& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of doubles.
   */
  void data(std::vector<double>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector<double>& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector<double>& v) { assert (reading()); data(v); return *this; };


#ifdef USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a vector of complex<float>.
   */
  void data(std::vector< std::complex<float> >& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector< std::complex<float> >& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector< std::complex<float> >& v) { assert (reading()); data(v); return *this; };

  /**
   * Inputs or outputs a vector of complex<double>.
   */
  void data(std::vector< std::complex<double> >& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::vector< std::complex<double> >& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::vector< std::complex<double> >& v) { assert (reading()); data(v); return *this; };

#endif


  /**
   * Inputs or outputs a single string.
   */
  void data(std::string& s, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::string& v) { assert (writing()); data(v); return *this; };

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::string& v) { assert (reading()); data(v); return *this; };


private:

  /**
   * The mode used for accessing the file.
   */ 
  const XdrMODE mode;

#ifdef HAVE_RPC_RPC_H
  
  /**
   * Pointer to the standard @p xdr
   * struct.  See the standard
   * header file rpc/rpc.h
   * for more information.
   */
  XDR* xdrs;
  
#endif

  /**
   * File pointer.
   */
  FILE* fp;

  /**
   * The output file stream.
   */
  std::ofstream out;

  /**
   * The input file stream.
   */
  std::ifstream in;

  /**
   * A buffer to put comment strings into.
   */
  const int comm_len;
  char comm[xdr_MAX_STRING_LENGTH];  
};



#endif
