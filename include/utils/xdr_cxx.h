// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_XDR_CXX_H
#define LIBMESH_XDR_CXX_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h"
#include "libmesh/enum_xdr_mode.h" // READ, WRITE, etc.

// C++ includes
#include <memory>
#include <cstdio> // FILE
#ifdef LIBMESH_HAVE_XDR
// I see a redundant declaration warning here on Ubuntu 20.10
#include "libmesh/ignore_warnings.h"
# include <rpc/rpc.h>
# include <rpc/xdr.h>
#include "libmesh/restore_warnings.h"
#endif

#include <iosfwd>
#include <vector>
#include <string>
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
# include <complex>
#endif

const unsigned int xdr_MAX_STRING_LENGTH=256;

#ifndef LIBMESH_DEFAULT_SINGLE_PRECISION
#define xdr_REAL xdr_double
#else
#define xdr_REAL xdr_float
#endif

namespace libMesh
{

/**
 * This class implements a C++ interface to the XDR
 * (eXternal Data Representation) format.  XDR is useful for
 * creating platform-independent binary files.  This class was
 * created to handle equation system output as a replacement for
 * XdrIO since that is somewhat limited.
 *
 * \author Benjamin Kirk
 * \date 2003
 * \brief C++ interface for the XDR (eXternal Data Representation) format.
 */
class Xdr
{

public:

  /**
   * Constructor.  Takes the filename and the mode.
   * Valid modes are ENCODE, DECODE, READ, and WRITE.
   */
  Xdr (const std::string & name="", const XdrMODE m=UNKNOWN);

  /**
   * Destructor.  Closes the file if it is open.
   */
  ~Xdr ();

  /**
   * Opens the file.
   */
  void open (const std::string & name);

  /**
   * Closes the file if it is open.
   */
  void close();

  /**
   * \returns \p true if the Xdr file is open, false
   * if it is closed.
   */
  bool is_open() const;

  /**
   * \returns \p true if the Xdr file being read is at End-Of-File.
   *
   * \note This is \e not a const method - the only portable way to
   * test for an impending EOF is to peek at the next byte of the file
   * first, which may set the eof flag on the istream.
   */
  bool is_eof();

  /**
   * \returns \p true if the file is opened in a reading
   * state, false otherwise.
   */
  bool reading() const { return ((mode == DECODE) || (mode == READ)); }

  /**
   * \returns \p true if the file is opened in a writing
   * state, false otherwise.
   */
  bool writing() const { return ((mode == ENCODE) || (mode == WRITE)); }

  /**
   * \returns The mode used to access the file.  Valid modes
   * are ENCODE, DECODE, READ, or WRITE.
   */
  XdrMODE access_mode () const { return mode; }

  // Data access methods

  /**
   * Inputs or outputs a single value.
   */
  template <typename T>
  void data(T & a, const char * comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template <typename T>
  Xdr & operator << (T & a) { libmesh_assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template <typename T>
  Xdr & operator >> (T & a) { libmesh_assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a raw data stream.
   */
  template <typename T>
  void data_stream (T * val, const unsigned int len, const unsigned int line_break=libMesh::invalid_uint);

  /**
   * Writes or reads (ignores) a comment line.
   */
  void comment (std::string &);

  /**
   * Sets the version of the file that is being read
   */
  void set_version(int ver) { version_number = ver; }

  /**
   * Gets the version of the file that is being read
   */
  int version() const { return version_number; }

private:

  /**
   * Helper method for reading different data types
   */
  template <typename T>
  void do_read(T & a);

  template <typename T>
  void do_read(std::complex<T> & a);

  template <typename T>
  void do_read(std::vector<T> & a);

  template <typename T>
  void do_read(std::vector<std::complex<T>> & a);

  /**
   * Helper method for writing different data types
   */
  template <typename T>
  void do_write(T & a);

  template <typename T>
  void do_write(std::complex<T> & a);

  template <typename T>
  void do_write(std::vector<T> & a);

  template <typename T>
  void do_write(std::vector<std::complex<T>> & a);

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
   * Pointer to the standard XDR struct.  See the standard header file
   * rpc/rpc.h for more information.
   */
  std::unique_ptr<XDR> xdrs;

  /**
   * File pointer.
   */
  FILE * fp;

#endif

  /**
   * The input file stream.
   */
  std::unique_ptr<std::istream> in;

  /**
   * The output file stream.
   */
  std::unique_ptr<std::ostream> out;

  /**
   * A buffer to put comment strings into.
   */
  const int comm_len;
  char comm[xdr_MAX_STRING_LENGTH];

  /**
   * Are we reading/writing zipped files?
   */
  bool gzipped_file, bzipped_file, xzipped_file;

  /**
   * Version of the file being read
   */
  int version_number;
};


} // namespace libMesh


#endif // LIBMESH_XDR_CXX_H
