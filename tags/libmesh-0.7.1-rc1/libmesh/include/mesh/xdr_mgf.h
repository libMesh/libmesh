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

#ifndef __xdr_mgf_h__
#define __xdr_mgf_h__

// C++ includes
#include <cstdio>  // for std::FILE
#include <string>
#include <cstring> // std::strlen, std::strcmp
#include <fstream> // for std::ifstream
#include <sstream>

// Local includes
#include "legacy_xdr_io.h"          // for LegacyXdrIO::FileFormat
#include "libmesh_config.h"  // for LIBMESH_HAVE_XDR
#include "o_f_stream.h"      // for OFStream

// Forward Declarations

#ifdef LIBMESH_HAVE_XDR
#  include <rpc/rpc.h>
#  ifndef LIBMESH_DEFAULT_SINGLE_PRECISION
#    ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
// #      define xdr_REAL xdr_quadruple
// For some reason my xdr implementation doesn't define
// xdr_quadruple... - RHS
#      define xdr_REAL xdr_double
#      define xdr_Real double
#    else
#      define xdr_REAL xdr_double
#      define xdr_Real Real
#    endif
#  else
#    define xdr_REAL xdr_float
#    define xdr_Real Real
#  endif
#else
#  define xdr_Real Real
#endif


namespace libMesh
{

/**
 * This class is taken directly from MGF.  It facilitates reading
 * and writing binary solution/mesh files using the \p xdr binary
 * format, which allows for portable binary files across various
 * platforms.  For more information on the \p xdr format, see the
 * standard C include file \p rpc/rpc.h.
 *
 * @author Bill Barth, Robert McLay.
 *
 * There are essentially two inheritance trees and six classes:
 *       XdrMGF                         XdrHEAD
 *     ^        ^                    ^          ^   
 *     |        |                    |          |
 *  XdrMESH  XdrSOLN             XdrMHEAD     XdrSHEAD
 *
 *
 *
 * XdrHEAD, XdrMHEAD, and XdrSHEAD just read the headers of solution and mesh files.
 *
 * XdrMGF, XdrMESH, and XdrSOLN handle the "meat" of the files: everything
 * other than the headers.
 */ 
class XdrMGF
{
public:
  /**
   * This enum specifies the access
   * permission which will be acquired
   * for the current \p xdr file.
   * Note that it is only possible
   * to read (\p DECODE) or write (\p ENCODE)
   * but not both.  For ASCII type files,
   * use WRITE or READ instead!
   */

  enum XdrIO_TYPE {UNKNOWN = -1, ENCODE=0, DECODE,
		   W_ASCII , R_ASCII};

     
  /**
   * Constructor.  Intializes
   * the access type, \p xdr file
   * handle, \p xdr file pointer,
   * and originator flag. Zero
   * is a good default value for
   * the flag, since that is the
   * DEAL identifier.
   * The \p xdr file handle 
   * is a struct defined in the
   * standard C header \p rpc/rpc.h.
   */
#ifdef LIBMESH_HAVE_XDR
  XdrMGF() : _num_levels(0), m_type(UNKNOWN), mp_xdr_handle(0), orig_flag(LegacyXdrIO::LIBM), mp_fp(0) {}
#else
  XdrMGF() : _num_levels(0), m_type(UNKNOWN), orig_flag(LegacyXdrIO::LIBM), mp_fp(0) {}
#endif
    
  /**
   * Initialization of the \p xdr file.
   * This function performs the following
   * operations:
   * @begin{itemize}
   * @item Closes the old \p xdr file if necessary.
   *
   * @item Creates a new \p xdr file name and opens this file.
   *
   * @item Opens the appropriate \p xdr file handle.
   *
   * @item Reads/Writes a signature to the file.
   *
   * @end{itemize}
   */
  void init(XdrIO_TYPE t, const char* fn, const char* type, int icnt);

  /**
   * Destructor. Frees the memory
   * which was allocated to contain
   * several strings.
   */
  virtual ~XdrMGF();

  /**
   * Finalizes operations on
   * the current \p xdr file handle,
   * and closes the \p xdr file.
   *
   * Uses \p xdr_destroy found in
   * \p rpc/rpc.h.
   */
  void fini();

  /**
   * Reads/Writes a block of \p ints
   * to/from the current \p xdr
   * file/file handle.
   * \param array Pointer to data to be read/written
   * \param numvar The total number of variables (size of the array)
   * \param size The size of each individual variable in the array
   */
  int dataBlk(int*  array, int numvar, int size);

  /**
   * Read/Writes a block of \p Reals
   * to/from the current \p xdr
   * file/file handle.
   */
  int dataBlk(Real* array, int numvar, int size);

  /**
   * Get the originator flag.
   */
  LegacyXdrIO::FileFormat get_orig_flag() const { return orig_flag; }

  /**
   * Set the originator flag.
   */
  void set_orig_flag(LegacyXdrIO::FileFormat in_orig_flag) { orig_flag = in_orig_flag; }


  /**
   * Set number of levels
   */
  void set_num_levels(unsigned int num_levels) { _num_levels = num_levels; }

  /**
   * Get number of levels
   */
  unsigned int get_num_levels() { return _num_levels; }
  
protected:

  /**
   * Number of levels of refinement in the mesh
   */
  unsigned int _num_levels;

  /**
   * Specifies the read/write
   * permission for the current
   * \p xdr file.  Possibilities
   * are:
   * @begin{itemize}
   * @item \p UNKNOWN = -1
   * @item \p ENCODE  = 0
   * @item \p DECODE  = 1
   * @end{itemize}
   */
  XdrIO_TYPE m_type;

#ifdef LIBMESH_HAVE_XDR
  
  /**
   * Pointer to the standard \p{xdr}
   * struct.  See the standard
   * header file \p rpc/rpc.h
   * for more information.
   */
  XDR*  mp_xdr_handle;
  
#endif
  
  /**
   * Flag indicating how much checking
   * we need to do.  We can read in
   * mgf meshes more quickly because
   * there is only one type of element
   * in these meshes.  Deal meshes
   * on the other hand will require
   * a check for each element to find
   * out what type it is.  Possible
   * values are:
   * @begin{itemize}
   * @item 0: It's an DEAL style mesh
   * @item 1: It's a MGF style mesh
   * @end{itemize}
   */
  LegacyXdrIO::FileFormat orig_flag;

  /**
   * An input file stream object
   */
  std::ifstream mp_in;

  /**
   * An output file stream object.
   * Use the customized class to enable
   * features also for compilers with broken
   * iostream
   */
  OFStream mp_out;
  
private:
  std::FILE* mp_fp;

  /**
   * This function allows us to set the number of levels in
   * the mesh when reading.
   */
  void tokenize_first_line(const char* p)
  {
    std::string buf_str(p);
    std::stringstream ss(buf_str);

    char token[256];
    ss >> token;
    if(std::strcmp(token,"LIBM") == 0)
      {
        ss >> token;
        _num_levels = std::atoi(token);
      }

  }
};


} // namespace libMesh

  
#endif // #ifndef __xdr_mgf_h__
