// $Id: xdrIO.h,v 1.14 2003-10-02 03:39:25 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __xdrIO_h__
#define __xdrIO_h__

#include "libmesh_config.h"

// C++ includes
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <string>

#ifdef HAVE_RPC_RPC_H
# include <rpc/rpc.h>
#endif


// Local includes
#include "libmesh_common.h"
#include "enum_elem_type.h"
#include "o_f_stream.h"

#ifndef SINGLE_PRECISION
#define xdr_REAL xdr_double
#else
#define xdr_REAL xdr_float
#endif


// forward declarations
class Originator; 

/**
 * This class is taken directly from MGF.
 * It facilitates reading and writing
 * binary solution/mesh
 * files using the \p xdr binary
 * format, which allows for portable
 * binary files across various
 * platforms.  For more information
 * on the \p xdr format, see
 * the standard C include file
 * \p rpc/rpc.h.
 *
 * @author Bill Barth, Robert McLay.
 */ 
class XdrIO
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
#ifdef HAVE_RPC_RPC_H
  XdrIO() : m_type(UNKNOWN), mp_xdr_handle(0), orig_flag(0), mp_fp(0) {}
#else
  XdrIO() : m_type(UNKNOWN), orig_flag(0), mp_fp(0) {}
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
  virtual ~XdrIO();

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
  int dataBlk(REAL* array, int numvar, int size);

  /**
   * Return an Originator data structure
   * for this code. We're using Deal 3.3
   * as our originator signature.
   */
  const Originator get_originator();

  /**
   * Get the originator flag.
   */
  int get_orig_flag() const { return orig_flag; }

  
 protected:

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

#ifdef HAVE_RPC_RPC_H
  
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
  int orig_flag;

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
  FILE* mp_fp;
};

/**
 * This is a helper class for XdrIO.
 * It defines a data structure which
 * contains information (name,
 * major version number, minor version
 * number) about the software which
 * originally wrote (or is writing)
 * the current XDR file.
 *
 * @author John W. Peterson, 2002
 */	     
class Originator
{
 public:
  
  /**
   * Initialization constructor.
   * Sets the name, major, and
   * minor version numbers of
   * the Originator data structure.
   * Use this one!
   */
  Originator(const char* theName,
	     int   theMajorVersion,
	     int   theMinorVersion)  : name(theName),
                                       major_version(theMajorVersion),
                                       minor_version(theMinorVersion)
    {}

  /**
   * Destructor.
   */
  ~Originator() {} 

  /**
   * Overloads equality operator.
   * Useful for checking to see
   * if the Originator read from
   * the XDR file is the same as
   * another Originator.
   * @return True when two Originators are equivalent
   */
  bool operator == (const Originator &rhs) const;

  /**
   * @return Name of the code, generally only four chars
   */
  const char* get_name() const { return name; }

  /**
   * @return Major version number of the code.
   */
  int get_major()  const { return major_version; }

  /**
   * @return Minor version number of the code.
   */
  int get_minor()  const { return minor_version; }
  
 private:

  /**
   * Default constructor -- you aren't allowed to call this.
   */
  Originator() {}
  
  const char* name;
  int   major_version;
  int   minor_version;
};

inline
bool Originator::operator == (const Originator &rhs) const
{
  bool b=false;

  //std::cout << "Values being compared are: \"" << this->get_name() <<
  //  "\" and \"" << rhs.get_name() << "\"" << std::endl;
  
  //int dbug = strncmp(this->get_name(), rhs.get_name(), 4);
  //std::cout << "Value of string comparison was: " << dbug << std::endl;
  
  if ( (!strncmp(this->get_name(), rhs.get_name(), 4)) &&
       (this->get_major() == rhs.get_major())      &&
       (this->get_minor() == rhs.get_minor()) )
    b=true;
      
  return b;
}

class XdrMHEAD; // forward declaration

/**
 * The \p XdrMESH class.
 * This class is responsible
 * for reading/writing
 * information about the mesh
 * to \p xdr style binary files.
 *
 * @author Bill Barth, Robert McLay.
 */
class XdrMESH: public XdrIO
{
 public:

  /**
   * Constructor.  Initializes
   * \p m_dim to -1.
   */
  XdrMESH() : m_dim(-1) {}

  /**
   * Calls the \p init method
   * in the parent class, \p XdrIO
   * with the appropriate parameters.
   *
   * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
   * \param fn const char pointer which points to the filename 
   * \param icnt Number to be appended to file e.g. \p name.mesh.0000
   * \param dim Problem dimension (always three in MGF)
   */
  void init(XdrIO_TYPE type, const char* fn, int icnt, int dim=3) 
    { XdrIO::init(type, fn, "mesh", icnt); m_dim = dim;}

  /**
   * Destructor.
   */
  ~XdrMESH() {}

  /**
   * Read/Write the mesh_base.header.
   * Uses \p xdr_int found
   * in \p rpc/rpc.h.
   *
   * \param hd Pointer to an \p xdr mesh_base.header object
   * @return 1 on success
   */
  int header(XdrMHEAD *hd);

  /**
   * Read/Write an integer connectivity array 
   *
   * \param array Pointer to an array of \p ints
   * \param numvar Total number of variables to be read/written
   * \param num Basically a dummy parameter
   * @return numvar*num
   */
  int Icon(int* array, int numvar, int num)   { return dataBlk(array, numvar, num);}

  /**
   * Read/Write a coord of appropriate size.
   *
   * \param array Pointer to an array of \p Reals
   * \param size Size of \p array (number of elements)
   * @return dim*size
   */
  int coord(Real* array, int dim, int size)  { return dataBlk(array, dim, size);}

  /**
   * Read/Write a BC of appropriate size
   *
   * \param array Pointer to an array of \p Reals
   * \param size Size of \p array (number of elements)
   * @return 3*size
   */
  int BC(int* array, int size)                { return dataBlk(array, 3, size);}

 private:

  /**
   * Dimension of the mesh
   */
  int m_dim;
};



class XdrSHEAD; // forward declaration

/**
 * The \p XdrSOLN class.
 * This class is responsible
 * for reading/writing
 * information about the solution
 * to \p xdr style binary files.
 *
 * @author Bill Barth, Robert McLay.
 */
class XdrSOLN: public XdrIO
{
 public:
  /**
   * Constructor.
   * Initializes \p m_wrtVar to -1.
   */
  XdrSOLN() : m_wrtVar(-1) {}

  /**
   * Calls the \p init method
   * in the parent class, \p XdrIO
   * with the appropriate parameters.
   *
   * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
   * \param fn const char pointer to a file name 
   * \param icnt Number to be appended to file e.g. \p name.soln.0000
   */
  void init(XdrIO_TYPE type, const char* fn, int icnt) 
    {XdrIO::init (type, fn, "soln",icnt);}

  /**
   * Destructor.
   */
  ~XdrSOLN() {}

  /**
   * Read/Write the solution header.
   * Uses \p xdr_int found
   * in \p rpc/rpc.h.
   *
   * \param hd Pointer to an \p xdr solution header object
   * @return 1 on success
   */
  int header(XdrSHEAD *hd);

  /**
   * Read/Write solution values.
   *
   * \param array Pointer to array of \p Reals to be read/written
   * \param size Size of individual variables to be written
   * @return m_wrtVar*size
   */
  int values(Real* array, int size) { return dataBlk(array, m_wrtVar, size);}

 private:
  int m_wrtVar;
};


/**
 * The \p XdrHEAD class.
 * This is a
 * base class for deriving
 * either solution (\p XdrSHEAD)
 * or mesh (\p XdrMHEAD)
 * header interface classes.
 *
 * @author Bill Barth, Robert McLay.
 */
class XdrHEAD
{
 public:
  /**
   * Constructor.
   */
  XdrHEAD();

  /**
   * Destructor.
   */
  virtual ~XdrHEAD();

  /**
   * Set the mesh/solution file id.
   */
  void setId(const char* id)           { delete [] mp_id; mp_id = cpyString(id); }

  /**
   * Get the mesh/solution file id.
   */
  const char* getId() const            { return mp_id; }

  /**
   * Set the mesh/solution file title.
   */
  void setTitle(const char* title)     { delete [] mp_title; mp_title = cpyString(title); }

  /**
   * Get the mesh/solution file title.
   */
  const char* getTitle() const         { return mp_title; }

  /**
   * Set the total number of
   * nodes in the mesh/solution file.
   */
  void setNumNodes(int numNodes)       { m_numNodes = numNodes; }

  /**
   * Get the total number of
   * nodes in the mesh/solution file. 
   */
  int  getNumNodes() const             { return m_numNodes; }

  /**
   * Set the number of
   * boundary conditions in the
   * mesh/solution file.
   */
  void setNumBCs(int numBCs)           { m_numBCs = numBCs; }

  /**
   * Get the number of
   * boundary conditions in
   * them mesh/solution file.
   */
  int  getNumBCs() const               { return m_numBCs; }

  /**
   * Set the string size of the
   * mesh/solution file. (?)
   */
  void setStrSize(int strSize)         { m_strSize = strSize; }

  /**
   * Set the string size of the
   * mesh /solutionfile. (?)
   */
  int  getStrSize() const              { return m_strSize; }
  
 protected:
       
  /**
   * Number of variables written
   * to output, e.g. u,v,w,p,T = 5
   */
  int m_wrtVar;

  /**
   * Total number of variables,
   * may differ from the total
   * number of variables actually
   * written.
   */
  int m_numvar;

  /**
   * The mesh file number
   * which corresponds to a given
   * solution file.
   */
  int m_meshCnt;

  /**
   * The internal solution number.
   */
  int m_kstep;

  /**
   * Number of elemetns in the
   * solution/mesh.
   */
  int m_numel;

  /**
   * Number of nodes in the
   * solution/mesh.
   */
  int m_numNodes;

  /**
   * Total mesh weighting i.e.
   * How many nodes are there
   * and where are they?
   */
  int m_sumWghts;

  /**
   * Number of boundary
   * conditions in the solution/mesh.
   */
  int m_numBCs;

  /**
   * String size (Not sure of what?)
   */
  int m_strSize;

  /**
   * An ID string for the file.
   */
  char* mp_id;

  /**
   * A title string for the file.
   */
  char* mp_title;

  /**
   * User's simulation title
   */
  char* mp_userTitle;

  /**
   * List of null-separated variable names.
   */
  char* mp_varTitle;

  /**
   * Current solution time.
   */
  Real m_time;

  /**
   * Uses memcpy to create an exact
   * copy of \p src, then returns
   * that copy.  Note: I don't know
   * where the memory allocated
   * for this copy gets deleted!
   *
   * @return Copy of \p src
   */
  char* cpyString(const char* src, int len = -1);
  
 private:
  XdrHEAD(const XdrHEAD&);
  const XdrHEAD& operator=(const XdrHEAD&);
};



/**
 * The \p XdrMHEAD class.
 * This class is responsible
 * for reading/writing \p xdr mesh file headers.
 *
 * @author Bill Barth, Robert McLay.  Modified: John W. Peterson
 */
class XdrMHEAD : public XdrHEAD
{
  friend class XdrMESH;
 public:
  /**
   * Constructor.
   */
  XdrMHEAD()                           {}

  /**
   * Destructor.
   */
  ~XdrMHEAD()                          {}

  /**
   * Set the number of
   * elements in the mesh.
   */
  void setNumEl(int numel)             { m_numel = numel; }

  /**
   * Get the number of
   * elements in the mesh.
   */
  int  getNumEl() const                { return m_numel; }

  /**
   * Set the mesh weighting.
   * You probably shouldn't
   * set this yourself ...
   */
  void setSumWghts(int sumWghts)       { m_sumWghts = sumWghts; }

  /**
   * Get the mesh weighting.
   *
   * @sect2{DEAL mesh specific get/set functions}
   */
  int  getSumWghts() const             { return m_sumWghts; }

  /**
   * A mesh block by definition contains
   * only a single type of element.
   *
   * @return The number of mesh blocks.
   */
  unsigned int get_n_blocks() const { return n_blocks; }

  /**
   * Sets the number of mesh blocks.
   */
  void set_n_blocks(const unsigned int nb) { n_blocks = nb; }

  /**
   * Element block types are defined in elem_type.h.
   * They may be for example TRI, TRI6, QUAD, etc.
   *
   * @return The vector of element block types.
   */
  std::vector<ElemType> get_block_elt_types() const { return block_elt_types; }

  /**
   * Set the vector of element block types
   */
  void set_block_elt_types(const std::vector<ElemType> bet) { block_elt_types = bet; }

  /**
   * The size of each element block is
   * the total number of a given type of
   * element in the mesh.
   *
   * @return The vector of block sizes
   */
  std::vector<unsigned int> get_num_elem_each_block() const { return num_elem_each_block; }

  /**
   * Set the vector of block sizes
   */
  void set_num_elem_each_block(const std::vector<unsigned int> neeb) { num_elem_each_block = neeb; }



  
 private:
  /**
   * DEAL mesh specific variables:
   *
   *
   * Tells the total number of element
   * blocks.  An element block is
   * contains only a single type of
   * element.
   */
  unsigned int n_blocks;

  /**
   * A vector of length n_blocks
   * which describes the elemnt type
   * in each block e.g. TRI, QUAD, etc.
   * Note: The element type uniquely
   * defines the number of nodes for
   * that element.
   * @see elem_type.h for more
   */
  std::vector<ElemType> block_elt_types;

  /**
   * A vector of length n_blocks
   * containing the number of elements
   * in each block.
   */
  std::vector<unsigned int> num_elem_each_block;
};


/**
 * The \p XdrSHEAD class.
 * This class is responsible
 * for reading/writing \p xdr solution file headers.
 *
 * @author Bill Barth, Robert McLay.
 */
class XdrSHEAD : public XdrHEAD
{
  friend class XdrSOLN;
 public:
  /**
   * Constructor.
   */
  XdrSHEAD()                                    {}

  /**
   * Destructor.
   */
  ~XdrSHEAD()                                   {}

  /**
   * Set the total number of
   * solution variables.
   */
  void setNumVar(int numvar)                    { m_numvar = numvar; }

  /**
   * Get the total number of
   * solution variables.
   */
  int  getNumVar() const                        { return m_numvar; }

  /**
   * Set the number of written
   * solution variables.
   */
  void setWrtVar(int wrtVar)                    { m_wrtVar = wrtVar; }

  /**
   * Get the number of written
   * solution variables.
   */
  int  getWrtVar() const                        { return m_wrtVar; }

  /**
   * Set the mesh file number.
   */
  void setMeshCnt(int meshCnt)                  { m_meshCnt = meshCnt; }

  /**
   * Get the mesh file number.
   */
  int  getMeshCnt() const                       { return m_meshCnt; }

  /**
   * Set the solution step
   * number.
   */
  void setKstep(int kstep)                      { m_kstep = kstep; }

  /**
   * Get the solution step
   * number.
   */
  int  getKstep() const                         { return m_kstep; }

  /**
   * Set the solution time.
   */
  void setTime(Real time)                       { m_time = time; }

  /**
   * Get the solution time.
   */
  Real getTime() const                          { return m_time; }

  /**
   * Set the user solution title.
   */
  void setUserTitle(const char* title)          { delete [] mp_userTitle; mp_userTitle = cpyString(title); }

  /**
   * Get the user solution title.
   */
  const char* getUserTitle() const              { return mp_userTitle; }

  /**
   * Set null-terminated list of
   * variable names.
   */
  void setVarTitle(const char* titles, int len) { delete [] mp_varTitle; mp_varTitle = cpyString(titles, len); }

  /**
   * Get null-terminated list of
   * variable names.
   */
  const char* getVarTitle() const               { return mp_varTitle; }

};



#endif
