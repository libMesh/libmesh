// $Id: xdr_io.C,v 1.5 2004-11-15 22:09:14 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// C++ includes
#include <iostream>
#include <iomanip>
#include <stdio.h>  // for FILE
#include <string.h> // for strncmp
#include <vector>
#include <string>

// Local includes
#include "mesh.h"
#include "mesh_data.h"
#include "mesh_tools.h"
#include "cell_hex27.h" // Needed for MGF-style Hex27 meshes
#include "xdr_io.h"
#include "enum_elem_type.h"
#include "o_f_stream.h"
#include "utility.h"

#ifdef HAVE_XDR
#  include <rpc/rpc.h>
#  ifndef SINGLE_PRECISION
#    define xdr_REAL xdr_double
#  else
#    define xdr_REAL xdr_float
#  endif
#endif



//-----------------------------------------------------------------------------
// anonymous namespace to hold helper classes
namespace
{
  //---------------------------------------------------------------------------
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
#ifdef HAVE_XDR
    XdrMGF() : m_type(UNKNOWN), mp_xdr_handle(0), orig_flag(0), mp_fp(0) {}
#else
    XdrMGF() : m_type(UNKNOWN), orig_flag(0), mp_fp(0) {}
#endif
    
    //-------------------------------------------------------------------------
    /**
     * This is a helper class for XdrMGF.
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
    int dataBlk(REAL* array, int numvar, int size);

    /**
     * Return an Originator data structure
     * for this code. We're using Deal 3.3
     * as our originator signature.
     */
    Originator get_originator();

    /**
     * Get the originator flag.
     */
    int get_orig_flag() const { return orig_flag; }

    /**
     * Set the originator flag.
     */
    void set_orig_flag(const unsigned int in_orig_flag) { orig_flag=in_orig_flag; }
  
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

#ifdef HAVE_XDR
  
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


  
  //---------------------------------------------------------------------------
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
  class XdrMESH: public XdrMGF
  {
  public:

    /**
     * Constructor.  Initializes
     * \p m_dim to -1.
     */
    XdrMESH() : m_dim(-1) {}

    /**
     * Calls the \p init method
     * in the parent class, \p XdrMGF
     * with the appropriate parameters.
     *
     * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
     * \param fn const char pointer which points to the filename 
     * \param icnt Number to be appended to file e.g. \p name.mesh.0000
     * \param dim Problem dimension (always three in MGF)
     */
    void init(XdrIO_TYPE type, const char* fn, int icnt, int dim=3) 
    { XdrMGF::init(type, fn, "mesh", icnt); m_dim = dim;}

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



  //---------------------------------------------------------------------------
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
  class XdrSOLN: public XdrMGF
  {
  public:
    /**
     * Constructor.
     * Initializes \p m_wrtVar to -1.
     */
    XdrSOLN() : m_wrtVar(-1) {}

    /**
     * Calls the \p init method
     * in the parent class, \p XdrMGF
     * with the appropriate parameters.
     *
     * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
     * \param fn const char pointer to a file name 
     * \param icnt Number to be appended to file e.g. \p name.soln.0000
     */
    void init(XdrIO_TYPE type, const char* fn, int icnt) 
    {XdrMGF::init (type, fn, "soln",icnt);}

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

  

  //---------------------------------------------------------------------------
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

//     /**
//      * Set the string size of the
//      * mesh /solutionfile. (?)
//      */
//     int  getStrSize() const              { return m_strSize; }
  
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



  //---------------------------------------------------------------------------
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
     * They may be for example TRI3, TRI6, QUAD4, etc.
     *
     * @return A writeable reference to the vector of element block types.
     */
    void get_block_elt_types(std::vector<ElemType>& bet) const { bet = block_elt_types; }

    /**
     * Set the vector of element block types
     */
    void set_block_elt_types(const std::vector<ElemType>& bet) { block_elt_types = bet; }

    /**
     * The size of each element block is
     * the total number of a given type of
     * element in the mesh.
     *
     * @return The vector of block sizes
     */
    void get_num_elem_each_block(std::vector<unsigned int>& neeb) const { neeb = num_elem_each_block; }

    /**
     * Set the vector of block sizes
     */
    void set_num_elem_each_block(const std::vector<unsigned int>& neeb) { num_elem_each_block = neeb; }

    
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



  //---------------------------------------------------------------------------
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

//     /**
//      * Get the total number of
//      * solution variables.
//      */
//     int  getNumVar() const                        { return m_numvar; }

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

//     /**
//      * Get the mesh file number.
//      */
//     int  getMeshCnt() const                       { return m_meshCnt; }

    /**
     * Set the solution step
     * number.
     */
    void setKstep(int kstep)                      { m_kstep = kstep; }

//     /**
//      * Get the solution step
//      * number.
//      */
//     int  getKstep() const                         { return m_kstep; }

    /**
     * Set the solution time.
     */
    void setTime(Real time)                       { m_time = time; }

//     /**
//      * Get the solution time.
//      */
//     Real getTime() const                          { return m_time; }

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



  // ------------------------------------------------------------
  // XdrMGF members
  XdrMGF::~XdrMGF()
  {
    this->fini();
  }



  void XdrMGF::fini()
  {
  
#ifdef HAVE_XDR
  
    if (mp_xdr_handle)
      {
	//std::cout << "Destroying XDR file handle." << std::endl;
	xdr_destroy(mp_xdr_handle);
      }
  
    //std::cout << "Deleting the file handle pointer." << std::endl;
    delete mp_xdr_handle;
  
    mp_xdr_handle = NULL;
  
#endif
  
    if (mp_fp)
      {
	//std::cout << "Closing file." << std::endl;
	fflush(mp_fp);
	fclose(mp_fp);
      }

    mp_fp = NULL;
  }



  XdrMGF::Originator XdrMGF::get_originator()
  {
    switch (orig_flag)
      {
      case 0:
	{
	  const Originator orig_deal("DEAL", 3, 3);
	  return orig_deal;
	}

      case 1:
	{
	  const Originator orig_mgf("MGF ", 2, 0);
	  return orig_mgf;
	}
      }

  
    // Huh?  
    std::cerr << "orig_flag set to unknown value: "
	      << orig_flag << std::endl;
    error();
    Originator never ("ERROR", 0, 0);
    return never;
  }



  void XdrMGF::init (XdrMGF::XdrIO_TYPE t, const char* fn, const char*, int)
  {
    m_type=t;

    // Close old file if necessary
    if (mp_fp) this->fini(); 

  
    // Open file 
    switch (m_type)
      {
      
#ifdef HAVE_XDR
      
      case (XdrMGF::ENCODE):
      case (XdrMGF::DECODE):
	{
	  mp_fp = fopen (fn, (m_type == ENCODE) ? "w" : "r");

	  // Make sure the file is ready for use
	  if (!mp_fp)
	    {
	      std::cerr << "XDR Error: Accessing file: "
			<< fn
			<< " failed."
			<< std::endl;
	      error();
	    }

	  // Create the XDR handle 
	  mp_xdr_handle = new XDR;
	  xdrstdio_create(mp_xdr_handle,
			  mp_fp,
			  ((m_type == ENCODE) ? XDR_ENCODE : XDR_DECODE));
	
	  break;
	}
      
#endif
      
      case (XdrMGF::R_ASCII):
	{
	  mp_in.open(fn, std::ios::in);

	  // Make sure the file is ready for use
	  if (!mp_in.good())
	    {
	      std::cerr << "XDR Error: Accessing file: "
			<< fn
			<< " failed."
			<< std::endl;
	      error();
	    }

	  break;
	}
      
      case (XdrMGF::W_ASCII):
	{
	  mp_out.open(fn, std::ios::out);

	  // Make sure the file is ready for use
	  if (!mp_out.good())
	    {
	      std::cerr << "XDR Error: Accessing file: "
			<< fn
			<< " failed."
			<< std::endl;
	      error();
	    }

	  break;
	}
      
      default:
	{
	  std::cout << "Unrecognized file access type!" << std::endl;
	  error();
	}
      }




  
    // Read/Write the file signature
    const int  bufLen = 12;
    char       buf[bufLen+1];

    switch (m_type)
      {
      
#ifdef HAVE_XDR
      
      case (XdrMGF::ENCODE):
	{
	  char* p = &buf[0];
	  const Originator orig = this->get_originator();
	
	  sprintf(&buf[0], "%s %03d:%03d",
		  orig.get_name(),
		  orig.get_major(),
		  orig.get_minor());
	  xdr_string(mp_xdr_handle, &p, bufLen);  // Writes binary signature

	  break;
	}
      
      case (XdrMGF::DECODE):
	{
	  char* p = &buf[0];
	  xdr_string(mp_xdr_handle, &p, bufLen); // Reads binary signature
	
	  break;
	}
      
#endif
      
      case (XdrMGF::W_ASCII):
	{
	  const Originator orig = this->get_originator();
	  sprintf(&buf[0], "%s %03d:%03d",
		  orig.get_name(),
		  orig.get_major(),
		  orig.get_minor());
	  mp_out << buf << '\n';
	
	  break;
	}
      
      case (XdrMGF::R_ASCII):
	{

#ifdef __HP_aCC
	  // weirdly, _only_ here aCC
	  // is not fond of mp_in.getline()
	  // however, using mp_in.getline()
	  // further below is ok...
	  std::string buf_buf;
	  std::getline (mp_in, buf_buf, '\n');
	  assert (buf_buf.size() <= bufLen);

	  buf_buf.copy (buf, std::string::npos);
#else
	  mp_in.getline(buf, bufLen+1);
#endif

	  break;
	}

      default:
	error();
      }



    // If you are reading or decoding, process the signature
    if ((m_type == R_ASCII) || (m_type == DECODE))
      {
	const Originator deal("DEAL", 3, 3);
	const Originator mgf("MGF ", 2, 0);
	int major =-1, minor = -1;
	char name[5];
	strncpy(name, &buf[0], 4);
	name[4] = '\0';
	sscanf(&buf[4], "%d:%d",&major,&minor);

	const Originator fromSig(name, major, minor);

	if (fromSig == deal)
	  {
	    orig_flag = 0; // 0 is the DEAL identifier by definition
	  }
	else if (fromSig == mgf)
	  {
	    orig_flag = 1; // 1 is the MGF identifier by definition
	  }
	else
	  {
	    std::cerr << "No originating software can be determined. Error." << std::endl;
	    std::cerr << "Values found in signature were: " << std::endl;
	    std::cerr << "Software name: " << fromSig.get_name() << std::endl;
	    std::cerr << "Major Verison: " << fromSig.get_major() << std::endl;
	    std::cerr << "Minor Version: " << fromSig.get_minor() << std::endl;
	    error();
	  }
      }
  
  }



  int XdrMGF::dataBlk(int* array, int numvar, int size)
  {
    int totalSize = numvar*size;

    switch (m_type)
      {

#ifdef HAVE_XDR
      
      case (XdrMGF::DECODE):
      case (XdrMGF::ENCODE):
	{
	  xdr_vector(mp_xdr_handle,
		     (char *) &array[0],
		     totalSize, 
		     sizeof(int),
		     (xdrproc_t) xdr_int);
	  break;
	}
      
#endif
      
      case (XdrMGF::W_ASCII):
	{	
	  for (int i=0; i<size; i++)
	    {
	      for (int j=0; j<numvar; j++)
		mp_out << array[i*numvar + j] << " ";
	  
	      mp_out << '\n';
	    }
	
	  mp_out.flush();
	  break;
	}

      case (XdrMGF::R_ASCII):
	{
	  assert (mp_in.good());
	
	  for (int i=0; i<size; i++)
	    {
	      for (int j=0; j<numvar; j++)
		mp_in >> array[i*numvar + j];
	  
	      mp_in.ignore(); // Read newline
	    }
	
	  break;
	}

      default:
	// Unknown access type
	error();
      }

    return totalSize;
  }



  int XdrMGF::dataBlk(REAL* array, int numvar, int size)
  {
    int totalSize = numvar*size;

    // If this function is called by coord(),
    // numvar is the problem dimension, and
    // size is the number of nodes in the problem.
  
    //std::cout << "Total amount of data to be written: " << totalSize << std::endl;
  
    switch (m_type)
      {
      
#ifdef HAVE_XDR
      
      case (XdrMGF::DECODE):
      case (XdrMGF::ENCODE):
	{ 
	  xdr_vector(mp_xdr_handle,
		     (char *) &array[0],
		     totalSize, 
		     sizeof(REAL),
		     (xdrproc_t) xdr_REAL);
	}
      
#endif
      
      case (XdrMGF::W_ASCII):
	{

	  for (int i=0; i<size; i++)
	    {
	      for (int j=0; j<numvar; j++)
		OFSRealscientific(mp_out,12,array[i*numvar + j]) << " \t";
	    
	      mp_out << '\n';
	    }
	
	  mp_out.flush();
	  break;
	}

      case (XdrMGF::R_ASCII):
	{
	  assert (mp_in.good());
	
	  for (int i=0; i<size; i++)
	    {
	      assert (mp_in.good());
	
	      for (int j=0; j<numvar; j++)
		mp_in >> array[i*numvar + j];
	  
	      mp_in.ignore(); // Read newline
	    }
	
	  break;
	}

      default:
	// Unknown access type
	error();
      }
      
    return totalSize;
  }



  // ------------------------------------------------------------
  // XdrMGF::Originator members
  bool XdrMGF::Originator::operator == (const XdrMGF::Originator &rhs) const
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



  // ------------------------------------------------------------
  // XdrMESH members
  int XdrMESH::header(XdrMHEAD *hd)
  {
    // Temporary variables to facilitate stream reading
    const int comm_len= 256;  
    char comment[comm_len];
  
    switch (m_type)
      {
      
#ifdef HAVE_XDR
      
      case (XdrMGF::DECODE):
      case (XdrMGF::ENCODE): 
	{
	  xdr_int(mp_xdr_handle, &(hd->m_numel));
	  xdr_int(mp_xdr_handle, &(hd->m_numNodes));
	  xdr_int(mp_xdr_handle, &(hd->m_sumWghts));
	  xdr_int(mp_xdr_handle, &(hd->m_numBCs));
	  xdr_int(mp_xdr_handle, &(hd->m_strSize));
	  break;
	}

#endif
      
      case (XdrMGF::W_ASCII):
	{
	  mp_out << hd->m_numel    << "\t # Num. Elements\n";
	  mp_out << hd->m_numNodes << "\t # Num. Nodes\n";
	  mp_out << hd->m_sumWghts << "\t # Sum of Element Weights\n";
	  mp_out << hd->m_numBCs   << "\t # Num. Boundary Conds.\n";
	  mp_out << hd->m_strSize  << "\t # String Size (ignore)\n";
	  break;
	}

      case (XdrMGF::R_ASCII):
	{
	  assert (mp_in.good());
	
	  mp_in >> hd->m_numel    ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_numNodes ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_sumWghts ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_numBCs   ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_strSize  ; mp_in.getline(comment, comm_len);

	  assert(mp_in.good());

	  break;
	}

      default:
	// Unknown access type
	error();
      
      }
  
    // Let's write the augmented header information
    // before we write the title and id strings.
  
    if (orig_flag == 0) // DEAL originator, has augmented header
      {

	switch (m_type)
	  {
	  
#ifdef HAVE_XDR
	  
	  case (XdrMGF::ENCODE):
	  case (XdrMGF::DECODE):
	    {
	      xdr_u_int(mp_xdr_handle, &(hd->n_blocks));
	      break;
	    }

#endif
	  
	  case (XdrMGF::W_ASCII):
	    {
	      mp_out << hd->n_blocks << "\t # Num. Element Blocks.\n";
	      break;
	    }

	  case (XdrMGF::R_ASCII):
	    {
	      assert (mp_in.good());
	
	      mp_in >> hd->n_blocks; mp_in.getline(comment, comm_len);
	      break;
	    }

	  default:
	    // Unknown access type
	    error();
	  }

      
	std::vector<ElemType> et;
	hd->get_block_elt_types(et);
      
	// Note:  If DECODING or READING, allocate space in the vector
	if ((m_type == DECODE) || (m_type == R_ASCII)) et.resize(hd->n_blocks);  

	switch (m_type)
	  {
	  
#ifdef HAVE_XDR
	  
	  case (XdrMGF::ENCODE):
	  case (XdrMGF::DECODE):
	    {
	      xdr_vector(mp_xdr_handle,
			 (char *) &et[0],
			 et.size(), 
			 sizeof(unsigned int),
			 (xdrproc_t) xdr_u_int);
	      break;
	    }

#endif

	  case (XdrMGF::W_ASCII):
	    {
	      for (unsigned int i=0; i<hd->n_blocks; i++)
		mp_out << et[i] << " ";
	      
	      mp_out << "\t # Element types in each block.\n";
	      break;
	    }

	  case (XdrMGF::R_ASCII):
	    {
	      assert (mp_in.good());
	
	      for (unsigned int i=0; i<hd->n_blocks; i++)
		{
		  // convoluted way of doing it to
		  // satisfy icc
		  unsigned int type;
		
		  mp_in >> type ;
		
		  et[i] = static_cast<ElemType>(type) ;
		}
	      mp_in.getline(comment, comm_len);
	      break;
	    }

	  default:
	    // Unknown access type
	    error();
	  }


      
	// Note:  If DECODING or READING, you need to set the value in the header data structure.
	if ((m_type == DECODE) || (m_type == R_ASCII)) hd->set_block_elt_types(et);


	std::vector<unsigned int> neeb;
	hd->get_num_elem_each_block(neeb);

	// If DECODING or READING, allocate space for the vector 
	if ((m_type == DECODE) || (m_type == R_ASCII)) neeb.resize(hd->n_blocks);

	switch (m_type)
	  {
	  
#ifdef HAVE_XDR
	  
	  case (XdrMGF::ENCODE):
	  case (XdrMGF::DECODE):
	    {
	      xdr_vector(mp_xdr_handle,
			 (char *) &neeb[0],
			 neeb.size(), 
			 sizeof(unsigned int),
			 (xdrproc_t) xdr_u_int);
	    }

#endif
	  
	  case (XdrMGF::W_ASCII):
	    {
	      for (unsigned int i=0; i<hd->n_blocks; i++)
		mp_out << neeb[i] << " ";
	      
	      mp_out << "\t # Num. of elements in each block.\n";
	      break;
	    }

	  case (XdrMGF::R_ASCII):
	    {
	      assert (mp_in.good());
	
	      for (unsigned int i=0; i<hd->n_blocks; i++)
		mp_in >> neeb[i] ;
	      
	      mp_in.getline(comment, comm_len);
	      break;
	    }

	  default:
	    // Unknown access type
	    error();
	  }
      
	if ((m_type == DECODE) || (m_type == R_ASCII)) hd->set_num_elem_each_block(neeb);      
      }


    else if (orig_flag == 1) // MGF originator
      {
      }
    else  // Unknown Originator!
      {
	error();
      }
  
  


    // Write the ID and TITLE strings (can be safely ignored)
    switch (m_type)
      {

#ifdef HAVE_XDR
      
      case (XdrMGF::ENCODE):
      case (XdrMGF::DECODE):
	{
	  char* temp = hd->cpyString(hd->getId());
	  xdr_string(mp_xdr_handle,&temp, ((m_type == XdrMGF::ENCODE) ? strlen(temp) : hd->m_strSize));
	  hd->setId(temp);
	  delete [] temp;

	  temp = hd->cpyString(hd->getTitle());

	  xdr_string(mp_xdr_handle,&temp, ((m_type == XdrMGF::ENCODE) ? strlen(temp) : hd->m_strSize));
	  hd->setTitle(temp);
	  if(temp) std::cout << strlen(hd->getTitle()) << std::endl;
	  delete [] temp;
	  break;
	}

#endif
      
      case (XdrMGF::W_ASCII):
	{
	  mp_out << hd->mp_id    << '\n';
	  mp_out << hd->mp_title << '\n';
	  break;
	}

      case (XdrMGF::R_ASCII):
	{
	  assert (mp_in.good());
	
	  mp_in.getline(comment, comm_len);
	  hd->setId(comment);

	  assert (mp_in.good());
	
	  mp_in.getline(comment, comm_len);
	  hd->setTitle(comment);

	  break;
	}

      default:
	// Unknown access type
	error();
      }
  
    return 1;
  }



  // ------------------------------------------------------------
  // XdrHEAD members
  XdrHEAD::XdrHEAD() 
  {
    m_wrtVar = 0;
    m_numvar = 0;
  
    m_meshCnt = 0;
    m_kstep = 0;
  
    m_numel = 0;
    m_numNodes = 0;
    m_sumWghts = 0;
    m_numBCs = 0;
    m_strSize = 0;
    mp_id = 0;
    mp_title = 0;
    mp_userTitle = 0;
    mp_varTitle = 0;
  
    m_time = 0;
  }



  XdrHEAD::~XdrHEAD()
  {
    delete [] mp_id;
    delete [] mp_title;
    delete [] mp_userTitle;
    delete [] mp_varTitle; 
  }



  char* XdrHEAD::cpyString(const char* src, int len)
  {
    char* temp = NULL;
    int myLen = len;
    if(src)
      {
	if (myLen == -1)
	  myLen = strlen(src)+1;
	temp = new char[myLen];
	temp = (char *) memcpy(temp, (char *) src, (myLen)*sizeof(char));
      }
    return temp;
  }



  // ------------------------------------------------------------
  // XdrSOLN members
  int XdrSOLN::header(XdrSHEAD *hd)
  {
    // Temporary variables to facilitate stream reading
    const int comm_len= 80;  
    char comment[comm_len];


  
    switch (m_type)
      {
      
#ifdef HAVE_XDR
      
      case (XdrMGF::ENCODE):
      case (XdrMGF::DECODE):
	{
  
	  xdr_int(mp_xdr_handle,  &(hd->m_wrtVar));
	  xdr_int(mp_xdr_handle,  &(hd->m_numvar));
	  xdr_int(mp_xdr_handle,  &(hd->m_numNodes));
	  xdr_int(mp_xdr_handle,  &(hd->m_meshCnt));
	  xdr_int(mp_xdr_handle,  &(hd->m_kstep));
	  xdr_int(mp_xdr_handle,  &(hd->m_strSize));
	  xdr_REAL(mp_xdr_handle, &(hd->m_time));
	
	  m_wrtVar=hd->m_wrtVar;

	  char* temp = const_cast<char *>(hd->getId());
	  xdr_string(mp_xdr_handle,&(temp),
		     ((m_type == XdrMGF::ENCODE) ? strlen(temp)    : hd->m_strSize));
	  hd->setId(temp);
	
	  temp = const_cast<char *>(hd->getTitle());
	  xdr_string(mp_xdr_handle,&(temp),
		     ((m_type == XdrMGF::ENCODE) ? strlen(temp) : hd->m_strSize));
	  hd->setTitle(temp);

	  temp = const_cast<char *>(hd->getUserTitle());
	  xdr_string(mp_xdr_handle,&(temp),
		     ((m_type == XdrMGF::ENCODE) ? strlen(temp) : hd->m_strSize));
	  hd->setUserTitle(temp);
		
	
	  char * tempTitle = new char[hd->m_strSize*m_wrtVar];
  
  
	  if (m_type == XdrMGF::DECODE)
	    {
	      int tempSize = 0;
	      xdr_string(mp_xdr_handle, &tempTitle, hd->m_strSize*m_wrtVar);
	      int olen= strlen(tempTitle);
	      char *p;
	      char *top = tempTitle;
	      for (int ivar = 0; ivar < m_wrtVar; ++ivar)
		{
		  p = strchr(tempTitle,' ');
		  *p = '\0';
		  tempSize = strlen(tempTitle) ;
		  tempTitle+=tempSize+1;
		}
	      tempTitle = top;
	      hd->mp_varTitle = new char[olen];
	      memcpy(hd->mp_varTitle,tempTitle,olen*sizeof(char));
	    }
	  else if (m_type == XdrMGF::ENCODE)
	    {
	      char *p = hd->mp_varTitle;
	      char *top = tempTitle;
	      for (int ivar = 0; ivar < m_wrtVar; ++ivar)
		{
		  int tempSize = strlen(p) + 1;
		  memcpy(tempTitle,p,tempSize*sizeof(char));
		  tempSize = strlen(tempTitle);
		  tempTitle[tempSize] = ' ';
		  tempTitle += tempSize+1;
		  p += tempSize+1;
		}
	      tempTitle = top;
	      xdr_string(mp_xdr_handle, &tempTitle, hd->m_strSize*m_wrtVar);
	    }
	  delete [] tempTitle;

	  return 0;
	}
#endif


      case (XdrMGF::R_ASCII):
	{
	  assert (mp_in.good());
	
	  mp_in >> hd->m_numNodes ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_wrtVar   ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_strSize  ; mp_in.getline(comment, comm_len);
	  mp_in >> hd->m_time     ; mp_in.getline(comment, comm_len);
	
	  mp_in.getline(comment, comm_len);
	  hd->setId(comment);

	  mp_in.getline(comment, comm_len);
	  hd->setTitle(comment);

	  mp_in.getline(comment, comm_len);
	  hd->setUserTitle(comment);

	  m_wrtVar = hd->m_wrtVar;

	  // Read the variable names
	  {
	    std::string var_name;
	    char* titles = new char[hd->m_wrtVar*hd->m_strSize];
	    unsigned int c=0;
	  
	    for (int var=0; var < hd->m_wrtVar; var++)
	      {
		mp_in >> var_name;

		for (unsigned int l=0; l<var_name.size(); l++)
		  titles[c++] = var_name[l];

		titles[c++] = '\0';
	      }

	    mp_in.getline(comment, comm_len);

	    hd->setVarTitle(titles, c);

	    delete [] titles;
	  }

	
	  return 0;
	}

      
      case (XdrMGF::W_ASCII):
	{
	  mp_out << hd->m_numNodes   << "\t # Num. Nodes\n";
	  mp_out << hd->m_wrtVar     << "\t # Num. of Vars\n";
	  mp_out << hd->m_strSize    << "\t # String Size (ignore)\n";
	  mp_out << hd->m_time       << "\t # Current Time\n";
	  mp_out << hd->mp_id        << '\n';
	  mp_out << hd->mp_title     << '\n';
	  mp_out << hd->mp_userTitle << '\n';

	  // write the variable names
	  {
	    const char* p = hd->getVarTitle();

	    for (int var=0; var<hd->m_wrtVar ; var++)
	      {
		mp_out << p << " ";
		p += strlen(p)+1;
	      }	  
	    mp_out << "\t # Variable Names\n";
	  }

	  m_wrtVar = hd->m_wrtVar;

	  return 0;
	}


      
      default:
	// Unknown access type
	error();

      }
  
    return 1;
  }
  
} // end anonymous namespace



// ------------------------------------------------------------
// XdrIO members
void XdrIO::read (const std::string& name)
{
  if (this->binary())
    this->read_binary (name);
  else
    this->read_ascii  (name);
}



void XdrIO::read_mgf (const std::string& name)
{
  if (this->binary())
    this->read_binary (name, 1);
  else
    this->read_ascii  (name, 1);
}



void XdrIO::write (const std::string& name)
{
  if (this->binary())
    this->write_binary (name);
  else
    this->write_ascii  (name);
}



void XdrIO::write_mgf (const std::string& name)
{
  if (this->binary())
    this->write_binary (name, 1);
  else
    this->write_ascii  (name, 1);
}



void XdrIO::read_mgf_soln (const std::string& name,
			   std::vector<Number>& soln,
			   std::vector<std::string>& var_names) const
{
  here();
  std::cerr << "WARNING: this method is deprecated and will disappear soon!"
	    << std::endl;
  
#ifdef USE_COMPLEX_NUMBERS
  
  // buffer for writing separately
  std::vector<Real> real_soln;
  std::vector<Real> imag_soln;
  
  Utility::prepare_complex_data (soln, real_soln, imag_soln);

  this->read_soln (Utility::complex_filename(name, 0), 
		   real_soln, 
		   var_names);
  
  this->read_soln (Utility::complex_filename(name, 1), 
		   imag_soln, 
		   var_names);
  
#else
  
  this->read_soln (name, soln, var_names);
      
#endif
}



void XdrIO::write_mgf_soln (const std::string& name,
			    std::vector<Number>& soln,
			    std::vector<std::string>& var_names) const
{
  here();
  std::cerr << "WARNING: this method is deprecated and will disappear soon!"
	    << std::endl;
  
#ifdef USE_COMPLEX_NUMBERS
  
  // buffer for writing separately
  std::vector<Real> real_soln;
  std::vector<Real> imag_soln;
  
  Utility::prepare_complex_data (soln, real_soln, imag_soln);
  
  this->write_soln (Utility::complex_filename(name, 0), 
		    real_soln, 
		    var_names);
  
  this->write_soln (Utility::complex_filename(name, 1), 
		    imag_soln, 
		    var_names);
  
#else
  
  this->write_soln (name, soln, var_names);
      
#endif
}



void XdrIO::read_ascii (const std::string& name, const unsigned int originator)
{
  // get a writeable reference to the underlying mesh
  Mesh& mesh = this->mesh();
  
  // clear any existing mesh data
  mesh.clear();
    
  // read the mesh
  this->read_mesh (name, originator);
}



void XdrIO::read_binary (const std::string& name, const unsigned int originator)
{
#ifndef HAVE_XDR

  std::cerr << "WARNING: Compiled without XDR binary support.\n"
	    << "Will try ASCII instead" << std::endl << std::endl;

  this->read_ascii (name);
  
#else
  
  // get a writeable reference to the underlying mesh
  Mesh& mesh = this->mesh();
  
  // clear any existing mesh data
  mesh.clear();

  // read the mesh
  this->read_mesh (name, originator);
  
#endif
}



void XdrIO::write_ascii (const std::string& name, const unsigned int originator)
{
  this->write_mesh (name, originator);
}



void XdrIO::write_binary (const std::string& name, const unsigned int originator)
{
#ifndef HAVE_XDR

  std::cerr << "WARNING: Compiled without XDR binary support.\n"
	    << "Will try ASCII instead" << std::endl << std::endl;

  this->write_ascii (name);

#else
  
  this->write_mesh (name, originator);  
  
#endif
}



void XdrIO::read_mesh (const std::string& name,
		       const unsigned int mgf_originator,
		       MeshData* mesh_data)
{
  // get a writeable reference to the mesh
  Mesh& mesh = this->mesh();

  // clear any data in the mesh
  mesh.clear();
  
  // Create an XdrMESH object.
  XdrMESH m;

  // Create a pointer
  // to an XdrMESH file
  // header.
  XdrMHEAD mh;

  // Open the XDR file for reading.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  m.set_orig_flag(mgf_originator);
  m.init((this->binary() ? XdrMGF::DECODE : XdrMGF::R_ASCII), name.c_str(), 0); // mesh files are always number 0 ...

  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  unsigned int              n_blocks = 0;
  std::vector<ElemType>     etypes;
  std::vector<unsigned int> neeb;
	
  // Get the information from
  // the header, and place it
  // in the header pointer.
  m.header(&mh);
	
  // Read information from the
  // file header.  This depends on
  // whether its a libMesh or MGF mesh.
  const int numElem     = mh.getNumEl();
  const int numNodes    = mh.getNumNodes();
  const int totalWeight = mh.getSumWghts();
  const int numBCs      = mh.getNumBCs();
  
  // If a libMesh-type mesh, read the augmented mesh information
  if (m.get_orig_flag() == 0)
    {
      // Read augmented header
      n_blocks = mh.get_n_blocks();
      
      etypes.resize(n_blocks);
      mh.get_block_elt_types(etypes);
      
      neeb.resize(n_blocks);
      mh.get_num_elem_each_block(neeb);
    }

  
  
  // Read the connectivity  
  std::vector<int> conn;
  
  // Now that we know the
  // number of nodes and elements,
  // we can resize the
  // appropriate vectors if we are
  // reading information in.
  mesh.reserve_nodes (numNodes);
  mesh.reserve_elem  (numElem);
  
  // Each element stores two extra
  // locations: one which tells
  // what type of element it is,
  // and one which tells how many
  // nodes it has. Therefore,
  // the total number of nodes
  // (totalWeight) must be augmented
  // by 2 times the number of elements
  // in order to read in the entire
  // connectivity array.
  
  // Note: This section now depends on
  // whether we are reading a libMesh or
  // MGF style mesh.  
  if (m.get_orig_flag() == 0) // libMesh 
    {
      conn.resize(totalWeight);
      m.Icon(&conn[0], 1, totalWeight);
    }
  
  else if (m.get_orig_flag() == 1) // MGF
    {
      conn.resize(totalWeight+(2*numElem));
      m.Icon(&conn[0], 1, totalWeight+(2*numElem));
    }
  
  else
    {
      // I don't know what type of mesh it is.
      error();
    }
    
  // read in the nodal coordinates and form points.
  {
    std::vector<Real> coords(numNodes*mesh.spatial_dimension()); // Always use three coords per node
    m.coord(&coords[0], mesh.spatial_dimension(), numNodes);
  
    // Form Nodes out of
    // the coordinates.  If the    
    // MeshData object is active,
    // add the nodes and ids also          
    // to its map.
    for (int innd=0; innd<numNodes; ++innd)      
      {
	Node* node = mesh.add_point (Point(coords[0+innd*3],  
					   coords[1+innd*3],
					   coords[2+innd*3]));
				       
	if (mesh_data != NULL)
	  if (mesh_data->active())
	    {
	      // add the id to the MeshData, so that
	      // it knows the foreign id, even when 
	      // the underlying mesh got re-numbered,
	      // refined, elements/nodes added...   
	      mesh_data->add_foreign_node_id(node, innd);
	    }
      }  
  }

  
  
  // Build the elements.
  // Note: If the originator was MGF, we don't
  // have to do much checking ...
  // all the elements are Hex27.
  // If the originator was
  // this code, we have to loop over
  // et and neeb to read in all the
  // elements correctly.
  //
  // (This used to be before the coords block, but it
  // had to change now that elements store pointers to
  // nodes.  The nodes must exist before we assign them to
  // the elements. BSK, 1/13/2003)
  if (m.get_orig_flag() == 0) // libMesh-style (0) hybrid mesh possible
    {
      unsigned int lastConnIndex = 0;
      unsigned int lastFaceIndex = 0;
      
      for (unsigned int idx=0; idx<etypes.size(); idx++)  
	{             
	  for (unsigned int e=lastFaceIndex; e<lastFaceIndex+neeb[idx]; e++)  
	    {
	      Elem* elem = mesh.add_elem(Elem::build(etypes[idx]).release());
	      
	      // Add elements with the same id as in libMesh.  
	      // Provided the data files that MeshData reads    
	      // were only written with MeshData, then this      
	      // should work properly.  This is an inline
	      // function, so that for disabled MeshData, this
	      // should not induce too much cost
	      if (mesh_data != NULL)
		mesh_data->add_foreign_elem_id (elem, e);
	      
	      for (unsigned int innd=0; innd < elem->n_nodes(); innd++)
		elem->set_node(innd) = mesh.node_ptr(conn[innd+lastConnIndex]);
	      
	      lastConnIndex += elem->n_nodes();
	    }
	  lastFaceIndex += neeb[idx];
	}
    }
  
  else if (m.get_orig_flag() == 1) // MGF-style (1) Hex27 mesh
    {
      
#ifdef DEBUG
      if (mesh_data != NULL)
	if (mesh_data->active())
	  {
	    std::cerr << "ERROR: MeshData not implemented for MGF-style mesh."
		      << std::endl;
	    error();
	  }
#endif
      
      for (int ielm=0; ielm < numElem; ++ielm)
	{
	  Elem* elem = mesh.add_elem(new Hex27);
	  
	  for (int innd=0; innd < 27; ++innd)
	    elem->set_node(innd) = mesh.node_ptr(conn[innd+2+(27+2)*ielm]);	
	}
    }

  
  // tell the MeshData object that we are finished 
  // reading data
  if (mesh_data != NULL)
    mesh_data->close_foreign_id_maps ();
  
  // Free memory used in
  // the connectivity
  // vector.
  conn.clear();


  // If we are reading,
  // read in the BCs
  // from the mesh file,
  // otherwise write the
  // boundary conditions
  // if the BoundaryInfo
  // object exists.
  if (numBCs > 0)
    {
      std::vector<int> bcs(numBCs*3);

      // Read the BCs from the XDR file
      m.BC(&bcs[0], numBCs);
  
      // Add to the boundary_info
      for (int ibc=0; ibc < numBCs; ibc++)
	mesh.boundary_info.add_side(bcs[0+ibc*3], bcs[1+ibc*3], bcs[2+ibc*3]);
    }
}



void XdrIO::write_mesh (const std::string& name,
			const unsigned int mgf_originator)
{
  // get a read-only reference to the mesh
  const Mesh& mesh = this->cmesh();
  
  // Create an XdrMESH object.
  XdrMESH m;

  // Create a pointer
  // to an XdrMESH file
  // header.
  XdrMHEAD mh;

  // Open the XDR file for writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  m.set_orig_flag(mgf_originator);
  m.init((this->binary() ? XdrMGF::ENCODE : XdrMGF::W_ASCII), name.c_str(), 0); // mesh files are always number 0 ...

  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  std::vector<unsigned int> neeb;

  const int                   numElem = mesh.n_elem();       
  const int                   numBCs  = mesh.boundary_info.n_boundary_conds();
  const std::vector<ElemType> etypes  = mesh.elem_types();
  neeb.resize(etypes.size());
	
  for (unsigned int i=0; i<etypes.size(); i++)
    neeb[i] = mesh.n_elem_of_type(etypes[i]);

  // Now we check to see if we're doing
  // MGF-style headers or libMesh-style
  // "augmented" headers.  An
  // augmented header contains 
  // information about mesh blocks,
  // allowing us to optimize storage
  // and minimize IO requirements
  // for these meshes.
  if (m.get_orig_flag() == 0)
    {
      mh.set_n_blocks(etypes.size());
      mh.set_block_elt_types(etypes);
      mh.set_num_elem_each_block(neeb);
    }
  else
    assert(etypes.size() == 1);
  
  mh.setNumEl(numElem);
  mh.setNumNodes(mesh.n_nodes());
  mh.setStrSize(65536);
  
  if (m.get_orig_flag() == 0)
    mh.setSumWghts(MeshTools::total_weight(mesh));
  else if (m.get_orig_flag() == 1)
    mh.setSumWghts(MeshTools::total_weight(mesh)+2*numElem);
  else
    error();
        
  mh.setNumBCs(numBCs);
  mh.setId("Id String");	   // Ignored
  mh.setTitle("Title String"); // Ignored
  
  // Put the information
  // in the XDR file.
  m.header(&mh); // Needs to work for both types of file
  
  
  // Write the connectivity  
  {
    std::vector<int> conn;
    int orig_type = m.get_orig_flag();
    int totalWeight =0;
    
    if (orig_type == 0)
      totalWeight = MeshTools::total_weight(mesh);
    else if(orig_type == 1)
      totalWeight = MeshTools::total_weight(mesh)+2*mh.getNumEl();
    else
      error();
    
    conn.resize(totalWeight);
    
    std::cout << orig_type << std::endl;
    	
    unsigned int lastConnIndex = 0;
    unsigned int nn = 0;
    
    for (unsigned int idx=0; idx<etypes.size(); idx++)
      {
	nn = lastConnIndex = 0;
	
	for (unsigned int e=0; e<mesh.n_elem(); e++)
	  if (mesh.elem(e)->type() == etypes[idx])
	    {
	      int nstart=0;
	      if (orig_type == 0)
		nn = mesh.elem(e)->n_nodes();
	      else if (orig_type == 1)
		{
		  nstart=2;
		  nn = mesh.elem(e)->n_nodes()+2;
		  conn[lastConnIndex + 0] = 27;
		  conn[lastConnIndex + 1] = 0;
		}
	      else
		error();
                  
	      //std::cout << " HERE " << std::endl;
	      //std::cout << e << std::endl;
        
	      for (unsigned int n=nstart; n<nn; n++)
		conn[lastConnIndex + n] = mesh.elem(e)->node(n-nstart);
		  
	      lastConnIndex += nn;
	    }
	    
	// Send conn to the XDR file
	// std::cout << nn << " " << lastConnIndex/nn << std::endl;
	// m.Icon(&conn[0], nn, totalWeight/nn);
	m.Icon(&conn[0], nn, lastConnIndex/nn);
      }
  }
    
  // create the vector of coords and send
  // it to the XDR file.
  {
    std::vector<Real> coords;
    
    coords.resize(mesh.spatial_dimension()*mesh.n_nodes()); 
    int lastIndex=0;
    for (unsigned int i=0; i<mesh.n_nodes(); i++)
      {
	const Point& p = mesh.node(i);
	
	coords[lastIndex+0] = p(0);
	coords[lastIndex+1] = p(1);
	coords[lastIndex+2] = p(2);
	
	lastIndex += 3;
      }
	  
    // Put the nodes in the XDR file
    m.coord(&coords[0], mesh.spatial_dimension(), mesh.n_nodes()); 
  }

  
  // write the
  // boundary conditions
  // if the BoundaryInfo
  // object exists.
  if (numBCs > 0)
    {
      std::vector<int> bcs(numBCs*3);
    
      std::cout << "numBCs=" << numBCs << std::endl;
    
      //std::cout << "Preparing to write boundary conditions." << std::endl;
      std::vector<unsigned int> elem_list;
      std::vector<unsigned short int> side_list;
      std::vector<short int> elem_id_list;
      
      mesh.boundary_info.build_side_list (elem_list, side_list, elem_id_list);
    
      for (int ibc=0;  ibc<numBCs; ibc++)
	{
	  bcs[0+ibc*3] = elem_list[ibc];
	  bcs[1+ibc*3] = side_list[ibc];
	  bcs[2+ibc*3] = elem_id_list[ibc];
	}
    
      // Put the BCs in the XDR file
      m.BC(&bcs[0], numBCs);
    }
}



void XdrIO::read_soln (const std::string& name,
		       std::vector<Real>& soln,
		       std::vector<std::string>& var_names) const
{
  // Create an XdrSOLN object.
  XdrSOLN s;
  
  // Create an XdrSHEAD object.
  XdrSHEAD sh;
  
  // Open the XDR file for
  // reading or writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  s.init((this->binary() ? XdrMGF::DECODE : XdrMGF::R_ASCII), name.c_str(), 0); // mesh files are always number 0 ...
  
  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  int         numVar      = 0;       
  int         numNodes    = 0;
  const char* varNames;
	
  // Get the information from
  // the header, and place it
  // in the header pointer.
  s.header(&sh);
	
  // Read information from the
  // file header.  This depends on
  // whether its a libMesh or MGF mesh.
  numVar   = sh.getWrtVar();
  numNodes = sh.getNumNodes();
  varNames = sh.getVarTitle();
	
  // Get the variable names
  {	  
    var_names.resize(numVar);
    
    const char* p = varNames;
    
    for (int i=0; i<numVar; i++)
      {
	var_names[i] = p;
	p += strlen(p) + 1;
      }
  }
  
  // Read the soln vector
  soln.resize(numVar*numNodes);
    
  s.values(&soln[0], numNodes);	
}  



void XdrIO::write_soln (const std::string& name,
			std::vector<Real>& soln,
			std::vector<std::string>& var_names) const
{
  // get a read-only reference to the mesh
  const Mesh& mesh = this->cmesh();
  
  // Create an XdrSOLN object.
  XdrSOLN s;
  
  // Create an XdrSHEAD object.
  XdrSHEAD sh;
  
  // Open the XDR file for
  // reading or writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  s.init((this->binary() ? XdrMGF::ENCODE : XdrMGF::W_ASCII), name.c_str(), 0); // mesh files are always number 0 ...

  // Build the header
  sh.setWrtVar(var_names.size());
  sh.setNumVar(var_names.size());
  sh.setNumNodes(mesh.n_nodes());
  sh.setNumBCs(mesh.boundary_info.n_boundary_conds());
  sh.setMeshCnt(0);
  sh.setKstep(0);
  sh.setTime(0.);
  sh.setStrSize(65536);
  sh.setId("Id String");	               // Ignored
  sh.setTitle("Title String");          // Ignored
  sh.setUserTitle("User Title String"); // Ignored
  
  // create the variable array
  {
    std::string var_title;
    
    for (unsigned int var=0; var<var_names.size(); var++)
      {
	for (unsigned int c=0; c<var_names[var].size(); c++)
	  var_title += var_names[var][c];
	
	var_title += '\0';
      }
    
    sh.setVarTitle(var_title.c_str(), var_title.size());
  }
  
  // Put the informationin the XDR file.
  s.header(&sh); // Needs to work for both types of file
  
  // Write the solution vector
  assert (soln.size() == var_names.size()*mesh.n_nodes());
  
  s.values(&soln[0], mesh.n_nodes());
}
