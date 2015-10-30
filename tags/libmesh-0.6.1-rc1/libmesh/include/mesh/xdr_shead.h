// $Id: xdr_shead.h,v 1.1 2007-01-22 17:40:45 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __xdr_shead_h__
#define __xdr_shead_h__

// Local includes
#include "xdr_head.h" // for base class
#include "xdr_soln.h" // for friend

/**
 * The \p XdrSHEAD class.  This class is responsible for
 * reading/writing \p xdr solution file headers.
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


#endif // #ifndef __xdr_shead_h__
