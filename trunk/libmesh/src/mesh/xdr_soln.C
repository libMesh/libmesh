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

// Local includes
#include "xdr_soln.h"
#include "xdr_shead.h"

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
		   ((m_type == XdrMGF::ENCODE) ? std::strlen(temp)    : hd->m_strSize));
	hd->setId(temp);
	
	temp = const_cast<char *>(hd->getTitle());
	xdr_string(mp_xdr_handle,&(temp),
		   ((m_type == XdrMGF::ENCODE) ? std::strlen(temp) : hd->m_strSize));
	hd->setTitle(temp);

	temp = const_cast<char *>(hd->getUserTitle());
	xdr_string(mp_xdr_handle,&(temp),
		   ((m_type == XdrMGF::ENCODE) ? std::strlen(temp) : hd->m_strSize));
	hd->setUserTitle(temp);
		
	
	char * tempTitle = new char[hd->m_strSize*m_wrtVar];
  
  
	if (m_type == XdrMGF::DECODE)
	  {
	    int tempSize = 0;
	    xdr_string(mp_xdr_handle, &tempTitle, hd->m_strSize*m_wrtVar);
	    int olen= std::strlen(tempTitle);
	    char *p;
	    char *top = tempTitle;
	    for (int ivar = 0; ivar < m_wrtVar; ++ivar)
	      {
		p = strchr(tempTitle,' ');
		*p = '\0';
		tempSize = std::strlen(tempTitle) ;
		tempTitle+=tempSize+1;
	      }
	    tempTitle = top;
	    hd->mp_varTitle = new char[olen];
	    std::memcpy(hd->mp_varTitle,tempTitle,olen*sizeof(char));
	  }
	else if (m_type == XdrMGF::ENCODE)
	  {
	    char *p = hd->mp_varTitle;
	    char *top = tempTitle;
	    for (int ivar = 0; ivar < m_wrtVar; ++ivar)
	      {
		int tempSize = std::strlen(p) + 1;
		std::memcpy(tempTitle,p,tempSize*sizeof(char));
		tempSize = std::strlen(tempTitle);
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
	      p += std::strlen(p)+1;
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
  

