// "$Id: xdrIO.C,v 1.22 2004-08-09 17:34:58 jwpeterson Exp $\n"

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
#include <stdio.h>    // for sprintf
#include <assert.h>

// Local includes
#include "xdrIO.h"




XdrIO::~XdrIO()
{
  this->fini();
}



void XdrIO::fini()
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



Originator XdrIO::get_originator()
{
  const Originator orig("DEAL", 3, 3);
  return orig;
}



void XdrIO::init (XdrIO::XdrIO_TYPE t, const char* fn, const char*, int)
{
  m_type=t;

  // Close old file if necessary
  if (mp_fp) this->fini(); 

  
  // Open file 
  switch (m_type)
    {
      
#ifdef HAVE_XDR
      
    case (XdrIO::ENCODE):
    case (XdrIO::DECODE):
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
      
    case (XdrIO::R_ASCII):
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
      
    case (XdrIO::W_ASCII):
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
      
    case (XdrIO::ENCODE):
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
      
    case (XdrIO::DECODE):
      {
	char* p = &buf[0];
	xdr_string(mp_xdr_handle, &p, bufLen); // Reads binary signature
	
	break;
      }
      
#endif
      
    case (XdrIO::W_ASCII):
      {
	const Originator orig = this->get_originator();
	sprintf(&buf[0], "%s %03d:%03d",
		orig.get_name(),
		orig.get_major(),
		orig.get_minor());
	mp_out << buf << std::endl;
	
	break;
      }
      
    case (XdrIO::R_ASCII):
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



int XdrIO::dataBlk(int* array, int numvar, int size)
{
  int totalSize = numvar*size;

  switch (m_type)
    {

#ifdef HAVE_XDR
      
    case (XdrIO::DECODE):
    case (XdrIO::ENCODE):
      {
	xdr_vector(mp_xdr_handle,
		   (char *) &array[0],
		   totalSize, 
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	break;
      }
      
#endif
      
    case (XdrIO::W_ASCII):
      {	
	for (int i=0; i<size; i++)
	  {
	    for (int j=0; j<numvar; j++)
	      mp_out << array[i*numvar + j] << " ";
	  
	    mp_out << std::endl;
	  }
	
	mp_out.flush();
	break;
      }

    case (XdrIO::R_ASCII):
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



int XdrIO::dataBlk(REAL* array, int numvar, int size)
{
  int totalSize = numvar*size;

  // If this function is called by coord(),
  // numvar is the problem dimension, and
  // size is the number of nodes in the problem.
  
  //std::cout << "Total amount of data to be written: " << totalSize << std::endl;
  
  switch (m_type)
    {
      
#ifdef HAVE_XDR
      
    case (XdrIO::DECODE):
    case (XdrIO::ENCODE):
      { 
	xdr_vector(mp_xdr_handle,
		   (char *) &array[0],
		   totalSize, 
		   sizeof(REAL),
		   (xdrproc_t) xdr_REAL);
      }
      
#endif
      
    case (XdrIO::W_ASCII):
      {

	for (int i=0; i<size; i++)
	  {
	    for (int j=0; j<numvar; j++)
	      OFSRealscientific(mp_out,12,array[i*numvar + j]) << " \t";
	    
	    mp_out << std::endl;
	  }
	
	mp_out.flush();
	break;
      }

    case (XdrIO::R_ASCII):
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



int XdrMESH::header(XdrMHEAD *hd)
{
  // Temporary variables to facilitate stream reading
  const int comm_len= 256;  
  char comment[comm_len];
  
  switch (m_type)
    {
      
#ifdef HAVE_XDR
      
    case (XdrIO::DECODE):
    case (XdrIO::ENCODE): 
      {
	xdr_int(mp_xdr_handle, &(hd->m_numel));
	xdr_int(mp_xdr_handle, &(hd->m_numNodes));
	xdr_int(mp_xdr_handle, &(hd->m_sumWghts));
	xdr_int(mp_xdr_handle, &(hd->m_numBCs));
	xdr_int(mp_xdr_handle, &(hd->m_strSize));
	break;
      }

#endif
      
    case (XdrIO::W_ASCII):
      {
	mp_out << hd->m_numel    << "\t # Num. Elements"          << std::endl;
	mp_out << hd->m_numNodes << "\t # Num. Nodes"             << std::endl;
	mp_out << hd->m_sumWghts << "\t # Sum of Element Weights" << std::endl;
	mp_out << hd->m_numBCs   << "\t # Num. Boundary Conds."   << std::endl;
	mp_out << hd->m_strSize  << "\t # String Size (ignore) "  << std::endl;
	break;
      }

    case (XdrIO::R_ASCII):
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
	  
	case (XdrIO::ENCODE):
	case (XdrIO::DECODE):
	  {
	    xdr_u_int(mp_xdr_handle, &(hd->n_blocks));
	    break;
	  }

#endif
	  
	case (XdrIO::W_ASCII):
	  {
	    mp_out << hd->n_blocks << "\t # Num. Element Blocks." << std::endl;
	    break;
	  }

	case (XdrIO::R_ASCII):
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
	  
	case (XdrIO::ENCODE):
	case (XdrIO::DECODE):
	  {
	    xdr_vector(mp_xdr_handle,
		       (char *) &et[0],
		       et.size(), 
		       sizeof(unsigned int),
		       (xdrproc_t) xdr_u_int);
	    break;
	  }

#endif

	case (XdrIO::W_ASCII):
	  {
	    for (unsigned int i=0; i<hd->n_blocks; i++)
	      mp_out << et[i] << " ";
	      
	    mp_out << "\t # Element types in each block." << std::endl;
	    break;
	  }

	case (XdrIO::R_ASCII):
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
	  
	case (XdrIO::ENCODE):
	case (XdrIO::DECODE):
	  {
	    xdr_vector(mp_xdr_handle,
		       (char *) &neeb[0],
		       neeb.size(), 
		       sizeof(unsigned int),
		       (xdrproc_t) xdr_u_int);
	  }

#endif
	  
	case (XdrIO::W_ASCII):
	  {
	    for (unsigned int i=0; i<hd->n_blocks; i++)
	      mp_out << neeb[i] << " ";
	      
	    mp_out << "\t # Num. of elements in each block." << std::endl;
	    break;
	  }

	case (XdrIO::R_ASCII):
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



  else if (orig_flag != 1) // Not MGF originator, unknown Originator!
    {
      error();
    }
  
  


  // Write the ID and TITLE strings (can be safely ignored)
  switch (m_type)
    {

#ifdef HAVE_XDR
      
    case (XdrIO::ENCODE):
    case (XdrIO::DECODE):
      {
	char* temp = const_cast<char *>(hd->getId());
	xdr_string(mp_xdr_handle,&temp,    ((m_type == ENCODE) ? strlen(temp)    : hd->m_strSize));
	hd->setId(temp);

	temp = const_cast<char *>(hd->getTitle());
	xdr_string(mp_xdr_handle,&temp, ((m_type == ENCODE) ? strlen(temp) : hd->m_strSize));
	hd->setTitle(temp);
	break;
      }

#endif
      
    case (XdrIO::W_ASCII):
      {
	mp_out << hd->mp_id    << std::endl;
	mp_out << hd->mp_title << std::endl;
	break;
      }

    case (XdrIO::R_ASCII):
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



int XdrSOLN::header(XdrSHEAD *hd)
{
  // Temporary variables to facilitate stream reading
  const int comm_len= 80;  
  char comment[comm_len];


  
  switch (m_type)
    {
      
#ifdef HAVE_XDR
      
    case (XdrIO::ENCODE):
    case (XdrIO::DECODE):
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
		   ((m_type == XdrIO::ENCODE) ? strlen(temp)    : hd->m_strSize));
	hd->setId(temp);
	
	temp = const_cast<char *>(hd->getTitle());
	xdr_string(mp_xdr_handle,&(temp),
		   ((m_type == XdrIO::ENCODE) ? strlen(temp) : hd->m_strSize));
	hd->setTitle(temp);

	temp = const_cast<char *>(hd->getUserTitle());
	xdr_string(mp_xdr_handle,&(temp),
		   ((m_type == XdrIO::ENCODE) ? strlen(temp) : hd->m_strSize));
	hd->setUserTitle(temp);
		
	
	char * tempTitle = new char[hd->m_strSize*m_wrtVar];
  
  
	if (m_type == XdrIO::DECODE)
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
	else if (m_type == XdrIO::ENCODE)
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


    case (XdrIO::R_ASCII):
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

      
    case (XdrIO::W_ASCII):
      {
	mp_out << hd->m_numNodes << "\t # Num. Nodes"             << std::endl;
	mp_out << hd->m_wrtVar   << "\t # Num. of Vars"           << std::endl;
	mp_out << hd->m_strSize  << "\t # String Size (ignore) "  << std::endl;
	mp_out << hd->m_time     << "\t # Current Time "          << std::endl;
	mp_out << hd->mp_id        << std::endl;
	mp_out << hd->mp_title     << std::endl;
	mp_out << hd->mp_userTitle << std::endl;

	// write the variable names
	{
	  const char* p = hd->getVarTitle();

	  for (int var=0; var<hd->m_wrtVar ; var++)
	    {
	      mp_out << p << " ";
	      p += strlen(p)+1;
	    }	  
	  mp_out << "\t # Variable Names " << std::endl;
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
  if (len == -1)
    len = strlen(src)+1;
  char* temp = new char[len];
  return (char *) memcpy(temp, (char *) src, (len)*sizeof(char));
}

#undef xdr_REAL
