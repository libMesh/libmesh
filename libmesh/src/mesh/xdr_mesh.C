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
#include "xdr_mesh.h"
#include "xdr_mhead.h"
#include "enum_elem_type.h" // for ElemType

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
  // before we write the title and id string

  // Both DEAL and LIBM style files have augmented headers.
  if ((orig_flag == 0) || (orig_flag == 2)) 
    {

      switch (m_type)
	{
	  
#ifdef HAVE_XDR
	  
	case (XdrMGF::ENCODE):
	case (XdrMGF::DECODE):
	  {
	    // this used to be 0.  How did that work?
	    unsigned int temp_n_blocks = hd->get_n_blocks(); 
	    xdr_u_int(mp_xdr_handle, &temp_n_blocks);
	    hd->set_n_blocks(temp_n_blocks);

	    // The number of blocks (i.e. the number of element types)
	    // for any mesh must always
	    // be at least 1.
	    assert(hd->get_n_blocks() != 0);
	    break;
	  }

#endif
	  
	case (XdrMGF::W_ASCII):
	  {
	    mp_out << hd->get_n_blocks() << "\t # Num. Element Blocks.\n";
	    break;
	  }

	case (XdrMGF::R_ASCII):
	  {
	    assert (mp_in.good());
	    unsigned int temp_n_blocks=0;
	    mp_in >> temp_n_blocks;
	    hd->set_n_blocks(temp_n_blocks);
	    mp_in.getline(comment, comm_len);
	    break;
	  }

	default:
	  // Unknown access type
	  error();
	}

      
      std::vector<ElemType> et;
      hd->get_block_elt_types(et);

      
      // Note:  If DECODING or READING, allocate space in the vector
      if ((m_type == DECODE) || (m_type == R_ASCII))
	et.resize(hd->get_n_blocks());  


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
	    for (unsigned int i=0; i<hd->get_n_blocks(); i++)
	      mp_out << et[i] << " ";
	      
	    mp_out << "\t # Element types in each block.\n";
	    break;
	  }

	case (XdrMGF::R_ASCII):
	  {
	    assert (mp_in.good());
	
	    for (unsigned int i=0; i<hd->get_n_blocks(); i++)
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


      
      // Note:  If DECODING or READING, you need to set the value 
      // in the header data structure.
      if ((m_type == DECODE) || (m_type == R_ASCII)) 
	hd->set_block_elt_types(et);


      std::vector<unsigned int> neeb;
      hd->get_num_elem_each_block(neeb);

      // If DECODING or READING, allocate space for the vector 
      if ((m_type == DECODE) || (m_type == R_ASCII))
	neeb.resize( hd->get_n_blocks()*(this->get_num_levels()+1) );

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
	    for (unsigned int i=0; i<neeb.size(); i++)
	      mp_out << neeb[i] << " ";
	      
	    mp_out << "\t # Num. of elements in each block at each level.\n";
	    break;
	  }

	case (XdrMGF::R_ASCII):
	  {

	    // We will treat this line as containing
	    // 1.) The number of elements in each block OR
	    // 2.) The number of elements at each level in each block
	    // Therefore, we don't know a-priori how many ints to read.

	    // Get the full line from the stream up to the newline
	    mp_in.getline(comment, comm_len);

	    // Construct a char buffer to hold the tokens as we
	    // process them, and construct a std::string object and
	    // a std::stringstream object for tokenizing this line.
	    char token[comm_len];
	    std::string s_temp(comment);
	    std::stringstream ss(s_temp);

	    // Resize the neeb vector to zero so we can push back
	    // values onto it.  Note that we are using a tokenizer
	    // scheme again here to read the line, but it's not entirely
	    // necessary since we know the size neeb should have.
	    neeb.resize(0);

	    // Process the tokens one at a time
	    while (ss >> token)
              {
                // If you reach the hash, the rest of the line is a comment,
                // so quit reading.
                if (token[0] == '#')
                  break;

                // If you reach an alphabetic character, this is an error
                if (!isdigit(token[0]))
		  {
		    std::cerr << "Error: Unrecognized character detected." 
			      << std::endl;
		    error();
		  }

                // Otherwise, add the value to the neeb vector
                neeb.push_back( std::atoi(token) );
              }
              
	    // Be sure you have the right number of entries in neeb
	    assert (neeb.size() == (hd->get_n_blocks() * (this->get_num_levels()+1)));

	    break;
	  }

	default:
	  // Unknown access type
	  error();
	}

      if ((m_type == DECODE) || (m_type == R_ASCII)) 
	hd->set_num_elem_each_block(neeb);      
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
	xdr_string(mp_xdr_handle,&temp, ((m_type == XdrMGF::ENCODE) ? std::strlen(temp) : hd->m_strSize));
	hd->setId(temp);
	delete [] temp;

	temp = hd->cpyString(hd->getTitle());

	xdr_string(mp_xdr_handle,&temp, ((m_type == XdrMGF::ENCODE) ? std::strlen(temp) : hd->m_strSize));
	hd->setTitle(temp);
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
