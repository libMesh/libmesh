// $Id: cell_tet10.C,v 1.5 2003-01-24 17:24:43 jwpeterson Exp $

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


// C++ includes

// Local includes
#include "mesh.h"
#include "cell_tet10.h"




// ------------------------------------------------------------
// Tet10 class member functions
AutoPtr<Elem> Tet10::build_side (const unsigned int i) const
{
  assert (i < n_sides());


  
  AutoPtr<Elem> face(Elem::build(TRI6));
  
  switch (i)
    {
    case 0:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(1);
	face->set_node(3) = get_node(6);
	face->set_node(4) = get_node(5);
	face->set_node(5) = get_node(4);

	return face;
      }
    case 1:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(3);
	face->set_node(3) = get_node(4);
	face->set_node(4) = get_node(8);
	face->set_node(5) = get_node(7);

	return face;
      }
    case 2:
      {
	face->set_node(0) = get_node(1);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(3);
	face->set_node(3) = get_node(5);
	face->set_node(4) = get_node(9);
	face->set_node(5) = get_node(8);

	return face;
      }
    case 3:
      {
	face->set_node(0) = get_node(2);
	face->set_node(1) = get_node(0);
	face->set_node(2) = get_node(3);
	face->set_node(3) = get_node(6);
	face->set_node(4) = get_node(7);
	face->set_node(5) = get_node(9);

	return face;
      }
    default:
      {
	error();
      }
    };

  // We'll never get here.
  error();
  return face;
};



const std::vector<unsigned int> Tet10::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());

  std::vector<unsigned int> conn(8);

  switch (sc)
    {
      
      
      // Linear sub-tet 0
    case 0:
      
      conn[0] = node(0)+1;
      conn[1] = node(4)+1;
      conn[2] = node(6)+1;
      conn[3] = node(6)+1;
      conn[4] = node(7)+1;
      conn[5] = node(7)+1;
      conn[6] = node(7)+1;
      conn[7] = node(7)+1;

      return conn;

      // Linear sub-tet 1
    case 1:
      
      conn[0] = node(4)+1;
      conn[1] = node(1)+1;
      conn[2] = node(5)+1;
      conn[3] = node(5)+1;
      conn[4] = node(8)+1;
      conn[5] = node(8)+1;
      conn[6] = node(8)+1;
      conn[7] = node(8)+1;

      return conn;

      // Linear sub-tet 2
    case 2:
      
      conn[0] = node(5)+1;
      conn[1] = node(2)+1;
      conn[2] = node(6)+1;
      conn[3] = node(6)+1;
      conn[4] = node(9)+1;
      conn[5] = node(9)+1;
      conn[6] = node(9)+1;
      conn[7] = node(9)+1;

      return conn;

      // Linear sub-tet 3
    case 3:
      
      conn[0] = node(7)+1;
      conn[1] = node(8)+1;
      conn[2] = node(9)+1;
      conn[3] = node(9)+1;
      conn[4] = node(3)+1;
      conn[5] = node(3)+1;
      conn[6] = node(3)+1;
      conn[7] = node(3)+1;

      return conn;

      // Linear sub-tet 4
    case 4:
      
      conn[0] = node(4)+1;
      conn[1] = node(8)+1;
      conn[2] = node(6)+1;
      conn[3] = node(6)+1;
      conn[4] = node(7)+1;
      conn[5] = node(7)+1;
      conn[6] = node(7)+1;
      conn[7] = node(7)+1;

      return conn;

      // Linear sub-tet 5
    case 5:
      
      conn[0] = node(4)+1;
      conn[1] = node(5)+1;
      conn[2] = node(6)+1;
      conn[3] = node(6)+1;
      conn[4] = node(8)+1;
      conn[5] = node(8)+1;
      conn[6] = node(8)+1;
      conn[7] = node(8)+1;

      return conn;

      // Linear sub-tet 6
    case 6:
      
      conn[0] = node(5)+1;
      conn[1] = node(9)+1;
      conn[2] = node(6)+1;
      conn[3] = node(6)+1;
      conn[4] = node(8)+1;
      conn[5] = node(8)+1;
      conn[6] = node(8)+1;
      conn[7] = node(8)+1;

      return conn;

      // Linear sub-tet 7
    case 7:
      
      conn[0] = node(7)+1;
      conn[1] = node(6)+1;
      conn[2] = node(9)+1;
      conn[3] = node(9)+1;
      conn[4] = node(8)+1;
      conn[5] = node(8)+1;
      conn[6] = node(8)+1;
      conn[7] = node(8)+1;

      return conn;


    default:

      error();
    };
  
  return conn;
};



void Tet10::vtk_connectivity(const unsigned int sc,
			     std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(4);


  switch (sc)
    {            
      // Linear sub-tet 0
    case 0:
      
      (*conn)[0] = node(0);
      (*conn)[1] = node(4);
      (*conn)[2] = node(6);
      (*conn)[3] = node(7);

      return;

      // Linear sub-tet 1
    case 1:
      
      (*conn)[0] = node(4);
      (*conn)[1] = node(1);
      (*conn)[2] = node(5);
      (*conn)[3] = node(8);

      return;

      // Linear sub-tet 2
    case 2:
      
      (*conn)[0] = node(5);
      (*conn)[1] = node(2);
      (*conn)[2] = node(6);
      (*conn)[3] = node(9);

      return;

      // Linear sub-tet 3
    case 3:
      
      (*conn)[0] = node(7);
      (*conn)[1] = node(8);
      (*conn)[2] = node(9);
      (*conn)[3] = node(3);

      return;

      // Linear sub-tet 4
    case 4:
      
      (*conn)[0] = node(4);
      (*conn)[1] = node(8);
      (*conn)[2] = node(6);
      (*conn)[3] = node(7);

      return;

      // Linear sub-tet 5
    case 5:
      
      (*conn)[0] = node(4);
      (*conn)[1] = node(5);
      (*conn)[2] = node(6);
      (*conn)[3] = node(8);

      return;

      // Linear sub-tet 6
    case 6:
      
      (*conn)[0] = node(5);
      (*conn)[1] = node(9);
      (*conn)[2] = node(6);
      (*conn)[3] = node(8);

      return;

      // Linear sub-tet 7
    case 7:
      
      (*conn)[0] = node(7);
      (*conn)[1] = node(6);
      (*conn)[2] = node(9);
      (*conn)[3] = node(8);

      return;


    default:

      error();
    };
  
  return;
};



#ifdef ENABLE_AMR



#endif
