// $Id: boundary_data_unv_support.C,v 1.1 2003-05-14 11:54:37 ddreyer Exp $

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
#include <fstream>

// Local includes
#include "boundary_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// BoundaryData UNV support functions
void BoundaryData::read_unv(const std::string& name)
{
  std::ifstream file (name.c_str());

  static const string label_unv_dataset   = "2414";

  // check file
  if ( !file.good() )
  {
    std::cerr << "ERROR: Input file not good!" << std::endl;
    error();
  }


  /*
   * locate & read the data set
   */
  {
    std::string olds,news;

    bool got_what_i_wanted = false;
    bool finished          = false;
    
    while (!finished)   
      {
        file >> olds >> news;

	// a "-1" followed by a number means the beginning of a dataset
	while( ((olds != "-1") || (news == "-1") )
	       &&
	       !file.eof() )
	  {
	    // proceed 
	    olds = news;  
	    file >> news;
	  }

	// end of file reached?
	finished = file.eof();

	/*
	 * we have the label of a dataset, is it the
	 * one we want?
	 */
	if (news == label_unv_dataset)
	  {
	    /* Determine from the header:
	     *
	     * - the entity for which data is stored (nodes/elements),
	     * - how many Real's are stored per entity (distinguish
	     *   between Number=Real and Number=Complex:  For complex,
	     *   read two Real's, store them
	     * - for how many entities is data stored?
	     */

//TODO:[SP] Continue...


	    // resize the local storage _node_data/_elem_data



	    // when reading nodal data, make sure the id map is ok
	    assert (_node_id_map_closed);


	    // when reading nodal data, make sure the id map is ok
	    assert (_elem_id_map_closed);



            /*
	     * finished reading the data, so 
	     * quit these loops
	     */
	    got_what_i_wanted = true;
	    finished = true;
	  }
      }

    if (!got_what_i_wanted)
      {
	std::cerr << "ERROR: Data set for boundary conditions not found"
		  << std::endl
		  << " in file: " << name
		  << std::endl;
	error();
      }
  }

  file.close();
  return;
}




void BoundaryData::write_unv (const std::string& /* name */)
{
  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}

