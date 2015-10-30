// $Id: mesh_data_unv_support.C,v 1.4 2003-05-22 16:20:52 benkirk Exp $

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
#include <stdio.h>

// Local includes
#include "mesh_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// MeshData UNV support functions
void MeshData::read_unv(const std::string& name)
{
  /*
   * When reading data, make sure the id maps are ok
   */
  assert (_node_id_map_closed);
  assert (_elem_id_map_closed);

  /**
   * clear the data, but keep the id maps
   */
  this->clear();

  std::ifstream in_file (name.c_str());

  if ( !in_file.good() )
    {
      std::cerr << "ERROR: Input file not good." 
		<< std::endl;
      error();
    }

  const std::string  _label_dataset_mesh_data    = "2414";

  /**
   * locate the beginning of data set
   * and read it.
   */
  {
    std::string olds, news;

    while (true)
    {
      in_file >> olds >> news;

      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while( ((olds != "-1") || (news == "-1") ) && !in_file.eof() )
	{
	  olds = news;
	  in_file >> news;
	}

      if(in_file.eof())
	break;

      /*
       * if beginning of dataset, buffer it in
       * temp_buffer, if desired
       */
      if (news == _label_dataset_mesh_data)
        {

	  /**
	   * Now read the data of interest.
	   */

	  /**
	   * Ignore the first lines that stand for
	   * analysis dataset label and name.
	   */ 
	  for(unsigned int i=0; i<3; i++)
	    in_file.ignore(256,'\n');

	  unsigned int Dataset_location;

	  unsigned int Model_type,          
	               Analysis_type,
	               Data_characteristic,
	               Result_type,
	               Data_type,
	               NVALDC;

	  Real dummy_Real;

	  /**
	   * Read the dataset location, where
	   * 1: Data at nodes
	   * 2: Data on elements
	   * other sets are currently not supported.
	   */
	  in_file >> Dataset_location;


	  if (Dataset_location != 1)
	    {
	      std::cerr << "ERROR: Currently only Data at nodes is supported." 
			<< std::endl;
	      error();
	    }

	  /**
	   * Ignore five ID lines.
	   */ 
	  for(unsigned int i=0; i<6; i++)
	    in_file.ignore(256,'\n');

	  /**
	   * Read record 9.
	   */
	  in_file >> Model_type           // not supported yet
		  >> Analysis_type        // not supported yet
		  >> Data_characteristic  // not supported yet
		  >> Result_type          // not supported yet
		  >> Data_type
		  >> NVALDC;


	  /**
	   * Ignore record 10 and 11
	   * (Integer analysis type specific data).
	   */ 
	  for (unsigned int i=0; i<3; i++)
	    in_file.ignore(256,'\n');

	  /**
	   * Ignore record 12 and record 13.
	   */
	  for (unsigned int i=0; i<12; i++)
	      in_file >> dummy_Real;

	  /**
	   * Now get the foreign node id number and the respective nodal data.
	   */
	  int f_n_id;
	  std::vector<Number> values;

	  while(true)
	    {
	      in_file >> f_n_id;
	  
	      /**
	       * if node_nr = -1 then we have reached the end of the dataset.
	       */
	      if (f_n_id==-1)
		  break;

	      /**
	       * Resize the values vector (usually data in three
	       * principal directions, i.e. NVALDC = 3).
	       */
	      values.resize(NVALDC);
	  
	      /**
	       * Read the meshdata for the respective node.
	       */	       
	      for (unsigned int data_cnt=0; data_cnt<NVALDC; data_cnt++)
		{
		  /**
		   * Check what data type we are reading.
		   * 2,4: Real
		   * 5,6: Complex
		   * other data types are not supported yet.
		   */
		  if (Data_type == 2 || Data_type == 4)
		    in_file >> values[data_cnt];

		  else if(Data_type == 5 || Data_type == 6)
		    {
#ifdef USE_COMPLEX_NUMBERS
		      
		      Real re_val, im_val;

		      in_file >> re_val
			      >> im_val;

		      values[data_cnt] = Complex(re_val,im_val);
#else

		      std::cerr << "ERROR: Complex data only supported" << std::endl
				<< "when libMesh is configured with --enable-complex!"
				<< std::endl;
		      error();
#endif
		    }

		  else
		    {
		      std::cerr << "ERROR: Data type not supported." 
				<< std::endl;
		      error();
		    }

		} // end loop data_cnt

	      /**
	       * Add the values vector to the MechData data structure.
	       */
	      const Node* node = foreign_id_to_node(f_n_id);
	      _node_data.insert (std::make_pair(node, values));

	    } // while(true)
	}


      else
        {
	  /**
	   * all other datasets are ignored
	   */
        }

    }
  }


  /*
   * finished reading.  Ready for use
   */
  this->_node_data_closed = true;
  this->_elem_data_closed = true;
}




void MeshData::write_unv (const std::string& name)
{
  /*
   * make sure the id maps are ready
   * and that we have data to write
   */
  assert (_node_id_map_closed);
  assert (_elem_id_map_closed);

  assert (_node_data_closed);
  assert (_elem_data_closed);

  std::ofstream out_file (name.c_str());

  const std::string  _label_dataset_mesh_data    = "2414";

  /**
   * Currently this function handles only nodal data.
   */
  assert (!_node_data.empty());

  /**
   * Write the beginning of the dataset.
   */
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_mesh_data
	   << "\n";

  /**
   * Record 1 (analysis dataset label, currently just a dummy).
   */
  char buf[78];
  sprintf(buf, "%10i\n",1);
  out_file << buf;

  /**
   * Record 2 (analysis dataset name, currently just a dummy).
   */
  std::string _dummy_string = "libMesh\n";
  out_file << _dummy_string;

  /**
   * Record 3 (dataset location).
   */
  sprintf(buf, "%10i\n",1);
  out_file << buf;

  /**
   * Record 4 trough 8 (ID lines).
   */
  for (unsigned int i=0; i<5;i++)
    out_file << _dummy_string;

  /**
   * Record 9 trough 11.
   */
  sprintf(buf, "%10i%10i%10i%10i%10i%10i\n", 0, 0, 0, 0, 0, 0);
  out_file << buf;

  sprintf(buf, "%10i%10i%10i%10i%10i%10i%10i%10i\n", 0, 0, 0, 0, 0, 0, 0, 0);
  out_file << buf;

  sprintf(buf, "%10i%10i\n", 0, 0);
  out_file << buf;

  /**
   * Record 12 and 13.
   */
  sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n", 0., 0., 0., 0., 0., 0.);
  out_file << buf;

  sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n", 0., 0., 0., 0., 0., 0.);
  out_file << buf;

  /**
   * Write the foreign node number and the respective data.
   */
  std::map<const Node*, 
    std::vector<Number> >::const_iterator nit = _node_data.begin();

  for (; nit != _node_data.end(); ++nit)
    {
      const Node* node = (*nit).first;

      unsigned int f_n_id = node_to_foreign_id (node);
      sprintf(buf, "%10i\n", f_n_id);
      out_file << buf;

      std::vector<Number> values;
      this->operator()(node, values);

      for (unsigned int v_cnt=0; v_cnt<values.size(); v_cnt++)
	{
#ifdef USE_COMPLEX_NUMBERS
	  sprintf(buf, "%13.5E%13.5E", values[v_cnt].real(),
		                       values[v_cnt].imag());
	  out_file << buf;
#else
	  sprintf(buf, "%13.5E", values[v_cnt]);
	  out_file << buf;
#endif
	}

      out_file << "\n";


    }

  /**
   * Write end of the dataset.
   */
  out_file << "    -1\n";


}

