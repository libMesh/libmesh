// $Id: mesh_data_unv_support.C,v 1.8 2003-08-04 17:23:51 ddreyer Exp $

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
#include <stdio.h>

// Local includes
#include "mesh_data.h"
#include "mesh_base.h"
#include "auto_ptr.h"



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

  const std::string _label_dataset_mesh_data = "2414";

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
       * if beginning of dataset
       */
      if (news == _label_dataset_mesh_data)
        {

	  /**
	   * Now read the data of interest.
	   * Start with the header.
	   */

	  unsigned int Dataset_location;
	      
	  unsigned int Model_type,          
	    Analysis_type,
	    Data_characteristic,
	    Result_type,
	    Data_type,
	    NVALDC;


	  /**
	   * If there is no MeshDataUnvHeader object
	   * attached
	   */
	  if (_unv_header==NULL)
	    {
	      /**
	       * Ignore the first lines that stand for
	       * analysis dataset label and name.
	       */ 
	      for(unsigned int i=0; i<3; i++)
		in_file.ignore(256,'\n');	      
	      
	      /**
	       * Read the dataset location, where
	       * 1: Data at nodes
	       * 2: Data on elements
	       * other sets are currently not supported.
	       */
	      in_file >> Dataset_location;

	      /**
	       * Ignore five ID lines.
	       */ 
	      for(unsigned int i=0; i<6; i++)
		in_file.ignore(256,'\n');
	      
	      /**
	       * Read record 9.
	       */
	      in_file >> Model_type           // not used here
		      >> Analysis_type        // not used here
		      >> Data_characteristic  // not used here
		      >> Result_type          // not used here
		      >> Data_type
		      >> NVALDC;
	      
	      
	      /**
	       * Ignore record 10 and 11
	       * (Integer analysis type specific data).
	       */ 
	      for (unsigned int i=0; i<3; i++)
		in_file.ignore(256,'\n');
	      
	      /**
	       * Ignore record 12 and record 13.  Since there
	       * exist UNV files with 'D' instead of 'e' as 
	       * 10th-power char, it is safer to use a string
	       * to read the dummy reals.
	       */
	      {
	        std::string dummy_Real;
		for (unsigned int i=0; i<12; i++)
		    in_file >> dummy_Real;
	      }

	    }
	  else
	    {
	      _unv_header->read (in_file);

	      Dataset_location = _unv_header->dataset_location;
	      NVALDC = _unv_header->nvaldc;
	      Data_type = _unv_header->data_type;
	    }

	  /**
	   * Check the location of the dataset.
	   */
	  if (Dataset_location != 1)
	    {
	      std::cerr << "ERROR: Currently only Data at nodes is supported." 
			<< std::endl;
	      error();
	    }


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
		   * As again, these floats may also be written
		   * using a 'D' instead of an 'e'.
		   */
		  if (Data_type == 2 || Data_type == 4)
		    {
		      std::string buf;
		      in_file >> buf;
		      MeshDataUnvHeader::need_D_to_e(buf);
		      values[data_cnt] = atof(buf.c_str());
		    }
		  else if(Data_type == 5 || Data_type == 6)
		    {
#ifdef USE_COMPLEX_NUMBERS
		      Real re_val, im_val;

		      std::string buf;
		      in_file >> buf;

		      if (MeshDataUnvHeader::need_D_to_e(buf))
		        {
			  re_val = atof(buf.c_str());
			  in_file >> buf;
			  MeshDataUnvHeader::need_D_to_e(buf);
			  im_val = atof(buf.c_str()); 
			}
		      else
		        {
			  re_val = atof(buf.c_str());
			  in_file >> im_val;
			}

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
	       * Add the values vector to the MeshData data structure.
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

  if (!_elem_data.empty())
      std::cerr << "WARNING: MeshData currently only supports nodal data for Universal files."
		<< std::endl
		<< "         Will proceed writing only nodal data, ignoring element data."
		<< std::endl;


  /**
   * Write the beginning of the dataset.
   */
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_mesh_data
	   << "\n";

  /**
   * Write the header
   */
  if (_unv_header==NULL)
    {
      /*
       * create a header that holds at
       * least sufficient data to specify
       * what this data set currently holds.
       * 
       * The empty constructor automatically
       * takes care of \p dataset_location
       * and \p data_type.
       */
      AutoPtr<MeshDataUnvHeader> my_header(new MeshDataUnvHeader);

      /*
       * It remains to set the correct nvaldc
       */
      my_header->nvaldc = this->n_val_per_node();

      /*
       * write this default header, then let
       * the AutoPtr go out of scope.  This 
       * will take care of memory management.
       */
      my_header->write (out_file);
    }
  else
      _unv_header->write (out_file);


  /**
   * Write the foreign node number and the respective data.
   */
  std::map<const Node*, 
    std::vector<Number> >::const_iterator nit = _node_data.begin();

  char buf[27];
  for (; nit != _node_data.end(); ++nit)
    {
      const Node* node = (*nit).first;

      unsigned int f_n_id = node_to_foreign_id (node);
      sprintf(buf, "%10i\n", f_n_id);
      out_file << buf;

      /* since we are iterating over our own map, this assert
       * should never break...
       */
      assert (this->has_data(node));

      // const reference to the nodal values
      const std::vector<Number>& values = this->get_data(node);

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





//------------------------------------------------------
// MeshDataUnvHeader functions
MeshDataUnvHeader::MeshDataUnvHeader() :
  dataset_label        (0),
  dataset_name         ("libMesh mesh data"),
  dataset_location     (1),  // default to nodal data
  id_line_1            ("libMesh default"),
  id_line_2            ("libMesh default"),
  id_line_3            ("libMesh default"),
  id_line_4            ("libMesh default"),
  id_line_5            ("libMesh default"),
  model_type           (0),          
  analysis_type        (0),
  data_characteristic  (0),
  result_type          (0),
#ifdef USE_COMPLEX_NUMBERS
  data_type            (5),  // default to single precision complex
#else
  data_type            (2),  // default to single precision real
#endif
  nvaldc               (3)   // default to 3 (principal directions)
{
  /**
   * resize analysis specific data.
   */
  record_10.resize(8);
  record_11.resize(2);
  record_12.resize(6);
  record_13.resize(8);
}





MeshDataUnvHeader::~MeshDataUnvHeader()
{
}




void MeshDataUnvHeader::read (std::ifstream& in_file)
{
  in_file >> this->dataset_label;

  in_file.ignore(256,'\n');
  std::getline(in_file, dataset_name, '\n');

  in_file >> this->dataset_location;

  // ID lines are ignored when reading
  for(unsigned int i=0; i<6; i++)
    in_file.ignore(256,'\n');

  in_file >> this->model_type     
	  >> this->analysis_type
	  >> this->data_characteristic
	  >> this->result_type
	  >> this->data_type
	  >> this->nvaldc;
    
  for (unsigned int i=0; i<8; i++)
    in_file >> this->record_10[i];

  for (unsigned int i=0; i<2; i++)
    in_file >> this->record_11[i];


  /*
   * There are UNV-files where floats are 
   * written with 'D' as the 10th-power 
   * character. Replace this 'D' by 'e',
   * so that atof() can work fine.
   */
  std::string buf;
  in_file >> buf;

  if (need_D_to_e(buf))
    {
      // have to convert _all_ 'D' to 'e'
      this->record_12[0] = atof(buf.c_str());

      for (unsigned int i=1; i<6; i++)
        {
	  in_file >> buf;
	  need_D_to_e(buf);
	  this->record_12[i] = atof(buf.c_str());
        }

      for (unsigned int i=0; i<6; i++)
        {
	  in_file >> buf;
	  need_D_to_e(buf);
	  this->record_13[i] = atof(buf.c_str());
        }
    }
  else
    {
      // no 'D', the stream will recognize the floats
      this->record_12[0] = atof(buf.c_str());

      for (unsigned int i=1; i<6; i++)
	  in_file >> this->record_12[i];

      for (unsigned int i=0; i<6; i++)
	  in_file >> this->record_13[i];
    }
}




void MeshDataUnvHeader::write (std::ofstream& out_file)
{
  char buf[81];

  sprintf(buf, "%6i\n",this->dataset_label);
  out_file << buf;
 
  out_file << this->dataset_name << "\n";

  sprintf(buf, "%6i\n",this->dataset_location);
  out_file << buf;

  out_file << this->id_line_1 << "\n";
  out_file << this->id_line_2 << "\n";
  out_file << this->id_line_3 << "\n";
  out_file << this->id_line_4 << "\n";
  out_file << this->id_line_5 << "\n";

  sprintf(buf, "%10i%10i%10i%10i%10i%10i\n",
	  model_type,  analysis_type, data_characteristic,
	  result_type, data_type,     nvaldc);
  out_file << buf;

  sprintf(buf, "%10i%10i%10i%10i%10i%10i%10i%10i\n",
	  record_10[0], record_10[1], record_10[2], record_10[3],
	  record_10[4], record_10[5], record_10[6], record_10[7]);
  out_file << buf;

  sprintf(buf, "%10i%10i\n", record_11[0], record_11[1]);
  out_file << buf;

  sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n",
	  record_12[0], record_12[1], record_12[2],
	  record_12[3], record_12[4], record_12[5]);
  out_file << buf;

  sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n",
	  record_13[0], record_13[1], record_13[2],
	  record_13[3], record_13[4], record_13[5]);
  out_file << buf;

}





bool MeshDataUnvHeader::need_D_to_e (std::string& number)
{
  // position of the "D" in the string
  unsigned int position;

  // find "D" in string, start looking at 6th element, to improve speed.
  // We dont expect a "D" earlier
  position = number.find("D",6);

  if(position!=std::string::npos)     // npos means no position
    {
      // replace "D" in string
      number.replace(position,1,"e"); 
      return true;
    }
  else
      // we assume that if this one number is written correctly, all numbers are
      return false;
}





