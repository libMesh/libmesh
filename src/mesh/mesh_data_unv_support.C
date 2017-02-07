// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cstdio> // for std::sprintf
#include <fstream>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh_data.h"
#include "libmesh/auto_ptr.h"

#ifdef  LIBMESH_HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif


namespace libMesh
{

//------------------------------------------------------
// MeshData UNV support functions
void MeshData::read_unv (const std::string & file_name)
{
  /*
   * We better be active or in compatibility mode
   */
  libmesh_assert (this->_active || this->_compatibility_mode);

  /*
   * When reading data, make sure the id maps are ok
   */
  libmesh_assert (this->_node_id_map_closed);
  libmesh_assert (this->_elem_id_map_closed);

  /*
   * clear the data, but keep the id maps
   */
  this->clear();

  /*
   * We can read either ".unv", or ".unv.gz"
   * files, provided zlib.h is there
   */
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef LIBMESH_HAVE_GZSTREAM
      igzstream in_stream(file_name.c_str());
      this->read_unv_implementation (in_stream);
#else
      libmesh_error_msg("ERROR:  You must have the zlib.h header files and libraries to read and write compressed streams.");
#endif
      return;
    }

  else
    {
      std::ifstream in_stream(file_name.c_str());
      this->read_unv_implementation (in_stream);
      return;
    }
}






void MeshData::read_unv_implementation (std::istream & in_file)
{
  /*
   * This is the actual implementation of
   * reading in UNV format.  This enables
   * to read either through the conventional
   * C++ stream, or through a stream that
   * allows to read .gz'ed files.
   */
  if ( !in_file.good() )
    libmesh_error_msg("ERROR: Input file not good.");

  const std::string _label_dataset_mesh_data = "2414";

  /*
   * locate the beginning of data set
   * and read it.
   */
  {
    std::string olds, news;

    while (true)
      {
        in_file >> olds >> news;

        /*
         * Yes, really dirty:
         *
         * When we found a dataset, and the user does
         * not want this dataset, we jump back here
         */
      go_and_find_the_next_dataset:

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

            /*
             * Now read the data of interest.
             * Start with the header.  For
             * explanation of the variable
             * dataset_location, see below.
             */
            unsigned int dataset_location;

            /*
             * the type of data (complex, real,
             * float, double etc, see below)
             */
            unsigned int data_type;

            /*
             * the number of floating-point values per entity
             */
            unsigned int NVALDC;


            /*
             * If there is no MeshDataUnvHeader object
             * attached
             */
            if (_unv_header == libmesh_nullptr)
              {
                /*
                 * Ignore the first lines that stand for
                 * analysis dataset label and name.
                 */
                for(unsigned int i=0; i<3; i++)
                  in_file.ignore(256,'\n');

                /*
                 * Read the dataset location, where
                 * 1: Data at nodes
                 * 2: Data on elements
                 * other sets are currently not supported.
                 */
                in_file >> dataset_location;

                /*
                 * Ignore five ID lines.
                 */
                for(unsigned int i=0; i<6; i++)
                  in_file.ignore(256,'\n');

                /*
                 * These data are all of no interest to us...
                 */
                unsigned int model_type,
                  analysis_type,
                  data_characteristic,
                  result_type;

                /*
                 * Read record 9.
                 */
                in_file >> model_type           // not used here
                        >> analysis_type        // not used here
                        >> data_characteristic  // not used here
                        >> result_type          // not used here
                        >> data_type
                        >> NVALDC;


                /*
                 * Ignore record 10 and 11
                 * (Integer analysis type specific data).
                 */
                for (unsigned int i=0; i<3; i++)
                  in_file.ignore(256,'\n');

                /*
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

                /*
                 * the read() method returns false when
                 * the user wanted a special header, and
                 * when the current header is _not_ the correct
                 * header
                 */
                if (_unv_header->read(in_file))
                  {
                    dataset_location = _unv_header->dataset_location;
                    NVALDC = _unv_header->nvaldc;
                    data_type = _unv_header->data_type;
                  }
                else
                  {
                    /*
                     * This is not the correct header.  Go
                     * and find the next.  For this to
                     * work correctly, shift to the
                     * next line, so that the "-1"
                     * disappears from olds
                     */
                    olds = news;
                    in_file >> news;

                    /*
                     * No good style, i know...
                     */
                    goto go_and_find_the_next_dataset;
                  }

              }

            /*
             * Check the location of the dataset.
             */
            if (dataset_location != 1)
              libmesh_error_msg("ERROR: Currently only Data at nodes is supported.");


            /*
             * Now get the foreign node id number and the respective nodal data.
             */
            int f_n_id;
            std::vector<Number> values;

            while(true)
              {
                in_file >> f_n_id;

                /*
                 * if node_nr = -1 then we have reached the end of the dataset.
                 */
                if (f_n_id==-1)
                  break;

                /*
                 * Resize the values vector (usually data in three
                 * principle directions, i.e. NVALDC = 3).
                 */
                values.resize(NVALDC);

                /*
                 * Read the meshdata for the respective node.
                 */
                for (unsigned int data_cnt=0; data_cnt<NVALDC; data_cnt++)
                  {
                    /*
                     * Check what data type we are reading.
                     * 2,4: Real
                     * 5,6: Complex
                     * other data types are not supported yet.
                     * As again, these floats may also be written
                     * using a 'D' instead of an 'e'.
                     */
                    if (data_type == 2 || data_type == 4)
                      {
                        std::string buf;
                        in_file >> buf;
                        MeshDataUnvHeader::need_D_to_e(buf);
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                        values[data_cnt] = Complex(std::atof(buf.c_str()), 0.);
#else
                        values[data_cnt] = std::atof(buf.c_str());
#endif
                      }

                    else if(data_type == 5 || data_type == 6)

                      {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                        Real re_val, im_val;

                        std::string buf;
                        in_file >> buf;

                        if (MeshDataUnvHeader::need_D_to_e(buf))
                          {
                            re_val = std::atof(buf.c_str());
                            in_file >> buf;
                            MeshDataUnvHeader::need_D_to_e(buf);
                            im_val = std::atof(buf.c_str());
                          }
                        else
                          {
                            re_val = std::atof(buf.c_str());
                            in_file >> im_val;
                          }

                        values[data_cnt] = Complex(re_val,im_val);
#else

                        libmesh_error_msg("ERROR: Complex data only supported when libMesh is configured with --enable-complex!");
#endif
                      }

                    else
                      libmesh_error_msg("ERROR: Data type not supported.");

                  } // end loop data_cnt

                /*
                 * Add the values vector to the MeshData data structure.
                 */
                const Node * node = foreign_id_to_node(f_n_id);
                _node_data.insert (std::make_pair(node, values));

              } // while(true)
          }


        else
          {
            /*
             * all other datasets are ignored
             */
          }

      }
  }


  /*
   * finished reading.  Ready for use, provided
   * there was any data contained in the file.
   */
  libmesh_assert ((this->_node_data.size() != 0) || (this->_elem_data.size() != 0));

  this->_node_data_closed = true;
  this->_elem_data_closed = true;
}






void MeshData::write_unv (const std::string & file_name)
{
  /*
   * We better be active or in compatibility mode
   */
  libmesh_assert (this->_active || this->_compatibility_mode);

  /*
   * make sure the id maps are ready
   * and that we have data to write
   */
  libmesh_assert (this->_node_id_map_closed);
  libmesh_assert (this->_elem_id_map_closed);

  libmesh_assert (this->_node_data_closed);
  libmesh_assert (this->_elem_data_closed);

  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef LIBMESH_HAVE_GZSTREAM
      ogzstream out_stream(file_name.c_str());
      this->write_unv_implementation (out_stream);
#else
      libmesh_error_msg("ERROR:  You must have the zlib.h header files and libraries to read and write compressed streams.");
#endif
      return;

    }

  else
    {
      std::ofstream out_stream(file_name.c_str());
      this->write_unv_implementation (out_stream);
      return;
    }
}






void MeshData::write_unv_implementation (std::ostream & out_file)
{
  /*
   * This is the actual implementation of writing
   * unv files, either as .unv or as .unv.gz file
   */
  if ( !out_file.good() )
    libmesh_error_msg("ERROR: Output file not good.");


  /*
   * the beginning marker of the dataset block for
   * nodal/element-associated data (not to be confused
   * with _desired_dataset_label!)
   */
  const std::string  _label_dataset_mesh_data    = "2414";

  /*
   * Currently this function handles only nodal data.
   */
  libmesh_assert (!_node_data.empty());

  if (!_elem_data.empty())
    libMesh::err << "WARNING: MeshData currently only supports nodal data for Universal files."
                 << std::endl
                 << "         Will proceed writing only nodal data, ignoring element data."
                 << std::endl;


  /*
   * Write the beginning of the dataset.
   */
  out_file << "    -1\n"
           << "  "
           << _label_dataset_mesh_data
           << "\n";

  /*
   * Write the header
   */
  if (_unv_header == libmesh_nullptr)
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
      MeshDataUnvHeader my_header;

      /*
       * It remains to set the correct nvaldc...
       */
      my_header.nvaldc = this->n_val_per_node();

      /*
       * and the correct data type.  By default
       * only distinguish complex or real data.
       */
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      my_header.data_type = 5;
#else
      my_header.data_type = 2;
#endif

      /*
       * write this default header, then let
       * the UniquePtr go out of scope.  This
       * will take care of memory management.
       */
      my_header.write (out_file);
    }

  else
    {
      /*
       * make sure our nvaldc coincide.
       */
      if (this->n_val_per_node() != _unv_header->nvaldc)
        {
          libMesh::err << "WARNING: nvaldc=" << _unv_header->nvaldc
                       << " of attached MeshDataUnvHeader object not valid!" << std::endl
                       << "         re-set nvaldc to " << this->n_val_per_node() << std::endl;
          _unv_header->nvaldc = this->n_val_per_node();
        }


      /*
       * only issue a warning when data_type does
       * not coincide.  Perhaps user provided some
       * other header in order to convert complex
       * to real...
       */
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      const unsigned int my_data_type = 5;
#else
      const unsigned int my_data_type = 2;
#endif
      if (my_data_type != _unv_header->data_type)
        {
          libMesh::err << "WARNING: data_type=" << _unv_header->data_type
                       << " of attached MeshDataUnvHeader differs from" << std::endl
                       << "         default value=" << my_data_type
                       << " Perhaps the user wanted this," << std::endl
                       << "         so I use the value from the MeshDataUnvHeader."
                       << std::endl;
        }
      _unv_header->write (out_file);
    }


  /*
   * Write the foreign node number and the respective data.
   */
  std::map<const Node *,
           std::vector<Number> >::const_iterator nit = _node_data.begin();

  char buf[27];
  for (; nit != _node_data.end(); ++nit)
    {
      const Node * node = (*nit).first;

      unsigned int f_n_id = node_to_foreign_id (node);
      std::sprintf(buf, "%10u\n", f_n_id);
      out_file << buf;

      /* since we are iterating over our own map, this libmesh_assert
       * should never break...
       */
      libmesh_assert (this->has_data(node));

      // const reference to the nodal values
      const std::vector<Number> & values = this->get_data(node);

      for (std::size_t v_cnt=0; v_cnt<values.size(); v_cnt++)
        {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
          std::sprintf(buf, "%13.5E%13.5E", values[v_cnt].real(),
                       values[v_cnt].imag());
          out_file << buf;
#else
          std::sprintf(buf, "%13.5E",
                       static_cast<double>(values[v_cnt]));
          out_file << buf;
#endif
        }

      out_file << "\n";


    }

  /*
   * Write end of the dataset.
   */
  out_file << "    -1\n";
}





//------------------------------------------------------
// MeshDataUnvHeader functions
MeshDataUnvHeader::MeshDataUnvHeader() :
  dataset_label          (0),
  dataset_name           ("libMesh mesh data"),
  dataset_location       (1),  // default to nodal data
  model_type             (0),
  analysis_type          (0),
  data_characteristic    (0),
  result_type            (0),
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  data_type              (5),  // default to single precision complex
#else
  data_type              (2),  // default to single precision real
#endif
  nvaldc                 (3),  // default to 3 (principle directions)
  _desired_dataset_label (libMesh::invalid_uint)
{
  id_lines_1_to_5.resize(5);
  std::fill (id_lines_1_to_5.begin(), id_lines_1_to_5.end(), std::string("libMesh default"));
  /*
   * resize analysis specific data.
   */
  record_10.resize(8);
  record_11.resize(2);
  record_12.resize(6);
  record_13.resize(6);
}





MeshDataUnvHeader::~MeshDataUnvHeader()
{
  // empty
}




bool MeshDataUnvHeader::read (std::istream & in_file)
{
  in_file >> this->dataset_label;

  /*
   * currently, we compare only the
   * dataset_label with the _desired_dataset_label,
   * but it may be easy to also compare the
   * dataset_name.
   *
   * When the user provided a dataset label, and
   * the current label does _not_ match, then just
   * return false.
   *
   * Otherwise: when the current label matches,
   * or when there is no desired dataset label,
   * simply proceed.
   */
  if ((this->_desired_dataset_label != libMesh::invalid_uint) &&
      (this->dataset_label != this->_desired_dataset_label))
    return false;


  in_file.ignore(256,'\n');
  std::getline(in_file, dataset_name, '\n');

  in_file >> this->dataset_location;
  in_file.ignore(256,'\n');


  for (unsigned int n=0; n<5; n++)
    std::getline(in_file, this->id_lines_1_to_5[n], '\n');


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
   * so that std::atof() can work fine.
   */
  std::string buf;
  in_file >> buf;

  if (need_D_to_e(buf))
    {
      // have to convert _all_ 'D' to 'e'
      this->record_12[0] = std::atof(buf.c_str());

      for (unsigned int i=1; i<6; i++)
        {
          in_file >> buf;
          need_D_to_e(buf);
          this->record_12[i] = std::atof(buf.c_str());
        }

      for (unsigned int i=0; i<6; i++)
        {
          in_file >> buf;
          need_D_to_e(buf);
          this->record_13[i] = std::atof(buf.c_str());
        }
    }
  else
    {
      // no 'D', the stream will recognize the floats
      this->record_12[0] = std::atof(buf.c_str());

      for (unsigned int i=1; i<6; i++)
        in_file >> this->record_12[i];

      for (unsigned int i=0; i<6; i++)
        in_file >> this->record_13[i];
    }

  /*
   * no matter whether the user provided a desired
   * dataset label or not: return true, b/c the
   * non-match was already caught before.
   */
  return true;
}




void MeshDataUnvHeader::write (std::ostream & out_file)
{


  char buf[82];

  std::sprintf(buf, "%6i\n",this->dataset_label);

  out_file << buf;

  out_file << this->dataset_name << "\n";

  std::sprintf(buf, "%6i\n",this->dataset_location);

  out_file << buf;

  for (unsigned int n=0; n<5; n++)
    out_file << this->id_lines_1_to_5[n] << "\n";

  std::sprintf(buf, "%10i%10i%10i%10i%10i%10i\n",
               model_type,  analysis_type, data_characteristic,
               result_type, data_type,     nvaldc);

  out_file << buf;

  std::sprintf(buf, "%10i%10i%10i%10i%10i%10i%10i%10i\n",
               record_10[0], record_10[1], record_10[2], record_10[3],
               record_10[4], record_10[5], record_10[6], record_10[7]);

  out_file << buf;

  std::sprintf(buf, "%10i%10i\n", record_11[0], record_11[1]);
  out_file << buf;

  std::sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n",
               static_cast<double>(record_12[0]),
               static_cast<double>(record_12[1]),
               static_cast<double>(record_12[2]),
               static_cast<double>(record_12[3]),
               static_cast<double>(record_12[4]),
               static_cast<double>(record_12[5]));

  out_file << buf;

  std::sprintf(buf, "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n",
               static_cast<double>(record_13[0]),
               static_cast<double>(record_13[1]),
               static_cast<double>(record_13[2]),
               static_cast<double>(record_13[3]),
               static_cast<double>(record_13[4]),
               static_cast<double>(record_13[5]));

  out_file << buf;
}





bool MeshDataUnvHeader::need_D_to_e (std::string & number)
{
  // find "D" in string, start looking at 6th element, to improve speed.
  // We dont expect a "D" earlier
  std::string::size_type position = number.find("D",6);

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



void MeshDataUnvHeader::which_dataset (const unsigned int ds_label)
{
  this->_desired_dataset_label = ds_label;
}



void MeshDataUnvHeader::operator = (const MeshDataUnvHeader & omduh)
{
  this->dataset_label          = omduh.dataset_label;
  this->dataset_name           = omduh.dataset_name;
  this->dataset_location       = omduh.dataset_location;
  this->id_lines_1_to_5        = omduh.id_lines_1_to_5;

  this->model_type             = omduh.model_type;
  this->analysis_type          = omduh.analysis_type;
  this->data_characteristic    = omduh.data_characteristic;
  this->result_type            = omduh.result_type;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  /*
   * in complex mode allow only
   * values 5 or 6 (complex) for data_type
   */
  if ((omduh.data_type == 5) ||
      (omduh.data_type == 6))
    this->data_type          = omduh.data_type;
  else
    {
#  ifdef DEBUG
      libMesh::err << "WARNING: MeshDataUnvHeader::operator=(): Other object has data_type for" << std::endl
                   << "         real values.  Will use default data_type=5 during assignment." << std::endl
                   << std::endl;
#  endif
      this->data_type          = 5;
    }

#else

  /*
   * in real mode allow only
   * values 2 or 4 (real) for data_type
   */
  if ((omduh.data_type == 2) ||
      (omduh.data_type == 4))
    this->data_type          = omduh.data_type;
  else
    {
#  ifdef DEBUG
      libMesh::err << "WARNING: Other MeshDataUnvHeader has data_type for complex values." << std::endl
                   << "         Data import will likely _not_ work and result in infinite loop," << std::endl
                   << "         provided the user forgot to re-size nvaldc to 2*nvaldc_old!" << std::endl
                   << std::endl;
#  endif
      this->data_type          = 2;
    }

#endif

  this->nvaldc                 = omduh.nvaldc;

  this->record_10              = omduh.record_10;
  this->record_11              = omduh.record_11;
  this->record_12              = omduh.record_12;
  this->record_13              = omduh.record_13;

  this->_desired_dataset_label = omduh._desired_dataset_label;
}




bool MeshDataUnvHeader::operator == (const MeshDataUnvHeader & omduh) const
{
  return (this->dataset_label          == omduh.dataset_label       &&
          this->dataset_name           == omduh.dataset_name        &&
          this->dataset_location       == omduh.dataset_location    &&
          this->id_lines_1_to_5        == omduh.id_lines_1_to_5     &&

          this->model_type             == omduh.model_type          &&
          this->analysis_type          == omduh.analysis_type       &&
          this->data_characteristic    == omduh.data_characteristic &&
          this->result_type            == omduh.result_type         &&

          this->data_type              == omduh.data_type           &&
          this->nvaldc                 == omduh.nvaldc              &&

          this->record_10              == omduh.record_10           &&
          this->record_11              == omduh.record_11           &&
          this->record_12              == omduh.record_12           &&
          this->record_13              == omduh.record_13           &&

          this->_desired_dataset_label == omduh._desired_dataset_label);
}

} // namespace libMesh
