// $Id: mesh.C,v 1.14 2003-05-15 23:34:35 benkirk Exp $

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
#include "libmesh.h"
#include "mesh_logging.h"


#ifdef HAVE_SFCURVES
// prototype for SFC code
namespace sfc {
  extern "C" {
#include "sfcurves.h"
  }
}
#endif



// ------------------------------------------------------------
// Mesh class member functions
Mesh::Mesh (unsigned int d,
	    unsigned int pid) :
  MeshBase           (d, pid),
  boundary_mesh      (d-1, pid)
{
  assert (libMesh::initialized());
}



Mesh::~Mesh ()
{
  clear ();
  
  assert (!libMesh::closed());
}



void Mesh::clear ()
{
  // Clear other data structures
  boundary_mesh.clear();

  MeshBase::clear();
}



void Mesh::read (const std::string& name)
{
  START_LOG("read()", "Mesh");

  
  // Read the file based on extension
  {
    if (name.rfind(".mat") < name.size())
      read_matlab (name);
    
    else if (name.rfind(".ucd") < name.size())
      read_ucd (name);

    else if (name.rfind(".exd") < name.size())
      read_exd (name);

    else if (name.rfind(".xda") < name.size())
      read_xdr (name);

    else if ((name.rfind(".off")  < name.size()) ||
	     (name.rfind(".ogl")  < name.size()) ||
	     (name.rfind(".oogl") < name.size()))
      read_off(name);

    else if ((name.rfind(".xdr")  < name.size()) ||
	     (name.rfind(".0000") < name.size()))
      read_xdr_binary (name);

    else if (name.rfind(".mesh") < name.size())
      read_shanee (name);

    else if (name.rfind(".unv") < name.size())
      read_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.mat  -- Matlab triangular ASCII file\n"
		  << "     *.ucd  -- AVS's ASCII UCD format\n"
		  << "     *.mesh -- Ben's \"shanee\" format\n"
		  << "     *.off  -- OOGL OFF surface format\n"
		  << "     *.exd  -- Sandia's ExodusII format\n"
		  << "     *.xda  -- Internal ASCII format\n"
		  << "     *.xdr  -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << "     *.unv  -- I-deas format\n"
		  << std::endl;
	error();

      }    
  }
  STOP_LOG("read()", "Mesh");


  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use();
}



void Mesh::write (const std::string& name)
{
  START_LOG("write()", "Mesh");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name);

    else if (name.rfind(".ucd") < name.size())
      write_ucd (name);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, NULL, NULL, true);
	else
	  write_gmv_binary(name);
      }


    else if (name.rfind(".ugrid") < name.size())
      write_diva (name);
    
    else if (name.rfind(".xda") < name.size())
      write_xdr (name);
    
    else if (name.rfind(".xdr") < name.size())
      write_xdr_binary (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat   -- Tecplot ASCII file\n"
		  << "     *.plt   -- Tecplot binary file\n"
		  << "     *.ucd   -- AVS's ASCII UCD format\n"
		  << "     *.ugrid -- Kelly's DIVA ASCII format\n"
		  << "     *.gmv   -- LANL's GMV (General Mesh Viewer) format\n"
		  << "     *.xda   -- Internal ASCII format\n"
		  << "     *.xdr   -- Internal binary format,\n"
		  << "               compatible with XdrMGF\n"
		  << std::endl
		  << "\n Exiting without writing output\n";
      }    
  }
  
  STOP_LOG("write()", "Mesh");
}



void Mesh::write (const std::string& name,
		  std::vector<Number>& v,
		  std::vector<std::string>& vn)
{
  START_LOG("write()", "Mesh");

  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name, &v, &vn);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name, &v, &vn);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, &v, &vn, true);
	else
	  write_gmv_binary(name, &v, &vn);
      }    
    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat  -- Tecplot ASCII file\n"
		  << "     *.plt  -- Tecplot binary file\n"
		  << "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
		  << "\n Exiting without writing output\n";
      }
  }

  STOP_LOG("write()", "Mesh");
}
