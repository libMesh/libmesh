// $Id: plt_loader_read.C,v 1.1 2004-09-30 21:21:49 benkirk Exp $

// Copyright (C) 2002-2004  Benjamin S. Kirk
  
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



#include <fstream>
#include "utility.h"
#include "plt_loader.h"



//-----------------------------------------------------------------------------
// PltLoader reading members
void PltLoader::read (const std::string& name)
{
  std::ifstream in (name.c_str(), std::ios::in|std::ios::binary);

  if (!in.good())
    {
      std::cerr << "Error reading input file " << name
		<< std::endl;

      abort();
    }
    

  if (this->verbose())
    std::cout << std::endl
	      << "Reading input file " << name
	      << std::endl
	      << "-------------------------------------------------------------------------"
	      << std::endl;

  this->read_header (in);
  this->read_data   (in);

  if (this->verbose())
    std::cout << std::endl
	      << "-------------------------------------------------------------------------"
	      << std::endl;

}



void PltLoader::read_header (std::istream& in)
{
  assert (in.good());
  
  //----------------------------------------------------
  // Read the TECPLOT header
  
  // Read the version number
  {
    in.read (buf, 8);

    this->version().clear();

    for (unsigned int i=0; i<8; i++)
      this->version() += buf[i];
    
    if (this->verbose())
      std::cout << "Tecplot Version: "
		<< this->version()
		<< std::endl;
  }


  //----------------------------------------------------
  // Read plt files written by older versions of Tecplot
  if (this->version().rfind("V7") < this->version().size())
    {
      if (this->verbose())
	std::cout << "Reading legacy .plt format (<= v9) ..."
		  << std::endl;
      
      // Read the value of 1 to determine byte ordering
      {
	int one = 0;	
	in.read (buf, SIZEOF_INT);
	memcpy  (&one, buf, SIZEOF_INT);
	
	if (one != 1)
	  {
	    if (this->verbose())
	      std::cout << "Tecplot data is Foreign!"
			<< std::endl;
	    this->is_foreign() = true;

	    // Make sure one reversed is one
	    Utility::ReverseBytes rb(this->is_foreign());
	    assert (rb(one) == 1);
	  }
      }

      // A byte-reverser in case the data is foreign
      Utility::ReverseBytes rb(this->is_foreign());

      // Read the title
      {
	int i=0;
	this->title().clear();
	
	do
	  {
	    in.read (buf, SIZEOF_INT);
	    memcpy  (&i, buf, SIZEOF_INT);
	    rb(i);

	    // Don't add trailing \0
	    if (i)
	      this->title() += static_cast<char>(i);
	  }
	while (i);
      }

      // Read the number of variables in the data set
      {
	int nv;
	in.read (buf, SIZEOF_INT);
	memcpy  (&nv, buf, SIZEOF_INT);
	rb(nv);
	
	this->set_n_vars (nv);
      }

      // Read the variable names
      for (unsigned int v=0; v<this->n_vars(); v++)
	{
	  int i=0;
	  this->var_name(v).clear();
      
	  do
	    {
	      in.read (buf, SIZEOF_INT);
	      memcpy  (&i, buf, SIZEOF_INT);
	      rb(i);
	      
	      // Don't add trailing \0
	      if (i)
		this->var_name(v) += static_cast<char>(i);
	    }
	  while (i);
	}

  
      
      // Read zones from the header.
      // Continue reading until the end-of-header
      // marker (357.) is found.
      int nz=0;
      std::vector<std::string> zname;
      std::vector<int>         ztype, zimax, zjmax, zkmax;
      
      {
	float f=0.;
    
	do
	  {
	    // find the next Zone marker	
	    do
	      {
		f = 0.;
		in.read (buf, SIZEOF_FLOAT);
		memcpy  (&f, buf, SIZEOF_FLOAT);
		rb(f);
	      }
	    while ((f != 299.) &&
		   (f != 357.) &&
		   in.good());
	  
	  
	    // Did we overrun the file?
	    if (!in.good())
	      {
		std::cerr << "ERROR: Unexpected end-of-file!"
			  << std::endl;
		abort();
	      }    
	  
	    // Found a Zone marker
	    else if (f == 299.) 
	      {
		// Incriment the Zone counter
		nz++;
	  
		// Read the zone name
		{
		  int i=0;
		  std::string name;
	    
		  do
		    {
		      in.read (buf, SIZEOF_INT);
		      memcpy  (&i, buf, SIZEOF_INT);
		      rb(i);
		      
		      // Don't add trailing \0
		      if (i)
			name += static_cast<char>(i);
		    }
		  while (i);
	    
		  zname.push_back(name);
		}
	      
		// Read the zone format
		{
		  int zt;
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&zt, buf, SIZEOF_INT);
		  rb(zt);
		  
		  ztype.push_back(zt);
		  //std::cout << "zone type=" << ztype.back() << std::endl;
		}
	      
		// Read the zone color
		{
		  int zc=0;
	    
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&zc, buf, SIZEOF_INT);
		  rb(zc);
		  
		  //std::cout << "zone color=" << zc << std::endl;
		}
	  
		// Read in the block dimensions
		{
		  int
		    imax=0,
		    jmax=0,
		    kmax=0;
	    
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&imax, buf, SIZEOF_INT);
		  rb(imax);
		  
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&jmax, buf, SIZEOF_INT);
		  rb(jmax);
		  
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&kmax, buf, SIZEOF_INT);
		  rb(kmax);
		  
		  zimax.push_back (imax);
		  zjmax.push_back (jmax);
		  zkmax.push_back (kmax);
		}
	      } // else if (f == 299.)
	  }
	while ((f != 357.) && in.good());
      }

      // Set the header data
      this->set_n_zones (nz);

      for (unsigned int z=0; z<this->n_zones(); z++)
	{
	  this->zone_type(z) = ztype[z];
	  this->zone_name(z) = zname[z];
	  this->imax(z)      = zimax[z];
	  this->jmax(z)      = zjmax[z];
	  this->kmax(z)      = zkmax[z];
	}
    }
  

  //----------------------------------------------------
  // Read plt files written by newer versions of Tecplot
  else if (this->version().rfind("V1") < this->version().size())
    {
      if (this->verbose())
	std::cout << "Reading new .plt format (>= v10)..."
		  << std::endl;
      
      // Read the value of 1 to determine byte ordering
      {
	int one = 0;
	
	in.read (buf, SIZEOF_INT);
	memcpy  (&one, buf, SIZEOF_INT);
	
	if (one != 1)
	  {
	    if (this->verbose())
	      std::cerr << "Tecplot data is Foreign!"
			<< std::endl;
	    this->is_foreign() = true;

	    // Make sure one reversed is one
	    Utility::ReverseBytes rb(this->is_foreign());
	    assert (rb(one) == 1);
	  }
      }

      // A byte-reverser in case the data is foreign
      Utility::ReverseBytes rb(this->is_foreign());

      // Read the title
      {
	int i=0;
    
	this->title().clear();
	do
	  {
	    in.read (buf, SIZEOF_INT);
	    memcpy  (&i, buf, SIZEOF_INT);
	    rb(i);
	    
	    // Don't add trailing \0
	    if (i)
	      this->title() += static_cast<char>(i);
	  }
	while (i);
      }

      // Read the number of variables in the data set
      {
	int nv;
	in.read (buf, SIZEOF_INT);
	memcpy  (&nv, buf, SIZEOF_INT);
	rb(nv);
	
	this->set_n_vars (nv);
      }

      // Read the variable names
      for (unsigned int v=0; v<this->n_vars(); v++)
	{
	  int i=0;
	  this->var_name(v).clear();
      
	  do
	    {
	      in.read (buf, SIZEOF_INT);
	      memcpy  (&i, buf, SIZEOF_INT);
	      rb(i);
	      
	      // Don't add trailing \0
	      if (i)
		this->var_name(v) += static_cast<char>(i);
	    }
	  while (i);
	}

  
      
      // Read zones from the header.
      // Continue reading until the end-of-header
      // marker (357.) is found.
      int nz=0;
      std::vector<std::string> zname;
      std::vector<int>         zpack, ztype, zimax, zjmax, zkmax;
      
      {
	float f=0.;
    
	do
	  {
	    // find the next Zone marker	
	    do
	      {
		f = 0.;
		in.read (buf, SIZEOF_FLOAT);
		memcpy  (&f, buf, SIZEOF_FLOAT);
		rb(f);
	      }
	    while ((f != 299.) &&
		   (f != 357.) &&
		   in.good());
	  
	  
	    // Did we overrun the file?
	    if (!in.good())
	      {
		std::cerr << "ERROR: Unexpected end-of-file!"
			  << std::endl;
		abort();
	      }    
	  
	    // Found a Zone marker
	    else if (f == 299.) 
	      {
		// Incriment the Zone counter
		nz++;
	  
		// Read the zone name
		{
		  int i=0;
		  std::string name;
	    
		  do
		    {
		      in.read (buf, SIZEOF_INT);
		      memcpy  (&i, buf, SIZEOF_INT);
		      rb(i);
		  
		      // Don't add trailing \0
		      if (i)
			name += static_cast<char>(i);
		    }
		  while (i);
	    
		  zname.push_back(name);
		}
	      
		// Read the zone color
		{
		  int zc=0;	    
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&zc, buf, SIZEOF_INT);
		  rb(zc);
		}
	      
		// Read the zone format
		{
		  int zt;
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&zt, buf, SIZEOF_INT);
		  rb(zt);
		  
		  ztype.push_back(zt);
		}

		// Read the data packing flag
		{
		  int dp=0;
		  in.read (buf, SIZEOF_INT);
		  memcpy (&dp, buf, SIZEOF_INT);
		  rb(dp);

		  zpack.push_back (dp);
		}

		// Will we specify the variable location?
		{
		  int svl=0;
		  int  vl=0;
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&svl, buf, SIZEOF_INT);
		  rb(svl);
		  
		  if (svl)
		    for (unsigned int v=0; v<this->n_vars(); v++)
		      {
			in.read (buf, SIZEOF_INT);
			memcpy  (&vl, buf, SIZEOF_INT);
			rb(vl);
			assert (vl == 0); // Only know about node-based data
			                  // right now
		      }
		    
		}

		// Get the number of user-defined face-neighbors
		{
		  int fn=0;
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&fn, buf, SIZEOF_INT);
		  rb(fn);
		}
	  
		// Read in the block dimensions
		{
		  if (ztype.back() != ORDERED)
		    {
		      int np=0, ne=0;

		      in.read (buf, SIZEOF_INT);
		      memcpy  (&np, buf, SIZEOF_INT);
		      rb(np);

		      in.read (buf, SIZEOF_INT);
		      memcpy  (&ne, buf, SIZEOF_INT);
		      rb(ne);
		    }

		  
		  int
		    imax=0,
		    jmax=0,
		    kmax=0;
		  
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&imax, buf, SIZEOF_INT);
		  rb(imax);
		  
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&jmax, buf, SIZEOF_INT);
		  rb(jmax);
		  
		  in.read (buf, SIZEOF_INT);
		  memcpy  (&kmax, buf, SIZEOF_INT);
		  rb(kmax);
		  
		  zimax.push_back (imax);
		  zjmax.push_back (jmax);
		  zkmax.push_back (kmax);
		}
	      } // else if (f == 299.)
	  }
	while ((f != 357.) && in.good());
      }

      // Set the header data
      this->set_n_zones (nz);

      for (unsigned int z=0; z<this->n_zones(); z++)
	{
	  this->zone_type(z) = ztype[z];
	  this->zone_name(z) = zname[z];
	  this->zone_pack(z) = zpack[z];
	  this->imax(z)      = zimax[z];
	  this->jmax(z)      = zjmax[z];
	  this->kmax(z)      = zkmax[z];
	}
    }

  

  //----------------------------------------------------
  // Unrecognized Tecplot Version!
  else
    {
      std::cerr << "ERROR:  This plot file was written by an unrecognized version of Tecplot!:"
		<< std::endl
		<< this->version()
		<< std::endl;
      abort();
    }
  


  




  // Print the data to the screen.
  if (this->verbose())
    {
      std::cout << "Tecplot Header: "
		<< this->title() << std::endl;

      std::cout << "Variables: ";      
      for (unsigned int v=0; v<this->n_vars(); v++)
	std::cout << "\"" << this->var_name (v) << "\"" << " ";  
      std::cout << std::endl;

      std::cout << "Variable Types: ";      
      for (unsigned int v=0; v<this->n_vars(); v++)
	std::cout << this->var_type (v) << " ";  
      std::cout << std::endl;

      std::cout << "Zone Names: ";
      for (unsigned int z=0; z<this->n_zones(); z++)
	std::cout << "\"" << this->zone_name (z) << "\"" << " ";
      std::cout << std::endl;

      std::cout << "Zone Types: ";
      for (unsigned int z=0; z<this->n_zones(); z++)
	std::cout << this->zone_type (z) << " ";
      std::cout << std::endl;

      std::cout << "Zone Dimensions: " << std::endl;
      for (unsigned int z=0; z<this->n_zones(); z++)
	std::cout << "  ("
		  << this->imax(z) << "," 
		  << this->jmax(z) << "," 
		  << this->kmax(z) << ")"
		  << std::endl;  
    }
}



void PltLoader::read_data (std::istream& in)
{
  assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());
  
  //----------------------------------------------------
  // Read the TECPLOT data for each zone
  if (this->verbose())
    {
      std::cout << "Reading Zones";
      std::cout.flush();
    }
      
  
  for (unsigned int zone=0; zone<this->n_zones(); zone++)
    {
      if (this->verbose())
	{
	  std::cout << ".";
	  std::cout.flush();
	}


      //----------------------------------------------------
      // Read plt files written by older versions of Tecplot
      if (this->version().rfind("V7") < this->version().size())
	{
	  float f = 0.;
	  
	  // Find the next Zone marker.
	  do
	    {
	      f = 0.;
	      in.read (buf, SIZEOF_FLOAT);
	      memcpy  (&f, buf, SIZEOF_FLOAT);
	      rb(f);
	    }
	  while ((f != 299.) && in.good());
	
	  // Did we overrun the file?
	  if (!in.good())
	    {
	      std::cerr << "ERROR: Unexpected end-of-file!"
			<< std::endl;
	      abort();
	    }

	  // Get the number of repeated vars.
	  unsigned int n_rep_vars=0;
	  std::vector<int> rep_vars;
	
	  {	
	    in.read (buf, SIZEOF_INT);
	    memcpy  (&n_rep_vars, buf, SIZEOF_INT);
	    rb(n_rep_vars);

	    rep_vars.resize (n_rep_vars);
	
	    // Get the repeated variables number.
	    for (unsigned int v=0; v<n_rep_vars; v++)
	      {
		std::cerr << "ERROR:  I don't understand repeated variables yet!"
			  << std::endl;
		abort();
	    
		in.read (buf, SIZEOF_INT);
		memcpy  (&rep_vars[v], buf, SIZEOF_INT);
		rb(rep_vars[v]);
	      }	
	  }

	  // Get the variable data type
	  //std::cout << "var_types=";
	  for (unsigned int v=0; v<this->n_vars(); v++)
	    {
	      in.read (buf, SIZEOF_INT);
	      memcpy  (&this->var_type(v), buf, SIZEOF_INT);
	      rb(this->var_type(v));

	      //std::cout << this->var_type(v) << " ";
	    }
	  //std::cout << std::endl;


      
	  // Read the data.
	  switch (this->zone_type(zone) )
	    {
	      // Block-based data.  Structured meshes.
	    case BLOCK:
	      {
		this->read_block_data (in, zone);
		break;
	      }
	  
	      // Point-based data.  Structured meshes.
	    case POINT:
	      {
		this->read_point_data (in, zone);
		break;
	      }
	  
	      // FE block data.  Unstructured meshes.
	    case FEBLOCK:
	      {
		this->read_feblock_data (in, zone);
	    
		if (this->verbose())
	      
		  std::cout << "Zone " << zone << ":" << std::endl
			    << "  nnodes   =" << this->imax(zone) << std::endl
			    << "  nelem    =" << this->jmax(zone) << std::endl
			    << "  elem_type=" << this->kmax(zone) << std::endl
			    << std::endl;
		break;
	      }

	      // FE point data.  Unstructured meshes.
	    case FEPOINT:
	      {
		this->read_fepoint_data (in, zone);
		break;
	      }
	  
	    default:
	      {
		std::cerr << "ERROR: Unsupported Zone type: "
			  << this->zone_type(zone)
			  << std::endl;
		abort();
	      }
	    } // end switch on zone type
	}
  

      //----------------------------------------------------
      // Read plt files written by newer versions of Tecplot
      else if (this->version().rfind("V1") < this->version().size())
	{
	  float f = 0.;
	  
	  // Find the next Zone marker.
	  do
	    {
	      f = 0.;
	      in.read (buf, SIZEOF_FLOAT);
	      memcpy  (&f, buf, SIZEOF_FLOAT);
	      rb(f);
	    }
	  while ((f != 299.) && in.good());
	
	  // Did we overrun the file?
	  if (!in.good())
	    {
	      std::cerr << "ERROR: Unexpected end-of-file!"
			<< std::endl;
	      abort();
	    }

	  // Get the variable data type
	  for (unsigned int v=0; v<this->n_vars(); v++)
	    {
	      in.read (buf, SIZEOF_INT);
	      memcpy  (&this->var_type(v), buf, SIZEOF_INT);
	      rb(this->var_type(v));

	      //std::cout << this->var_type(v) << " ";
	    }

	  // Get the variable sharing flag
	  {
	    int vs=0;
	    int sv=0;

	    in.read (buf, SIZEOF_INT);
	    memcpy  (&vs, buf, SIZEOF_INT);
	    rb(vs);

	    if (vs)
	      {
		for (unsigned int v=0; v<this->n_vars(); v++)
		  {
		    in.read (buf, SIZEOF_INT);
		    memcpy  (&sv, buf, SIZEOF_INT);
		    rb(sv);
		    
		    if (sv != -1)
		      {
			std::cerr << "ERROR:  I don't understand variable sharing!"
				  << std::endl;
			abort();
		      }
		  }
	      }
	  }

	  // Get zone to share connectivity with
	  {
	    int sc=0;
	    in.read (buf, SIZEOF_INT);
	    memcpy  (&sc, buf, SIZEOF_INT);
	    rb(sc);
	    
	    assert (sc == -1);
	  }

      
	  // Read the data.
	  if (this->zone_type(zone) == ORDERED)
	    {
		// Block-based data.  Structured meshes.
	      if (this->zone_pack(zone) == 0)
		this->read_block_data (in, zone);
	  
	      // Point-based data.  Structured meshes.
	      else if (this->zone_pack(zone) == 1)
		this->read_point_data (in, zone);

	      else
		abort();
	    }
	  else
	    {
	      // Block-based data.  Unstructured meshes.
	      if (this->zone_pack(zone) == 0)
		this->read_feblock_data (in, zone);

	      // Point-based data.  Unstructured meshes.
	      else if (this->zone_pack(zone) == 1)
		this->read_fepoint_data (in, zone);

	      else
		abort();
	    }
	}

  

      //----------------------------------------------------
      // Unrecognized Tecplot Version!
      else
	{
	  std::cerr << "ERROR:  This plot file was written by an unrecognized version of Tecplot!:"
		    << std::endl
		    << this->version()
		    << std::endl;
	  abort();
	}
      
    } // end loop on zones
}



void PltLoader::read_block_data (std::istream& in, const unsigned int zone)
{
  assert (in.good());

  
  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  
  for (unsigned int var=0; var<this->n_vars(); var++) 
    {
      
      switch (this->var_type(var))
	{
	  
	  // Read a single-precision variable
	case FLOAT:
	  {
	    std::vector<float> & data = _data[zone][var];
	    
	    data.clear();
	    data.resize (this->imax(zone)*
			 this->jmax(zone)*
			 this->kmax(zone));
	    
	    in.read ((char*) &data[0], SIZEOF_FLOAT*data.size());

	    for (unsigned int i=0; i<data.size(); i++)
	      rb(data[i]);
	    
	    break;
	  }
	  
	  // Read a double-precision variable
	case DOUBLE:
	  {
	    std::vector<double> ddata;
	    std::vector<float> & data = _data[zone][var];
	    
	    data.clear();
	    data.resize (this->imax(zone)*
			 this->jmax(zone)*
			 this->kmax(zone));
	    
	    ddata.resize (this->imax(zone)*
			  this->jmax(zone)*
			  this->kmax(zone));
	    
	    in.read ((char*) &ddata[0], SIZEOF_DOUBLE*ddata.size());

	    for (unsigned int i=0; i<data.size(); i++)
	      data[i] = rb(ddata[i]);
	    
	    break;
	  }
		      
	default:
	  {
	    std::cerr << "ERROR: Unsupported data type: "
		      << this->var_type(var)
		      << std::endl;
	    abort();
	  }
	}
    }
}



void PltLoader::read_point_data (std::istream& in, const unsigned int zone)
{
  assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // First allocate space
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      std::vector<float> & data = _data[zone][var];
      
      data.clear();
      data.reserve (this->imax(zone)*
		    this->jmax(zone)*
		    this->kmax(zone));
    }

  
  for (unsigned int k=0; k<this->kmax(zone); k++)
    for (unsigned int j=0; j<this->jmax(zone); j++)
      for (unsigned int i=0; i<this->imax(zone); i++)
	for (unsigned int var=0; var<this->n_vars(); var++)
	  if (this->var_type(var) == FLOAT)
	    {
	      float f = 0.;
	      
	      assert (in.good());
	      
	      in.read (buf, SIZEOF_FLOAT);
	      memcpy  (&f, buf, SIZEOF_FLOAT);
	      rb(f);
		    
	      _data[zone][var].push_back(f);			
	    }
	  else if (this->var_type(var) == DOUBLE)
	    {
	      double d = 0.;
	      
	      assert (in.good());
	      
	      in.read (buf, SIZEOF_DOUBLE);
	      memcpy  (&d, buf, SIZEOF_DOUBLE);
	      rb(d);
	      
	      _data[zone][var].push_back(d);			
	    }
	  else
	    {
	      std::cerr << "ERROR: unsupported data type: "
			<< this->var_type(var)
			<< std::endl;
	      abort();
	    }
}



void PltLoader::read_feblock_data (std::istream& in, const unsigned int zone)
{
  assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // Read the variable values at each node.
  for (unsigned int var=0; var<this->n_vars(); var++) 
    {      
      switch (this->var_type(var))
	{
	  
	  // Read a single-precision variable
	case FLOAT:
	  {
	    std::vector<float> & data = _data[zone][var];
	    
	    data.clear();
	    data.resize (this->imax(zone));
	    
	    in.read ((char*) &data[0], SIZEOF_FLOAT*data.size());

	    for (unsigned int i=0; i<data.size(); i++)
	      rb(data[i]);
	    
	    break;
	  }
	  
	  // Read a double-precision variable
	case DOUBLE:
	  {
	    std::vector<double> ddata;
	    std::vector<float> & data = _data[zone][var];
	    
	    data.clear();
	    data.resize (this->imax(zone));
	    ddata.resize (this->imax(zone));
	    
	    in.read ((char*) &ddata[0], SIZEOF_DOUBLE*ddata.size());

	    for (unsigned int i=0; i<data.size(); i++)
	      data[i] = rb(ddata[i]);
	    
	    break;
	  }
		      
	default:
	  {
	    std::cerr << "ERROR: Unsupported data type: "
		      << this->var_type(var)
		      << std::endl;
	    abort();
	  }
	}
    }

  // Read the connectivity
  {
    // Get the connectivity repetition flag
    int rep=0;
    in.read ((char*) &rep, SIZEOF_INT);
    rb(rep);

    if (rep == 1)
      {
	std::cerr << "ERROR:  Repeated connectivity not supported!"
		  << std::endl;
	abort();
      }

    // Read the connectivity
    else
      {
	assert (zone < _conn.size());
	assert (this->kmax(zone) < 4);
	
	_conn[zone].resize (this->jmax(zone)*NNodes[this->kmax(zone)]);
	
	in.read ((char*) &_conn[zone][0], SIZEOF_INT*_conn[zone].size());

	for (unsigned int i=0; i<_conn[zone].size(); i++)
	  rb(_conn[zone][i]);
      }
  }  
}



void PltLoader::read_fepoint_data (std::istream& in, const unsigned int zone)
{
  assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // First allocate space
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      std::vector<float> & data = _data[zone][var];
      
      data.clear();
      data.reserve (this->imax(zone));
    }

  
  for (unsigned int i=0; i<this->imax(zone); i++)
    for (unsigned int var=0; var<this->n_vars(); var++)
      if (this->var_type(var) == FLOAT)
	{
	  float f = 0.;
	  
	  assert (in.good());
	  
	  in.read (buf, SIZEOF_FLOAT);
	  memcpy  (&f, buf, SIZEOF_FLOAT);
	  rb(f);
	  
	  _data[zone][var].push_back(f);			
	}
      else if (this->var_type(var) == DOUBLE)
	{
	  double d = 0.;
	  
	  assert (in.good());
	  
	  in.read (buf, SIZEOF_DOUBLE);
	  memcpy  (&d, buf, SIZEOF_DOUBLE);
	  rb(d);
	  
	  _data[zone][var].push_back(d);			
	}
      else
	{
	  std::cerr << "ERROR: unsupported data type: "
		    << this->var_type(var)
		    << std::endl;
	  abort();
	}

  // Read the connectivity
  {
    // Get the connectivity repetition flag
    int rep=0;

    in.read ((char*) &rep, SIZEOF_INT);
    rb(rep);

    if (rep == 1)
      {
	std::cerr << "ERROR:  Repeated connectivity not supported!"
		  << std::endl;
	abort();
      }

    // Read the connectivity
    else
      {
	assert (zone < _conn.size());
	assert (this->kmax(zone) < 4);
	
	_conn[zone].resize (this->jmax(zone)*NNodes[this->kmax(zone)]);
	
	in.read ((char*) &_conn[zone][0], SIZEOF_INT*_conn[zone].size());

	for (unsigned int i=0; i<_conn[zone].size(); i++)
	  rb(_conn[zone][i]);
      }
  }
}
