#include "mesh_config.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#ifdef HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  Hower IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>

// Local Includes
#include "libmesh.h"
#include "mesh.h"
#include "dof_map.h"
#include "perfmon.h"
#include "point.h"
#include "elem_quality.h"
#include "statistics.h"

// Conditionally include Petsc stuff
#ifdef HAVE_PETSC
#include "petsc_interface.h"
#include "petsc_matrix.h"
#endif


void usage(char *progName)
{
  std::string baseName;
  static std::string helpList =
    "usage:\n"
    "	%s [options] ...\n"
    "\n"
    "options:\n"
    "    -i <string>                   Input file name\n"
    "    -o <string>                   Output file name\n"
    "    -s <string>                   Solution file name\n"
    "    -d <dim>                      <dim>-dimensional mesh\n"
    "    -D <factor>                   Randomly move interior nodes by D*hmin\n"
#ifdef ENABLE_AMR
    "    -r <count>                    Globally refine <count> times\n"
#endif
    "    -p <count>                    Partition into <count> subdomains\n"
    "    -b                            Write the boundary conditions\n"
#ifdef ENABLE_INFINITE_ELEMENTS
    "    -a                            Add infinite elemens\n"
    "    -x <coord>                    Specify infinite element origin,\n"
    "    -y <coord>                    defaults to 0.,0.,0. if none given\n"
    "    -z <coord>                    \n"
    "    -X                            When building infinite elements \n"
    "    -Y                            treats mesh as x/y/z-symmetric\n"
    "    -Z                            \n"
#endif
    "    -l                            Build the L connectivity matrix \n"
    "    -L                            Build the script L connectivity matrix \n"
    "    -v                            Verbose\n"
    "    -h                            Print help menu\n"
    "\n"
    "\n"
    " This program is used to convert and partions from/to a variety of\n"
    " formats.  File types are inferred from file extensions.  For example,\n"
    " the command:\n"
    "\n"
    "  ./meshtool -d 2 -i in.exd -o out.plt\n"
    "\n"
    " will read a 2D mesh in the ExodusII format (from Cubit, for example)\n"
    " from the file in.exd.  It will then write the mesh in the Tecplot\n"
    " binary format to out.plt.\n"
    "\n"
    " and\n"
    "\n"
    "  ./meshtool -d 3 -i bench12.mesh.0000 -o out.gmv -s bench12.soln.0137\n"
    "\n"
    " will read a 3D MGF mesh from the file bench12.mesh.0000, read a\n"
    " solution from bench12.soln.0137, and write the output in GMV format\n"
    " to out.gmv\n"
#ifdef ENABLE_INFINITE_ELEMENTS
    "\n"
    " and\n"
    "\n"
    "  ./meshtool -d 3 -i dry.unv -o packed.gmv -a -x 30.5 -y -10.5 -Z\n"
    "\n"
    " will read a 3D Universal file, build infinite elements with the\n"
    " origin (30.5, -10.5, 0.0) on top of volume elements, while preserving\n"
    " a symmetry plane through (30.5, -10.5, 0.0) perpendicular to z.\n"
#endif
    "\n"
    " Currently this program supports the following formats:\n"
    "\n"
    "INPUT:\n"
    "     .ucd -- AVS unstructured ASCII grid format\n"
    "     .exd -- Sandia's ExodusII binary grid format\n"
    "     .unv -- SDRC I-Deas Universal File ASCII format\n"
    "     .xdr -- Internal binary mesh format\n"
    "     .xda -- Same format, but ASCII and human-readable\n"
    "\n"
    "OUTPUT:\n"
    "     .plt   -- Tecplot binary format\n"
    "     .dat   -- Tecplot ASCII format\n"
    "     .gmv   -- LANL's General Mesh Viewer (~benkirk/work/GMV/linuxogl)\n"
    "     .ugrid -- Kelly's DIVA ASCII format (3D only)\n"
    "     .xdr   -- Internal binary mesh format\n"
    "     .xda   -- Same format, but ASCII and human-readable\n"
    "\n"
    " Direct questions to:\n"
    " benkirk@cfdlab.ae.utexas.edu\n";

    
  if (progName == NULL)
    baseName = "UNKNOWN";
  else 
    baseName = progName;

  
  fprintf(stderr, helpList.c_str(), baseName.c_str());
  fflush(stderr);

  abort();
};



void process_cmd_line(int argc, char **argv,
		      std::vector<std::string>& names,
		      unsigned int& n_subdomains,
		      unsigned int& n_rsteps,
		      unsigned int& dim,
		      double& dist_fact,
		      bool& verbose,
		      bool& write_bndry,
		      bool& addinfelems,
		      double& origin_x,
		      double& origin_y,
		      double& origin_z,
		      bool& x_sym,
		      bool& y_sym,
		      bool& z_sym,
		      bool& build_l,
		      bool& build_script_l
		      )
{

#ifndef ENABLE_INFINITE_ELEMENTS

  addinfelems = false;
  origin_x = origin_y = origin_z = 0.;
  x_sym    = y_sym    = z_sym    = false;

  char optionStr[] =
    "i:o:s:d:D:r:p:bvlLm?h";

#else

  char optionStr[] =
    "i:o:s:d:D:r:p:ba::x:y:z:XYZvlLm?h";

#endif

  int opt;  
  
  if (argc < 3)
    usage(argv[0]);


  
  while ((opt = getopt(argc, argv, optionStr)) != -1)
    {
      switch (opt)
	{
	  
	  /**
	   * Get input file name
	   */
	case 'i':
	  {
	    if (names.empty())
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Input name must preceed output name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }
	  
	  /**
	   * Get output file name
	   */
	case 'o':
	  {
	    if (!names.empty())
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Input name must preceed output name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }

	  /**
	   * Get solution file name
	   */
	case 's':
	  {
	    if (names.size() == 2)
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Input and output names must preceed solution name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }

	  /**
	   * Get the mesh dimension
	   */
	case 'd':
	  {
	    dim = atoi(optarg);
	    break;
	  }

	  /**
	   * Get the mesh distortion factor
	   */
	case 'D':
	  {
	    dist_fact = atof(optarg);
	    break;
	  }

	  /**
	   * Get the number of refinements to do
	   */
	case 'r':
	  {
	    n_rsteps = atoi(optarg);
	    break;
	  }

	  /**
	   * Get the number of subdomains for partitioning
	   */
	case 'p':
	  {
	    n_subdomains = atoi(optarg);
	    break;
	  }
	  
	  /**
	   * Be verbose
	   */
	case 'v':
	  {
	    verbose = true;
	    break;
	  }
	  
	  /**
	   * Try to write the boundary
	   */
	case 'b':
	  {
	    write_bndry = true;
	    break;
	  }
	  

#ifdef ENABLE_INFINITE_ELEMENTS

	  /**
	   * Add infinite elements
	   */
	case 'a':
	  {
	    addinfelems = true;
	    break;
	  }
	  
	  /**
	   * Specify origin coordinates
	   */
	  
	case 'x':
	  {
	    origin_x = atof(optarg);
	    break;
	  }
	  
	case 'y':
	  {
	    origin_y = atof(optarg);
	    break;
	  }
	  
	case 'z':
	  {
	    origin_z = atof(optarg);
	    break;
	  } 
	  
	  /**
	   *  Symmetries
	   */	  
	  
	case 'X':
	  {
	    x_sym = true;
	    break;
	  }
  	case 'Y':
	  {
	    y_sym = true;
	    break;
	  }
  	case 'Z':
	  {
	    z_sym = true;
	    break;
	  }

#endif //ifdef ENABLE_INFINITE_ELEMENTS

	case 'l':
	  {
	    build_l = true;
	    break;
	  }
	case 'L':
	  {
	    build_script_l = true;
	    break;
	  }


	case 'h':
	case '?':
	  usage(argv[0]);
	  
	default:
	  return;
	};
    };

};



int main (int argc, char** argv)
{
  libMesh::init (argc, argv);
  
  {
    PerfMon perfmon(argv[0]);
    
    unsigned int n_subdomains = 1;
    unsigned int n_rsteps = 0;
    unsigned int dim = static_cast<unsigned int>(-1); // invalid dimension
    double dist_fact = 0.;
    bool verbose = false;
    bool write_bndry = false;
    bool addinfelems = false;
    double origin_x=0.;
    double origin_y=0.;
    double origin_z=0.;
    bool x_sym=false;
    bool y_sym=false;
    bool z_sym=false;
    bool build_l=false;
    bool build_script_l=false;
    
      
    std::vector<std::string> names;
    std::vector<std::string> var_names;
    std::vector<Number>      soln;

    process_cmd_line(argc, argv, names,
		     n_subdomains, n_rsteps,
		     dim, dist_fact, verbose, write_bndry, 
		     addinfelems, origin_x, origin_y, origin_z, x_sym, y_sym, z_sym,
		     build_l,
		     build_script_l);

    if (dim == static_cast<unsigned int>(-1))
      {
	std::cout << "ERROR:  you must specify the dimension on "
		  << "the command line!\n\n"
		  << argv[0] << " -d 3 ... for example\n\n";
	error();
      }

    Mesh mesh(dim);
       
    /**
     * Read the input mesh
     */      
    if (!names.empty())
      {
	mesh.read(names[0]);

	if (verbose)
	  mesh.print_info();
      }
    
    else
      {
	std::cout << "No input specified." << std::endl;
	return 1;
      }
       	


#ifdef ENABLE_INFINITE_ELEMENTS

    if(addinfelems)
      {
	if (names.size() == 3)
	  {
	    std::cout << "ERROR: Invalid combination: Building infinite elements " << std::endl
		      << "not compatible with solution import." << std::endl;
	    exit(1);
	  }
	
	if (write_bndry)
	  {
	    std::cout << "ERROR: Invalid combination: Building infinite elements " << std::endl
		      << "not compatible with writing boundary conditions." << std::endl;
	    exit(1);
	  }
	
	mesh.build_inf_elem(Point(origin_x, origin_y, origin_z),
			    x_sym, y_sym, z_sym, 
			    verbose);
	
	if (verbose)
	  mesh.print_info();
      }
    
#endif


    /**
     * Possibly partition the mesh
     */
    if (n_subdomains > 1)
      mesh.metis_partition(n_subdomains);
    
    
    /**
     * Possibly read the solution
     */
    if (names.size() == 3)
      {
	mesh.read_xdr_soln_binary(names[2],
				  soln,
				  var_names);
      }


    /**
     * Compute Shape quality metrics
     */
    const bool do_quality = false;

    if (do_quality)
      {
	StatisticsVector<Real> sv;
	sv.resize(mesh.n_elem());
	
	const ElemQuality q = DIAGONAL;

	std::cout << "Quality type is: " << Quality::name(q) << std::endl;
	
	// What are the quality bounds for this element?
	std::pair<Real, Real> bounds = mesh.elem(0)->qual_bounds(q);
	std::cout << "Quality bounds for this element type are: (" << bounds.first
		  << ", " << bounds.second << ") "
		  << std::endl;

	for (unsigned int e=0; e<mesh.n_elem(); e++)
	  {
	    sv[e] = mesh.elem(e)->quality(q);
	  }

	const unsigned int n_bins = 10;
	std::cout << "Avg. shape quality: " << sv.mean() << std::endl;

	// Find element indices below the specified cutoff.
	// These might be considered "bad" elements which need refinement.
	std::vector<unsigned int> bad_elts  = sv.cut_below(0.8);
	std::cout << "Found " << bad_elts.size()
		  << " of " << mesh.n_elem()
		  << " elements below the cutoff." << std::endl;
	
	/*
	for (unsigned int i=0; i<bad_elts.size(); i++)
	  std::cout << bad_elts[i] << " ";
	std::cout << std::endl;
	*/
	
	// Compute the histogram for this distribution
	std::vector<unsigned int> histogram = sv.histogram(n_bins);
	
	/*
	for (unsigned int i=0; i<n_bins; i++)
	  histogram[i] = histogram[i] / mesh.n_elem();
	*/
	
	const bool do_matlab = true;

	if (do_matlab)
	  {
	    std::ofstream out ("histo.m");
	    
	    out << "% This is a sample histogram plot for Matlab." << std::endl;
	    out << "bin_members = [" << std::endl;
	    for (unsigned int i=0; i<n_bins; i++)
	      out << static_cast<Real>(histogram[i]) / static_cast<Real>(mesh.n_elem())
		  << std::endl;
	    out << "];" << std::endl;
	    
	    std::vector<Real> bin_coords(n_bins);
	    const Real max   = *(std::max_element(sv.begin(), sv.end()));
	    const Real min   = *(std::min_element(sv.begin(), sv.end()));
	    const Real delta = (max - min) / static_cast<Real>(n_bins);
	    for (unsigned int i=0; i<n_bins; i++)
	      bin_coords[i] = min + (i * delta) + delta / 2.0 ;
	    
	    out << "bin_coords = [" << std::endl;
	    for (unsigned int i=0; i<n_bins; i++)
	      out << bin_coords[i] << std::endl;
	    out << "];" << std::endl;
	    
	    out << "bar(bin_coords, bin_members, 1);" << std::endl;
	    out << "hold on" << std::endl;
	    out << "plot (bin_coords, 0, 'kx');" << std::endl;
	    out << "xlabel('Quality (0=Worst, 1=Best)');" << std::endl;
	    out << "ylabel('Percentage of elements in each bin');" << std::endl;
	    out << "axis([" << min << "," << max << ",0, max(bin_members)]);" << std::endl;

	    out << "title('" << Quality::name(q) << "');" << std::endl;
	    
	  }
      }
    
#ifdef ENABLE_AMR
    
    /**
     * Possibly refine the mesh
     */
    if (n_rsteps > 0)
      {
	if (verbose)
	  std::cout << "Refining the mesh "
		    << n_rsteps << " times"
		    << std::endl;
	
	mesh.mesh_refinement.uniformly_refine(n_rsteps);
	
	if (verbose)
	  mesh.print_info();
      };

    
    /**
     * Possibly distort the mesh
     */
    if (dist_fact > 0.)
      {
	std::cout << "Distoring the mesh by a factor of "
		  << dist_fact
		  << std::endl;
	
	mesh.distort(dist_fact);
      };


    /*
    char filechar[81];
    sprintf(filechar,"%s-%04d.plt", "out", 0);
    std::string oname(filechar);
    
    mesh.write(oname);
    
    for (unsigned int step=0; step<100; step++)
      {
//	const Real x = .5 + .25*cos((((Real) step)/100.)*3.1415927); 
//	const Real y = .5 + .25*sin((((Real) step)/100.)*3.1415927);
	const Real x = 2.5*cos((((Real) step)/100.)*3.1415927); 
	const Real y = 2.5*sin((((Real) step)/100.)*3.1415927);

	const Point p(x,y);
	
	for (unsigned int e=0; e<mesh.n_elem(); e++)
	  if (mesh.elem(e)->active())
	    mesh.elem(e)->set_refinement_flag() = -1;
	
	
	
	for (unsigned int e=0; e<mesh.n_elem(); e++)
	  if (mesh.elem(e)->active())
	    {
	      const Point diff = mesh.elem(e)->centroid(mesh) - p;

	      if (diff.size() < .5)
		{
		  if (mesh.elem(e)->level() < 4)
		    mesh.elem(e)->set_refinement_flag() = 1;
		  else if (mesh.elem(e)->level() == 4)
		    mesh.elem(e)->set_refinement_flag() = 0;
		}
	    }
	

	mesh.mesh_refinement.refine_and_coarsen_elements();
	
	char filechar[81];
	sprintf(filechar,"%s-%04d.plt", "out", step+1);
	std::string oname(filechar);

	mesh.write(oname);
      }
    */

#endif

    
    /**
     * Maybe create Damien's connectivity
     * graph matrices.
     */
    {
      if (build_l)
	{
#ifdef HAVE_PETSC	  
	  PetscMatrix<Number> conn;
	  mesh.build_L_graph (conn);
	  conn.print_matlab();
#else
	  std::cerr << "This functionality requires PETSC support!"
		    << std::endl;
	  error();
#endif
	};
      
      
      if (build_script_l)
	{
#ifdef HAVE_PETSC	  
	  PetscMatrix<Number> conn;
	  mesh.build_script_L_graph (conn);
	  conn.print_matlab();
#else
	  std::cerr << "This functionality requires PETSC support!"
		    << std::endl;
	  error();
#endif
	};
    };


    
    /**
     * Possibly write the mesh
     */
    {
      if (names.size() >= 2)
	{
	  if (names.size() == 2)
	    mesh.write(names[1]);
	  else if (names.size() == 3)
	    mesh.write(names[1], soln, var_names);
	  else
	    error();

	  
	  /**
	   * Possibly write the BCs
	   */
	  if (write_bndry)
	    {
	      std::string boundary_name = "bndry_";
	      boundary_name += names[1];
	      
	      mesh.boundary_info.sync(mesh.boundary_mesh);
	      
	      if (names.size() == 2)
		mesh.boundary_mesh.write(boundary_name);
	      else if (names.size() == 3)
		mesh.boundary_mesh.write(boundary_name,
					 soln, var_names);
	    }
	}
    };
          
    /*
    std::cout << "Infinite loop, look at memory footprint" << std::endl;
    for (;;)
      ;
    */
  };

  
  return libMesh::close();
};
  
