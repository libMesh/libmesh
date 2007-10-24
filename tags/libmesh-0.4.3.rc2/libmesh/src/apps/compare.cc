#include "libmesh_config.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#ifdef HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  Hower IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
#include <stdio.h>
#include <fstream>

// Local Includes
#include "libmesh.h"
#include "equation_systems.h"
#include "mesh.h"
#include "perfmon.h"
#include "enum_xdr_mode.h"




/** 
 * how to use this, and command line processor
 */
void usage(char *progName)
{
  std::string baseName;
  static std::string helpList =
    "usage:\n"
    "	%s [options] ...\n"
    "\n"
    "options:\n"
    "    -d <dim>                      <dim>-dimensional mesh\n"
    "    -m <string>                   Mesh file name\n"
    "    -l <string>                   Left Equation Systems file name\n"
    "    -r <string>                   Right Equation Systems file name\n"
    "    -t <float>                    threshold\n"
    "    -a                            ASCII format (default)\n"
    "    -b                            binary format\n"
    "    -v                            Verbose\n"
    "    -q                            really quiet\n"
    "    -h                            Print help menu\n"
    "\n"
    "\n"
    " This program is used to compare equation systems to a user-specified\n"
    " threshold.  Equation systems are imported in the libMesh format\n"
    " provided through the read and write methods in class EquationSystems.\n"
    " \n"
    "  ./compare -d 3 -m grid.xda -l leftfile.dat -r rightfile.dat -b -t 1.e-8\n"
    "\n"
    " will read in the mesh grid.xda, the equation systems leftfile.dat and\n"
    " rightfile.dat in binary format and compare systems, and especially the\n"
    " floats stored in vectors.  The comparison is said to be passed when the\n"
    " floating point values agree up to the given threshold.  When no threshold\n"
    " is set the default libMesh tolerance is used.  If neither -a or -b are set,\n"
    " ASCII format is assumed.\n"
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
		      unsigned int& dim,
		      double& threshold,
		      libMeshEnums::XdrMODE& format,
		      bool& verbose,
		      bool& quiet)
{
  char optionStr[] =
    "d:m:l:r:t:abvq?h";

  int opt;

  bool format_set    = false;
  bool left_name_set = false;

  if (argc < 3)
    usage(argv[0]);


  while ((opt = getopt(argc, argv, optionStr)) != -1)
    {
      switch (opt)
	{
	  
	  /**
	   * Get mesh file name
	   */
	case 'm':
	  {
	    if (names.empty())
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Mesh file name must preceed left file name!"
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
	   * Get left file name
	   */
	case 'l':
	  {
	    if (!left_name_set)
	      {
	        names.push_back(optarg);
		left_name_set = true;
	      }
	    else
	      {
		std::cout << "ERROR: Mesh file name must preceed right file name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }
	  
	  /**
	   * Get right file name
	   */
	case 'r':
	  {
	    if ((!names.empty()) && (left_name_set))
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Mesh file name and left file name must preceed "
			  << "right file name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }

	  /**
	   * Get the comparison threshold
	   */
	case 't':
	  {
	    threshold = atof(optarg);
	    break;
	  }

	  /**
	   * Use ascii format
	   */
	case 'a':
	  {
	    if (format_set)
	      {
		std::cout << "ERROR: Equation system file format already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
		format = libMeshEnums::READ;
		format_set = true;
	      }
	    break;
	  }
	  
	  /**
	   * Use binary format
	   */
	case 'b':
	  {
	    if (format_set)
	      {
		std::cout << "ERROR: Equation system file format already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        format = libMeshEnums::DECODE;
		format_set = true;
	      }
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
	   * Be totally quiet, no matter what -v says
	   */
	case 'q':
	  {
	    quiet = true;
	    break;
	  }
	  
	case 'h':
	case '?':
	  usage(argv[0]);
	  
	default:
	  return;
	}
    }

}





/**
 * everything that is identical for the systems, and
 * should _not_ go into EquationSystems::compare(),
 * can go in this do_compare().
 */
bool do_compare (EquationSystems& les,
		 EquationSystems& res,
		 double threshold,
		 bool verbose)
{

  if (verbose)
    {
      std::cout	<< "*********   LEFT SYSTEM    *********" << std::endl;
      les.print_info  ();
      std::cout << "*********   RIGHT SYSTEM   *********" << std::endl;
      res.print_info ();
      std::cout << "********* COMPARISON PHASE *********" << std::endl
		<< std::endl;
    }
 
  /**
   * start comparing
   */
  bool result = les.compare(res, threshold, verbose);
  if (verbose)
    {
      std::cout	<< "*********     FINISHED     *********" << std::endl;
    }
  return result;
}










int main (int argc, char** argv)
{
  libMesh::init (argc, argv);
  
  // these should better be not contained in the following braces
  bool quiet = false;
  bool are_equal;

  {
    PerfMon perfmon(argv[0]);
    
    // default values
    std::vector<std::string> names;
    unsigned int dim                = static_cast<unsigned int>(-1);
    double threshold                = TOLERANCE;
    libMeshEnums::XdrMODE format    = libMeshEnums::READ;
    bool verbose                    = false;
 
    // get commands
    process_cmd_line(argc, argv, 
		     names,
		     dim,
		     threshold,
		     format,
		     verbose,
		     quiet);


    if (dim == static_cast<unsigned int>(-1))
      {
	std::cout << "ERROR:  you must specify the dimension on "
		  << "the command line!\n\n"
		  << argv[0] << " -d 3 ... for example\n\n";
	error();
      }

    if (quiet)
      verbose = false;

    if (verbose)
      {
	  std::cout << "Settings:" << std::endl
		    << " dimensionality = " << dim << std::endl
		    << " mesh           = " << names[0] << std::endl
		    << " left system    = " << names[1] << std::endl
		    << " right system   = " << names[2] << std::endl
		    << " threshold      = " << threshold << std::endl
		    << " read format    = " << format << std::endl 
		    << std::endl;
      }	  


    /**
     * build the left and right mesh for left, inut them
     */
    Mesh left_mesh  (dim);
    Mesh right_mesh (dim);


    if (!names.empty())
      {
	left_mesh.read  (names[0]);
	right_mesh.read (names[0]);

	if (verbose)
	  left_mesh.print_info();
      }
    
    else
      {
	std::cout << "No input specified." << std::endl;
	return 1;
      }


    /**
     * build EquationSystems objects, read them
     */
    EquationSystems left_system  (left_mesh);
    EquationSystems right_system (right_mesh);

    if (names.size() == 3)
      {
	left_system.read  (names[1], format);
	right_system.read (names[2], format);	  
      }    
    else
      {
	std::cout << "Bad input specified." << std::endl;
	error();
      }

    are_equal = do_compare (left_system, right_system, threshold, verbose);

  }

  /**
   * let's see what do_compare found out
   */
  unsigned int our_result;

  if (are_equal)
    {
      if (!quiet)
	  std::cout << std::endl
		    << " Congrat's, up to the defined threshold, the two"  
		    << std::endl
		    << " are identical." 
		    << std::endl;
      our_result=0;
    }
  else
    {
      if (!quiet)
	  std::cout << std::endl
		    << " Oops, differences occured!"  
		    << std::endl
		    << " Use -v to obtain more information where differences occured."
		    << std::endl;
      our_result=1;
    }

//  return libMesh::close();
  return our_result;
}
