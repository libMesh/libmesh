#include "mesh_config.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
//#include <unistd.h>
#ifdef HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  Hower IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
//#include <math.h>
//#include <string.h>
#include <stdio.h>
#include <fstream>
//#include <algorithm>

// Local Includes
#include "libmesh.h"
#include "equation_systems.h"
#include "general_system.h"
#include "frequency_system.h"
#include "mesh.h"
// #include "dof_map.h"
#include "perfmon.h"
// #include "point.h"
// #include "elem_quality.h"
// #include "statistics.h"

/** 
 * convenient enums
 */
enum EqnSystemsFormat {EQNSYS_ASCII=0, 
		       EQNSYS_BINARY, 
		       EQNSYS_INVALID_IO};

enum SystemTypeHandled {GENERAL_SYSTEM=0,
			FREQUENCY_SYSTEM, 
			INVALID_SYSTEM};



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
    "    -m <string>                   Mesh file name\n"
    "    -d <dim>                      <dim>-dimensional mesh\n"
    "    -l <string>                   Left Equationn Systems file name\n"
    "    -r <string>                   Right Equation Systems file name\n"
    "    -t <float>                    threshold\n"
    "    -a                            ASCII format (default)\n"
    "    -b                            binary (XDR) format\n"
    "    -g                            GeneralSystem system type (default)\n"
    "    -f                            FrequencySystem system type\n"
    "    -v                            Verbose\n"
    "    -h                            Print help menu\n"
    "\n"
    "\n"
    " This program is used to compare equation systems to a user-specified\n"
    " threshold.  Equation systems are imported in the libMesh format\n"
    " provided through the read and write methods in class EquationSystems.\n"
    " \n"
    "  ./compare -d 3 -m grid.xda -l leftfile.dat -r rightfile.dat -b -t 1.e-8\n"
    "\n"
    " will read in the mesh grid.xdr, the equation systems leftfile.dat and\n"
    " rightfile.dat in binary format and compare string expressions, if given,\n"
    " and all numerical values.  The comparison is said to be passed when the\n"
    " floating point values agree up to the given threshold.  When no threshold\n"
    " is set the default libMesh tolerance is used.  If neither -a or -b are set,\n"
    " ASCII format is assumed.  The same holds for -g and -f, which defaults to -g.\n"
    " For mesh import formats, please consult the meshtool help.\n"
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
		      EqnSystemsFormat& format,
		      SystemTypeHandled& type,
		      bool& verbose
		      )
{
  char optionStr[] =
    "d:m:l:r:t:abgfv?h";

  int opt;

  bool format_set    = false;
  bool type_set      = false;
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
	   * Get output file name
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
		format = EQNSYS_ASCII;
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
	        format = EQNSYS_BINARY;
		format_set = true;
	      }
	    break;
	  }

	  
	  /**
	   * Load GeneralSystem
	   */
	case 'g':
	  {
	    if (type_set)
	      {
		std::cout << "ERROR: Equation system type already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        type = GENERAL_SYSTEM;
		type_set = true;
	      }
	    break;
	  }
	  
	  /**
	   * Load FrequencySystem
	   */
	case 'f':
	  {
	    if (type_set)
	      {
		std::cout << "ERROR: Equation system type already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        type = FREQUENCY_SYSTEM;
		type_set = true;
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
	  
	case 'h':
	case '?':
	  usage(argv[0]);
	  
	default:
	  return;
	}
    }

}






template <typename T>
bool do_compare (T&, T&, double, bool)
{
    std::cout << "Do nothing" << std::endl;
  error();
}






/**
 * specialization for GeneralSystem
 */
template <>
bool do_compare (EquationSystems<GeneralSystem>& les,
		 EquationSystems<GeneralSystem>& res,
		 double threshold,
		 bool verbose)
{

  if (verbose)
    {
      std::cout << "Comparing GeneralSystems" << std::endl << std::endl;
      std::cout << std::endl << " *** LEFT SYSTEM ***" << std::endl;
      les.print_info  ();
      std::cout << std::endl << " *** RIGHT SYSTEM ***" << std::endl;
      res.print_info ();
    }
 
  /**
   * ----------------------------------------------------------------------------
   *
   * start comparing values, depending on type
   */

  return true;
}





/**
 * specialization for FrequencySystem
 */
template <>
bool do_compare (EquationSystems<FrequencySystem>& les,
		 EquationSystems<FrequencySystem>& res,
		 double threshold,
		 bool verbose)
{

  if (verbose)
    {
      std::cout << "Comparing FrequencySystems" << std::endl << std::endl;
      std::cout << std::endl << " *** LEFT SYSTEM ***" << std::endl;
      les.print_info  ();
      std::cout << std::endl << " *** RIGHT SYSTEM ***" << std::endl;
      res.print_info ();
    }
      

  /**
   * ----------------------------------------------------------------------------
   *
   * start comparing values, depending on type
   */

  return false;
}






int main (int argc, char** argv)
{
  libMesh::init (argc, argv);
  
  {
    PerfMon perfmon(argv[0]);
    
    // default values
    std::vector<std::string> names;
    unsigned int dim                = static_cast<unsigned int>(-1);
    double threshold                = TOLERANCE;
    EqnSystemsFormat format         = EQNSYS_ASCII;
    SystemTypeHandled type          = GENERAL_SYSTEM;
    bool verbose                    = false;

    bool are_equal;
 
    // get commands
    process_cmd_line(argc, argv, 
		     names,
		     dim,
		     threshold,
		     format,
		     type,
		     verbose);


    if (dim == static_cast<unsigned int>(-1))
      {
	std::cout << "ERROR:  you must specify the dimension on "
		  << "the command line!\n\n"
		  << argv[0] << " -d 3 ... for example\n\n";
	error();
      }

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
     * build an appropriate EquationSystems object, read them
     */
    switch (type)
      {
	  
      case GENERAL_SYSTEM:
        {
	  EquationSystems<GeneralSystem> left_system  (left_mesh, 
						       INVALID_SOLVER_PACKAGE);
	  EquationSystems<GeneralSystem> right_system (right_mesh, 
						       INVALID_SOLVER_PACKAGE);

	  if (names.size() == 3)
	    {
	      switch (format)
	        {	  
		  case EQNSYS_ASCII:
		  {
		      left_system.read  (names[1], Xdr::READ);
		      right_system.read (names[2], Xdr::READ);
		      break;
		  }
		  
		  case EQNSYS_BINARY:
		  {
		      left_system.read  (names[1], Xdr::DECODE);
		      right_system.read (names[2], Xdr::DECODE);
		      break;
		  }

		  default:
		  {
		      std::cout << "Bad system type." << std::endl;
		      error();
		  }
		}
	  
	    }    
	  else
	    {
	      std::cout << "Bad input specified." << std::endl;
	      error();
	    }

	  are_equal = do_compare (left_system, right_system, threshold, verbose);

	  break;
	}

      case FREQUENCY_SYSTEM:
        {
	  EquationSystems<FrequencySystem> left_system  (left_mesh, 
							 INVALID_SOLVER_PACKAGE);;
	  EquationSystems<FrequencySystem> right_system (right_mesh, 
							 INVALID_SOLVER_PACKAGE);;

	  if (names.size() == 3)
	    {
	      switch (format)
	        {	  
		  case EQNSYS_ASCII:
		  {
		      left_system.read  (names[1], Xdr::READ);
		      right_system.read (names[2], Xdr::READ);
		      break;
		  }
		  
		  case EQNSYS_BINARY:
		  {
		      left_system.read  (names[1], Xdr::DECODE);
		      right_system.read (names[2], Xdr::DECODE);
		      break;
		  }

		  default:
		  {
		      std::cout << "Bad system type." << std::endl;
		      error();
		  }
		}
	  
	    }    
	  else
	    {
	      std::cout << "Bad input specified." << std::endl;
	      error();
	    }

	  are_equal = do_compare (left_system, right_system, threshold, verbose);

	  break;
	}
      default:
        {
	  std::cout << "Bad system type." << std::endl;
	  error();
	}
      }

    /**
     * let's see what do_compare found out
     */

    if (are_equal)
    {
	std::cout << " Congrat's, up to the defined threshold, the two"  << std::endl
		  << " are identical." << std::endl;
    }
    else
    {
	std::cout << " Oops, differences occured!"  << std::endl
		  << " Use -v to obtain more information where differences occured."
		  << std::endl;
    }
 
  }

  return libMesh::close();
}
