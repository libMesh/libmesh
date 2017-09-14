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


#include "libmesh/libmesh_config.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#ifdef LIBMESH_HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  However IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
#include <stdio.h>
#include <fstream>

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/perfmon.h"
#include "libmesh/enum_xdr_mode.h"


using namespace libMesh;


/**
 * how to use this, and command line processor
 */
void usage(const std::string & progName)
{
  std::ostringstream helpList;
  helpList << "usage:\n"
           << "        "
           << progName
           << " [options] ...\n"
           << "\n"
           << "options:\n"
           << "    -d <dim>                      <dim>-dimensional mesh\n"
           << "    -m <string>                   Mesh file name\n"
           << "    -l <string>                   Left Equation Systems file name\n"
           << "    -r <string>                   Right Equation Systems file name\n"
           << "    -t <float>                    threshold\n"
           << "    -a                            ASCII format (default)\n"
           << "    -b                            binary format\n"
           << "    -v                            Verbose\n"
           << "    -q                            really quiet\n"
           << "    -h                            Print help menu\n"
           << "\n"
           << "\n"
           << " This program is used to compare equation systems to a user-specified\n"
           << " threshold.  Equation systems are imported in the libMesh format\n"
           << " provided through the read and write methods in class EquationSystems.\n"
           << " \n"
           << "  ./compare -d 3 -m grid.xda -l leftfile.dat -r rightfile.dat -b -t 1.e-8\n"
           << "\n"
           << " will read in the mesh grid.xda, the equation systems leftfile.dat and\n"
           << " rightfile.dat in binary format and compare systems, and especially the\n"
           << " floats stored in vectors.  The comparison is said to be passed when the\n"
           << " floating point values agree up to the given threshold.  When no threshold\n"
           << " is set the default libMesh tolerance is used.  If neither -a or -b are set,\n"
           << " ASCII format is assumed.\n"
           << "\n";

  libmesh_error_msg(helpList.str());
}



void process_cmd_line(int argc,
                      char ** argv,
                      std::vector<std::string> & names,
                      unsigned char & dim,
                      double & threshold,
                      XdrMODE & format,
                      bool & verbose,
                      bool & quiet)
{
  char optionStr[] =
    "d:m:l:r:t:abvq?h";

  int opt;

  bool format_set    = false;
  bool left_name_set = false;

  if (argc < 3)
    usage(std::string(argv[0]));


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
              libmesh_error_msg("ERROR: Mesh file name must precede left file name!");
            break;
          }

          /**
           * Get the mesh dimension
           */
        case 'd':
          {
            dim = cast_int<unsigned char>(atoi(optarg));
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
              libmesh_error_msg("ERROR: Mesh file name must precede right file name!");
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
              libmesh_error_msg("ERROR: Mesh file name and left file name must precede right file name!");
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
              libmesh_error_msg("ERROR: Equation system file format already set!");
            else
              {
                format = READ;
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
              libmesh_error_msg("ERROR: Equation system file format already set!");
            else
              {
                format = DECODE;
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
bool do_compare (EquationSystems & les,
                 EquationSystems & res,
                 double threshold,
                 bool verbose)
{

  if (verbose)
    {
      libMesh::out        << "*********   LEFT SYSTEM    *********" << std::endl;
      les.print_info  ();
      libMesh::out << "*********   RIGHT SYSTEM   *********" << std::endl;
      res.print_info ();
      libMesh::out << "********* COMPARISON PHASE *********" << std::endl
                   << std::endl;
    }

  /**
   * start comparing
   */
  bool result = les.compare(res, threshold, verbose);
  if (verbose)
    {
      libMesh::out        << "*********     FINISHED     *********" << std::endl;
    }
  return result;
}










int main (int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  // these should not be contained in the following braces
  bool quiet = false;
  bool are_equal;

  PerfMon perfmon(argv[0]);

  // default values
  std::vector<std::string> names;
  unsigned char dim               = static_cast<unsigned char>(-1);
  double threshold                = TOLERANCE;
  XdrMODE format                  = READ;
  bool verbose                    = false;

  // get commands
  process_cmd_line(argc, argv,
                   names,
                   dim,
                   threshold,
                   format,
                   verbose,
                   quiet);


  if (dim == static_cast<unsigned char>(-1))
    libmesh_error_msg("ERROR:  you must specify the dimension on "      \
                      << "the command line!\n\n"                        \
                      << argv[0]                                        \
                      << " -d 3 ... for example");

  if (quiet)
    verbose = false;

  if (verbose)
    {
      libMesh::out << "Settings:" << std::endl
                   << " dimensionality = " << +dim << std::endl
                   << " mesh           = " << names[0] << std::endl
                   << " left system    = " << names[1] << std::endl
                   << " right system   = " << names[2] << std::endl
                   << " threshold      = " << threshold << std::endl
                   << " read format    = " << format << std::endl
                   << std::endl;
    }


  /**
   * build the left and right mesh for left, init them
   */
  Mesh left_mesh  (init.comm(), dim);
  Mesh right_mesh (init.comm(), dim);


  if (!names.empty())
    {
      left_mesh.read  (names[0]);
      right_mesh.read (names[0]);

      if (verbose)
        left_mesh.print_info();
    }
  else
    {
      libMesh::out << "No input specified." << std::endl;
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
    libmesh_error_msg("Bad input specified.");

  are_equal = do_compare (left_system, right_system, threshold, verbose);


  /**
   * let's see what do_compare found out
   */
  unsigned int our_result;

  if (are_equal)
    {
      if (!quiet)
        libMesh::out << std::endl
                     << " Congrat's, up to the defined threshold, the two"
                     << std::endl
                     << " are identical."
                     << std::endl;
      our_result=0;
    }
  else
    {
      if (!quiet)
        libMesh::out << std::endl
                     << " Oops, differences occurred!"
                     << std::endl
                     << " Use -v to obtain more information where differences occurred."
                     << std::endl;
      our_result=1;
    }

  //  return libMesh::close();
  return our_result;
}
