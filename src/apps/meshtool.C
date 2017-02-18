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
#include <algorithm>
#include <fstream>
#ifdef LIBMESH_HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  Hower IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <vector>

// Local Includes
#include "libmesh/boundary_info.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/elem_quality.h"
#include "libmesh/gmv_io.h"
#include "libmesh/inf_elem_builder.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/perfmon.h"
#include "libmesh/statistics.h"
#include "libmesh/string_to_enum.h"


using namespace libMesh;


/*
 * convenient enum for the mode in which the boundary mesh
 * should be written
 */
enum BoundaryMeshWriteMode {BM_DISABLED=0, BM_MESH_ONLY};


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
           << "    -i <string>                   Input file name\n"
           << "    -o <string>                   Output file name\n"
           << "    -s <string>                   Solution file name\n"
           << "\n    -b                            Write the boundary conditions\n"
           << "    -D <factor>                   Randomly move interior nodes by D*hmin\n"
           << "    -h                            Print help menu\n"
           << "    -p <count>                    Partition into <count> subdomains\n"
#ifdef LIBMESH_ENABLE_AMR
           << "    -r <count>                    Globally refine <count> times\n"
#endif
           << "    -t (-d 2 only)                Convert to triangles first\n"
           << "                                  (allows to write .unv file of the\n"
           << "                                  boundary with the correct node ids)\n"
           << "    -v                            Verbose\n"
           << "    -q <metric>                   Evaluates the named element quality metric\n"
           << "    -1                            Converts a mesh of higher order elements\n"
           << "                                  to their first-order counterparts:\n"
           << "                                  Quad8 -> Quad4, Tet10 -> Tet4 etc\n"
           << "    -2                            Converts a mesh of linear elements\n"
           << "                                  to their second-order counterparts:\n"
           << "                                  Quad4 -> Quad8, Tet4 -> Tet10 etc\n"
           << "    -3                            Same, but to the highest possible:\n"
           << "                                  Quad4 -> Quad9, Hex8 -> Hex27 etc\n"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
           << "\n    -a                            Add infinite elements\n"
           << "    -x <coord>                    Specify infinite element origin\n"
           << "    -y <coord>                    coordinates. If none given, origin\n"
           << "    -z <coord>                    is determined automatically.\n"
           << "    -X                            When building infinite elements \n"
           << "    -Y                            treat mesh as x/y/z-symmetric.\n"
           << "    -Z                            When -X is given, -x <coord> also\n"
           << "                                  has to be given.  Similar for y,z.\n"
#endif
    // << "\n    -l                            Build the L connectivity matrix \n"
    // <<   "    -L                            Build the script L connectivity matrix \n"
           << "\n"
           << "\n"
           << " This program is used to convert and partions from/to a variety of\n"
           << " formats.  File types are inferred from file extensions.  For example,\n"
           << " the command:\n"
           << "\n"
           << "  ./meshtool -d 2 -i in.e -o out.plt\n"
           << "\n"
           << " will read a 2D mesh in the ExodusII format (from Cubit, for example)\n"
           << " from the file in.e.  It will then write the mesh in the Tecplot\n"
           << " binary format to out.plt.\n"
           << "\n"
           << " and\n"
           << "\n"
           << "  ./meshtool -d 3 -i bench12.mesh.0000 -o out.gmv -s bench12.soln.0137\n"
           << "\n"
           << " will read a 3D MGF mesh from the file bench12.mesh.0000, read a\n"
           << " solution from bench12.soln.0137, and write the output in GMV format\n"
           << " to out.gmv\n"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
           << "\n"
           << " and\n"
           << "\n"
           << "  ./meshtool -d 3 -i dry.unv -o packed.gmv -a -x 30.5 -y -10.5 -X\n"
           << "\n"
           << " will read a 3D Universal file, determine the z-coordinate of the origin\n"
           << " automatically, e.g. z_origin = 3., build infinite elements with the\n"
           << " origin (30.5, -10.5, 3.) on top of volume elements, while preserving\n"
           << " a symmetry plane through (30.5, -10.5, 3.) perpendicular to x.\n"
           << " It is imperative that the origin lies _inside_ the given volume mesh.\n"
           << " If not, infinite elements are not correctly built!\n"
#endif
           << "\n"
           << " Currently this program supports the following formats:\n"
           << "\n"
           << "INPUT:\n"
           << "     .e   -- Sandia's ExodusII binary grid format\n"
           << "     .ucd -- AVS unstructured ASCII grid format\n"
           << "     .unv -- SDRC I-Deas Universal File ASCII format\n"
           << "     .xda -- libMesh human-readable ASCII format\n"
           << "     .xdr -- libMesh binary format\n"
           << "\n"
           << "OUTPUT:\n"
           << "     .dat   -- Tecplot ASCII format\n"
           << "     .e     -- Sandia's ExodusII format\n"
           << "     .exd   -- Sandia's ExodusII format\n"
           << "     .fro   -- ACDL's .fro format\n"
           << "     .gmv   -- LANL's General Mesh Viewer format\n"
           << "     .mesh  -- MEdit mesh format\n"
           << "     .msh   -- GMSH ASCII file\n"
           << "     .plt   -- Tecplot binary format\n"
           << "     .poly  -- TetGen ASCII file\n"
           << "     .pvtu  -- Paraview VTK format\n"
           << "     .ucd   -- AVS's ASCII UCD format\n"
           << "     .unv   -- I-deas Universal format\n"
           << "     .xda   -- libMesh ASCII format\n"
           << "     .xdr   -- libMesh binary format\n"
           << "     .gz    -- any above format gzipped\n"
           << "     .bz2   -- any above format bzip2'ed\n"
           << "\n";

  libmesh_error_msg(helpList.str());
}



void process_cmd_line(int argc,
                      char ** argv,
                      std::vector<std::string> & names,
                      unsigned int & n_subdomains,
                      unsigned int & n_rsteps,
                      unsigned char & dim,
                      double & dist_fact,
                      bool & verbose,
                      BoundaryMeshWriteMode & write_bndry,
                      bool & convert_first_order,
                      unsigned int & convert_second_order,
                      bool & triangulate,
                      bool & do_quality,
                      ElemQuality & quality_type,
                      bool & addinfelems,

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                      InfElemBuilder::InfElemOriginValue & origin_x,
                      InfElemBuilder::InfElemOriginValue & origin_y,
                      InfElemBuilder::InfElemOriginValue & origin_z,
#endif

                      bool & x_sym,
                      bool & y_sym,
                      bool & z_sym
                      )
{

#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /*
   * initialize these to some values,
   * so that the compiler does not complain
   */
  triangulate     = false;
  do_quality      = false;
  quality_type    = DIAGONAL;
  addinfelems     = false;
  x_sym           = y_sym           = z_sym           = false;

  char optionStr[] =
    "i:o:s:d:D:r:p:tbB123vlLm?h";

#else

  char optionStr[] =
    "i:o:q:s:d:D:r:p:tbB123a::x:y:z:XYZvlLm?h";

#endif

  bool b_mesh_b_given = false;
  bool b_mesh_B_given = false;

  int opt;

  if (argc < 3)
    usage(std::string(argv[0]));



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
              libmesh_error_msg("ERROR: Input name must precede output name!");
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
              libmesh_error_msg("ERROR: Input name must precede output name!");
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
              libmesh_error_msg("ERROR: Input and output names must precede solution name!");
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
           * Triangulate in 2D
           */
        case 't':
          {
            triangulate = true;
            break;
          }

          /**
           * Calculate element qualities
           */
        case 'q':
          {
            do_quality = true;
            quality_type = Utility::string_to_enum<ElemQuality>(optarg);
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
            if (b_mesh_B_given)
              libmesh_error_msg("ERROR: Do not use -b and -B concurrently!");

            b_mesh_b_given = true;
            write_bndry = BM_MESH_ONLY;
            break;
          }

          /**
           * Convert elements to first-order
           * counterparts
           */
        case '1':
          {
            convert_first_order = true;
            break;
          }


          /**
           * Convert elements to second-order
           * counterparts
           */
        case '2':
          {
            convert_second_order = 2;
            break;
          }

          /**
           * Convert elements to second-order
           * counterparts, highest order possible
           */
        case '3':
          {
            convert_second_order = 22;
            break;
          }


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

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
            origin_x.first  = true;
            origin_x.second = atof(optarg);
            break;
          }

        case 'y':
          {
            origin_y.first  = true;
            origin_y.second = atof(optarg);
            break;
          }

        case 'z':
          {
            origin_z.first  = true;
            origin_z.second = atof(optarg);
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

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS



        case 'h':
        case '?':
          usage(argv[0]);

        default:
          return;
        };
    };

}



// A helper function for creating a submesh of active elements.  Why do we need this?
// Well, the create submesh function is expecting const_element_iterators,
// and this is just one way to get them...
void construct_mesh_of_active_elements(Mesh & new_mesh, const Mesh & mesh)
{
  MeshBase::const_element_iterator       it     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator it_end = mesh.active_elements_end();
  mesh.create_submesh(new_mesh, it, it_end);
}






int main (int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  PerfMon perfmon(argv[0]);

  unsigned int n_subdomains = 1;
  unsigned int n_rsteps = 0;
  unsigned char dim = static_cast<unsigned char>(-1); // invalid dimension
  double dist_fact = 0.;
  bool verbose = false;
  BoundaryMeshWriteMode write_bndry = BM_DISABLED;
  bool convert_first_order = false;
  unsigned int convert_second_order = 0;
  bool addinfelems = false;
  bool triangulate = false;
  bool do_quality = false;
  ElemQuality quality_type = DIAGONAL;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  InfElemBuilder::InfElemOriginValue origin_x(false, 0.);
  InfElemBuilder::InfElemOriginValue origin_y(false, 0.);
  InfElemBuilder::InfElemOriginValue origin_z(false, 0.);
#endif

  bool x_sym=false;
  bool y_sym=false;
  bool z_sym=false;


  std::vector<std::string> names;
  std::vector<std::string> var_names;
  std::vector<Number>      soln;

  process_cmd_line(argc, argv, names,
                   n_subdomains, n_rsteps, dim,
                   dist_fact, verbose, write_bndry,
                   convert_first_order,
                   convert_second_order,

                   triangulate,

                   do_quality,
                   quality_type,

                   addinfelems,

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                   origin_x, origin_y, origin_z,
#endif

                   x_sym, y_sym, z_sym);

  UniquePtr<Mesh> mesh_ptr;
  if (dim == static_cast<unsigned char>(-1))
    {
      mesh_ptr.reset(new Mesh(init.comm()));
    }
  else
    {
      mesh_ptr.reset(new Mesh(init.comm(),dim));
    }

  Mesh & mesh = *mesh_ptr;

  /**
   * Read the input mesh
   */
  if (!names.empty())
    {
      mesh.read(names[0]);

      if (verbose)
        {
          mesh.print_info();
          mesh.get_boundary_info().print_summary();
        }

    }

  else
    {
      libMesh::out << "No input specified." << std::endl;
      return 1;
    }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if(addinfelems)
    {
      if (names.size() == 3)
        libmesh_error_msg("ERROR: Invalid combination: Building infinite elements\n"
                          << "not compatible with solution import.");

      if (write_bndry != BM_DISABLED)
        libmesh_error_msg("ERROR: Invalid combination: Building infinite elements\n"
                          << "not compatible with writing boundary conditions.");

      /*
       * Sanity checks: -X/Y/Z can only be used, when the
       * corresponding coordinate is also given (using -x/y/z)
       */
      if ((x_sym && !origin_x.first) ||     // claim x-symmetry, but x-coordinate of origin not given!
          (y_sym && !origin_y.first) ||     // the same for y
          (z_sym && !origin_z.first))       // the same for z
        libmesh_error_msg("ERROR: When x-symmetry is requested using -X, then\n"
                          << "the option -x <coord> also has to be given.\n"
                          << "This holds obviously for y and z, too.");

      // build infinite elements
      InfElemBuilder(mesh).build_inf_elem(origin_x, origin_y, origin_z,
                                          x_sym, y_sym, z_sym,
                                          verbose);

      if (verbose)
        {
          mesh.print_info();
          mesh.get_boundary_info().print_summary();
        }

    }

  // sanity check
  else if ((origin_x.first ||  origin_y.first || origin_z.first) ||
           (x_sym          ||  y_sym          || z_sym))
    libmesh_error_msg("ERROR:  -x/-y/-z/-X/-Y/-Z is only to be used when\n"
                      << "the option -a is also specified!");

#endif


  /**
   * Possibly read the solution (-s option).
   */
  if (names.size() == 3)
    {
      // TODO: Read XDR/A mesh file, contstruct an EquationSystems
      // object, read XDR/A solution file by calling
      // es.read(file, READ_HEADER|READ_DATA|READ_ADDITIONAL_DATA);
      // then store a localized copy of the solution vector into 'soln'.
      libmesh_error_msg("Importing an XDA solution file with -s is not supported.");
    }



  /**
   * Maybe Triangulate
   */
  //  if (dim == 2 && triangulate)
  if (triangulate)
    {
      if (verbose)
        libMesh::out << "...Converting to all simplices...\n";

      MeshTools::Modification::all_tri(mesh);
    }

  /**
   * Compute Shape quality metrics
   */
  if (do_quality)
    {
      StatisticsVector<Real> sv;
      sv.reserve(mesh.n_elem());

      libMesh::out << "Quality type is: " << Quality::name(quality_type) << std::endl;

      // What are the quality bounds for this element?
      std::pair<Real, Real> bounds = mesh.elem_ref(0).qual_bounds(quality_type);
      libMesh::out << "Quality bounds for this element type are: ("
                   << bounds.first
                   << ", "
                   << bounds.second
                   << ") "
                   << std::endl;

      MeshBase::const_element_iterator it  = mesh.active_elements_begin(),
        end = mesh.active_elements_end();
      for (; it != end; ++it)
        {
          Elem * e = *it;
          sv.push_back(e->quality(quality_type));
        }

      const unsigned int n_bins = 10;
      libMesh::out << "Avg. shape quality: " << sv.mean() << std::endl;

      // Find element indices below the specified cutoff.
      // These might be considered "bad" elements which need refinement.
      std::vector<dof_id_type> bad_elts  = sv.cut_below(0.8);
      libMesh::out << "Found " << bad_elts.size()
                   << " of " << mesh.n_elem()
                   << " elements below the cutoff." << std::endl;

      // Compute the histogram for this distribution
      std::vector<dof_id_type> histogram;
      sv.histogram(histogram, n_bins);

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

          out << "title('" << Quality::name(quality_type) << "');" << std::endl;

        }
    }


  /**
   * Possibly convert all linear
   * elements to second-order counterparts
   */
  if (convert_first_order)
    {
      if (verbose)
        libMesh::out << "Converting elements to first order counterparts\n";

      mesh.all_first_order();

      if (verbose)
        {
          mesh.print_info();
          mesh.get_boundary_info().print_summary();
        }
    }

  /**
   * Possibly convert all linear
   * elements to second-order counterparts
   */
  if (convert_second_order > 0)
    {
      bool second_order_mode = true;
      std:: string message = "Converting elements to second order counterparts";
      if (convert_second_order == 2)
        {
          second_order_mode = false;
          message += ", lower version: Quad4 -> Quad8, not Quad9";
        }

      else if (convert_second_order == 22)
        {
          second_order_mode = true;
          message += ", highest version: Quad4 -> Quad9";
        }

      else
        libmesh_error_msg("Invalid value, convert_second_order = " << convert_second_order);

      if (verbose)
        libMesh::out << message << std::endl;

      mesh.all_second_order(second_order_mode);

      if (verbose)
        {
          mesh.print_info();
          mesh.get_boundary_info().print_summary();
        }
    }


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Possibly refine the mesh
   */
  if (n_rsteps > 0)
    {
      if (verbose)
        libMesh::out << "Refining the mesh "
                     << n_rsteps << " times"
                     << std::endl;

      MeshRefinement mesh_refinement (mesh);
      mesh_refinement.uniformly_refine(n_rsteps);

      if (verbose)
        {
          mesh.print_info();
          mesh.get_boundary_info().print_summary();
        }
    }


  /**
   * Possibly distort the mesh
   */
  if (dist_fact > 0.)
    {
      libMesh::out << "Distoring the mesh by a factor of "
                   << dist_fact
                   << std::endl;

      MeshTools::Modification::distort(mesh,dist_fact);
    };


  /*
    char filechar[81];
    sprintf(filechar,"%s-%04d.plt", "out", 0);
    std::string oname(filechar);

    mesh.write(oname);

    for (unsigned int step=0; step<100; step++)
    {
    //        const Real x = .5 + .25*cos((((Real) step)/100.)*3.1415927);
    //        const Real y = .5 + .25*sin((((Real) step)/100.)*3.1415927);
    const Real x = 2.5*cos((((Real) step)/100.)*3.1415927);
    const Real y = 2.5*sin((((Real) step)/100.)*3.1415927);

    const Point p(x,y);

    for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem_ref(e).active())
    mesh.elem_ref(e).set_refinement_flag() = -1;



    for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem_ref(e).active())
    {
    const Point diff = mesh.elem_ref(e).centroid(mesh) - p;

    if (diff.size() < .5)
    {
    if (mesh.elem_ref(e).level() < 4)
    mesh.elem_ref(e).set_refinement_flag() = 1;
    else if (mesh.elem_ref(e).level() == 4)
    mesh.elem_ref(e).set_refinement_flag() = 0;
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



  //     /**
  //      * Possibly partition the mesh
  //      */
  if (n_subdomains > 1)
    mesh.partition(n_subdomains);


  /**
   * Possibly write the mesh
   */
  {
    if (names.size() >= 2)
      {
        /*
         * When the mesh got refined, it is likely that
         * the user does _not_ want to write also
         * the coarse elements, but only the active ones.
         * Use Mesh::create_submesh() to create a mesh
         * of only active elements, and then write _this_
         * new mesh.
         */
        if (n_rsteps > 0)
          {
            if (verbose)
              libMesh::out << " Mesh got refined, will write only _active_ elements." << std::endl;

            Mesh new_mesh (init.comm(), mesh.mesh_dimension());

            construct_mesh_of_active_elements(new_mesh, mesh);

            // now write the new_mesh
            if (names.size() == 2)
              new_mesh.write(names[1]);
            else if (names.size() == 3)
              new_mesh.write(names[1], soln, var_names);
            else
              libmesh_error_msg("Invalid names.size() = " << names.size());
          }
        else
          {
            if (names.size() == 2)
              mesh.write(names[1]);
            else if (names.size() == 3)
              mesh.write(names[1], soln, var_names);
            else
              libmesh_error_msg("Invalid names.size() = " << names.size());
          }



        /**
         * Possibly write the BCs
         */
        if (write_bndry != BM_DISABLED)
          {
            BoundaryMesh boundary_mesh
              (mesh.comm(), cast_int<unsigned char>(mesh.mesh_dimension()-1));

            std::string boundary_name = "bndry_";
            boundary_name += names[1];

            if (write_bndry == BM_MESH_ONLY)
              mesh.get_boundary_info().sync(boundary_mesh);

            else
              libmesh_error_msg("Invalid value write_bndry = " << write_bndry);

            if (names.size() == 2)
              boundary_mesh.write(boundary_name);
            else if (names.size() == 3)
              boundary_mesh.write(boundary_name, soln, var_names);
          }
      }
  };

  return 0;
}
