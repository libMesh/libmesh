// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/namebased_io.h"
#include "libmesh/dyna_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/tetgen_io.h"
#include "libmesh/ucd_io.h"
#include "libmesh/unv_io.h"
#include "libmesh/utility.h"
#include "libmesh/matlab_io.h"
#include "libmesh/off_io.h"
#include "libmesh/medit_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/fro_io.h"
#include "libmesh/stl_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/abaqus_io.h"
#include "libmesh/checkpoint_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/parallel.h" // broadcast

// C++ includes
#include <iomanip>
#include <fstream>
#include <vector>

#ifdef LIBMESH_HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>  // for getpid() on Unix
#endif

#ifdef LIBMESH_HAVE_PROCESS_H
#include <process.h> // for getpid() on Windows
#endif


namespace libMesh
{



// ------------------------------------------------------------
// NameBasedIO members
void NameBasedIO::read (const std::string & name)
{
  MeshBase & mymesh = MeshInput<MeshBase>::mesh();

  const std::string_view basename = Utility::basename_of(name);

  // See if the file exists.  Perform this check on all processors
  // so that the code is terminated properly in the case that the
  // file does not exist.

  // For Nemesis files, the name we try to read will have suffixes
  // identifying processor rank
  if (basename.rfind(".nem") == basename.size() - 4 ||
      basename.rfind(".n") == basename.size() - 2)
    {
      std::ostringstream full_name;

      // Find the length of a string which represents the highest processor ID
      full_name << (mymesh.n_processors());
      int field_width = cast_int<int>(full_name.str().size());

      // reset the string stream
      full_name.str("");

      // And build up the full filename
      full_name << name
                << '.' << mymesh.n_processors()
                << '.' << std::setfill('0') << std::setw(field_width) << mymesh.processor_id();

      std::ifstream in (full_name.str().c_str());
      libmesh_error_msg_if(!in.good(), "ERROR: cannot locate specified file:\n\t" << full_name.str());
    }
  else if (basename.rfind(".cp")) {} // Do error checking in the reader
  else
    {
      std::ifstream in (name.c_str());
      libmesh_error_msg_if(!in.good(), "ERROR: cannot locate specified file:\n\t" << name);
    }

  // Look for parallel formats first
  if (is_parallel_file_format(basename))
    {
      // no need to handle bz2 files here -- the Xdr class does that.
      if ((basename.rfind(".xda") < basename.size()) ||
          (basename.rfind(".xdr") < basename.size()))
        {
          XdrIO xdr_io(mymesh);

          // .xda* ==> bzip2/gzip/ASCII flavors
          if (basename.rfind(".xda") < basename.size())
            {
              xdr_io.binary() = false;
              xdr_io.read (name);
            }
          else // .xdr* ==> true binary XDR file
            {
              xdr_io.binary() = true;
              xdr_io.read (name);
            }

          // The xdr_io object gets constructed with legacy() == false.
          // if legacy() == true then it means that a legacy file was detected and
          // thus processor 0 performed the read. We therefore need to broadcast the
          // mesh.  Further, for this flavor of mesh solution data ordering is tied
          // to the node ordering, so we better not reorder the nodes!
          if (xdr_io.legacy())
            {
              mymesh.allow_renumbering(false);
              MeshCommunication().broadcast(mymesh);
            }

          // libHilbert-enabled libMesh builds should construct files
          // with a canonical node ordering, which libHilbert-enabled
          // builds will be able to read in again regardless of any
          // renumbering.  So in that case we're free to renumber.
          // However, if either the writer or the reader of this file
          // don't have libHilbert, then we'll have to skip
          // renumbering because we need the numbering to remain
          // consistent with any solution file we read in next.
#ifdef LIBMESH_HAVE_LIBHILBERT
          // if (!xdr_io.libhilbert_ordering())
          //   skip_renumber_nodes_and_elements = true;
#else
          mymesh.allow_renumbering(false);
#endif
        }
      else if (basename.rfind(".nem") < basename.size() ||
               basename.rfind(".n")   < basename.size())
        Nemesis_IO(mymesh).read (name);
      else if (basename.rfind(".cp") < basename.size())
        {
          if (basename.rfind(".cpa") < basename.size())
            CheckpointIO(mymesh, false).read(name);
          else
            CheckpointIO(mymesh, true).read(name);
        }
    }

  // Serial mesh formats
  else
    {
      // Read the file based on extension.  Only processor 0
      // needs to read the mesh.  It will then broadcast it and
      // the other processors will pick it up
      if (mymesh.processor_id() == 0)
        {
          LOG_SCOPE("read()", "NameBasedIO");

          std::ostringstream pid_suffix;
          pid_suffix << '_' << getpid();
          // Nasty hack for reading/writing zipped files
          std::string new_name = name;
          if (name.rfind(".bz2") == name.size() - 4)
            {
#ifdef LIBMESH_HAVE_BZIP
              new_name.erase(new_name.end() - 4, new_name.end());
              new_name += pid_suffix.str();
              std::string system_string = "bunzip2 -f -k -c ";
              system_string += name + " > " + new_name;
              LOG_SCOPE("system(bunzip2)", "NameBasedIO");
              if (std::system(system_string.c_str()))
                libmesh_file_error(system_string);
#else
              libmesh_error_msg("ERROR: need bzip2/bunzip2 to open .bz2 file " << name);
#endif
            }
          else if (name.rfind(".xz") == name.size() - 3)
            {
#ifdef LIBMESH_HAVE_XZ
              new_name.erase(new_name.end() - 3, new_name.end());
              new_name += pid_suffix.str();
              std::string system_string = "xz -f -d -k -c ";
              system_string += name + " > " + new_name;
              LOG_SCOPE("system(xz -d)", "XdrIO");
              if (std::system(system_string.c_str()))
                libmesh_file_error(system_string);
#else
              libmesh_error_msg("ERROR: need xz to open .xz file " << name);
#endif
            }

          if (basename.rfind(".mat") < basename.size())
            MatlabIO(mymesh).read(new_name);

          else if (basename.rfind(".ucd") < basename.size())
            UCDIO(mymesh).read (new_name);

          else if ((basename.rfind(".off")  < basename.size()) ||
                   (basename.rfind(".ogl")  < basename.size()) ||
                   (basename.rfind(".oogl") < basename.size()))
            OFFIO(mymesh).read (new_name);

          else if (basename.rfind(".unv") < basename.size())
            UNVIO(mymesh).read (new_name);

          else if ((basename.rfind(".node")  < basename.size()) ||
                   (basename.rfind(".ele")   < basename.size()))
            TetGenIO(mymesh).read (new_name);

          else if (basename.rfind(".exd") < basename.size() ||
                   basename.rfind(".e") < basename.size())
            ExodusII_IO(mymesh).read (new_name);

          else if (basename.rfind(".msh") < basename.size())
            GmshIO(mymesh).read (new_name);

          else if (basename.rfind(".gmv") < basename.size())
            GMVIO(mymesh).read (new_name);

          else if (basename.rfind(".stl") < basename.size())
            STLIO(mymesh).read (new_name);

          else if (basename.rfind(".pvtu") < basename.size() ||
                   basename.rfind(".vtu") < basename.size())
            VTKIO(mymesh).read(new_name);

          else if (basename.rfind(".inp") < basename.size())
            AbaqusIO(mymesh).read(new_name);

          else if ((basename.rfind(".bext")  < basename.size()) ||
                   (basename.rfind(".bxt")   < basename.size()))
            DynaIO(mymesh).read (new_name);

          else if (basename.rfind(".bez")  < basename.size())
            DynaIO(mymesh, false).read (new_name);

          else
            {
              libmesh_error_msg(" ERROR: Unrecognized file extension: " \
                                << name                                 \
                                << "\n   I understand the following:\n\n" \
                                << "     *.bext -- Bezier files in DYNA format\n" \
                                << "     *.bez  -- Bezier DYNA files, omit spline nodes\n" \
                                << "     *.bxt  -- Bezier files in DYNA format\n" \
                                << "     *.cpa  -- libMesh Checkpoint ASCII format\n" \
                                << "     *.cpr  -- libMesh Checkpoint binary format\n" \
                                << "     *.e    -- Sandia's ExodusII format\n" \
                                << "     *.exd  -- Sandia's ExodusII format\n" \
                                << "     *.gmv  -- LANL's General Mesh Viewer format\n" \
                                << "     *.inp  -- Abaqus .inp format\n" \
                                << "     *.mat  -- Matlab triangular ASCII file\n" \
                                << "     *.n    -- Sandia's Nemesis format\n" \
                                << "     *.nem  -- Sandia's Nemesis format\n" \
                                << "     *.off  -- OOGL OFF surface format\n" \
                                << "     *.ogl  -- OOGL OFF surface format\n" \
                                << "     *.oogl -- OOGL OFF surface format\n" \
                                << "     *.pvtu -- Paraview VTK format\n" \
                                << "     *.stl  -- STereoLithography triangulation format\n" \
                                << "     *.ucd  -- AVS's ASCII UCD format\n" \
                                << "     *.unv  -- I-deas Universal format\n" \
                                << "     *.vtu  -- Paraview VTK format\n" \
                                << "     *.xda  -- libMesh ASCII format\n" \
                                << "     *.xdr  -- libMesh binary format\n" \
                                << "     *.gz   -- any above format gzipped\n" \
                                << "     *.bz2  -- any above format bzip2'ed\n" \
                                << "     *.xz   -- any above format xzipped\n" \
                                );
            }

          // If we temporarily decompressed a file, remove the
          // uncompressed version
          if (name.rfind(".bz2") == name.size() - 4)
            std::remove(new_name.c_str());
          if (name.rfind(".xz") == name.size() - 3)
            std::remove(new_name.c_str());
        }

      // Send the mesh & bcs (which are now only on processor 0) to the other
      // processors
      MeshCommunication().broadcast (mymesh);
    }
}


void NameBasedIO::write (const std::string & name)
{
  const MeshBase & mymesh = MeshOutput<MeshBase>::mesh();

  const std::string_view basename = Utility::basename_of(name);

  // parallel formats are special -- they may choose to write
  // separate files, let's not try to handle the zipping here.
  if (is_parallel_file_format(basename))
    {
      // no need to handle bz2 files here -- the Xdr class does that.
      if (basename.rfind(".xda") < basename.size())
        XdrIO(mymesh).write(name);

      else if (basename.rfind(".xdr") < basename.size())
        XdrIO(mymesh,true).write(name);

      else if (basename.rfind(".nem") < basename.size() ||
               basename.rfind(".n")   < basename.size())
        Nemesis_IO(mymesh).write(name);

      else if (basename.rfind(".cpa") < basename.size())
        CheckpointIO(mymesh,false).write(name);

      else if (basename.rfind(".cpr") < basename.size())
        CheckpointIO(mymesh,true).write(name);

      else
        libmesh_error_msg("Couldn't deduce filetype for " << name);
    }

  // serial file formats
  else
    {
      // Nasty hack for reading/writing zipped files
      std::string new_name = name;
      int pid_0 = 0;
      if (mymesh.processor_id() == 0)
        pid_0 = getpid();
      mymesh.comm().broadcast(pid_0);
      std::ostringstream pid_suffix;
      pid_suffix << '_' << pid_0;

      if (name.rfind(".bz2") == name.size() - 4)
        {
          new_name.erase(new_name.end() - 4, new_name.end());
          new_name += pid_suffix.str();
        }
      else if (name.rfind(".xz") == name.size() - 3)
        {
          new_name.erase(new_name.end() - 3, new_name.end());
          new_name += pid_suffix.str();
        }

      // New scope so that io will close before we try to zip the file
      {
        // Write the file based on extension
        if (basename.rfind(".dat") < basename.size())
          TecplotIO(mymesh).write (new_name);

        else if (basename.rfind(".plt") < basename.size())
          TecplotIO(mymesh,true).write (new_name);

        else if (basename.rfind(".ucd") < basename.size())
          UCDIO (mymesh).write (new_name);

        else if (basename.rfind(".gmv") < basename.size())
          if (mymesh.n_partitions() > 1)
            GMVIO(mymesh).write (new_name);
          else
            {
              GMVIO io(mymesh);
              io.partitioning() = false;
              io.write (new_name);
            }

        else if (basename.rfind(".exd") < basename.size() ||
                 basename.rfind(".e") < basename.size())
          ExodusII_IO(mymesh).write(new_name);

        else if (basename.rfind(".unv") < basename.size())
          UNVIO(mymesh).write (new_name);

        else if (basename.rfind(".mesh") < basename.size())
          MEDITIO(mymesh).write (new_name);

        else if (basename.rfind(".poly") < basename.size())
          TetGenIO(mymesh).write (new_name);

        else if (basename.rfind(".msh") < basename.size())
          GmshIO(mymesh).write (new_name);

        else if (basename.rfind(".fro") < basename.size())
          FroIO(mymesh).write (new_name);

        else if (basename.rfind(".pvtu") < basename.size())
          VTKIO(mymesh).write (name);

        else if (basename.rfind(".stl") < basename.size())
          STLIO(mymesh).write (new_name);

        else
          {
            libMesh::err
              << " ERROR: Unrecognized file extension: " << name
              << "\n   I understand the following:\n\n"
              << "     *.cpa   -- libMesh ASCII checkpoint format\n"
              << "     *.cpr   -- libMesh binary checkpoint format,\n"
              << "     *.dat   -- Tecplot ASCII file\n"
              << "     *.e     -- Sandia's ExodusII format\n"
              << "     *.exd   -- Sandia's ExodusII format\n"
              << "     *.fro   -- ACDL's surface triangulation file\n"
              << "     *.gmv   -- LANL's GMV (General Mesh Viewer) format\n"
              << "     *.mesh  -- MEdit mesh format\n"
              << "     *.msh   -- GMSH ASCII file\n"
              << "     *.n     -- Sandia's Nemesis format\n"
              << "     *.nem   -- Sandia's Nemesis format\n"
              << "     *.plt   -- Tecplot binary file\n"
              << "     *.poly  -- TetGen ASCII file\n"
              << "     *.pvtu  -- VTK (paraview-readable) format\n"
              << "     *.stl   -- STereoLithography triangulation format\n" \
              << "     *.ucd   -- AVS's ASCII UCD format\n"
              << "     *.unv   -- I-deas Universal format\n"
              << "     *.xda   -- libMesh ASCII format\n"
              << "     *.xdr   -- libMesh binary format,\n"
              << std::endl
              << "\n Exiting without writing output\n";
          }
      }

      // Nasty hack for reading/writing zipped files
      if (name.rfind(".bz2") == name.size() - 4)
        {
          LOG_SCOPE("system(bzip2)", "NameBasedIO");
          if (mymesh.processor_id() == 0)
            {
              std::string system_string = "bzip2 -f -c ";
              system_string += new_name + " > " + name;
              if (std::system(system_string.c_str()))
                libmesh_file_error(system_string);
              std::remove(new_name.c_str());
            }
          mymesh.comm().barrier();
        }
      if (name.rfind(".xz") == name.size() - 3)
        {
          LOG_SCOPE("system(xz)", "NameBasedIO");
          if (mymesh.processor_id() == 0)
            {
              std::string system_string = "xz -f -c ";
              system_string += new_name + " > " + name;
              if (std::system(system_string.c_str()))
                libmesh_file_error(system_string);
              std::remove(new_name.c_str());
            }
          mymesh.comm().barrier();
        }
    }
}


void NameBasedIO::write_nodal_data (const std::string & name,
                                    const std::vector<Number> & v,
                                    const std::vector<std::string> & vn)
{
  const MeshBase & mymesh = MeshOutput<MeshBase>::mesh();

  // Write the file based on extension
  if (name.rfind(".dat") < name.size())
    TecplotIO(mymesh).write_nodal_data (name, v, vn);

  else if (name.rfind(".exd") < name.size() ||
           name.rfind(".e") < name.size())
    ExodusII_IO(mymesh).write_nodal_data(name, v, vn);

  else if (name.rfind(".gmv") < name.size())
    {
      if (mymesh.n_subdomains() > 1)
        GMVIO(mymesh).write_nodal_data (name, v, vn);
      else
        {
          GMVIO io(mymesh);
          io.partitioning() = false;
          io.write_nodal_data (name, v, vn);
        }
    }

  else if (name.rfind(".mesh") < name.size())
    MEDITIO(mymesh).write_nodal_data (name, v, vn);

  else if (name.rfind(".msh") < name.size())
    GmshIO(mymesh).write_nodal_data (name, v, vn);

  else if (name.rfind(".nem") < name.size() ||
           name.rfind(".n")   < name.size())
    Nemesis_IO(mymesh).write_nodal_data(name, v, vn);

  else if (name.rfind(".plt") < name.size())
    TecplotIO(mymesh,true).write_nodal_data (name, v, vn);

  else if (name.rfind(".pvtu") < name.size())
    VTKIO(mymesh).write_nodal_data (name, v, vn);

  else if (name.rfind(".ucd") < name.size())
    UCDIO (mymesh).write_nodal_data (name, v, vn);

  else
    {
      libMesh::err
        << " ERROR: Unrecognized file extension: " << name
        << "\n   I understand the following:\n\n"
        << "     *.dat  -- Tecplot ASCII file\n"
        << "     *.e    -- Sandia's ExodusII format\n"
        << "     *.exd  -- Sandia's ExodusII format\n"
        << "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
        << "     *.mesh -- MEdit mesh format\n"
        << "     *.msh  -- GMSH ASCII file\n"
        << "     *.n    -- Sandia's Nemesis format\n"
        << "     *.nem  -- Sandia's Nemesis format\n"
        << "     *.plt  -- Tecplot binary file\n"
        << "     *.pvtu -- Paraview VTK file\n"
        << "     *.ucd  -- AVS's ASCII UCD format\n"
        << "\n Exiting without writing output\n";
    }
}


void NameBasedIO::write_equation_systems (const std::string & filename,
                                          const EquationSystems & es,
                                          const std::set<std::string> * system_names)
{
  // XDA/XDR require a separate code path, and currently only support
  // writing complete restarts
  if (!system_names)
    {
      const std::string_view basename =
        Utility::basename_of(filename);

      if (basename.rfind(".xda") < basename.size())
        {
          es.write(filename,WRITE,
                   EquationSystems::WRITE_DATA |
                   EquationSystems::WRITE_ADDITIONAL_DATA);
          return;
        }
      else if (basename.rfind(".xdr") < basename.size())
        {
          es.write(filename,ENCODE,
                   EquationSystems::WRITE_DATA |
                   EquationSystems::WRITE_ADDITIONAL_DATA);
          return;
        }
    }

  // Other formats just use the default "write nodal values" path
  MeshOutput<MeshBase>::write_equation_systems
    (filename, es, system_names);
}



} // namespace libMesh
