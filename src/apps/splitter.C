// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Read in a Mesh file and write out partitionings of it that are suitable
// for reading into a DistributedMesh
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/default_coupling.h"
#include "libmesh/checkpoint_io.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/getpot.h"

using namespace libMesh;

// From: http://stackoverflow.com/a/6417908/2042320
std::string remove_extension (const std::string & filename)
{
  size_t lastdot = filename.find_last_of(".");

  if (lastdot == std::string::npos)
    return filename;

  return filename.substr(0, lastdot);
}

int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  if (libMesh::on_command_line("--help") || argc < 3)
    {
      libMesh::out << "Example: " << argv[0] << " --mesh=filename.e --n-procs='4 8 16' "
                                                "[--num-ghost-layers <n>] [--dry-run] [--ascii]\n\n"
                   << "--mesh             Full name of the mesh file to read in. \n"
                   << "--n-procs          Vector of number of processors.\n"
                   << "--num-ghost-layers Number of layers to ghost when partitioning (Default: 1).\n"
                   << "--dry-run          Only test the partitioning, don't write any files.\n"
                   << "--ascii            Write ASCII cpa files rather than binary cpr files.\n"
                   << std::endl;

      return 0;
    }

  std::string filename = libMesh::command_line_value("--mesh", std::string());

  std::vector<int> all_n_procs;
  libMesh::command_line_vector("--n-procs", all_n_procs);

  unsigned int num_ghost_layers = libMesh::command_line_value("--num-ghost-layers", 1);

  ReplicatedMesh mesh(init.comm());

  // If the user has requested additional ghosted layers, we need to add a ghosting functor.
  DefaultCoupling default_coupling;
  if (num_ghost_layers > 1)
    {
      default_coupling.set_n_levels(num_ghost_layers);
      mesh.add_ghosting_functor(default_coupling);
    }

  libMesh::out << "Reading " << filename << std::endl;

  mesh.read(filename);

  for (std::size_t i = 0; i < all_n_procs.size(); i++)
    {
      processor_id_type n_procs = all_n_procs[i];
      libMesh::out << "splitting " << n_procs << " ways..." << std::endl;

      auto cpr = split_mesh(mesh, n_procs);

      if (!libMesh::on_command_line("--dry-run"))
        {
          libMesh::out << "    * writing " << cpr->current_processor_ids().size() << " files per process..." << std::endl;

          const bool binary = !libMesh::on_command_line("--ascii");

          cpr->binary() = binary;
          std::ostringstream outputname;
          outputname << remove_extension(filename) << (binary ? ".cpr" : ".cpa");
          cpr->write(outputname.str());
        }
    }

  return 0;
}
