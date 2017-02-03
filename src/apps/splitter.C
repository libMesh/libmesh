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


// Read in a Mesh file and write out partitionings of it that are suitable
// for reading into a DistributedMesh
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
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
      libMesh::out << "Example: ./splitter-opt --mesh=filename.e --n-procs='4 8 16' --dry-run\n\n"
                   << "--mesh             Full name of the mesh file to read in. \n"
                   << "--n-procs          Vector of number of processors.\n"
                   << "--dry-run          Only test the partitioning, don't write any files.\n"
                   << std::endl;

      return 0;
    }

  std::string filename = libMesh::command_line_value("--mesh", std::string());

  std::vector<int> all_n_procs;
  libMesh::command_line_vector("--n-procs", all_n_procs);

  Parallel::Communicator & comm = init.comm();

  Mesh mesh(init.comm());

  libMesh::out << "Reading " << filename << std::endl;

  mesh.read(filename);

  MetisPartitioner partitioner;

  for (unsigned int i = 0; i < all_n_procs.size(); i++)
    {
      processor_id_type n_procs = all_n_procs[i];

      libMesh::out << "\nWriting out files for " << n_procs << " processors...\n\n" << std::endl;

      // Reset the partitioning each time after the first one
      if (i > 0)
        {
          libMesh::out << "Resetting Partitioning" << std::endl;
          partitioner.partition(mesh, 1);
        }

      libMesh::out << "Partitioning" << std::endl;

      // Partition it to how we want it
      partitioner.partition(mesh, n_procs);

      mesh.print_info();

      // When running in parallel each processor will write out a portion of the mesh files

      processor_id_type num_chunks = n_procs / comm.size();
      processor_id_type remaining_chunks = n_procs % comm.size();

      processor_id_type my_num_chunks = num_chunks;

      processor_id_type my_first_chunk = 0;

      processor_id_type rank = comm.rank();
      processor_id_type comm_size = comm.size();

      if (n_procs >= comm_size) // More partitions than processors
        {
          if (remaining_chunks) // Means that it doesn't split up evenly
            {
              // Spread the remainder over the first few processors
              // There will be "remaining_chunks" number of processors that will each
              // get one extra chunk
              if (rank < remaining_chunks)
                {
                  my_num_chunks += 1;
                  my_first_chunk = my_num_chunks * rank;
                }
              else // The processors beyond the "first" set that don't get an extra chunk
                {
                  // The number of chunks dealt with by the first processors
                  // num chunks         // num procs
                  processor_id_type num_chunks_in_first_procs = (my_num_chunks + 1) * remaining_chunks;
                  processor_id_type distance_to_first_procs = rank - remaining_chunks;

                  my_first_chunk = num_chunks_in_first_procs + (my_num_chunks * distance_to_first_procs);
                }
            }
          else // Splits evenly
            my_first_chunk = my_num_chunks * rank;
        }
      else // More processors than partitions
        {
          if (rank < n_procs)
            {
              my_num_chunks = 1;
              my_first_chunk = rank;
            }
          else
            {
              my_num_chunks = 0;
              my_first_chunk = std::numeric_limits<processor_id_type>::max();
            }
        }

      if (!libMesh::on_command_line("--dry-run"))
        {
          libMesh::out << "Writing " << my_num_chunks << " Files" << std::endl;

          for (unsigned int i = my_first_chunk; i < my_first_chunk + my_num_chunks; i++)
            {
              libMesh::out << " " << 100 * (static_cast<double>(i) / static_cast<double>(my_num_chunks)) << "% Complete" << std::endl;

              CheckpointIO cpr(mesh);
              cpr.current_processor_id() = i;
              cpr.current_n_processors() = n_procs;
              cpr.binary() = true;
              cpr.parallel() = true;
              cpr.write(remove_extension(filename) + ".cpr");
            }

          libMesh::out << " 100% Complete" << std::endl;
        }
    }

  return 0;
}
