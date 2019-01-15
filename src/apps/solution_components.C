// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Open the mesh and solution file given, create a new solution file,
// and copy all listed variables from the old solution to the new.

#include "libmesh/libmesh.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/id_types.h"
#include "libmesh/elem.h"

unsigned char dim = 2; // This gets overridden by most mesh formats

int main(int argc, char ** argv)
{
  using namespace libMesh;

  LibMeshInit init(argc, argv);

  Mesh mesh1(init.comm(), dim);
  EquationSystems es1(mesh1);

  libMesh::out << "Usage: " << argv[0]
               << " mesh oldsolution newsolution system1 variable1 [sys2 var2...]" << std::endl;

  // We should have one system name for each variable name, and those
  // get preceded by an even number of arguments.
  libmesh_assert (!(argc % 2));

  // We should have at least one system/variable pair following the
  // initial arguments
  libmesh_assert_greater_equal (argc, 6);

  mesh1.read(argv[1]);
  libMesh::out << "Loaded mesh " << argv[1] << std::endl;
  Mesh mesh2(mesh1);
  EquationSystems es2(mesh2);

  es1.read(argv[2],
           EquationSystems::READ_HEADER |
           EquationSystems::READ_DATA |
           EquationSystems::READ_ADDITIONAL_DATA |
           EquationSystems::READ_BASIC_ONLY);
  libMesh::out << "Loaded solution " << argv[2] << std::endl;

  std::vector<unsigned int> old_sys_num((argc-4)/2),
    new_sys_num((argc-4)/2),
    old_var_num((argc-4)/2),
    new_var_num((argc-4)/2);

  std::vector<const System *> old_system((argc-4)/2);
  std::vector<System *> new_system((argc-4)/2);

  for (int argi = 4; argi < argc; argi += 2)
    {
      const char * sysname = argv[argi];
      const char * varname = argv[argi+1];

      const unsigned int pairnum = (argi-4)/2;

      libmesh_assert(es1.has_system(sysname));

      const System & old_sys = es1.get_system(sysname);
      old_system[pairnum] = &old_sys;
      old_sys_num[pairnum] = old_sys.number();

      libmesh_assert(old_sys.has_variable(varname));

      old_var_num[pairnum] = old_sys.variable_number(varname);

      const Variable & variable = old_sys.variable(old_var_num[pairnum]);

      std::string systype = old_sys.system_type();

      System & new_sys = es2.add_system(systype, sysname);
      new_system[pairnum] = &new_sys;
      new_sys_num[pairnum] = new_sys.number();

      new_var_num[pairnum] =
        new_sys.add_variable(varname, variable.type(),
                             &variable.active_subdomains());
    }

  es2.init();

  // A future version of this app should copy variables for
  // non-solution vectors too.

  // Copy over any nodal degree of freedom coefficients
  MeshBase::const_node_iterator new_nit = mesh2.local_nodes_begin();

  for (const auto & old_node : mesh1.local_node_ptr_range())
    {
      const Node * new_node = *new_nit++;

      // Mesh::operator= hopefully preserved elem/node orderings...
      libmesh_assert (*old_node == *new_node);

      for (int argi = 4; argi < argc; argi += 2)
        {
          const unsigned int pairnum = (argi-4)/2;

          const System & old_sys = *old_system[pairnum];
          System & new_sys = *new_system[pairnum];

          const unsigned int n_comp =
            old_node->n_comp(old_sys_num[pairnum],old_var_num[pairnum]);
          libmesh_assert_equal_to(n_comp,
                                  new_node->n_comp(new_sys_num[pairnum],new_var_num[pairnum]));

          for (unsigned int i=0; i<n_comp; i++)
            {
              const dof_id_type
                old_index = old_node->dof_number
                (old_sys_num[pairnum], old_var_num[pairnum], i),
                new_index = new_node->dof_number
                (new_sys_num[pairnum], new_var_num[pairnum], i);
              new_sys.solution->set(new_index,(*old_sys.solution)(old_index));
            }
        }
    }


  // Copy over any element degree of freedom coefficients
  MeshBase::const_element_iterator new_eit = mesh2.active_local_elements_begin();

  for (const auto & old_elem : mesh1.active_local_element_ptr_range())
    {
      const Elem * new_elem = *new_eit++;

      // Mesh::operator= hopefully preserved elem/node orderings...
      libmesh_assert (*old_elem == *new_elem);

      for (int argi = 4; argi < argc; argi += 2)
        {
          const unsigned int pairnum = (argi-4)/2;

          const System & old_sys = *old_system[pairnum];
          System & new_sys = *new_system[pairnum];

          const unsigned int n_comp =
            old_elem->n_comp(old_sys_num[pairnum],old_var_num[pairnum]);
          libmesh_assert_equal_to(n_comp,
                                  new_elem->n_comp(new_sys_num[pairnum],new_var_num[pairnum]));

          for (unsigned int i=0; i<n_comp; i++)
            {
              const dof_id_type
                old_index = old_elem->dof_number
                (old_sys_num[pairnum], old_var_num[pairnum], i),
                new_index = new_elem->dof_number
                (new_sys_num[pairnum], new_var_num[pairnum], i);
              new_sys.solution->set(new_index,(*old_sys.solution)(old_index));
            }
        }
    }

  es2.write(argv[3], EquationSystems::WRITE_DATA);

  return 0;
}
