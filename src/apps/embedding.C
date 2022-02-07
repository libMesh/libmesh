// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Open the input mesh and corresponding solution file named in command line
// arguments, open the output mesh, project that solution onto the
// output mesh, and write a corresponding output solution file.

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/getpot.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh.h"
#include "libmesh/reference_elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/point.h"

// C++ includes


using namespace libMesh;


// If there's a missing input argument, then print a help message
void usage_error(const char * progname)
{
  libMesh::out << "Options: " << progname << '\n'
               << " --dim d               elem dimension\n"
               << " --elem type           elem type (e.g. TET14)\n"
               << " --childnum num        child number\n"
               << " --denominator num     denominator to use\n"
               << std::endl;

  exit(1);
}

// Get an input argument, or print a help message if it's missing
template <typename T>
T assert_argument (GetPot & cl,
                   const std::string & argname,
                   const char * progname,
                   const T & defaultarg)
{
  if (!cl.search(argname))
    {
      libMesh::err << ("No " + argname + " argument found!") << std::endl;
      usage_error(progname);
    }
  return cl.next(defaultarg);
}


int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  GetPot cl(argc, argv);

  const int dim =
    assert_argument(cl, "--dim", argv[0], 0);

  const std::string elem_type_string =
    assert_argument(cl, "--elem", argv[0], std::string(""));

  const int childnum =
    assert_argument(cl, "--childnum", argv[0], 0);

  const int denominator =
    assert_argument(cl, "--denominator", argv[0], 0);

  const ElemType elem_type =
    Utility::string_to_enum<ElemType>(elem_type_string);

  // Getting an embedding matrix isn't a static function, thanks to
  // situations like Tet diagonal selection
  std::unique_ptr<Elem> elem = Elem::build(elem_type);

  const Elem & ref = ReferenceElem::get(elem_type);

  // Lagrange FE for nodal calculations
  FEType fe_type(elem->default_order());

  const unsigned int n_nodes = FEInterface::n_dofs(fe_type, elem.get());
  libmesh_error_msg_if(n_nodes != elem->n_nodes(), "Bad FEInterface value?");

  std::vector<Node> nodes(n_nodes);

  // Get the child vertex positions from childnum; those are easy.
  for (auto v : make_range(elem->n_vertices()))
    {
      for (auto n : make_range(n_nodes))
        {
          const Real embed = ref.embedding_matrix(childnum, v, n);
          if (embed == 1.0)
            {
              nodes[v] = ref.point(n);
              elem->set_node(v) = &nodes[v];
            }
          else if (embed != 0.0)
            libmesh_error_msg("Found fractional embedding on vertex!?");
        }
    }

  // Now figure out the others
  for (auto n : make_range(elem->n_vertices(), n_nodes))
    {
      const auto & pbns = ref.parent_bracketing_nodes(childnum, n);
      if (pbns.empty())
        libmesh_error();

      for (auto pbn : pbns)
        nodes[n] += (ref.point(pbn.first) + ref.point(pbn.second))/2;
      nodes[n] /= pbns.size();
      elem->set_node(n) = &nodes[n];
    }

  for (auto i : elem->node_index_range())
    {
      const Point & pt = elem->point(i);
      std::cout << '{';
      for (auto j : make_range(n_nodes))
        {
          Real shape = FEInterface::shape(dim, fe_type, elem_type, j, pt);
          const Real embed = ref.embedding_matrix(childnum, i, j);
          if (std::abs(shape - embed) < TOLERANCE*TOLERANCE)
            std::cout << "++++++,";
          else
            {
              std::cout << (shape*denominator);
              if (shape != 0.0)
                std::cout << "/r" << denominator;
              std::cout << ", ";
            }
        }
      std::cout << '}' << std::endl;
    }

  return 0;
}
