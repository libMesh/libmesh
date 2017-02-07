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
#ifdef LIBMESH_ENABLE_VSMOOTHER

// C++ includes
#include <time.h> // for clock_t, clock()
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <iomanip>

// Local includes
#include "libmesh/mesh_smoother_vsmoother.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"

namespace libMesh
{

// Optimization at -O2 or greater seem to break Intel's icc. So if we are
// being compiled with icc let's dumb-down the optimizations for this file
#ifdef __INTEL_COMPILER
#  pragma optimize ( "", off )
#endif

// Member functions for the Variational Smoother
VariationalMeshSmoother::VariationalMeshSmoother(UnstructuredMesh & mesh,
                                                 double theta,
                                                 unsigned miniter,
                                                 unsigned maxiter,
                                                 unsigned miniterBC) :
  MeshSmoother(mesh),
  _percent_to_move(1),
  _dist_norm(0.),
  _adapt_data(libmesh_nullptr),
  _dim(mesh.mesh_dimension()),
  _miniter(miniter),
  _maxiter(maxiter),
  _miniterBC(miniterBC),
  _metric(UNIFORM),
  _adaptive_func(NONE),
  _theta(theta),
  _generate_data(false),
  _n_nodes(0),
  _n_cells(0),
  _n_hanging_edges(0),
  _area_of_interest(libmesh_nullptr)
{}




VariationalMeshSmoother::VariationalMeshSmoother(UnstructuredMesh & mesh,
                                                 std::vector<float> * adapt_data,
                                                 double theta,
                                                 unsigned miniter,
                                                 unsigned maxiter,
                                                 unsigned miniterBC,
                                                 double percent_to_move) :
  MeshSmoother(mesh),
  _percent_to_move(percent_to_move),
  _dist_norm(0.),
  _adapt_data(adapt_data),
  _dim(mesh.mesh_dimension()),
  _miniter(miniter),
  _maxiter(maxiter),
  _miniterBC(miniterBC),
  _metric(UNIFORM),
  _adaptive_func(CELL),
  _theta(theta),
  _generate_data(false),
  _n_nodes(0),
  _n_cells(0),
  _n_hanging_edges(0),
  _area_of_interest(libmesh_nullptr)
{}



VariationalMeshSmoother::VariationalMeshSmoother(UnstructuredMesh & mesh,
                                                 const UnstructuredMesh * area_of_interest,
                                                 std::vector<float> * adapt_data,
                                                 double theta,
                                                 unsigned miniter,
                                                 unsigned maxiter,
                                                 unsigned miniterBC,
                                                 double percent_to_move) :
  MeshSmoother(mesh),
  _percent_to_move(percent_to_move),
  _dist_norm(0.),
  _adapt_data(adapt_data),
  _dim(mesh.mesh_dimension()),
  _miniter(miniter),
  _maxiter(maxiter),
  _miniterBC(miniterBC),
  _metric(UNIFORM),
  _adaptive_func(CELL),
  _theta(theta),
  _generate_data(false),
  _n_nodes(0),
  _n_cells(0),
  _n_hanging_edges(0),
  _area_of_interest(area_of_interest)
{}



double VariationalMeshSmoother::smooth(unsigned int)
{
  // If the log file is already open, for example on subsequent calls
  // to smooth() on the same object, we'll just keep writing to it,
  // otherwise we'll open it...
  if (!_logfile.is_open())
    _logfile.open("smoother.out");

  int
    me = _metric,
    gr = _generate_data ? 0 : 1,
    adp = _adaptive_func,
    miniter = _miniter,
    maxiter = _maxiter,
    miniterBC = _miniterBC;

  double theta = _theta;

  // Metric file name
  std::string metric_filename = "smoother.metric";
  if (gr == 0 && me > 1)
    {
      // grid filename
      std::string grid_filename = "smoother.grid";

      // generate metric from initial mesh (me = 2,3)
      metr_data_gen(grid_filename, metric_filename, me);
    }

  // Initialize the _n_nodes and _n_cells member variables
  this->_n_nodes = _mesh.n_nodes();
  this->_n_cells = _mesh.n_active_elem();

  // Initialize the _n_hanging_edges member variable
  MeshTools::find_hanging_nodes_and_parents(_mesh, _hanging_nodes);
  this->_n_hanging_edges =
    cast_int<dof_id_type>(_hanging_nodes.size());

  std::vector<int>
    mask(_n_nodes),
    edges(2*_n_hanging_edges),
    mcells(_n_cells),
    hnodes(_n_hanging_edges);

  Array2D<double> R(_n_nodes, _dim);
  Array2D<int> cells(_n_cells, 3*_dim + _dim%2);
  Array3D<double> H(_n_cells, _dim, _dim);

  // initial grid
  int vms_err = readgr(R, mask, cells, mcells, edges, hnodes);
  if (vms_err < 0)
    {
      _logfile << "Error reading input mesh file" << std::endl;
      return _dist_norm;
    }

  if (me > 1)
    vms_err = readmetr(metric_filename, H);

  if (vms_err < 0)
    {
      _logfile << "Error reading metric file" << std::endl;
      return _dist_norm;
    }

  std::vector<int> iter(4);
  iter[0] = miniter;
  iter[1] = maxiter;
  iter[2] = miniterBC;

  // grid optimization
  _logfile << "Starting Grid Optimization" << std::endl;
  clock_t ticks1 = clock();
  full_smooth(R, mask, cells, mcells, edges, hnodes, theta, iter, me, H, adp, gr);
  clock_t ticks2 = clock();
  _logfile << "full_smooth took ("
           << ticks2
           << "-"
           << ticks1
           << ")/"
           << CLOCKS_PER_SEC
           << " = "
           << static_cast<double>(ticks2-ticks1)/static_cast<double>(CLOCKS_PER_SEC)
           << " seconds"
           << std::endl;

  // save result
  _logfile << "Saving Result" << std::endl;
  writegr(R);

  libmesh_assert_greater (_dist_norm, 0);
  return _dist_norm;
}



// save grid
int VariationalMeshSmoother::writegr(const Array2D<double> & R)
{
  libMesh::out << "Starting writegr" << std::endl;

  // Adjust nodal coordinates to new positions
  {
    MeshBase::node_iterator       it  = _mesh.nodes_begin();
    const MeshBase::node_iterator end = _mesh.nodes_end();

    libmesh_assert_equal_to(_dist_norm, 0.);
    _dist_norm = 0;
    for (int i=0; it!=end; ++it)
      {
        double total_dist = 0.;

        // For each node set its X Y [Z] coordinates
        for (unsigned int j=0; j<_dim; j++)
          {
            // Get a reference to the node
            Node & node = *(*it);

            double distance = R[i][j] - node(j);

            // Save the squares of the distance
            total_dist += Utility::pow<2>(distance);

            node(j) += distance*_percent_to_move;
          }

        libmesh_assert_greater_equal (total_dist, 0.);

        // Add the distance this node moved to the global distance
        _dist_norm += total_dist;

        i++;
      }

    // Relative "error"
    _dist_norm = std::sqrt(_dist_norm/_mesh.n_nodes());
  }

  libMesh::out << "Finished writegr" << std::endl;
  return 0;
}



// reading grid from input file
int VariationalMeshSmoother::readgr(Array2D<double> & R,
                                    std::vector<int> & mask,
                                    Array2D<int> & cells,
                                    std::vector<int> & mcells,
                                    std::vector<int> & edges,
                                    std::vector<int> & hnodes)
{
  libMesh::out << "Sarting readgr" << std::endl;
  // add error messages where format can be inconsistent

  // Find the boundary nodes
  std::vector<bool> on_boundary;
  MeshTools::find_boundary_nodes(_mesh, on_boundary);

  // Grab node coordinates and set mask
  {
    MeshBase::const_node_iterator       it  = _mesh.nodes_begin();
    const MeshBase::const_node_iterator end = _mesh.nodes_end();

    // Only compute the node to elem map once
    std::vector<std::vector<const Elem *> > nodes_to_elem_map;
    MeshTools::build_nodes_to_elem_map(_mesh, nodes_to_elem_map);

    for (int i=0; it != end; ++it)
      {
        // Get a reference to the node
        Node & node = *(*it);

        // For each node grab its X Y [Z] coordinates
        for (unsigned int j=0; j<_dim; j++)
          R[i][j] = node(j);

        // Set the Proper Mask
        // Internal nodes are 0
        // Immovable boundary nodes are 1
        // Movable boundary nodes are 2
        if (on_boundary[i])
          {
            // Only look for sliding edge nodes in 2D
            if (_dim == 2)
              {
                // Find all the nodal neighbors... that is the nodes directly connected
                // to this node through one edge
                std::vector<const Node *> neighbors;
                MeshTools::find_nodal_neighbors(_mesh, node, nodes_to_elem_map, neighbors);

                std::vector<const Node *>::const_iterator ne = neighbors.begin();
                std::vector<const Node *>::const_iterator ne_end = neighbors.end();

                // Grab the x,y coordinates
                Real x = node(0);
                Real y = node(1);

                // Theta will represent the atan2 angle (meaning with the proper quadrant in mind)
                // of the neighbor node in a system where the current node is at the origin
                Real theta = 0;
                std::vector<Real> thetas;

                // Calculate the thetas
                for (; ne != ne_end; ne++)
                  {
                    const Node & neighbor = *(*ne);

                    // Note that the x and y values of this node are subtracted off
                    // this centers the system around this node
                    theta = atan2(neighbor(1)-y, neighbor(0)-x);

                    // Save it for later
                    thetas.push_back(theta);
                  }

                // Assume the node is immovable... then prove otherwise
                mask[i] = 1;

                // Search through neighbor nodes looking for two that form a straight line with this node
                for (std::size_t a=0; a<thetas.size()-1; a++)
                  {
                    // Only try each pairing once
                    for (std::size_t b=a+1; b<thetas.size(); b++)
                      {
                        // Find if the two neighbor nodes angles are 180 degrees (pi) off of eachother (withing a tolerance)
                        // In order to make this a true movable boundary node... the two that forma  straight line with
                        // it must also be on the boundary
                        if (on_boundary[neighbors[a]->id()] &&
                            on_boundary[neighbors[b]->id()] &&
                            (
                             (std::abs(thetas[a] - (thetas[b] + (libMesh::pi))) < .001) ||
                             (std::abs(thetas[a] - (thetas[b] - (libMesh::pi))) < .001)
                             )
                            )
                          {
                            // if ((*(*it))(1) > 0.25 || (*(*it))(0) > .7 || (*(*it))(0) < .01)
                            mask[i] = 2;
                          }

                      }
                  }
              }
            else // In 3D set all boundary nodes to be fixed
              mask[i] = 1;
          }
        else
          mask[i] = 0;  // Internal Node

        // libMesh::out << "Node: " << i << "  Mask: " << mask[i] << std::endl;
        i++;
      }
  }

  // Grab the connectivity
  // FIXME: Generalize this!
  {
    MeshBase::const_element_iterator it  = _mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = _mesh.active_elements_end();

    for (int i=0; it != end; ++it)
      {
        const Elem * elem = *it;

        // Keep track of the number of nodes
        // there must be 6 for 2D
        // 10 for 3D
        // If there are actually less than that -1 must be used
        // to fill out the rest
        int num = 0;
        /*
          int num_necessary = 0;

          if (_dim == 2)
          num_necessary = 6;
          else
          num_necessary = 10;
        */

        if (_dim == 2)
          {
            switch (elem->n_vertices())
              {
                // Grab nodes that do exist
              case 3:  // Tri
                for (unsigned int k=0; k<elem->n_vertices(); k++)
                  cells[i][k] = elem->node_id(k);

                num = elem->n_vertices();
                break;

              case 4:  // Quad 4
                cells[i][0] = elem->node_id(0);
                cells[i][1] = elem->node_id(1);
                cells[i][2] = elem->node_id(3); // Note that 2 and 3 are switched!
                cells[i][3] = elem->node_id(2);
                num = 4;
                break;

              default:
                libmesh_error_msg("Unknown number of vertices = " << elem->n_vertices());
              }
          }
        else
          {
            // Grab nodes that do exist
            switch (elem->n_vertices())
              {
                // Tet 4
              case 4:
                for (unsigned int k=0; k<elem->n_vertices(); k++)
                  cells[i][k] = elem->node_id(k);
                num = elem->n_vertices();
                break;

                // Hex 8
              case 8:
                cells[i][0] = elem->node_id(0);
                cells[i][1] = elem->node_id(1);
                cells[i][2] = elem->node_id(3); // Note that 2 and 3 are switched!
                cells[i][3] = elem->node_id(2);

                cells[i][4] = elem->node_id(4);
                cells[i][5] = elem->node_id(5);
                cells[i][6] = elem->node_id(7); // Note that 6 and 7 are switched!
                cells[i][7] = elem->node_id(6);
                num=8;
                break;

              default:
                libmesh_error_msg("Unknown 3D element with = " << elem->n_vertices() << " vertices.");
              }
          }

        // Fill in the rest with -1
        for (int j=num; j<static_cast<int>(cells[i].size()); j++)
          cells[i][j] = -1;

        // Mask it with 0 to state that this is an active element
        // FIXME: Could be something other than zero
        mcells[i] = 0;
        i++;
      }
  }

  // Grab hanging node connectivity
  {
    std::map<dof_id_type, std::vector<dof_id_type> >::iterator
      it = _hanging_nodes.begin(),
      end = _hanging_nodes.end();

    for (int i=0; it!=end; it++)
      {
        libMesh::out << "Parent 1: " << (it->second)[0] << std::endl;
        libMesh::out << "Parent 2: " << (it->second)[1] << std::endl;
        libMesh::out << "Hanging Node: " << it->first << std::endl << std::endl;

        // First Parent
        edges[2*i] = (it->second)[1];

        // Second Parent
        edges[2*i+1] = (it->second)[0];

        // Hanging Node
        hnodes[i] = it->first;

        i++;
      }
  }
  libMesh::out << "Finished readgr" << std::endl;

  return 0;
}



// Read Metrics
int VariationalMeshSmoother::readmetr(std::string name,
                                      Array3D<double> & H)
{
  std::ifstream infile(name.c_str());
  std::string dummy;

  for (dof_id_type i=0; i<_n_cells; i++)
    for (unsigned j=0; j<_dim; j++)
      {
        for (unsigned k=0; k<_dim; k++)
          infile >> H[i][j][k];

        // Read to end of line and discard
        std::getline(infile, dummy);
      }

  return 0;
}



// Stolen from ErrorVector!
float VariationalMeshSmoother::adapt_minimum() const
{
  float min = 1.e30;

  for (std::size_t i=0; i<_adapt_data->size(); i++)
    {
      // Only positive (or zero) values in the error vector
      libmesh_assert_greater_equal ((*_adapt_data)[i], 0.);
      min = std::min (min, (*_adapt_data)[i]);
    }

  // ErrorVectors are for positive values
  libmesh_assert_greater_equal (min, 0.);

  return min;
}



void VariationalMeshSmoother::adjust_adapt_data()
{
  // For convenience
  const UnstructuredMesh & aoe_mesh = *_area_of_interest;
  std::vector<float> & adapt_data = *_adapt_data;

  float min = adapt_minimum();

  MeshBase::const_element_iterator       el     = _mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.elements_end();

  MeshBase::const_element_iterator       aoe_el     = aoe_mesh.elements_begin();
  const MeshBase::const_element_iterator aoe_end_el = aoe_mesh.elements_end();

  // Counter to keep track of which spot in adapt_data we're looking at
  for (int i=0; el!=end_el; el++)
    {
      // Only do this for active elements
      if (adapt_data[i])
        {
          Point centroid = (*el)->centroid();
          bool in_aoe = false;

          // See if the elements centroid lies in the aoe mesh
          // for (aoe_el=aoe_mesh.elements_begin(); aoe_el != aoe_end_el; ++aoe_el)
          // {
          if ((*aoe_el)->contains_point(centroid))
            in_aoe = true;
          // }

          // If the element is not in the area of interest... then set
          // it's adaptivity value to the minimum.
          if (!in_aoe)
            adapt_data[i] = min;
        }
      i++;
    }
}



// Read Adaptivity
int VariationalMeshSmoother::read_adp(std::vector<double> & afun)
{
  std::vector<float> & adapt_data = *_adapt_data;

  if (_area_of_interest)
    adjust_adapt_data();

  std::size_t m = adapt_data.size();

  std::size_t j =0 ;
  for (std::size_t i=0; i<m; i++)
    if (adapt_data[i] != 0)
      {
        afun[j] = adapt_data[i];
        j++;
      }

  return 0;
}



double VariationalMeshSmoother::jac3(double x1,
                                     double y1,
                                     double z1,
                                     double x2,
                                     double y2,
                                     double z2,
                                     double x3,
                                     double y3,
                                     double z3)
{
  return x1*(y2*z3 - y3*z2) + y1*(z2*x3 - z3*x2) + z1*(x2*y3 - x3*y2);
}



double VariationalMeshSmoother::jac2(double x1,
                                     double y1,
                                     double x2,
                                     double y2)
{
  return x1*y2 - x2*y1;
}



// BasisA determines matrix H^(-T)Q on one Jacobian matrix
int VariationalMeshSmoother::basisA(Array2D<double> & Q,
                                    int nvert,
                                    const std::vector<double> & K,
                                    const Array2D<double> & H,
                                    int me)
{
  Array2D<double> U(_dim, nvert);

  // Some useful constants
  const double
    sqrt3 = std::sqrt(3.),
    sqrt6 = std::sqrt(6.);

  if (_dim == 2)
    {
      if (nvert == 4)
        {
          // quad
          U[0][0] = -(1-K[1]);
          U[0][1] =  (1-K[1]);
          U[0][2] = -K[1];
          U[0][3] =  K[1];

          U[1][0] = -(1-K[0]);
          U[1][1] = -K[0];
          U[1][2] =  (1-K[0]);
          U[1][3] =  K[0];
        }
      else if (nvert == 3)
        {
          // tri
          // for right target triangle
          // U[0][0] = -1.; U[1][0] = -1.;
          // U[0][1] = 1.;  U[1][1] = 0.;
          // U[0][2] = 0.;  U[1][2] = 1.;

          // for regular triangle
          U[0][0] = -1.;
          U[0][1] =  1.;
          U[0][2] =    0;

          U[1][0] = -1./sqrt3;
          U[1][1] = -1./sqrt3;
          U[1][2] =  2./sqrt3;
        }
      else if (nvert == 6)
        {
          // curved triangle
          U[0][0] = -3*(1-K[0]-K[1]*0.5)+(K[0]-0.5*K[1])+K[1];
          U[0][1] = -(1-K[0]-K[1]*0.5)+(K[0]-0.5*K[1])*3-K[1];
          U[0][2] = 0;
          U[0][3] = K[1]*4;
          U[0][4] = -4*K[1];
          U[0][5] = 4*(1-K[0]-K[1]*0.5)-4*(K[0]-0.5*K[1]);

          U[1][0] = (1-K[2]-K[3]*0.5)*(-1.5)+(K[2]-0.5*K[3])*(0.5)+K[3]*(0.5);
          U[1][1] = (1-K[2]-K[3]*0.5)*(0.5)+(K[2]-0.5*K[3])*(-1.5)+K[3]*(0.5);
          U[1][2] = (1-K[2]-K[3]*0.5)*(-1)+(K[2]-0.5*K[3])*(-1)+K[3]*(3);
          U[1][3] = (K[2]-0.5*K[3])*(4.)+K[3]*(-2.);
          U[1][4] = (1-K[2]-K[3]*0.5)*(4.)+K[3]*(-2.);
          U[1][5] = (1-K[2]-K[3]*0.5)*(-2.)+(K[2]-0.5*K[3])*(-2.);
        }
      else
        libmesh_error_msg("Invalid nvert = " << nvert);
    }

  if (_dim == 3)
    {
      if (nvert == 8)
        {
          // hex
          U[0][0] = -(1-K[0])*(1-K[1]);
          U[0][1] =  (1-K[0])*(1-K[1]);
          U[0][2] = -K[0]*(1-K[1]);
          U[0][3] =  K[0]*(1-K[1]);
          U[0][4] = -(1-K[0])*K[1];
          U[0][5] =  (1-K[0])*K[1];
          U[0][6] = -K[0]*K[1];
          U[0][7] =  K[0]*K[1];

          U[1][0] = -(1-K[2])*(1-K[3]);
          U[1][1] = -K[2]*(1-K[3]);
          U[1][2] =  (1-K[2])*(1-K[3]);
          U[1][3] =  K[2]*(1-K[3]);
          U[1][4] = -(1-K[2])*K[3];
          U[1][5] = -K[2]*K[3];
          U[1][6] =  (1-K[2])*K[3];
          U[1][7] =  K[2]*K[3];

          U[2][0] = -(1-K[4])*(1-K[5]);
          U[2][1] = -K[4]*(1-K[5]);
          U[2][2] = -(1-K[4])*K[5];
          U[2][3] = -K[4]*K[5];
          U[2][4] =  (1-K[4])*(1-K[5]);
          U[2][5] =  K[4]*(1-K[5]);
          U[2][6] =  (1-K[4])*K[5];
          U[2][7] =  K[4]*K[5];
        }
      else if (nvert == 4)
        {
          // linear tetr
          // for regular reference tetrahedron
          U[0][0] = -1;
          U[0][1] = 1;
          U[0][2] = 0;
          U[0][3] = 0;

          U[1][0] = -1./sqrt3;
          U[1][1] = -1./sqrt3;
          U[1][2] = 2./sqrt3;
          U[1][3] = 0;

          U[2][0] = -1./sqrt6;
          U[2][1] = -1./sqrt6;
          U[2][2] = -1./sqrt6;
          U[2][3] = 3./sqrt6;

          // for right corner reference tetrahedron
          // U[0][0] = -1; U[1][0] = -1; U[2][0] = -1;
          // U[0][1] = 1;  U[1][1] = 0;  U[2][1] = 0;
          // U[0][2] = 0;  U[1][2] = 1;  U[2][2] = 0;
          // U[0][3] = 0;  U[1][3] = 0;  U[2][3] = 1;
        }
      else if (nvert == 6)
        {
          // prism
          // for regular triangle in the prism base
          U[0][0] = -1+K[0];
          U[0][1] = 1-K[0];
          U[0][2] = 0;
          U[0][3] = -K[0];
          U[0][4] = K[0];
          U[0][5] = 0;

          U[1][0] = 0.5*(-1+K[1]);
          U[1][1] = 0.5*(-1+K[1]);
          U[1][2] = (1-K[1]);
          U[1][3] = -0.5*K[1];
          U[1][4] = -0.5*K[1];
          U[1][5] = K[1];

          U[2][0] = -1. + K[2] + 0.5*K[3];
          U[2][1] = -K[2] + 0.5*K[3];
          U[2][2] = -K[3];
          U[2][3] = 1 - K[2] - 0.5*K[3];
          U[2][4] = K[2] - 0.5*K[3];
          U[2][5] = K[3];
        }
      else if (nvert == 10)
        {
          // quad tetr
          U[0][0] = -3.*(1 - K[0] - 0.5*K[1] - K[2]/3.) + (K[0] - 0.5*K[1] - K[2]/3.)    + (K[1] - K[2]/3.) + K[2];
          U[0][1] = -1.*(1 - K[0] - 0.5*K[1] - K[2]/3.) + 3.*(K[0] - 0.5*K[1] - K[2]/3.) - (K[1] - K[2]/3.) - K[2];
          U[0][2] = 0.;
          U[0][3] = 0.;
          U[0][4] = 4.*(K[1] - K[2]/3.);
          U[0][5] = -4.*(K[1] - K[2]/3.);
          U[0][6] = 4.*(1. - K[0] - 0.5*K[1] - K[2]/3.) - 4.*(K[0] - 0.5*K[1] - K[2]/3.);
          U[0][7] = 4.*K[2];
          U[0][8] = 0.;
          U[0][9] = -4.*K[2];

          U[1][0] = -1.5*(1. - K[3] - K[4]*0.5 - K[5]/3.) + 0.5*(K[3] - 0.5*K[4] - K[5]/3.) + 0.5*(K[4] - K[5]/3.) + 0.5*K[5];
          U[1][1] =  0.5*(1. - K[3] - K[4]*0.5 - K[5]/3.) - 1.5*(K[3] - 0.5*K[4] - K[5]/3.) + 0.5*(K[4] - K[5]/3.) + 0.5*K[5];
          U[1][2] = -1.*(1. - K[3] - K[4]*0.5 - K[5]/3.) - (K[3] - 0.5*K[4] - K[5]/3.) + 3.*(K[4] - K[5]/3.) - K[5];
          U[1][3] = 0.;
          U[1][4] =  4.*(     K[3] - 0.5*K[4] - K[5]/3.) - 2.*(K[4] - K[5]/3.);
          U[1][5] =  4.*(1. - K[3] - 0.5*K[4] - K[5]/3.) - 2.*(K[4] - K[5]/3.);
          U[1][6] = -2.*(1. - K[3] - 0.5*K[4] - K[5]/3.) - 2.*(K[3] - 0.5*K[4] - K[5]/3.);
          U[1][7] = -2.*K[5];
          U[1][8] = 4.*K[5];
          U[1][9] = -2.*K[5];

          U[2][0] = -(1. - K[6] - 0.5*K[7] - K[8]/3.)    + (K[6] - 0.5*K[7] - K[8]/3.)/3. + (K[7] - K[8]/3.)/3. + K[8]/3.;
          U[2][1] =  (1. - K[6] - 0.5*K[7] - K[8]/3.)/3. - (K[6] - 0.5*K[7] - K[8]/3.)    + (K[7] - K[8]/3.)/3. + K[8]/3.;
          U[2][2] =  (1. - K[6] - 0.5*K[7] - K[8]/3.)/3. + (K[6] - 0.5*K[7] - K[8]/3.)/3. - (K[7] - K[8]/3.)    + K[8]/3.;
          U[2][3] = -(1. - K[6] - 0.5*K[7] - K[8]/3.)    - (K[6] - 0.5*K[7] - K[8]/3.)    - (K[7] - K[8]/3.)    + 3.*K[8];
          U[2][4] = -4.*(K[6] - K[7]*0.5 - K[8]/3.)/3. - 4.*(K[7] - K[8]/3.)/3.;
          U[2][5] = -4.*(1. - K[6] - K[7]*0.5 - K[8]/3.)/3. - 4.*(           K[7] - K[8]/3.)/3.;
          U[2][6] = -4.*(1. - K[6] - K[7]*0.5 - K[8]/3.)/3. - 4.*(K[6] - 0.5*K[7] - K[8]/3.)/3.;
          U[2][7] =  4.*(K[6] - K[7]*0.5 - K[8]/3.) - 4.*K[8]/3.;
          U[2][8] =  4.*(K[7] - K[8]/3.) - 4.*K[8]/3.;
          U[2][9] =  4.*(1. - K[6] - K[7]*0.5 - K[8]/3.) - 4.*K[8]/3.;
        }
      else
        libmesh_error_msg("Invalid nvert = " << nvert);
    }

  if (me == 1)
    Q = U;

  else
    {
      for (unsigned i=0; i<_dim; i++)
        for (int j=0; j<nvert; j++)
          {
            Q[i][j] = 0;
            for (unsigned k=0; k<_dim; k++)
              Q[i][j] += H[k][i]*U[k][j];
          }
    }

  return 0;
}



// Specify adaptive function
void VariationalMeshSmoother::adp_renew(const Array2D<double> & R,
                                        const Array2D<int> & cells,
                                        std::vector<double> & afun,
                                        int adp)
{
  // evaluates given adaptive function on the provided mesh
  if (adp < 0)
    {
      for (dof_id_type i=0; i<_n_cells; i++)
        {
          double
            x = 0.,
            y = 0.,
            z = 0.;
          int nvert = 0;

          while (cells[i][nvert] >= 0)
            {
              x += R[cells[i][nvert]][0];
              y += R[cells[i][nvert]][1];
              if (_dim == 3)
                z += R[cells[i][nvert]][2];
              nvert++;
            }

          // adaptive function, cell based
          afun[i] = 5.*(x+y+z);
        }
    }

  else if (adp > 0)
    {
      for (dof_id_type i=0; i<_n_nodes; i++)
        {
          double z = 0;
          for (unsigned j=0; j<_dim; j++)
            z += R[i][j];

          // adaptive function, node based
          afun[i] = 5*std::sin(R[i][0]);
        }
    }
}



// Preprocess mesh data and control smoothing/untangling iterations
void VariationalMeshSmoother::full_smooth(Array2D<double> & R,
                                          const std::vector<int> & mask,
                                          const Array2D<int> & cells,
                                          const std::vector<int> & mcells,
                                          const std::vector<int> & edges,
                                          const std::vector<int> & hnodes,
                                          double w,
                                          const std::vector<int> & iter,
                                          int me,
                                          const Array3D<double> & H,
                                          int adp,
                                          int gr)
{
  // Control the amount of print statements in this funcion
  int msglev = 1;

  dof_id_type afun_size = 0;

  // Adaptive function is on cells
  if (adp < 0)
    afun_size = _n_cells;

  // Adaptive function is on nodes
  else if (adp > 0)
    afun_size = _n_nodes;

  std::vector<double> afun(afun_size);
  std::vector<int> maskf(_n_nodes);
  std::vector<double> Gamma(_n_cells);

  if (msglev >= 1)
    _logfile << "N=" << _n_nodes
             << " ncells=" << _n_cells
             << " nedges=" << _n_hanging_edges
             << std::endl;


  // Boundary node counting
  int NBN=0;
  for (dof_id_type i=0; i<_n_nodes; i++)
    if (mask[i] == 2 || mask[i] == 1)
      NBN++;

  if (NBN > 0)
    {
      if (msglev >= 1)
        _logfile << "# of Boundary Nodes=" << NBN << std::endl;

      NBN=0;
      for (dof_id_type i=0; i<_n_nodes; i++)
        if (mask[i] == 2)
          NBN++;

      if (msglev >= 1)
        _logfile << "# of moving Boundary Nodes=" << NBN << std::endl;
    }

  for (dof_id_type i=0; i<_n_nodes; i++)
    {
      if ((mask[i]==1) || (mask[i]==2))
        maskf[i] = 1;
      else
        maskf[i] = 0;
    }

  // determination of min jacobian
  double vol, Vmin;
  double qmin = minq(R, cells, mcells, me, H, vol, Vmin);

  if (me > 1)
    vol = 1.;

  if (msglev >= 1)
    _logfile << "vol=" << vol
             << " qmin=" << qmin
             << " min volume = " << Vmin
             << std::endl;

  // compute max distortion measure over all cells
  double epsilon = 1.e-9;
  double eps = qmin < 0 ? std::sqrt(epsilon*epsilon+0.004*qmin*qmin*vol*vol) : epsilon;
  double emax = maxE(R, cells, mcells, me, H, vol, eps, w, Gamma, qmin);

  if (msglev >= 1)
    _logfile << " emax=" << emax << std::endl;

  // unfolding/smoothing

  // iteration tolerance
  double Enm1 = 1.;

  // read adaptive function from file
  if (adp*gr != 0)
    read_adp(afun);

  {
    int
      counter = 0,
      ii = 0;

    while (((qmin <= 0) || (counter < iter[0]) || (std::abs(emax-Enm1) > 1e-3)) &&
           (ii < iter[1]) &&
           (counter < iter[1]))
      {
        libmesh_assert_less (counter, iter[1]);

        Enm1 = emax;

        if ((ii >= 0) && (ii % 2 == 0))
          {
            if (qmin < 0)
              eps = std::sqrt(epsilon*epsilon + 0.004*qmin*qmin*vol*vol);
            else
              eps = epsilon;
          }

        int ladp = adp;

        if ((qmin <= 0) || (counter < ii))
          ladp = 0;

        // update adaptation function before each iteration
        if ((ladp != 0) && (gr == 0))
          adp_renew(R, cells, afun, adp);

        double Jk = minJ(R, maskf, cells, mcells, eps, w, me, H, vol, edges, hnodes,
                         msglev, Vmin, emax, qmin, ladp, afun);

        if (qmin > 0)
          counter++;
        else
          ii++;

        if (msglev >= 1)
          _logfile << "niter=" << counter
                   << ", qmin*G/vol=" << qmin
                   << ", Vmin=" << Vmin
                   << ", emax=" << emax
                   << ", Jk=" << Jk
                   << ", Enm1=" << Enm1
                   << std::endl;
      }
  }

  // BN correction - 2D only!
  epsilon = 1.e-9;
  if (NBN > 0)
    for (int counter=0; counter<iter[2]; counter++)
      {
        // update adaptation function before each iteration
        if ((adp != 0) && (gr == 0))
          adp_renew(R, cells, afun, adp);

        double Jk = minJ_BC(R, mask, cells, mcells, eps, w, me, H, vol, msglev, Vmin, emax, qmin, adp, afun, NBN);

        if (msglev >= 1)
          _logfile << "NBC niter=" << counter
                   << ", qmin*G/vol=" << qmin
                   << ", Vmin=" << Vmin
                   << ", emax=" << emax
                   << std::endl;

        // Outrageous Enm1 to make sure we hit this at least once
        Enm1 = 99999;

        // Now that we've moved the boundary nodes (or not) we need to resmoooth
        for (int j=0; j<iter[1]; j++)
          {
            if (std::abs(emax-Enm1) < 1e-2)
              break;

            // Save off the error from the previous smoothing step
            Enm1 = emax;

            // update adaptation function before each iteration
            if ((adp != 0) && (gr == 0))
              adp_renew(R, cells, afun, adp);

            Jk = minJ(R, maskf, cells, mcells, eps, w, me, H, vol, edges, hnodes, msglev, Vmin, emax, qmin, adp, afun);

            if (msglev >= 1)
              _logfile << "  Re-smooth: niter=" << j
                       << ", qmin*G/vol=" << qmin
                       << ", Vmin=" << Vmin
                       << ", emax=" << emax
                       << ", Jk=" << Jk
                       << std::endl;
          }

        if (msglev >= 1)
          _logfile << "NBC smoothed niter=" << counter
                   << ", qmin*G/vol=" << qmin
                   << ", Vmin=" << Vmin
                   << ", emax=" << emax
                   << std::endl;
      }
}



// Determines the values of maxE_theta
double VariationalMeshSmoother::maxE(Array2D<double> & R,
                                     const Array2D<int> & cells,
                                     const std::vector<int> & mcells,
                                     int me,
                                     const Array3D<double> & H,
                                     double v,
                                     double epsilon,
                                     double w,
                                     std::vector<double> & Gamma,
                                     double & qmin)
{
  Array2D<double> Q(3, 3*_dim + _dim%2);
  std::vector<double> K(9);

  double
    gemax = -1.e32,
    vmin = 1.e32;

  for (dof_id_type ii=0; ii<_n_cells; ii++)
    if (mcells[ii] >= 0)
      {
        // Final value of E will be saved in Gamma at the end of this loop
        double E = 0.;

        if (_dim == 2)
          {
            if (cells[ii][3] == -1)
              {
                // tri
                basisA(Q, 3, K, H[ii], me);

                std::vector<double> a1(3), a2(3);
                for (int k=0; k<2; k++)
                  for (int l=0; l<3; l++)
                    {
                      a1[k] += Q[k][l]*R[cells[ii][l]][0];
                      a2[k] += Q[k][l]*R[cells[ii][l]][1];
                    }

                double det = jac2(a1[0], a1[1], a2[0], a2[1]);
                double tr = 0.5*(a1[0]*a1[0] + a2[0]*a2[0] + a1[1]*a1[1] + a2[1]*a2[1]);
                double chi = 0.5*(det+std::sqrt(det*det+epsilon*epsilon));
                E = (1-w)*tr/chi + 0.5*w*(v + det*det/v)/chi;

                if (E > gemax)
                  gemax = E;
                if (vmin > det)
                  vmin = det;

              }
            else if (cells[ii][4] == -1)
              {
                // quad
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;
                        basisA(Q, 4, K, H[ii], me);

                        std::vector<double> a1(3), a2(3);
                        for (int k=0; k<2; k++)
                          for (int l=0; l<4; l++)
                            {
                              a1[k] += Q[k][l]*R[cells[ii][l]][0];
                              a2[k] += Q[k][l]*R[cells[ii][l]][1];
                            }

                        double det = jac2(a1[0],a1[1],a2[0],a2[1]);
                        double tr = 0.5*(a1[0]*a1[0] + a2[0]*a2[0] + a1[1]*a1[1] + a2[1]*a2[1]);
                        double chi = 0.5*(det+std::sqrt(det*det+epsilon*epsilon));
                        E += 0.25*((1-w)*tr/chi + 0.5*w*(v + det*det/v)/chi);

                        if (vmin > det)
                          vmin = det;
                      }
                  }

                if (E > gemax)
                  gemax = E;
              }
            else
              {
                // quad tri
                for (int i=0; i<3; i++)
                  {
                    K[0] = i*0.5;
                    int k = i/2;
                    K[1] = static_cast<double>(k);

                    for (int j=0; j<3; j++)
                      {
                        K[2] = j*0.5;
                        k = j/2;
                        K[3] = static_cast<double>(k);

                        basisA(Q, 6, K, H[ii], me);

                        std::vector<double> a1(3), a2(3);
                        for (int k=0; k<2; k++)
                          for (int l=0; l<6; l++)
                            {
                              a1[k] += Q[k][l]*R[cells[ii][l]][0];
                              a2[k] += Q[k][l]*R[cells[ii][l]][1];
                            }

                        double det = jac2(a1[0],a1[1],a2[0],a2[1]);
                        double sigma = 1./24.;

                        if (i==j)
                          sigma = 1./12.;

                        double tr = 0.5*(a1[0]*a1[0] + a2[0]*a2[0] + a1[1]*a1[1] + a2[1]*a2[1]);
                        double chi = 0.5*(det + std::sqrt(det*det + epsilon*epsilon));
                        E += sigma*((1-w)*tr/chi + 0.5*w*(v + det*det/v)/chi);
                        if (vmin > det)
                          vmin = det;
                      }
                  }

                if (E > gemax)
                  gemax = E;
              }
          }

        if (_dim == 3)
          {
            if (cells[ii][4] == -1)
              {
                // tetr
                basisA(Q, 4, K, H[ii], me);

                std::vector<double> a1(3), a2(3), a3(3);
                for (int k=0; k<3; k++)
                  for (int l=0; l<4; l++)
                    {
                      a1[k] += Q[k][l]*R[cells[ii][l]][0];
                      a2[k] += Q[k][l]*R[cells[ii][l]][1];
                      a3[k] += Q[k][l]*R[cells[ii][l]][2];
                    }

                double det = jac3(a1[0], a1[1], a1[2],
                                  a2[0], a2[1], a2[2],
                                  a3[0], a3[1], a3[2]);
                double tr = 0.;
                for (int k=0; k<3; k++)
                  tr += (a1[k]*a1[k] + a2[k]*a2[k] + a3[k]*a3[k])/3.;

                double chi = 0.5*(det+std::sqrt(det*det+epsilon*epsilon));
                E = (1-w)*pow(tr,1.5)/chi + 0.5*w*(v + det*det/v)/chi;

                if (E > gemax)
                  gemax = E;

                if (vmin > det)
                  vmin = det;
              }
            else if (cells[ii][6] == -1)
              {
                // prism
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;
                        for (int k=0; k<3; k++)
                          {
                            K[2] = 0.5*static_cast<double>(k);
                            K[3] = static_cast<double>(k % 2);
                            basisA(Q, 6, K, H[ii], me);

                            std::vector<double> a1(3), a2(3), a3(3);
                            for (int kk=0; kk<3; kk++)
                              for (int ll=0; ll<6; ll++)
                                {
                                  a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                  a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                  a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                }

                            double det = jac3(a1[0], a1[1], a1[2],
                                              a2[0], a2[1], a2[2],
                                              a3[0], a3[1], a3[2]);
                            double tr = 0;
                            for (int kk=0; kk<3; kk++)
                              tr += (a1[kk]*a1[kk] + a2[kk]*a2[kk] + a3[kk]*a3[kk])/3.;

                            double chi = 0.5*(det+std::sqrt(det*det+epsilon*epsilon));
                            E += ((1-w)*pow(tr,1.5)/chi + 0.5*w*(v + det*det/v)/chi)/12.;
                            if (vmin > det)
                              vmin = det;
                          }
                      }
                  }

                if (E > gemax)
                  gemax = E;
              }
            else if (cells[ii][8] == -1)
              {
                // hex
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;
                        for (int k=0; k<2; k++)
                          {
                            K[2] = k;
                            for (int l=0; l<2; l++)
                              {
                                K[3] = l;
                                for (int m=0; m<2; m++)
                                  {
                                    K[4] = m;
                                    for (int nn=0; nn<2; nn++)
                                      {
                                        K[5] = nn;
                                        basisA(Q, 8, K, H[ii], me);

                                        std::vector<double> a1(3), a2(3), a3(3);
                                        for (int kk=0; kk<3; kk++)
                                          for (int ll=0; ll<8; ll++)
                                            {
                                              a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                              a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                              a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                            }

                                        double det = jac3(a1[0], a1[1], a1[2],
                                                          a2[0], a2[1], a2[2],
                                                          a3[0], a3[1], a3[2]);
                                        double sigma = 0.;

                                        if ((i==nn) && (j==l) && (k==m))
                                          sigma = 1./27.;

                                        if (((i==nn) && (j==l) && (k!=m)) ||
                                            ((i==nn) && (j!=l) && (k==m)) ||
                                            ((i!=nn) && (j==l) && (k==m)))
                                          sigma = 1./54.;

                                        if (((i==nn) && (j!=l) && (k!=m)) ||
                                            ((i!=nn) && (j!=l) && (k==m)) ||
                                            ((i!=nn) && (j==l) && (k!=m)))
                                          sigma = 1./108.;

                                        if ((i!=nn) && (j!=l) && (k!=m))
                                          sigma = 1./216.;

                                        double tr = 0;
                                        for (int kk=0; kk<3; kk++)
                                          tr += (a1[kk]*a1[kk] + a2[kk]*a2[kk] + a3[kk]*a3[kk])/3.;

                                        double chi = 0.5*(det+std::sqrt(det*det + epsilon*epsilon));
                                        E += ((1-w)*pow(tr,1.5)/chi + 0.5*w*(v + det*det/v)/chi)*sigma;
                                        if (vmin > det)
                                          vmin = det;
                                      }
                                  }
                              }
                          }
                      }
                  }
                if (E > gemax)
                  gemax = E;
              }
            else
              {
                // quad tetr
                for (int i=0; i<4; i++)
                  {
                    for (int j=0; j<4; j++)
                      {
                        for (int k=0; k<4; k++)
                          {
                            switch (i)
                              {
                              case 0:
                                K[0] = 0;
                                K[1] = 0;
                                K[2] = 0;
                                break;
                              case 1:
                                K[0] = 1;
                                K[1] = 0;
                                K[2] = 0;
                                break;
                              case 2:
                                K[0] = 0.5;
                                K[1] = 1;
                                K[2] = 0;
                                break;
                              case 3:
                                K[0] = 0.5;
                                K[1] = 1./3.;
                                K[2] = 1;
                                break;
                              }

                            switch (j)
                              {
                              case 0:
                                K[3] = 0;
                                K[4] = 0;
                                K[5] = 0;
                                break;
                              case 1:
                                K[3] = 1;
                                K[4] = 0;
                                K[5] = 0;
                                break;
                              case 2:
                                K[3] = 0.5;
                                K[4] = 1;
                                K[5] = 0;
                                break;
                              case 3:
                                K[3] = 0.5;
                                K[4] = 1./3.;
                                K[5] = 1;
                                break;
                              }

                            switch (k)
                              {
                              case 0:
                                K[6] = 0;
                                K[7] = 0;
                                K[8] = 0;
                                break;
                              case 1:
                                K[6] = 1;
                                K[7] = 0;
                                K[8] = 0;
                                break;
                              case 2:
                                K[6] = 0.5;
                                K[7] = 1;
                                K[8] = 0;
                                break;
                              case 3:
                                K[6] = 0.5;
                                K[7] = 1./3.;
                                K[8] = 1;
                                break;
                              }

                            basisA(Q, 10, K, H[ii], me);

                            std::vector<double> a1(3), a2(3), a3(3);
                            for (int kk=0; kk<3; kk++)
                              for (int ll=0; ll<10; ll++)
                                {
                                  a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                  a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                  a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                }

                            double det = jac3(a1[0], a1[1], a1[2],
                                              a2[0], a2[1], a2[2],
                                              a3[0], a3[1], a3[2]);
                            double sigma = 0.;

                            if ((i==j) && (j==k))
                              sigma = 1./120.;
                            else if ((i==j) || (j==k) || (i==k))
                              sigma = 1./360.;
                            else
                              sigma = 1./720.;

                            double tr = 0;
                            for (int kk=0; kk<3; kk++)
                              tr += (a1[kk]*a1[kk] + a2[kk]*a2[kk] + a3[kk]*a3[kk])/3.;

                            double chi = 0.5*(det+std::sqrt(det*det+epsilon*epsilon));
                            E += ((1-w)*pow(tr,1.5)/chi + 0.5*w*(v+det*det/v)/chi)*sigma;
                            if (vmin > det)
                              vmin = det;
                          }
                      }
                  }

                if (E > gemax)
                  gemax = E;
              }
          }
        Gamma[ii] = E;
      }

  qmin = vmin;

  return gemax;
}



// Compute min Jacobian determinant (minq), min cell volume (Vmin), and average cell volume (vol).
double VariationalMeshSmoother::minq(const Array2D<double> & R,
                                     const Array2D<int> & cells,
                                     const std::vector<int> & mcells,
                                     int me,
                                     const Array3D<double> & H,
                                     double & vol,
                                     double & Vmin)
{
  std::vector<double> K(9);
  Array2D<double> Q(3, 3*_dim + _dim%2);

  double v = 0;
  double vmin = 1.e32;
  double gqmin = 1.e32;

  for (dof_id_type ii=0; ii<_n_cells; ii++)
    if (mcells[ii] >= 0)
      {
        if (_dim == 2)
          {
            // 2D
            if (cells[ii][3] == -1)
              {
                // tri
                basisA(Q, 3, K, H[ii], me);

                std::vector<double> a1(3), a2(3);
                for (int k=0; k<2; k++)
                  for (int l=0; l<3; l++)
                    {
                      a1[k] += Q[k][l]*R[cells[ii][l]][0];
                      a2[k] += Q[k][l]*R[cells[ii][l]][1];
                    }

                double det = jac2(a1[0],a1[1],a2[0],a2[1]);
                if (gqmin > det)
                  gqmin = det;

                if (vmin > det)
                  vmin = det;

                v += det;
              }
            else if (cells[ii][4] == -1)
              {
                // quad
                double vcell = 0.;
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;
                        basisA(Q, 4, K, H[ii], me);

                        std::vector<double> a1(3), a2(3);
                        for (int k=0; k<2; k++)
                          for (int l=0; l<4; l++)
                            {
                              a1[k] += Q[k][l]*R[cells[ii][l]][0];
                              a2[k] += Q[k][l]*R[cells[ii][l]][1];
                            }

                        double det = jac2(a1[0],a1[1],a2[0],a2[1]);
                        if (gqmin > det)
                          gqmin = det;

                        v += 0.25*det;
                        vcell += 0.25*det;
                      }
                  }
                if (vmin > vcell)
                  vmin = vcell;
              }
            else
              {
                // quad tri
                double vcell = 0.;
                for (int i=0; i<3; i++)
                  {
                    K[0] = i*0.5;
                    int k = i/2;
                    K[1] = static_cast<double>(k);

                    for (int j=0; j<3; j++)
                      {
                        K[2] = j*0.5;
                        k = j/2;
                        K[3] = static_cast<double>(k);
                        basisA(Q, 6, K, H[ii], me);

                        std::vector<double> a1(3), a2(3);
                        for (int k=0; k<2; k++)
                          for (int l=0; l<6; l++)
                            {
                              a1[k] += Q[k][l]*R[cells[ii][l]][0];
                              a2[k] += Q[k][l]*R[cells[ii][l]][1];
                            }

                        double det = jac2(a1[0], a1[1], a2[0], a2[1]);
                        if (gqmin > det)
                          gqmin = det;

                        double sigma = 1./24.;
                        if (i == j)
                          sigma = 1./12.;

                        v += sigma*det;
                        vcell += sigma*det;
                      }
                  }
                if (vmin > vcell)
                  vmin = vcell;
              }
          }
        if (_dim == 3)
          {
            // 3D
            if (cells[ii][4] == -1)
              {
                // tetr
                basisA(Q, 4, K, H[ii], me);

                std::vector<double> a1(3), a2(3), a3(3);
                for (int k=0; k<3; k++)
                  for (int l=0; l<4; l++)
                    {
                      a1[k] += Q[k][l]*R[cells[ii][l]][0];
                      a2[k] += Q[k][l]*R[cells[ii][l]][1];
                      a3[k] += Q[k][l]*R[cells[ii][l]][2];
                    }

                double det = jac3(a1[0], a1[1], a1[2],
                                  a2[0], a2[1], a2[2],
                                  a3[0], a3[1], a3[2]);

                if (gqmin > det)
                  gqmin = det;

                if (vmin > det)
                  vmin = det;
                v += det;
              }
            else if (cells[ii][6] == -1)
              {
                // prism
                double vcell = 0.;
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;

                        for (int k=0; k<3; k++)
                          {
                            K[2] = 0.5*k;
                            K[3] = static_cast<double>(k%2);
                            basisA(Q, 6, K, H[ii], me);

                            std::vector<double> a1(3), a2(3), a3(3);
                            for (int kk=0; kk<3; kk++)
                              for (int ll=0; ll<6; ll++)
                                {
                                  a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                  a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                  a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                }

                            double det = jac3(a1[0], a1[1], a1[2],
                                              a2[0], a2[1], a2[2],
                                              a3[0], a3[1], a3[2]);
                            if (gqmin > det)
                              gqmin = det;

                            double sigma = 1./12.;
                            v += sigma*det;
                            vcell += sigma*det;
                          }
                      }
                  }
                if (vmin > vcell)
                  vmin = vcell;
              }
            else if (cells[ii][8] == -1)
              {
                // hex
                double vcell = 0.;
                for (int i=0; i<2; i++)
                  {
                    K[0] = i;
                    for (int j=0; j<2; j++)
                      {
                        K[1] = j;
                        for (int k=0; k<2; k++)
                          {
                            K[2] = k;
                            for (int l=0; l<2; l++)
                              {
                                K[3] = l;
                                for (int m=0; m<2; m++)
                                  {
                                    K[4] = m;
                                    for (int nn=0; nn<2; nn++)
                                      {
                                        K[5] = nn;
                                        basisA(Q, 8, K, H[ii], me);

                                        std::vector<double> a1(3), a2(3), a3(3);
                                        for (int kk=0; kk<3; kk++)
                                          for (int ll=0; ll<8; ll++)
                                            {
                                              a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                              a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                              a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                            }

                                        double det = jac3(a1[0], a1[1], a1[2],
                                                          a2[0], a2[1], a2[2],
                                                          a3[0], a3[1], a3[2]);

                                        if (gqmin > det)
                                          gqmin = det;

                                        double sigma = 0.;

                                        if ((i==nn) && (j==l) && (k==m))
                                          sigma = 1./27.;

                                        if (((i==nn) && (j==l) && (k!=m)) ||
                                            ((i==nn) && (j!=l) && (k==m)) ||
                                            ((i!=nn) && (j==l) && (k==m)))
                                          sigma = 1./54.;

                                        if (((i==nn) && (j!=l) && (k!=m)) ||
                                            ((i!=nn) && (j!=l) && (k==m)) ||
                                            ((i!=nn) && (j==l) && (k!=m)))
                                          sigma = 1./108.;

                                        if ((i!=nn) && (j!=l) && (k!=m))
                                          sigma = 1./216.;

                                        v += sigma*det;
                                        vcell += sigma*det;
                                      }
                                  }
                              }
                          }
                      }
                  }

                if (vmin > vcell)
                  vmin = vcell;
              }
            else
              {
                // quad tetr
                double vcell = 0.;
                for (int i=0; i<4; i++)
                  {
                    for (int j=0; j<4; j++)
                      {
                        for (int k=0; k<4; k++)
                          {
                            switch (i)
                              {
                              case 0:
                                K[0] = 0;
                                K[1] = 0;
                                K[2] = 0;
                                break;

                              case 1:
                                K[0] = 1;
                                K[1] = 0;
                                K[2] = 0;
                                break;

                              case 2:
                                K[0] = 0.5;
                                K[1] = 1;
                                K[2] = 0;
                                break;

                              case 3:
                                K[0] = 0.5;
                                K[1] = 1./3.;
                                K[2] = 1;
                                break;

                              }
                            switch (j)
                              {
                              case 0:
                                K[3] = 0;
                                K[4] = 0;
                                K[5] = 0;
                                break;

                              case 1:
                                K[3] = 1;
                                K[4] = 0;
                                K[5] = 0;
                                break;

                              case 2:
                                K[3] = 0.5;
                                K[4] = 1;
                                K[5] = 0;
                                break;

                              case 3:
                                K[3] = 0.5;
                                K[4] = 1./3.;
                                K[5] = 1;
                                break;

                              }
                            switch (k)
                              {
                              case 0:
                                K[6] = 0;
                                K[7] = 0;
                                K[8] = 0;
                                break;

                              case 1:
                                K[6] = 1;
                                K[7] = 0;
                                K[8] = 0;
                                break;

                              case 2:
                                K[6] = 0.5;
                                K[7] = 1;
                                K[8] = 0;
                                break;

                              case 3:
                                K[6] = 0.5;
                                K[7] = 1./3.;
                                K[8] = 1;
                                break;
                              }

                            basisA(Q, 10, K, H[ii], me);

                            std::vector<double> a1(3), a2(3), a3(3);
                            for (int kk=0; kk<3; kk++)
                              for (int ll=0; ll<10; ll++)
                                {
                                  a1[kk] += Q[kk][ll]*R[cells[ii][ll]][0];
                                  a2[kk] += Q[kk][ll]*R[cells[ii][ll]][1];
                                  a3[kk] += Q[kk][ll]*R[cells[ii][ll]][2];
                                }

                            double det = jac3(a1[0], a1[1], a1[2],
                                              a2[0], a2[1], a2[2],
                                              a3[0], a3[1], a3[2]);
                            if (gqmin > det)
                              gqmin = det;

                            double sigma = 0.;

                            if ((i==j) && (j==k))
                              sigma = 1./120.;

                            else if ((i==j) || (j==k) || (i==k))
                              sigma = 1./360.;

                            else
                              sigma = 1./720.;

                            v += sigma*det;
                            vcell += sigma*det;
                          }
                      }
                  }
                if (vmin > vcell)
                  vmin = vcell;
              }
          }
      }

  // Fill in return value references
  vol = v/static_cast<double>(_n_cells);
  Vmin = vmin;

  return gqmin;
}



// Executes one step of minimization algorithm:
// finds minimization direction (P=H^{-1} \grad J) and solves approximately
// local minimization problem for optimal step in this minimization direction (tau=min J(R+tau P))
double VariationalMeshSmoother::minJ(Array2D<double> & R,
                                     const std::vector<int> & mask,
                                     const Array2D<int> & cells,
                                     const std::vector<int> & mcells,
                                     double epsilon,
                                     double w,
                                     int me,
                                     const Array3D<double> & H,
                                     double vol,
                                     const std::vector<int> & edges,
                                     const std::vector<int> & hnodes,
                                     int msglev,
                                     double & Vmin,
                                     double & emax,
                                     double & qmin,
                                     int adp,
                                     const std::vector<double> & afun)
{
  // columns - max number of nonzero entries in every row of global matrix
  int columns = _dim*_dim*10;

  // local Hessian matrix
  Array3D<double> W(_dim, 3*_dim + _dim%2, 3*_dim + _dim%2);

  // F - local gradient
  Array2D<double> F(_dim, 3*_dim + _dim%2);

  Array2D<double> Rpr(_n_nodes, _dim);

  // P - minimization direction
  Array2D<double> P(_n_nodes, _dim);

  // A, JA - internal form of global matrix
  Array2D<int> JA(_dim*_n_nodes, columns);
  Array2D<double> A(_dim*_n_nodes, columns);

  // G - adaptation metric
  Array2D<double> G(_n_cells, _dim);

  // rhs for solver
  std::vector<double> b(_dim*_n_nodes);

  // u - solution vector
  std::vector<double> u(_dim*_n_nodes);

  // matrix
  std::vector<double> a(_dim*_n_nodes*columns);
  std::vector<int> ia(_dim*_n_nodes + 1);
  std::vector<int> ja(_dim*_n_nodes*columns);

  // nonzero - norm of gradient
  double nonzero = 0.;

  // Jpr - value of functional
  double Jpr = 0.;

  // find minimization direction P
  for (dof_id_type i=0; i<_n_cells; i++)
    {
      int nvert = 0;
      while (cells[i][nvert] >= 0)
        nvert++;

      // determination of local matrices on each cell
      for (unsigned j=0; j<_dim; j++)
        {
          G[i][j] = 0;  // adaptation metric G is held constant throughout minJ run
          if (adp < 0)
            {
              for (int k=0; k<std::abs(adp); k++)
                G[i][j] += afun[i*(-adp)+k];  // cell-based adaptivity is computed here
            }
        }
      for (unsigned index=0; index<_dim; index++)
        {
          // initialise local matrices
          for (unsigned k=0; k<3*_dim + _dim%2; k++)
            {
              F[index][k] = 0;

              for (unsigned j=0; j<3*_dim + _dim%2; j++)
                W[index][k][j] = 0;
            }
        }
      if (mcells[i] >= 0)
        {
          // if cell is not excluded
          double lVmin, lqmin;
          Jpr += localP(W, F, R, cells[i], mask, epsilon, w, nvert, H[i],
                        me, vol, 0, lVmin, lqmin, adp, afun, G[i]);
        }
      else
        {
          for (unsigned index=0; index<_dim; index++)
            for (int j=0; j<nvert; j++)
              W[index][j][j] = 1;
        }

      // assembly of an upper triangular part of a global matrix A
      for (unsigned index=0; index<_dim; index++)
        {
          for (int l=0; l<nvert; l++)
            {
              for (int m=0; m<nvert; m++)
                {
                  if ((W[index][l][m] != 0) &&
                      (cells[i][m] >= cells[i][l]))
                    {
                      int sch = 0;
                      int ind = 1;
                      while (ind != 0)
                        {
                          if (A[cells[i][l] + index*_n_nodes][sch] != 0)
                            {
                              if (JA[cells[i][l] + index*_n_nodes][sch] == static_cast<int>(cells[i][m] + index*_n_nodes))
                                {
                                  A[cells[i][l] + index*_n_nodes][sch] = A[cells[i][l] + index*_n_nodes][sch] + W[index][l][m];
                                  ind=0;
                                }
                              else
                                sch++;
                            }
                          else
                            {
                              A[cells[i][l] + index*_n_nodes][sch] = W[index][l][m];
                              JA[cells[i][l] + index*_n_nodes][sch] = cells[i][m] + index*_n_nodes;
                              ind = 0;
                            }

                          if (sch > columns-1)
                            _logfile << "error: # of nonzero entries in the "
                                     << cells[i][l]
                                     << " row of Hessian ="
                                     << sch
                                     << ">= columns="
                                     << columns
                                     << std::endl;
                        }
                    }
                }
              b[cells[i][l] + index*_n_nodes] = b[cells[i][l] + index*_n_nodes] - F[index][l];
            }
        }
      // end of matrix A
    }

  // HN correction

  // tolerance for HN being the mid-edge node for its parents
  double Tau_hn = pow(vol, 1./static_cast<double>(_dim))*1e-10;
  for (dof_id_type i=0; i<_n_hanging_edges; i++)
    {
      int ind_i = hnodes[i];
      int ind_j = edges[2*i];
      int ind_k = edges[2*i+1];

      for (unsigned j=0; j<_dim; j++)
        {
          int g_i = R[ind_i][j] - 0.5*(R[ind_j][j]+R[ind_k][j]);
          Jpr += g_i*g_i/(2*Tau_hn);
          b[ind_i + j*_n_nodes] -= g_i/Tau_hn;
          b[ind_j + j*_n_nodes] += 0.5*g_i/Tau_hn;
          b[ind_k + j*_n_nodes] += 0.5*g_i/Tau_hn;
        }

      for (int ind=0; ind<columns; ind++)
        {
          if (JA[ind_i][ind] == ind_i)
            A[ind_i][ind] += 1./Tau_hn;

          if (JA[ind_i+_n_nodes][ind] == static_cast<int>(ind_i+_n_nodes))
            A[ind_i+_n_nodes][ind] += 1./Tau_hn;

          if (_dim == 3)
            if (JA[ind_i+2*_n_nodes][ind] == static_cast<int>(ind_i+2*_n_nodes))
              A[ind_i+2*_n_nodes][ind] += 1./Tau_hn;

          if ((JA[ind_i][ind] == ind_j) ||
              (JA[ind_i][ind] == ind_k))
            A[ind_i][ind] -= 0.5/Tau_hn;

          if ((JA[ind_i+_n_nodes][ind] == static_cast<int>(ind_j+_n_nodes)) ||
              (JA[ind_i+_n_nodes][ind] == static_cast<int>(ind_k+_n_nodes)))
            A[ind_i+_n_nodes][ind] -= 0.5/Tau_hn;

          if (_dim == 3)
            if ((JA[ind_i+2*_n_nodes][ind] == static_cast<int>(ind_j+2*_n_nodes)) ||
                (JA[ind_i+2*_n_nodes][ind] == static_cast<int>(ind_k+2*_n_nodes)))
              A[ind_i+2*_n_nodes][ind] -= 0.5/Tau_hn;

          if (JA[ind_j][ind] == ind_i)
            A[ind_j][ind] -= 0.5/Tau_hn;

          if (JA[ind_k][ind] == ind_i)
            A[ind_k][ind] -= 0.5/Tau_hn;

          if (JA[ind_j+_n_nodes][ind] == static_cast<int>(ind_i+_n_nodes))
            A[ind_j+_n_nodes][ind] -= 0.5/Tau_hn;

          if (JA[ind_k+_n_nodes][ind] == static_cast<int>(ind_i+_n_nodes))
            A[ind_k+_n_nodes][ind] -= 0.5/Tau_hn;

          if (_dim == 3)
            if (JA[ind_j+2*_n_nodes][ind] == static_cast<int>(ind_i+2*_n_nodes))
              A[ind_j+2*_n_nodes][ind] -= 0.5/Tau_hn;

          if (_dim == 3)
            if (JA[ind_k+2*_n_nodes][ind] == static_cast<int>(ind_i+2*_n_nodes))
              A[ind_k+2*_n_nodes][ind] -= 0.5/Tau_hn;

          if ((JA[ind_j][ind] == ind_j) ||
              (JA[ind_j][ind] == ind_k))
            A[ind_j][ind] += 0.25/Tau_hn;

          if ((JA[ind_k][ind] == ind_j) ||
              (JA[ind_k][ind] == ind_k))
            A[ind_k][ind] += 0.25/Tau_hn;

          if ((JA[ind_j+_n_nodes][ind] == static_cast<int>(ind_j+_n_nodes)) ||
              (JA[ind_j+_n_nodes][ind] == static_cast<int>(ind_k+_n_nodes)))
            A[ind_j+_n_nodes][ind] += 0.25/Tau_hn;

          if ((JA[ind_k+_n_nodes][ind] == static_cast<int>(ind_j+_n_nodes)) ||
              (JA[ind_k+_n_nodes][ind] == static_cast<int>(ind_k+_n_nodes)))
            A[ind_k+_n_nodes][ind] += 0.25/Tau_hn;

          if (_dim == 3)
            if ((JA[ind_j+2*_n_nodes][ind] == static_cast<int>(ind_j+2*_n_nodes)) ||
                (JA[ind_j+2*_n_nodes][ind] == static_cast<int>(ind_k+2*_n_nodes)))
              A[ind_j+2*_n_nodes][ind] += 0.25/Tau_hn;

          if (_dim==3)
            if ((JA[ind_k+2*_n_nodes][ind] == static_cast<int>(ind_j+2*_n_nodes)) ||
                (JA[ind_k+2*_n_nodes][ind] == static_cast<int>(ind_k+2*_n_nodes)))
              A[ind_k+2*_n_nodes][ind] += 0.25/Tau_hn;
        }
    }

  // ||\grad J||_2
  for (dof_id_type i=0; i<_dim*_n_nodes; i++)
    nonzero += b[i]*b[i];

  // sort matrix A
  for (dof_id_type i=0; i<_dim*_n_nodes; i++)
    {
      for (int j=columns-1; j>1; j--)
        {
          for (int k=0; k<j; k++)
            {
              if (JA[i][k] > JA[i][k+1])
                {
                  int sch = JA[i][k];
                  JA[i][k] = JA[i][k+1];
                  JA[i][k+1] = sch;
                  double tmp = A[i][k];
                  A[i][k] = A[i][k+1];
                  A[i][k+1] = tmp;
                }
            }
        }
    }

  double eps = std::sqrt(vol)*1e-9;

  // solver for P (unconstrained)
  ia[0] = 0;
  {
    int j = 0;
    for (dof_id_type i=0; i<_dim*_n_nodes; i++)
      {
        u[i] = 0;
        int nz = 0;
        for (int k=0; k<columns; k++)
          {
            if (A[i][k] != 0)
              {
                nz++;
                ja[j] = JA[i][k]+1;
                a[j] = A[i][k];
                j++;
              }
          }
        ia[i+1] = ia[i] + nz;
      }
  }

  dof_id_type m = _dim*_n_nodes;
  int sch = (msglev >= 3) ? 1 : 0;

  solver(m, ia, ja, a, u, b, eps, 100, sch);
  // sol_pcg_pj(m, ia, ja, a, u, b, eps, 100, sch);

  for (dof_id_type i=0; i<_n_nodes; i++)
    {
      //ensure fixed nodes are not moved
      for (unsigned index=0; index<_dim; index++)
        if (mask[i] == 1)
          P[i][index] = 0;
        else
          P[i][index] = u[i+index*_n_nodes];
    }

  // P is determined
  if (msglev >= 4)
    {
      for (dof_id_type i=0; i<_n_nodes; i++)
        {
          for (unsigned j=0; j<_dim; j++)
            if (P[i][j] != 0)
              _logfile << "P[" << i << "][" << j << "]=" << P[i][j];
        }
    }

  // local minimization problem, determination of tau
  if (msglev >= 3)
    _logfile << "dJ=" << std::sqrt(nonzero) << " J0=" << Jpr << std::endl;

  double
    J = 1.e32,
    tau = 0.,
    gVmin = 0.,
    gemax = 0.,
    gqmin = 0.;

  int j = 1;

  while ((Jpr <= J) && (j > -30))
    {
      j = j-1;

      // search through finite # of values for tau
      tau = pow(2., j);
      for (dof_id_type i=0; i<_n_nodes; i++)
        for (unsigned k=0; k<_dim; k++)
          Rpr[i][k] = R[i][k] + tau*P[i][k];

      J = 0;
      gVmin = 1e32;
      gemax = -1e32;
      gqmin = 1e32;
      for (dof_id_type i=0; i<_n_cells; i++)
        {
          if (mcells[i] >= 0)
            {
              int nvert = 0;
              while (cells[i][nvert] >= 0)
                nvert++;

              double lVmin, lqmin;
              double lemax = localP(W, F, Rpr, cells[i], mask, epsilon, w, nvert, H[i], me, vol, 1, lVmin, lqmin, adp, afun, G[i]);

              J += lemax;
              if (gVmin > lVmin)
                gVmin = lVmin;

              if (gemax < lemax)
                gemax = lemax;

              if (gqmin > lqmin)
                gqmin = lqmin;

              // HN correction
              for (dof_id_type ii=0; ii<_n_hanging_edges; ii++)
                {
                  int ind_i = hnodes[ii];
                  int ind_j = edges[2*ii];
                  int ind_k = edges[2*ii+1];
                  for (unsigned jj=0; jj<_dim; jj++)
                    {
                      int g_i = Rpr[ind_i][jj] - 0.5*(Rpr[ind_j][jj]+Rpr[ind_k][jj]);
                      J += g_i*g_i/(2*Tau_hn);
                    }
                }
            }
        }
      if (msglev >= 3)
        _logfile << "tau=" << tau << " J=" << J << std::endl;
    }


  double
    T = 0.,
    gtmin0 = 0.,
    gtmax0 = 0.,
    gqmin0 = 0.;

  if (j != -30)
    {
      Jpr = J;
      for (unsigned i=0; i<_n_nodes; i++)
        for (unsigned k=0; k<_dim; k++)
          Rpr[i][k] = R[i][k] + tau*0.5*P[i][k];

      J = 0;
      gtmin0 = 1e32;
      gtmax0 = -1e32;
      gqmin0 = 1e32;
      for (dof_id_type i=0; i<_n_cells; i++)
        {
          if (mcells[i] >= 0)
            {
              int nvert = 0;
              while (cells[i][nvert] >= 0)
                nvert++;

              double lVmin, lqmin;
              double lemax = localP(W, F, Rpr, cells[i], mask, epsilon, w, nvert, H[i], me, vol, 1, lVmin, lqmin, adp, afun, G[i]);
              J += lemax;

              if (gtmin0 > lVmin)
                gtmin0 = lVmin;

              if (gtmax0 < lemax)
                gtmax0 = lemax;

              if (gqmin0 > lqmin)
                gqmin0 = lqmin;

              // HN correction
              for (dof_id_type ii=0; ii<_n_hanging_edges; ii++)
                {
                  int ind_i = hnodes[ii];
                  int ind_j = edges[2*ii];
                  int ind_k = edges[2*ii+1];

                  for (unsigned jj=0; jj<_dim; jj++)
                    {
                      int g_i = Rpr[ind_i][jj] - 0.5*(Rpr[ind_j][jj] + Rpr[ind_k][jj]);
                      J += g_i*g_i/(2*Tau_hn);
                    }
                }
            }
        }
    }

  if (Jpr > J)
    {
      T = 0.5*tau;
      // Set up return values passed by reference
      Vmin = gtmin0;
      emax = gtmax0;
      qmin = gqmin0;
    }
  else
    {
      T = tau;
      J = Jpr;
      // Set up return values passed by reference
      Vmin = gVmin;
      emax = gemax;
      qmin = gqmin;
    }

  nonzero = 0;
  for (dof_id_type j=0; j<_n_nodes; j++)
    for (unsigned k=0; k<_dim; k++)
      {
        R[j][k] = R[j][k] + T*P[j][k];
        nonzero += T*P[j][k]*T*P[j][k];
      }

  if (msglev >= 2)
    _logfile << "tau=" << T << ", J=" << J << std::endl;

  return std::sqrt(nonzero);
}



// minJ() with sliding Boundary Nodes constraints and no account for HN,
// using Lagrange multiplier formulation: minimize L=J+\sum lam*g;
// only works in 2D
double VariationalMeshSmoother::minJ_BC(Array2D<double> & R,
                                        const std::vector<int> & mask,
                                        const Array2D<int> & cells,
                                        const std::vector<int> & mcells,
                                        double epsilon,
                                        double w,
                                        int me,
                                        const Array3D<double> & H,
                                        double vol,
                                        int msglev,
                                        double & Vmin,
                                        double & emax,
                                        double & qmin,
                                        int adp,
                                        const std::vector<double> & afun,
                                        int NCN)
{
  // new form of matrices, 5 iterations for minL
  double tau = 0., J = 0., T, Jpr, L, lVmin, lqmin, gVmin = 0., gqmin = 0., gVmin0 = 0.,
    gqmin0 = 0., lemax, gemax = 0., gemax0 = 0.;

  // array of sliding BN
  std::vector<int> Bind(NCN);
  std::vector<double> lam(NCN);
  std::vector<double> hm(2*_n_nodes);
  std::vector<double> Plam(NCN);

  // holds constraints = local approximation to the boundary
  std::vector<double> constr(4*NCN);

  Array2D<double> F(2, 6);
  Array3D<double> W(2, 6, 6);
  Array2D<double> Rpr(_n_nodes, 2);
  Array2D<double> P(_n_nodes, 2);

  std::vector<double> b(2*_n_nodes);

  Array2D<double> G(_n_cells, 6);

  // assembler of constraints
  const double eps = std::sqrt(vol)*1e-9;

  for (int i=0; i<4*NCN; i++)
    constr[i] = 1./eps;

  {
    int I = 0;
    for (dof_id_type i=0; i<_n_nodes; i++)
      if (mask[i] == 2)
        {
          Bind[I] = i;
          I++;
        }
  }

  for (int I=0; I<NCN; I++)
    {
      // The boundary connectivity loop sets up the j and k indices
      // but I am not sure about the logic of this: j and k are overwritten
      // at every iteration of the boundary connectivity loop and only used
      // *after* the boundary connectivity loop -- this seems like a bug but
      // I maintained the original behavior of the algorithm [JWP].
      int
        i = Bind[I],
        j = 0,
        k = 0,
        ind = 0;

      // boundary connectivity
      for (dof_id_type l=0; l<_n_cells; l++)
        {
          int nvert = 0;

          while (cells[l][nvert] >= 0)
            nvert++;

          switch (nvert)
            {
            case 3: // tri
              if (i == cells[l][0])
                {
                  if (mask[cells[l][1]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][1];
                      else
                        k = cells[l][1];
                      ind++;
                    }

                  if (mask[cells[l][2]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][2];
                      else
                        k = cells[l][2];
                      ind++;
                    }
                }

              if (i == cells[l][1])
                {
                  if (mask[cells[l][0]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][0];
                      else
                        k = cells[l][0];
                      ind++;
                    }

                  if (mask[cells[l][2]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][2];
                      else
                        k = cells[l][2];
                      ind++;
                    }
                }

              if (i == cells[l][2])
                {
                  if (mask[cells[l][1]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][1];
                      else
                        k = cells[l][1];
                      ind++;
                    }

                  if (mask[cells[l][0]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][0];
                      else
                        k = cells[l][0];
                      ind++;
                    }
                }
              break;

            case 4: // quad
              if ((i == cells[l][0]) ||
                  (i == cells[l][3]))
                {
                  if (mask[cells[l][1]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][1];
                      else
                        k = cells[l][1];
                      ind++;
                    }

                  if (mask[cells[l][2]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][2];
                      else
                        k = cells[l][2];
                      ind++;
                    }
                }

              if ((i == cells[l][1]) ||
                  (i == cells[l][2]))
                {
                  if (mask[cells[l][0]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][0];
                      else
                        k = cells[l][0];
                      ind++;
                    }

                  if (mask[cells[l][3]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][3];
                      else
                        k = cells[l][3];
                      ind++;
                    }
                }
              break;

            case 6: //quad tri
              if (i == cells[l][0])
                {
                  if (mask[cells[l][1]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][5];
                      else
                        k = cells[l][5];
                      ind++;
                    }

                  if (mask[cells[l][2]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][4];
                      else
                        k = cells[l][4];
                      ind++;
                    }
                }

              if (i == cells[l][1])
                {
                  if (mask[cells[l][0]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][5];
                      else
                        k = cells[l][5];
                      ind++;
                    }

                  if (mask[cells[l][2]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][3];
                      else
                        k = cells[l][3];
                      ind++;
                    }
                }

              if (i == cells[l][2])
                {
                  if (mask[cells[l][1]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][3];
                      else
                        k = cells[l][3];
                      ind++;
                    }

                  if (mask[cells[l][0]] > 0)
                    {
                      if (ind == 0)
                        j = cells[l][4];
                      else
                        k = cells[l][4];
                      ind++;
                    }
                }

              if (i == cells[l][3])
                {
                  j = cells[l][1];
                  k = cells[l][2];
                }

              if (i == cells[l][4])
                {
                  j = cells[l][2];
                  k = cells[l][0];
                }

              if (i == cells[l][5])
                {
                  j = cells[l][0];
                  k = cells[l][1];
                }
              break;

            default:
              libmesh_error_msg("Unrecognized nvert = " << nvert);
            }
        } // end boundary connectivity

      // lines
      double del1 = R[j][0] - R[i][0];
      double del2 = R[i][0] - R[k][0];

      if ((std::abs(del1) < eps) &&
          (std::abs(del2) < eps))
        {
          constr[I*4] = 1;
          constr[I*4+1] = 0;
          constr[I*4+2] = R[i][0];
          constr[I*4+3] = R[i][1];
        }
      else
        {
          del1 = R[j][1] - R[i][1];
          del2 = R[i][1] - R[k][1];
          if ((std::abs(del1) < eps) &&
              (std::abs(del2) < eps))
            {
              constr[I*4] = 0;
              constr[I*4+1] = 1;
              constr[I*4+2] = R[i][0];
              constr[I*4+3] = R[i][1];
            }
          else
            {
              J =
                (R[i][0]-R[j][0])*(R[k][1]-R[j][1]) -
                (R[k][0]-R[j][0])*(R[i][1]-R[j][1]);

              if (std::abs(J) < eps)
                {
                  constr[I*4] = 1./(R[k][0]-R[j][0]);
                  constr[I*4+1] = -1./(R[k][1]-R[j][1]);
                  constr[I*4+2] = R[i][0];
                  constr[I*4+3] = R[i][1];
                }
              else
                {
                  // circle
                  double x0 = ((R[k][1]-R[j][1])*(R[i][0]*R[i][0]-R[j][0]*R[j][0]+R[i][1]*R[i][1]-R[j][1]*R[j][1]) +
                               (R[j][1]-R[i][1])*(R[k][0]*R[k][0]-R[j][0]*R[j][0]+R[k][1]*R[k][1]-R[j][1]*R[j][1]))/(2*J);
                  double y0 = ((R[j][0]-R[k][0])*(R[i][0]*R[i][0]-R[j][0]*R[j][0]+R[i][1]*R[i][1]-R[j][1]*R[j][1]) +
                               (R[i][0]-R[j][0])*(R[k][0]*R[k][0]-R[j][0]*R[j][0]+R[k][1]*R[k][1]-R[j][1]*R[j][1]))/(2*J);
                  constr[I*4] = x0;
                  constr[I*4+1] = y0;
                  constr[I*4+2] = (R[i][0]-x0)*(R[i][0]-x0) + (R[i][1]-y0)*(R[i][1]-y0);
                }
            }
        }
    }

  // for (int i=0; i<NCN; i++){
  // for (int j=0; j<4; j++) fprintf(stdout," %e ",constr[4*i+j]);
  // fprintf(stdout, "constr %d node %d \n",i,Bind[i]);
  // }

  // initial values
  for (int i=0; i<NCN; i++)
    lam[i] = 0;

  // Eventual return value
  double nonzero = 0.;

  // Temporary result variable
  double qq = 0.;

  for (int nz=0; nz<5; nz++)
    {
      // find H and -grad J
      nonzero = 0.;
      Jpr = 0;
      for (dof_id_type i=0; i<2*_n_nodes; i++)
        {
          b[i] = 0;
          hm[i] = 0;
        }

      for (dof_id_type i=0; i<_n_cells; i++)
        {
          int nvert = 0;

          while (cells[i][nvert] >= 0)
            nvert++;

          for (int j=0; j<nvert; j++)
            {
              G[i][j] = 0;
              if (adp < 0)
                for (int k=0; k<std::abs(adp); k++)
                  G[i][j] += afun[i*(-adp) + k];
            }

          for (int index=0; index<2; index++)
            for (int k=0; k<nvert; k++)
              {
                F[index][k] = 0;
                for (int j=0; j<nvert; j++)
                  W[index][k][j] = 0;
              }

          if (mcells[i] >= 0)
            Jpr += localP(W, F, R, cells[i], mask, epsilon, w, nvert, H[i],
                          me, vol, 0, lVmin, lqmin, adp, afun, G[i]);

          else
            {
              for (unsigned index=0; index<2; index++)
                for (int j=0; j<nvert; j++)
                  W[index][j][j] = 1;
            }

          for (unsigned index=0; index<2; index++)
            for (int l=0; l<nvert; l++)
              {
                // diagonal Hessian
                hm[cells[i][l] + index*_n_nodes] += W[index][l][l];
                b[cells[i][l] + index*_n_nodes] -= F[index][l];
              }
        }

      // ||grad J||_2
      for (dof_id_type i=0; i<2*_n_nodes; i++)
        nonzero += b[i]*b[i];

      // solve for Plam
      for (int I=0; I<NCN; I++)
        {
          int i = Bind[I];
          double
            Bx = 0.,
            By = 0.,
            g = 0.;

          if (constr[4*I+3] < 0.5/eps)
            {
              Bx = constr[4*I];
              By = constr[4*I+1];
              g = (R[i][0]-constr[4*I+2])*constr[4*I] + (R[i][1]-constr[4*I+3])*constr[4*I+1];
            }
          else
            {
              Bx = 2*(R[i][0] - constr[4*I]);
              By = 2*(R[i][1] - constr[4*I+1]);
              hm[i] += 2*lam[I];
              hm[i+_n_nodes] += 2*lam[I];
              g =
                (R[i][0] - constr[4*I])*(R[i][0]-constr[4*I]) +
                (R[i][1] - constr[4*I+1])*(R[i][1]-constr[4*I+1]) - constr[4*I+2];
            }

          Jpr += lam[I]*g;
          qq = Bx*b[i]/hm[i] + By*b[i+_n_nodes]/hm[i+_n_nodes] - g;
          double a = Bx*Bx/hm[i] + By*By/hm[i+_n_nodes];

          if (a != 0)
            Plam[I] = qq/a;
          else
            _logfile << "error: B^TH-1B is degenerate" << std::endl;

          b[i] -= Plam[I]*Bx;
          b[i+_n_nodes] -= Plam[I]*By;
          Plam[I] -= lam[I];
        }

      // solve for P
      for (dof_id_type i=0; i<_n_nodes; i++)
        {
          P[i][0] = b[i]/hm[i];
          P[i][1] = b[i+_n_nodes]/hm[i+_n_nodes];
        }

      // correct solution
      for (dof_id_type i=0; i<_n_nodes; i++)
        for (unsigned j=0; j<2; j++)
          if ((std::abs(P[i][j]) < eps) || (mask[i] == 1))
            P[i][j] = 0;

      // P is determined
      if (msglev >= 3)
        {
          for (dof_id_type i=0; i<_n_nodes; i++)
            for (unsigned j=0; j<2; j++)
              if (P[i][j] != 0)
                _logfile << "P[" << i << "][" << j << "]=" << P[i][j] << "  ";
        }

      // local minimization problem, determination of tau
      if (msglev >= 3)
        _logfile << "dJ=" << std::sqrt(nonzero) << " L0=" << Jpr << std::endl;

      L = 1.e32;
      int j = 1;

      while ((Jpr <= L) && (j > -30))
        {
          j = j-1;
          tau = pow(2.,j);
          for (dof_id_type i=0; i<_n_nodes; i++)
            for (unsigned k=0; k<2; k++)
              Rpr[i][k] = R[i][k] + tau*P[i][k];

          J = 0;
          gVmin = 1.e32;
          gemax = -1.e32;
          gqmin = 1.e32;
          for (dof_id_type i=0; i<_n_cells; i++)
            if (mcells[i] >= 0)
              {
                int nvert = 0;
                while (cells[i][nvert] >= 0)
                  nvert++;

                lemax = localP(W, F, Rpr, cells[i], mask, epsilon, w, nvert, H[i], me, vol, 1, lVmin,
                               lqmin, adp, afun, G[i]);
                J += lemax;

                if (gVmin > lVmin)
                  gVmin = lVmin;

                if (gemax < lemax)
                  gemax = lemax;

                if (gqmin > lqmin)
                  gqmin = lqmin;
              }

          L = J;

          // constraints contribution
          for (int I=0; I<NCN; I++)
            {
              int i = Bind[I];
              double g = 0.;

              if (constr[4*I+3] < 0.5/eps)
                g = (Rpr[i][0] - constr[4*I+2])*constr[4*I] + (Rpr[i][1]-constr[4*I+3])*constr[4*I+1];

              else
                g = (Rpr[i][0]-constr[4*I])*(Rpr[i][0]-constr[4*I]) +
                  (Rpr[i][1]-constr[4*I+1])*(Rpr[i][1]-constr[4*I+1]) - constr[4*I+2];

              L += (lam[I] + tau*Plam[I])*g;
            }

          // end of constraints
          if (msglev >= 3)
            _logfile << " tau=" << tau << " J=" << J << std::endl;
        } // end while

      if (j == -30)
        T = 0;
      else
        {
          Jpr = L;
          qq = J;
          for (dof_id_type i=0; i<_n_nodes; i++)
            for (unsigned k=0; k<2; k++)
              Rpr[i][k] = R[i][k] + tau*0.5*P[i][k];

          J = 0;
          gVmin0 = 1.e32;
          gemax0 = -1.e32;
          gqmin0 = 1.e32;

          for (dof_id_type i=0; i<_n_cells; i++)
            if (mcells[i] >= 0)
              {
                int nvert = 0;
                while (cells[i][nvert] >= 0)
                  nvert++;

                lemax = localP(W, F, Rpr, cells[i], mask, epsilon, w, nvert, H[i], me, vol, 1, lVmin,
                               lqmin, adp, afun, G[i]);
                J += lemax;

                if (gVmin0 > lVmin)
                  gVmin0 = lVmin;

                if (gemax0 < lemax)
                  gemax0 = lemax;

                if (gqmin0 > lqmin)
                  gqmin0 = lqmin;
              }

          L = J;

          // constraints contribution
          for (int I=0; I<NCN; I++)
            {
              int i = Bind[I];
              double g = 0.;

              if (constr[4*I+3] < 0.5/eps)
                g = (Rpr[i][0]-constr[4*I+2])*constr[4*I] + (Rpr[i][1]-constr[4*I+3])*constr[4*I+1];

              else
                g = (Rpr[i][0]-constr[4*I])*(Rpr[i][0]-constr[4*I]) +
                  (Rpr[i][1]-constr[4*I+1])*(Rpr[i][1]-constr[4*I+1]) - constr[4*I+2];

              L += (lam[I] + tau*0.5*Plam[I])*g;
            }
          // end of constraints
        }

      if (Jpr > L)
        {
          T = 0.5*tau;
          // Set return values by reference
          Vmin = gVmin0;
          emax = gemax0;
          qmin = gqmin0;
        }
      else
        {
          T = tau;
          L = Jpr;
          J = qq;
          // Set return values by reference
          Vmin = gVmin;
          emax = gemax;
          qmin = gqmin;
        }

      for (dof_id_type i=0; i<_n_nodes; i++)
        for (unsigned k=0; k<2; k++)
          R[i][k] += T*P[i][k];

      for (int i=0; i<NCN; i++)
        lam[i] += T*Plam[i];

    } // end Lagrangian iter

  if (msglev >= 2)
    _logfile << "tau=" << T << ", J=" << J << ", L=" << L << std::endl;

  return std::sqrt(nonzero);
}



// composes local matrix W and right side F from all quadrature nodes of one cell
double VariationalMeshSmoother::localP(Array3D<double> & W,
                                       Array2D<double> & F,
                                       Array2D<double> & R,
                                       const std::vector<int> & cell_in,
                                       const std::vector<int> & mask,
                                       double epsilon,
                                       double w,
                                       int nvert,
                                       const Array2D<double> & H,
                                       int me,
                                       double vol,
                                       int f,
                                       double & Vmin,
                                       double & qmin,
                                       int adp,
                                       const std::vector<double> & afun,
                                       std::vector<double> & Gloc)
{
  // K - determines approximation rule for local integral over the cell
  std::vector<double> K(9);

  // f - flag, f=0 for determination of Hessian and gradient of the functional,
  // f=1 for determination of the functional value only;
  // adaptivity is determined on the first step for adp>0 (nodal based)
  if (f == 0)
    {
      if (adp > 0)
        avertex(afun, Gloc, R, cell_in, nvert, adp);
      if (adp == 0)
        {
          for (unsigned i=0; i<_dim; i++)
            Gloc[i] = 1.;
        }
    }

  double
    sigma = 0.,
    fun = 0,
    gqmin = 1e32,
    g = 0,  // Vmin
    lqmin = 0.;

  // cell integration depending on cell type
  if (_dim == 2)
    {
      // 2D
      if (nvert == 3)
        {
          // tri
          sigma = 1.;
          fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me, vol, f, lqmin, adp, Gloc, sigma);
          g += sigma*lqmin;
          if (gqmin > lqmin)
            gqmin = lqmin;
        }
      if (nvert == 4)
        {
          //quad
          for (unsigned i=0; i<2; i++)
            {
              K[0] = i;
              for (unsigned j=0; j<2; j++)
                {
                  K[1] = j;
                  sigma = 0.25;
                  fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                                vol, f, lqmin, adp, Gloc, sigma);
                  g += sigma*lqmin;
                  if (gqmin > lqmin)
                    gqmin = lqmin;
                }
            }
        }
      else
        {
          // quad tri
          for (unsigned i=0; i<3; i++)
            {
              K[0] = i*0.5;
              unsigned k = i/2;
              K[1] = static_cast<double>(k);

              for (unsigned j=0; j<3; j++)
                {
                  K[2] = j*0.5;
                  k = j/2;
                  K[3] = static_cast<double>(k);
                  if (i == j)
                    sigma = 1./12.;
                  else
                    sigma = 1./24.;

                  fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                                vol, f, lqmin, adp, Gloc, sigma);
                  g += sigma*lqmin;
                  if (gqmin > lqmin)
                    gqmin = lqmin;
                }
            }
        }
    }
  if (_dim == 3)
    {
      // 3D
      if (nvert == 4)
        {
          // tetr
          sigma = 1.;
          fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                        vol, f, lqmin, adp, Gloc, sigma);
          g += sigma*lqmin;
          if (gqmin > lqmin)
            gqmin = lqmin;
        }
      if (nvert == 6)
        {
          //prism
          for (unsigned i=0; i<2; i++)
            {
              K[0] = i;
              for (unsigned j=0; j<2; j++)
                {
                  K[1] = j;
                  for (unsigned k=0; k<3; k++)
                    {
                      K[2] = 0.5*static_cast<double>(k);
                      K[3] = static_cast<double>(k % 2);
                      sigma = 1./12.;
                      fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                                    vol, f, lqmin, adp, Gloc, sigma);
                      g += sigma*lqmin;
                      if (gqmin > lqmin)
                        gqmin = lqmin;
                    }
                }
            }
        }
      if (nvert == 8)
        {
          // hex
          for (unsigned i=0; i<2; i++)
            {
              K[0] = i;
              for (unsigned j=0; j<2; j++)
                {
                  K[1] = j;
                  for (unsigned k=0; k<2; k++)
                    {
                      K[2] = k;
                      for (unsigned l=0; l<2; l++)
                        {
                          K[3] = l;
                          for (unsigned m=0; m<2; m++)
                            {
                              K[4] = m;
                              for (unsigned nn=0; nn<2; nn++)
                                {
                                  K[5] = nn;

                                  if ((i==nn) && (j==l) && (k==m))
                                    sigma = 1./27.;

                                  if (((i==nn) && (j==l) && (k!=m)) ||
                                      ((i==nn) && (j!=l) && (k==m)) ||
                                      ((i!=nn) && (j==l) && (k==m)))
                                    sigma = 1./54.;

                                  if (((i==nn) && (j!=l) && (k!=m)) ||
                                      ((i!=nn) && (j!=l) && (k==m)) ||
                                      ((i!=nn) && (j==l) && (k!=m)))
                                    sigma = 1./108.;

                                  if ((i!=nn) && (j!=l) && (k!=m))
                                    sigma=1./216.;

                                  fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                                                vol, f, lqmin, adp, Gloc, sigma);
                                  g += sigma*lqmin;

                                  if (gqmin > lqmin)
                                    gqmin = lqmin;
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          // quad tetr
          for (unsigned i=0; i<4; i++)
            {
              for (unsigned j=0; j<4; j++)
                {
                  for (unsigned k=0; k<4; k++)
                    {
                      switch (i)
                        {
                        case 0:
                          K[0] = 0;
                          K[1] = 0;
                          K[2] = 0;
                          break;

                        case 1:
                          K[0] = 1;
                          K[1] = 0;
                          K[2] = 0;
                          break;

                        case 2:
                          K[0] = 0.5;
                          K[1] = 1;
                          K[2] = 0;
                          break;

                        case 3:
                          K[0] = 0.5;
                          K[1] = 1./3.;
                          K[2] = 1;
                          break;
                        }

                      switch (j)
                        {
                        case 0:
                          K[3] = 0;
                          K[4] = 0;
                          K[5] = 0;
                          break;

                        case 1:
                          K[3] = 1;
                          K[4] = 0;
                          K[5] = 0;
                          break;

                        case 2:
                          K[3] = 0.5;
                          K[4] = 1;
                          K[5] = 0;
                          break;

                        case 3:
                          K[3] = 0.5;
                          K[4] = 1./3.;
                          K[5] = 1;
                          break;

                        }

                      switch (k)
                        {
                        case 0:
                          K[6] = 0;
                          K[7] = 0;
                          K[8] = 0;
                          break;

                        case 1:
                          K[6] = 1;
                          K[7] = 0;
                          K[8] = 0;
                          break;

                        case 2:
                          K[6] = 0.5;
                          K[7] = 1;
                          K[8] = 0;
                          break;

                        case 3:
                          K[6] = 0.5;
                          K[7] = 1./3.;
                          K[8] = 1;
                          break;

                        }

                      if ((i==j) && (j==k))
                        sigma = 1./120.;

                      else if ((i==j) || (j==k) || (i==k))
                        sigma = 1./360.;

                      else
                        sigma = 1./720.;

                      fun += vertex(W, F, R, cell_in, epsilon, w, nvert, K, H, me,
                                    vol, f, lqmin, adp, Gloc, sigma);

                      g += sigma*lqmin;
                      if (gqmin > lqmin)
                        gqmin = lqmin;
                    }
                }
            }
        }
    }

  // fixed nodes correction
  for (int ii=0; ii<nvert; ii++)
    {
      if (mask[cell_in[ii]] == 1)
        {
          for (unsigned kk=0; kk<_dim; kk++)
            {
              for (int jj=0; jj<nvert; jj++)
                {
                  W[kk][ii][jj] = 0;
                  W[kk][jj][ii] = 0;
                }

              W[kk][ii][ii] = 1;
              F[kk][ii] = 0;
            }
        }
    }
  // end of fixed nodes correction

  // Set up return value references
  Vmin = g;
  qmin = gqmin/vol;

  return fun;
}



// avertex - assembly of adaptivity metric on a cell
double VariationalMeshSmoother::avertex(const std::vector<double> & afun,
                                        std::vector<double> & G,
                                        const Array2D<double> & R,
                                        const std::vector<int> & cell_in,
                                        int nvert,
                                        int adp)
{
  std::vector<double> a1(3), a2(3), a3(3), qu(3), K(8);
  Array2D<double> Q(3, nvert);

  for (int i=0; i<8; i++)
    K[i] = 0.5;  // cell center

  basisA(Q, nvert, K, Q, 1);

  for (unsigned i=0; i<_dim; i++)
    for (int j=0; j<nvert; j++)
      {
        a1[i] += Q[i][j]*R[cell_in[j]][0];
        a2[i] += Q[i][j]*R[cell_in[j]][1];
        if (_dim == 3)
          a3[i] += Q[i][j]*R[cell_in[j]][2];

        qu[i] += Q[i][j]*afun[cell_in[j]];
      }

  double det = 0.;

  if (_dim == 2)
    det = jac2(a1[0], a1[1],
               a2[0], a2[1]);
  else
    det = jac3(a1[0], a1[1], a1[2],
               a2[0], a2[1], a2[2],
               a3[0], a3[1], a3[2]);

  // Return val
  double g = 0.;

  if (det != 0)
    {
      if (_dim == 2)
        {
          double df0 = jac2(qu[0], qu[1], a2[0], a2[1])/det;
          double df1 = jac2(a1[0], a1[1], qu[0], qu[1])/det;
          g = (1. + df0*df0 + df1*df1);

          if (adp == 2)
            {
              // directional adaptation G=diag(g_i)
              G[0] = 1. + df0*df0;
              G[1] = 1. + df1*df1;
            }
          else
            {
              for (unsigned i=0; i<_dim; i++)
                G[i] = g;  // simple adaptation G=gI
            }
        }
      else
        {
          double df0 = (qu[0]*(a2[1]*a3[2]-a2[2]*a3[1]) + qu[1]*(a2[2]*a3[0]-a2[0]*a3[2]) + qu[2]*(a2[0]*a3[1]-a2[1]*a3[0]))/det;
          double df1 = (qu[0]*(a3[1]*a1[2]-a3[2]*a1[1]) + qu[1]*(a3[2]*a1[0]-a3[0]*a1[2]) + qu[2]*(a3[0]*a1[1]-a3[1]*a1[0]))/det;
          double df2 = (qu[0]*(a1[1]*a2[2]-a1[2]*a2[1]) + qu[1]*(a1[2]*a2[0]-a1[0]*a2[2]) + qu[2]*(a1[0]*a2[1]-a1[1]*a2[0]))/det;
          g = 1. + df0*df0 + df1*df1 + df2*df2;
          if (adp == 2)
            {
              // directional adaptation G=diag(g_i)
              G[0] = 1. + df0*df0;
              G[1] = 1. + df1*df1;
              G[2] = 1. + df2*df2;
            }
          else
            {
              for (unsigned i=0; i<_dim; i++)
                G[i] = g;  // simple adaptation G=gI
            }
        }
    }
  else
    {
      g = 1.;
      for (unsigned i=0; i<_dim; i++)
        G[i] = g;
    }

  return g;
}



// Computes local matrics W and local rhs F on one basis
double VariationalMeshSmoother::vertex(Array3D<double> & W,
                                       Array2D<double> & F,
                                       const Array2D<double> & R,
                                       const std::vector<int> & cell_in,
                                       double epsilon,
                                       double w,
                                       int nvert,
                                       const std::vector<double> & K,
                                       const Array2D<double> & H,
                                       int me,
                                       double vol,
                                       int f,
                                       double & qmin,
                                       int adp,
                                       const std::vector<double> & g,
                                       double sigma)
{
  // hessian, function, gradient
  Array2D<double> Q(3, nvert);
  basisA(Q, nvert, K, H, me);

  std::vector<double> a1(3), a2(3), a3(3);
  for (unsigned i=0; i<_dim; i++)
    for (int j=0; j<nvert; j++)
      {
        a1[i] += Q[i][j]*R[cell_in[j]][0];
        a2[i] += Q[i][j]*R[cell_in[j]][1];
        if (_dim == 3)
          a3[i] += Q[i][j]*R[cell_in[j]][2];
      }

  // account for adaptation
  double G = 0.;
  if (adp < 2)
    G = g[0];
  else
    {
      G = 1.;
      for (unsigned i=0; i<_dim; i++)
        {
          a1[i] *= std::sqrt(g[0]);
          a2[i] *= std::sqrt(g[1]);
          if (_dim == 3)
            a3[i] *= std::sqrt(g[2]);
        }
    }

  double
    det = 0.,
    tr = 0.,
    phit = 0.;

  std::vector<double> av1(3), av2(3), av3(3);

  if (_dim == 2)
    {
      av1[0] =  a2[1];
      av1[1] = -a2[0];
      av2[0] = -a1[1];
      av2[1] =  a1[0];
      det = jac2(a1[0], a1[1], a2[0], a2[1]);
      for (int i=0; i<2; i++)
        tr += 0.5*(a1[i]*a1[i] + a2[i]*a2[i]);

      phit = (1-w)*tr + w*0.5*(vol/G + det*det*G/vol);
    }

  if (_dim == 3)
    {
      av1[0] = a2[1]*a3[2] - a2[2]*a3[1];
      av1[1] = a2[2]*a3[0] - a2[0]*a3[2];
      av1[2] = a2[0]*a3[1] - a2[1]*a3[0];

      av2[0] = a3[1]*a1[2] - a3[2]*a1[1];
      av2[1] = a3[2]*a1[0] - a3[0]*a1[2];
      av2[2] = a3[0]*a1[1] - a3[1]*a1[0];

      av3[0] = a1[1]*a2[2] - a1[2]*a2[1];
      av3[1] = a1[2]*a2[0] - a1[0]*a2[2];
      av3[2] = a1[0]*a2[1] - a1[1]*a2[0];

      det = jac3(a1[0], a1[1], a1[2],
                 a2[0], a2[1], a2[2],
                 a3[0], a3[1], a3[2]);

      for (int i=0; i<3; i++)
        tr += (1./3.)*(a1[i]*a1[i] + a2[i]*a2[i] + a3[i]*a3[i]);

      phit = (1-w)*pow(tr,1.5) + w*0.5*(vol/G + det*det*G/vol);
    }

  double dchi = 0.5 + det*0.5/std::sqrt(epsilon*epsilon + det*det);
  double chi = 0.5*(det + std::sqrt(epsilon*epsilon + det*det));
  double fet = (chi != 0.) ? phit/chi : 1.e32;

  // Set return value reference
  qmin = det*G;

  // gradient and Hessian
  if (f == 0)
    {
      Array3D<double> P(3, 3, 3);
      Array3D<double> d2phi(3, 3, 3);
      Array2D<double> dphi(3, 3);
      Array2D<double> dfe(3, 3);

      if (_dim == 2)
        {
          for (int i=0; i<2; i++)
            {
              dphi[0][i] = (1-w)*a1[i] + w*G*det*av1[i]/vol;
              dphi[1][i] = (1-w)*a2[i] + w*G*det*av2[i]/vol;
            }

          for (int i=0; i<2; i++)
            for (int j=0; j<2; j++)
              {
                d2phi[0][i][j] = w*G*av1[i]*av1[j]/vol;
                d2phi[1][i][j] = w*G*av2[i]*av2[j]/vol;

                if (i == j)
                  for (int k=0; k<2; k++)
                    d2phi[k][i][j] += 1.-w;
              }

          for (int i=0; i<2; i++)
            {
              dfe[0][i] = (dphi[0][i] - fet*dchi*av1[i])/chi;
              dfe[1][i] = (dphi[1][i] - fet*dchi*av2[i])/chi;
            }

          for (int i=0; i<2; i++)
            for (int j=0; j<2; j++)
              {
                P[0][i][j] = (d2phi[0][i][j] - dchi*(av1[i]*dfe[0][j] + dfe[0][i]*av1[j]))/chi;
                P[1][i][j] = (d2phi[1][i][j] - dchi*(av2[i]*dfe[1][j] + dfe[1][i]*av2[j]))/chi;
              }
        }

      if (_dim == 3)
        {
          for (int i=0; i<3; i++)
            {
              dphi[0][i] = (1-w)*std::sqrt(tr)*a1[i] + w*G*det*av1[i]/vol;
              dphi[1][i] = (1-w)*std::sqrt(tr)*a2[i] + w*G*det*av2[i]/vol;
              dphi[2][i] = (1-w)*std::sqrt(tr)*a3[i] + w*G*det*av3[i]/vol;
            }

          for (int i=0; i<3; i++)
            {
              if (tr != 0)
                {
                  d2phi[0][i][i] = (1-w)/(3.*std::sqrt(tr))*a1[i]*a1[i] + w*G*av1[i]*av1[i]/vol;
                  d2phi[1][i][i] = (1-w)/(3.*std::sqrt(tr))*a2[i]*a2[i] + w*G*av2[i]*av2[i]/vol;
                  d2phi[2][i][i] = (1-w)/(3.*std::sqrt(tr))*a3[i]*a3[i] + w*G*av3[i]*av3[i]/vol;
                }
              else
                {
                  for (int k=0; k<3; k++)
                    d2phi[k][i][i] = 0.;
                }

              for (int k=0; k<3; k++)
                d2phi[k][i][i] += (1-w)*std::sqrt(tr);
            }

          const double con = 100.;

          for (int i=0; i<3; i++)
            {
              dfe[0][i] = (dphi[0][i] - fet*dchi*av1[i] + con*epsilon*a1[i])/chi;
              dfe[1][i] = (dphi[1][i] - fet*dchi*av2[i] + con*epsilon*a2[i])/chi;
              dfe[2][i] = (dphi[2][i] - fet*dchi*av3[i] + con*epsilon*a3[i])/chi;
            }

          for (int i=0; i<3; i++)
            {
              P[0][i][i] = (d2phi[0][i][i] - dchi*(av1[i]*dfe[0][i] + dfe[0][i]*av1[i]) + con*epsilon)/chi;
              P[1][i][i] = (d2phi[1][i][i] - dchi*(av2[i]*dfe[1][i] + dfe[1][i]*av2[i]) + con*epsilon)/chi;
              P[2][i][i] = (d2phi[2][i][i] - dchi*(av3[i]*dfe[2][i] + dfe[2][i]*av3[i]) + con*epsilon)/chi;
            }

          for (int k=0; k<3; k++)
            for (int i=0; i<3; i++)
              for (int j=0; j<3; j++)
                if (i != j)
                  P[k][i][j] = 0.;
        }

      /*--------matrix W, right side F---------------------*/
      Array2D<double> gpr(3, nvert);

      for (unsigned i1=0; i1<_dim; i1++)
        {
          for (unsigned i=0; i<_dim; i++)
            for (int j=0; j<nvert; j++)
              for (unsigned k=0; k<_dim; k++)
                gpr[i][j] += P[i1][i][k]*Q[k][j];

          for (int i=0; i<nvert; i++)
            for (int j=0; j<nvert; j++)
              for (unsigned k=0; k<_dim; k++)
                W[i1][i][j] += Q[k][i]*gpr[k][j]*sigma;

          for (int i=0; i<nvert; i++)
            for (int k=0; k<2; k++)
              F[i1][i] += Q[k][i]*dfe[i1][k]*sigma;
        }
    }

  return fet*sigma;
}



// Solve Symmetric Positive Definite (SPD) linear system
// by Conjugate Gradient (CG) method preconditioned by
// Point Jacobi (diagonal scaling)
int VariationalMeshSmoother::solver(int n,
                                    const std::vector<int> & ia,
                                    const std::vector<int> & ja,
                                    const std::vector<double> & a,
                                    std::vector<double> & x,
                                    const std::vector<double> & b,
                                    double eps,
                                    int maxite,
                                    int msglev)
{
  _logfile << "Beginning Solve " << n << std::endl;

  // Check parameters
  int info = pcg_par_check(n, ia, ja, a, eps, maxite, msglev);
  if (info != 0)
    return info;

  // PJ preconditioner construction
  std::vector<double> u(n);
  for (int i=0; i<n; i++)
    u[i] = 1./a[ia[i]];

  // PCG iterations
  std::vector<double> r(n), p(n), z(n);
  info = pcg_ic0(n, ia, ja, a, u, x, b, r, p, z, eps, maxite, msglev);

  // Mat sparse_global;
  // MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,n,n,ia,ja,a,&sparse_global);

  _logfile << "Finished Solve " << n << std::endl;

  return info;
}



// Input parameter check
int VariationalMeshSmoother::pcg_par_check(int n,
                                           const std::vector<int> & ia,
                                           const std::vector<int> & ja,
                                           const std::vector<double> & a,
                                           double eps,
                                           int maxite,
                                           int msglev)
{
  int i, j, jj, k;

  if (n <= 0)
    {
      if (msglev > 0)
        _logfile << "sol_pcg: incorrect input parameter: n =" << n << std::endl;
      return -1;
    }

  if (ia[0] != 0)
    {
      if (msglev > 0)
        _logfile << "sol_pcg: incorrect input parameter: ia(1) =" << ia[0] << std::endl;
      return -2;
    }

  for (i=0; i<n; i++)
    {
      if (ia[i+1] < ia[i])
        {
          if (msglev >= 1)
            _logfile << "sol_pcg: incorrect input parameter: i ="
                     << i+1
                     << "; ia(i) ="
                     << ia[i]
                     << " must be <= a(i+1) ="
                     << ia[i+1]
                     << std::endl;
          return -2;
        }
    }

  for (i=0; i<n; i++)
    {
      if (ja[ia[i]] != (i+1))
        {
          if (msglev >= 1)
            _logfile << "sol_pcg: incorrect input parameter: i ="
                     << i+1
                     << " ; ia(i) ="
                     << ia[i]
                     << " ; absence of diagonal entry; k ="
                     << ia[i]+1
                     << " ; the value ja(k) ="
                     << ja[ia[i]]
                     << " must be equal to i"
                     << std::endl;

          return -3;
        }

      jj = 0;
      for (k=ia[i]; k<ia[i+1]; k++)
        {
          j=ja[k];
          if ((j<(i+1)) || (j>n))
            {
              if (msglev >= 1)
                _logfile << "sol_pcg: incorrect input parameter: i ="
                         << i+1
                         << " ; ia(i) ="
                         << ia[i]
                         << " ; a(i+1) ="
                         << ia[i+1]
                         << " ; k ="
                         << k+1
                         << " ; the value ja(k) ="
                         << ja[k]
                         << " is out of range"
                         << std::endl;

              return -3;
            }
          if (j <= jj)
            {
              if (msglev >= 1)
                _logfile << "sol_pcg: incorrect input parameter: i ="
                         << i+1
                         << " ; ia(i) ="
                         << ia[i]
                         << " ; a(i+1) ="
                         << ia[i+1]
                         << " ; k ="
                         << k+1
                         << " ; the value ja(k) ="
                         << ja[k]
                         << " must be less than ja(k+1) ="
                         << ja[k+1]
                         << std::endl;

              return -3;
            }
          jj = j;
        }
    }

  for (i=0; i<n; i++)
    {
      if (a[ia[i]] <= 0.)
        {
          if (msglev >= 1)
            _logfile << "sol_pcg: incorrect input parameter: i ="
                     << i+1
                     << " ; ia(i) ="
                     << ia[i]
                     << " ; the diagonal entry a(ia(i)) ="
                     << a[ia[i]]
                     << " must be > 0"
                     << std::endl;

          return -4;
        }
    }

  if (eps < 0.)
    {
      if (msglev > 0)
        _logfile << "sol_pcg: incorrect input parameter: eps =" << eps << std::endl;
      return -7;
    }

  if (maxite < 1)
    {
      if (msglev > 0)
        _logfile << "sol_pcg: incorrect input parameter: maxite =" << maxite << std::endl;
      return -8;
    }

  return 0;
}



// Solve the SPD linear system by PCG method
int VariationalMeshSmoother::pcg_ic0(int n,
                                     const std::vector<int> & ia,
                                     const std::vector<int> & ja,
                                     const std::vector<double> & a,
                                     const std::vector<double> & u,
                                     std::vector<double> & x,
                                     const std::vector<double> & b,
                                     std::vector<double> & r,
                                     std::vector<double> & p,
                                     std::vector<double> & z,
                                     double eps,
                                     int maxite,
                                     int msglev)
{
  // Return value
  int i = 0;

  double
    rhr = 0.,
    rhr0 = 0.,
    rhrold = 0.,
    rr0 = 0.;

  for (i=0; i<=maxite; i++)
    {
      if (i == 0)
        p = x;

      // z=Ap
      for (int ii=0; ii<n; ii++)
        z[ii] = 0.;

      for (int ii=0; ii<n; ii++)
        {
          z[ii] += a[ia[ii]]*p[ii];

          for (int j=ia[ii]+1; j<ia[ii+1]; j++)
            {
              z[ii] += a[j]*p[ja[j]-1];
              z[ja[j]-1] += a[j]*p[ii];
            }
        }

      if (i == 0)
        for (int k=0; k<n; k++)
          r[k] = b[k] - z[k]; // r(0) = b - Ax(0)

      if (i > 0)
        {
          double pap = 0.;
          for (int k=0; k<n; k++)
            pap += p[k]*z[k];

          double alpha = rhr/pap;
          for (int k=0; k<n; k++)
            {
              x[k] += p[k]*alpha;
              r[k] -= z[k]*alpha;
            }
          rhrold = rhr;
        }

      // z = H * r
      for (int ii=0; ii<n; ii++)
        z[ii] = r[ii]*u[ii];

      if (i == 0)
        p = z;

      rhr = 0.;
      double rr = 0.;
      for (int k=0; k<n; k++)
        {
          rhr += r[k]*z[k];
          rr += r[k]*r[k];
        }

      if (i == 0)
        {
          rhr0 = rhr;
          rr0 = rr;
        }

      if (msglev > 0)
        _logfile << i
                 << " )  rHr ="
                 << std::sqrt(rhr/rhr0)
                 << " rr ="
                 << std::sqrt(rr/rr0)
                 << std::endl;

      if (rr <= eps*eps*rr0)
        return i;

      if (i >= maxite)
        return i;

      if (i > 0)
        {
          double beta = rhr/rhrold;
          for (int k=0; k<n; k++)
            p[k] = z[k] + p[k]*beta;
        }
    } // end for i

  return i;
}




// Sample mesh file generation
void VariationalMeshSmoother::gener(char grid[],
                                    int n)
{
  const int n1 = 3;

  int
    N = 1,
    ncells = 1;

  for (int i=0; i<n; i++)
    {
      N *= n1;
      ncells *= (n1-1);
    }

  const double x = 1./static_cast<double>(n1-1);

  std::ofstream outfile(grid);

  outfile << n << "\n" << N << "\n" << ncells << "\n0\n" << std::endl;

  for (int i=0; i<N; i++)
    {
      // node coordinates
      int k = i;
      int mask = 0;
      for (int j=0; j<n; j++)
        {
          int ii = k % n1;
          if ((ii == 0) || (ii == n1-1))
            mask = 1;

          outfile << static_cast<double>(ii)*x << " ";
          k /= n1;
        }
      outfile << mask << std::endl;
    }

  for (int i=0; i<ncells; i++)
    {
      // cell connectivity
      int nc = i;
      int ii = nc%(n1-1);
      nc /= (n1-1);
      int jj = nc%(n1-1);
      int kk = nc/(n1-1);

      if (n == 2)
        outfile << ii+n1*jj << " "
                << ii+1+jj*n1 << " "
                << ii+(jj+1)*n1 << " "
                << ii+1+(jj+1)*n1 << " ";

      if (n == 3)
        outfile << ii+n1*jj+n1*n1*kk << " "
                << ii+1+jj*n1+n1*n1*kk << " "
                << ii+(jj+1)*n1+n1*n1*kk << " "
                << ii+1+(jj+1)*n1+n1*n1*kk << " "
                << ii+n1*jj+n1*n1*(kk+1) << " "
                << ii+1+jj*n1+n1*n1*(kk+1) << " "
                << ii+(jj+1)*n1+n1*n1*(kk+1) << " "
                << ii+1+(jj+1)*n1+n1*n1*(kk+1) << " ";

      outfile << "-1 -1 0 \n";
    }
}



// Metric Generation
void VariationalMeshSmoother::metr_data_gen(std::string grid,
                                            std::string metr,
                                            int me)
{
  double det, g1, g2, g3, det_o, g1_o, g2_o, g3_o, eps=1e-3;

  std::vector<double> K(9);
  Array2D<double> Q(3, 3*_dim + _dim%2);

  // Use _mesh to determine N and ncells
  this->_n_nodes = _mesh.n_nodes();
  this->_n_cells = _mesh.n_active_elem();

  std::vector<int>
    mask(_n_nodes),
    mcells(_n_cells);

  Array2D<int> cells(_n_cells, 3*_dim + _dim%2);
  Array2D<double> R(_n_nodes,_dim);

  readgr(R, mask, cells, mcells, mcells, mcells);

  // genetrate metric file
  std::ofstream metric_file(metr.c_str());
  if (!metric_file.good())
    libmesh_error_msg("Error opening metric output file.");

  // Use scientific notation with 6 digits
  metric_file << std::scientific << std::setprecision(6);

  int Ncells = 0;
  det_o = 1.;
  g1_o = 1.;
  g2_o = 1.;
  g3_o = 1.;
  for (unsigned i=0; i<_n_cells; i++)
    if (mcells[i] >= 0)
      {
        int nvert=0;
        while (cells[i][nvert] >= 0)
          nvert++;

        if (_dim == 2)
          {
            // 2D - tri and quad
            if (nvert == 3)
              {
                // tri
                basisA(Q, 3, K, Q, 1);

                Point a1, a2;
                for (int k=0; k<2; k++)
                  for (int l=0; l<3; l++)
                    {
                      a1(k) += Q[k][l]*R[cells[i][l]][0];
                      a2(k) += Q[k][l]*R[cells[i][l]][1];
                    }

                det = jac2(a1(0), a1(1), a2(0), a2(1));
                g1 = std::sqrt(a1(0)*a1(0) + a2(0)*a2(0));
                g2 = std::sqrt(a1(1)*a1(1) + a2(1)*a2(1));

                // need to keep data from previous cell
                if ((std::abs(det) < eps*eps*det_o) || (det < 0))
                  det = det_o;

                if ((std::abs(g1) < eps*g1_o) || (g1<0))
                  g1 = g1_o;

                if ((std::abs(g2) < eps*g2_o) || (g2<0))
                  g2 = g2_o;

                // write to file
                if (me == 2)
                  metric_file << 1./std::sqrt(det)
                              << " 0.000000e+00 \n0.000000e+00 "
                              << 1./std::sqrt(det)
                              << std::endl;

                if (me == 3)
                  metric_file << 1./g1
                              << " 0.000000e+00 \n0.000000e+00 "
                              << 1./g2
                              << std::endl;

                det_o = det;
                g1_o = g1;
                g2_o = g2;
                Ncells++;
              }

            if (nvert == 4)
              {
                // quad

                // Set up the node edge neighbor indices
                const unsigned node_indices[4] = {0, 1, 3, 2};
                const unsigned first_neighbor_indices[4] = {1, 3, 2, 0};
                const unsigned second_neighbor_indices[4] = {2, 0, 1, 3};

                // Loop over each node, compute some quantities associated
                // with its edge neighbors, and write them to file.
                for (unsigned ni=0; ni<4; ++ni)
                  {
                    unsigned
                      node_index = node_indices[ni],
                      first_neighbor_index = first_neighbor_indices[ni],
                      second_neighbor_index = second_neighbor_indices[ni];

                    Real
                      node_x = R[cells[i][node_index]][0],
                      node_y = R[cells[i][node_index]][1],
                      first_neighbor_x = R[cells[i][first_neighbor_index]][0],
                      first_neighbor_y = R[cells[i][first_neighbor_index]][1],
                      second_neighbor_x = R[cells[i][second_neighbor_index]][0],
                      second_neighbor_y = R[cells[i][second_neighbor_index]][1];


                    // "dx"
                    Point a1(first_neighbor_x - node_x,
                             second_neighbor_x - node_x);

                    // "dy"
                    Point a2(first_neighbor_y - node_y,
                             second_neighbor_y - node_y);

                    det = jac2(a1(0), a1(1), a2(0), a2(1));
                    g1 = std::sqrt(a1(0)*a1(0) + a2(0)*a2(0));
                    g2 = std::sqrt(a1(1)*a1(1) + a2(1)*a2(1));

                    // need to keep data from previous cell
                    if ((std::abs(det) < eps*eps*det_o) || (det < 0))
                      det = det_o;

                    if ((std::abs(g1) < eps*g1_o) || (g1 < 0))
                      g1 = g1_o;

                    if ((std::abs(g2) < eps*g2_o) || (g2 < 0))
                      g2 = g2_o;

                    // write to file
                    if (me == 2)
                      metric_file << 1./std::sqrt(det) << " "
                                  << 0.5/std::sqrt(det) << " \n0.000000e+00 "
                                  << 0.5*std::sqrt(3./det)
                                  << std::endl;

                    if (me == 3)
                      metric_file << 1./g1 << " "
                                  << 0.5/g2 << "\n0.000000e+00 "
                                  << 0.5*std::sqrt(3.)/g2
                                  << std::endl;

                    det_o = det;
                    g1_o = g1;
                    g2_o = g2;
                    Ncells++;
                  }
              } // end QUAD case
          } // end _dim==2

        if (_dim == 3)
          {
            // 3D - tetr and hex

            if (nvert == 4)
              {
                // tetr
                basisA(Q, 4, K, Q, 1);

                Point a1, a2, a3;

                for (int k=0; k<3; k++)
                  for (int l=0; l<4; l++)
                    {
                      a1(k) += Q[k][l]*R[cells[i][l]][0];
                      a2(k) += Q[k][l]*R[cells[i][l]][1];
                      a3(k) += Q[k][l]*R[cells[i][l]][2];
                    }

                det = jac3(a1(0), a1(1), a1(2),
                           a2(0), a2(1), a2(2),
                           a3(0), a3(1), a3(2));
                g1 = std::sqrt(a1(0)*a1(0) + a2(0)*a2(0) + a3(0)*a3(0));
                g2 = std::sqrt(a1(1)*a1(1) + a2(1)*a2(1) + a3(1)*a3(1));
                g3 = std::sqrt(a1(2)*a1(2) + a2(2)*a2(2) + a3(2)*a3(2));

                // need to keep data from previous cell
                if ((std::abs(det) < eps*eps*eps*det_o) || (det < 0))
                  det = det_o;

                if ((std::abs(g1) < eps*g1_o) || (g1 < 0))
                  g1 = g1_o;

                if ((std::abs(g2) < eps*g2_o) || (g2 < 0))
                  g2 = g2_o;

                if ((std::abs(g3) < eps*g3_o) || (g3 < 0))
                  g3 = g3_o;

                // write to file
                if (me == 2)
                  metric_file << 1./pow(det, 1./3.)
                              << " 0.000000e+00  0.000000e+00\n0.000000e+00 "
                              << 1./pow(det, 1./3.)
                              << " 0.000000e+00\n0.000000e+00 0.000000e+00 "
                              << 1./pow(det, 1./3.)
                              << std::endl;

                if (me == 3)
                  metric_file << 1./g1
                              << " 0.000000e+00  0.000000e+00\n0.000000e+00 "
                              << 1./g2
                              << " 0.000000e+00\n0.000000e+00 0.000000e+00 "
                              << 1./g3
                              << std::endl;

                det_o = det;
                g1_o = g1;
                g2_o = g2;
                g3_o = g3;
                Ncells++;
              }

            if (nvert == 8)
              {
                // hex

                // Set up the node edge neighbor indices
                const unsigned node_indices[8] = {0, 1, 3, 2, 4, 5, 7, 6};
                const unsigned first_neighbor_indices[8] = {1, 3, 2, 0, 6, 4, 5, 7};
                const unsigned second_neighbor_indices[8] = {2, 0, 1, 3, 5, 7, 6, 4};
                const unsigned third_neighbor_indices[8] = {4, 5, 7, 6, 0, 1, 3, 2};

                // Loop over each node, compute some quantities associated
                // with its edge neighbors, and write them to file.
                for (unsigned ni=0; ni<8; ++ni)
                  {
                    unsigned
                      node_index = node_indices[ni],
                      first_neighbor_index = first_neighbor_indices[ni],
                      second_neighbor_index = second_neighbor_indices[ni],
                      third_neighbor_index = third_neighbor_indices[ni];

                    Real
                      node_x = R[cells[i][node_index]][0],
                      node_y = R[cells[i][node_index]][1],
                      node_z = R[cells[i][node_index]][2],
                      first_neighbor_x = R[cells[i][first_neighbor_index]][0],
                      first_neighbor_y = R[cells[i][first_neighbor_index]][1],
                      first_neighbor_z = R[cells[i][first_neighbor_index]][2],
                      second_neighbor_x = R[cells[i][second_neighbor_index]][0],
                      second_neighbor_y = R[cells[i][second_neighbor_index]][1],
                      second_neighbor_z = R[cells[i][second_neighbor_index]][2],
                      third_neighbor_x = R[cells[i][third_neighbor_index]][0],
                      third_neighbor_y = R[cells[i][third_neighbor_index]][1],
                      third_neighbor_z = R[cells[i][third_neighbor_index]][2];

                    // "dx"
                    Point a1(first_neighbor_x - node_x,
                             second_neighbor_x - node_x,
                             third_neighbor_x - node_x);

                    // "dy"
                    Point a2(first_neighbor_y - node_y,
                             second_neighbor_y - node_y,
                             third_neighbor_y - node_y);

                    // "dz"
                    Point a3(first_neighbor_z - node_z,
                             second_neighbor_z - node_z,
                             third_neighbor_z - node_z);

                    det = jac3(a1(0), a1(1), a1(2),
                               a2(0), a2(1), a2(2),
                               a3(0), a3(1), a3(2));
                    g1 = std::sqrt(a1(0)*a1(0) + a2(0)*a2(0) + a3(0)*a3(0));
                    g2 = std::sqrt(a1(1)*a1(1) + a2(1)*a2(1) + a3(1)*a3(1));
                    g3 = std::sqrt(a1(2)*a1(2) + a2(2)*a2(2) + a3(2)*a3(2));

                    // need to keep data from previous cell
                    if ((std::abs(det) < eps*eps*eps*det_o) || (det < 0))
                      det = det_o;

                    if ((std::abs(g1) < eps*g1_o) || (g1 < 0))
                      g1 = g1_o;

                    if ((std::abs(g2) < eps*g2_o) || (g2 < 0))
                      g2 = g2_o;

                    if ((std::abs(g3) < eps*g3_o) || (g3 < 0))
                      g3 = g3_o;

                    // write to file
                    if (me == 2)
                      metric_file << 1./pow(det, 1./3.) << " "
                                  << 0.5/pow(det, 1./3.) << " "
                                  << 0.5/pow(det, 1./3.) << "\n0.000000e+00 "
                                  << 0.5*std::sqrt(3.)/pow(det, 1./3.) << " "
                                  << 0.5/(std::sqrt(3.)*pow(det, 1./3.)) << "\n0.000000e+00 0.000000e+00 "
                                  << std::sqrt(2/3.)/pow(det, 1./3.)
                                  << std::endl;

                    if (me == 3)
                      metric_file << 1./g1 << " "
                                  << 0.5/g2 << " "
                                  << 0.5/g3 << "\n0.000000e+00 "
                                  << 0.5*std::sqrt(3.)/g2 << " "
                                  << 0.5/(std::sqrt(3.)*g3) << "\n0.000000e+00 0.000000e+00 "
                                  << std::sqrt(2./3.)/g3
                                  << std::endl;

                    det_o = det;
                    g1_o = g1;
                    g2_o = g2;
                    g3_o = g3;
                    Ncells++;
                  } // end for ni
              } // end hex
          } // end dim==3
      }

  // Done with the metric file
  metric_file.close();

  std::ofstream grid_file(grid.c_str());
  if (!grid_file.good())
    libmesh_error_msg("Error opening file: " << grid);

  grid_file << _dim << "\n" << _n_nodes << "\n" << Ncells << "\n0" << std::endl;

  // Use scientific notation with 6 digits
  grid_file << std::scientific << std::setprecision(6);

  for (dof_id_type i=0; i<_n_nodes; i++)
    {
      // node coordinates
      for (unsigned j=0; j<_dim; j++)
        grid_file << R[i][j] << " ";

      grid_file << mask[i] << std::endl;
    }

  for (dof_id_type i=0; i<_n_cells; i++)
    if (mcells[i] >= 0)
      {
        // cell connectivity
        int nvert = 0;
        while (cells[i][nvert] >= 0)
          nvert++;

        if ((nvert == 3) || ((_dim == 3) && (nvert == 4)))
          {
            // tri & tetr
            for (int j=0; j<nvert; j++)
              grid_file << cells[i][j] << " ";

            for (unsigned j=nvert; j<3*_dim + _dim%2; j++)
              grid_file << "-1 ";

            grid_file << mcells[i] << std::endl;
          }

        if ((_dim == 2) && (nvert == 4))
          {
            // quad
            grid_file << cells[i][0] << " " << cells[i][1] << " " << cells[i][2] << " -1 -1 -1 0" << std::endl;
            grid_file << cells[i][1] << " " << cells[i][3] << " " << cells[i][0] << " -1 -1 -1 0" << std::endl;
            grid_file << cells[i][3] << " " << cells[i][2] << " " << cells[i][1] << " -1 -1 -1 0" << std::endl;
            grid_file << cells[i][2] << " " << cells[i][0] << " " << cells[i][3] << " -1 -1 -1 0" << std::endl;
          }

        if (nvert == 8)
          {
            // hex
            grid_file << cells[i][0] << " " << cells[i][1] << " " << cells[i][2] << " " << cells[i][4] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][1] << " " << cells[i][3] << " " << cells[i][0] << " " << cells[i][5] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][3] << " " << cells[i][2] << " " << cells[i][1] << " " << cells[i][7] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][2] << " " << cells[i][0] << " " << cells[i][3] << " " << cells[i][6] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][4] << " " << cells[i][6] << " " << cells[i][5] << " " << cells[i][0] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][5] << " " << cells[i][4] << " " << cells[i][7] << " " << cells[i][1] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][7] << " " << cells[i][5] << " " << cells[i][6] << " " << cells[i][3] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
            grid_file << cells[i][6] << " " << cells[i][7] << " " << cells[i][4] << " " << cells[i][2] << " -1 -1 -1 -1 -1 -1 0" << std::endl;
          }
      }
}

} // namespace libMesh

#endif // LIBMESH_ENABLE_VSMOOTHER
