// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes
#include <vector>

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"
#include "libmesh/wrapped_function.h"

namespace libMesh
{

// ------------------------------------------------------------
// Helper class definitions

  /**
   * This class implements projecting a vector from
   * an old mesh to the newly refined mesh.  This
   * may be executed in parallel on multiple threads.
   */
  class ProjectVector
  {
  private:
    const System                &system;
    const NumericVector<Number> &old_vector;
    NumericVector<Number>       &new_vector;

  public:
    ProjectVector (const System &system_in,
		   const NumericVector<Number> &old_v_in,
		   NumericVector<Number> &new_v_in) :
    system(system_in),
    old_vector(old_v_in),
    new_vector(new_v_in)
    {}

    void operator()(const ConstElemRange &range) const;
  };


  /**
   * This class builds the send_list of old dof indices
   * whose coefficients are needed to perform a projection.
   * This may be executed in parallel on multiple threads.
   * The end result is a \p send_list vector which is
   * unsorted and may contain duplicate elements.
   * The \p unique() method can be used to sort and
   * create a unique list.
   */
  class BuildProjectionList
  {
  private:
    const System              &system;

  public:
    BuildProjectionList (const System &system_in) :
      system(system_in),
      send_list()
    {}

    BuildProjectionList (BuildProjectionList &other, Threads::split) :
      system(other.system),
      send_list()
    {}

    void unique();
    void operator()(const ConstElemRange &range);
    void join (const BuildProjectionList &other);
    std::vector<dof_id_type> send_list;
  };



  /**
   * This class implements projecting an arbitrary
   * function to the current mesh.  This
   * may be exectued in parallel on multiple threads.
   */
  class ProjectSolution
  {
  private:
    const System                &system;

    AutoPtr<FunctionBase<Number> > f;
    AutoPtr<FunctionBase<Gradient> > g;
    const Parameters &parameters;
    NumericVector<Number>       &new_vector;

  public:
    ProjectSolution (const System &system_in,
		     FunctionBase<Number>* f_in,
		     FunctionBase<Gradient>* g_in,
		     const Parameters &parameters_in,
		     NumericVector<Number> &new_v_in) :
    system(system_in),
    f(f_in ? f_in->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(g_in ? g_in->clone() : AutoPtr<FunctionBase<Gradient> >(NULL)),
    parameters(parameters_in),
    new_vector(new_v_in)
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

    ProjectSolution (const ProjectSolution &in) :
    system(in.system),
    f(in.f.get() ? in.f->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(in.g.get() ? in.g->clone() : AutoPtr<FunctionBase<Gradient> >(NULL)),
    parameters(in.parameters),
    new_vector(in.new_vector)
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

    void operator()(const ConstElemRange &range) const;
  };


  /**
   * This class implements projecting an arbitrary
   * function to the current mesh.  This
   * may be exectued in parallel on multiple threads.
   */
  class ProjectFEMSolution
  {
  private:
    const System                &system;

    AutoPtr<FEMFunctionBase<Number> > f;
    AutoPtr<FEMFunctionBase<Gradient> > g;
    NumericVector<Number>       &new_vector;

  public:
    ProjectFEMSolution (const System &system_in,
			FEMFunctionBase<Number>* f_in,
			FEMFunctionBase<Gradient>* g_in,
			NumericVector<Number> &new_v_in) :
    system(system_in),
    f(f_in ? f_in->clone() : AutoPtr<FEMFunctionBase<Number> >(NULL)),
    g(g_in ? g_in->clone() : AutoPtr<FEMFunctionBase<Gradient> >(NULL)),
    new_vector(new_v_in)
    {
      libmesh_assert(f.get());
    }

    ProjectFEMSolution (const ProjectFEMSolution &in) :
    system(in.system),
    f(in.f.get() ? in.f->clone() : AutoPtr<FEMFunctionBase<Number> >(NULL)),
    g(in.g.get() ? in.g->clone() : AutoPtr<FEMFunctionBase<Gradient> >(NULL)),
    new_vector(in.new_vector)
    {
      libmesh_assert(f.get());
    }

    void operator()(const ConstElemRange &range) const;
  };


  /**
   * This class implements projecting an arbitrary
   * boundary function to the current mesh.  This
   * may be exectued in parallel on multiple threads.
   */
  class BoundaryProjectSolution
  {
  private:
    const std::set<boundary_id_type> &b;
    const std::vector<unsigned int>  &variables;
    const System                     &system;
    AutoPtr<FunctionBase<Number> >    f;
    AutoPtr<FunctionBase<Gradient> >  g;
    const Parameters                 &parameters;
    NumericVector<Number>            &new_vector;

  public:
    BoundaryProjectSolution (const std::set<boundary_id_type> &b_in,
                             const std::vector<unsigned int> &variables_in,
                             const System &system_in,
		             FunctionBase<Number>* f_in,
		             FunctionBase<Gradient>* g_in,
		             const Parameters &parameters_in,
		             NumericVector<Number> &new_v_in) :
    b(b_in),
    variables(variables_in),
    system(system_in),
    f(f_in ? f_in->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(g_in ? g_in->clone() : AutoPtr<FunctionBase<Gradient> >(NULL)),
    parameters(parameters_in),
    new_vector(new_v_in)
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

    BoundaryProjectSolution (const BoundaryProjectSolution &in) :
    b(in.b),
    variables(in.variables),
    system(in.system),
    f(in.f.get() ? in.f->clone() : AutoPtr<FunctionBase<Number> >(NULL)),
    g(in.g.get() ? in.g->clone() : AutoPtr<FunctionBase<Gradient> >(NULL)),
    parameters(in.parameters),
    new_vector(in.new_vector)
    {
      libmesh_assert(f.get());
      f->init();
      if (g.get())
        g->init();
    }

    void operator()(const ConstElemRange &range) const;
  };



// ------------------------------------------------------------
// System implementation
void System::project_vector (NumericVector<Number>& vector) const
{
  // Create a copy of the vector, which currently
  // contains the old data.
  AutoPtr<NumericVector<Number> >
    old_vector (vector.clone());

  // Project the old vector to the new vector
  this->project_vector (*old_vector, vector);
}


/**
 * This method projects the vector
 * via L2 projections or nodal
 * interpolations on each element.
 */
void System::project_vector (const NumericVector<Number>& old_v,
			     NumericVector<Number>& new_v) const
{
  START_LOG ("project_vector()", "System");

  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_v gives the solution on the
   * old mesh, while the \p new_v gives the solution (to be computed)
   * on the new mesh.
   */
  new_v.clear();

#ifdef LIBMESH_ENABLE_AMR

  // Resize the new vector and get a serial version.
  NumericVector<Number> *new_vector_ptr = NULL;
  AutoPtr<NumericVector<Number> > new_vector_built;
  NumericVector<Number> *local_old_vector;
  AutoPtr<NumericVector<Number> > local_old_vector_built;
  const NumericVector<Number> *old_vector_ptr = NULL;

  ConstElemRange active_local_elem_range
    (this->get_mesh().active_local_elements_begin(),
     this->get_mesh().active_local_elements_end());

  // If the old vector was uniprocessor, make the new
  // vector uniprocessor
  if (old_v.type() == SERIAL)
    {
      new_v.init (this->n_dofs(), false, SERIAL);
      new_vector_ptr = &new_v;
      old_vector_ptr = &old_v;
    }

  // Otherwise it is a parallel, distributed vector, which
  // we need to localize.
  else if (old_v.type() == PARALLEL)
    {
      // Build a send list for efficient localization
      BuildProjectionList projection_list(*this);
      Threads::parallel_reduce (active_local_elem_range,
				projection_list);

      // Create a sorted, unique send_list
      projection_list.unique();

      new_v.init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      new_vector_built = NumericVector<Number>::build(this->comm());
      local_old_vector_built = NumericVector<Number>::build(this->comm());
      new_vector_ptr = new_vector_built.get();
      local_old_vector = local_old_vector_built.get();
      new_vector_ptr->init(this->n_dofs(), false, SERIAL);
      local_old_vector->init(old_v.size(), false, SERIAL);
      old_v.localize(*local_old_vector, projection_list.send_list);
      local_old_vector->close();
      old_vector_ptr = local_old_vector;
    }
  else if (old_v.type() == GHOSTED)
    {
      // Build a send list for efficient localization
      BuildProjectionList projection_list(*this);
      Threads::parallel_reduce (active_local_elem_range,
				projection_list);

      // Create a sorted, unique send_list
      projection_list.unique();

      new_v.init (this->n_dofs(), this->n_local_dofs(),
                  this->get_dof_map().get_send_list(), false, GHOSTED);

      local_old_vector_built = NumericVector<Number>::build(this->comm());
      new_vector_ptr = &new_v;
      local_old_vector = local_old_vector_built.get();
      local_old_vector->init(old_v.size(), old_v.local_size(),
                             projection_list.send_list, false, GHOSTED);
      old_v.localize(*local_old_vector, projection_list.send_list);
      local_old_vector->close();
      old_vector_ptr = local_old_vector;
    }
  else // unknown old_v.type()
    {
      libMesh::err << "ERROR: Unknown old_v.type() == " << old_v.type()
		    << std::endl;
      libmesh_error();
    }

  // Note that the above will have zeroed the new_vector.
  // Just to be sure, assert that new_vector_ptr and old_vector_ptr
  // were successfully set before trying to deref them.
  libmesh_assert(new_vector_ptr);
  libmesh_assert(old_vector_ptr);

  NumericVector<Number> &new_vector = *new_vector_ptr;
  const NumericVector<Number> &old_vector = *old_vector_ptr;

  Threads::parallel_for (active_local_elem_range,
			 ProjectVector(*this,
				       old_vector,
				       new_vector)
			 );

  // Copy the SCALAR dofs from old_vector to new_vector
  // Note: We assume that all SCALAR dofs are on the
  // processor with highest ID
  if(this->processor_id() == (this->n_processors()-1))
  {
    const DofMap& dof_map = this->get_dof_map();
    for (unsigned int var=0; var<this->n_vars(); var++)
      if(this->variable(var).type().family == SCALAR)
        {
          // We can just map SCALAR dofs directly across
          std::vector<dof_id_type> new_SCALAR_indices, old_SCALAR_indices;
          dof_map.SCALAR_dof_indices (new_SCALAR_indices, var, false);
          dof_map.SCALAR_dof_indices (old_SCALAR_indices, var, true);
          const unsigned int new_n_dofs =
	    libmesh_cast_int<unsigned int>(new_SCALAR_indices.size());

          for (unsigned int i=0; i<new_n_dofs; i++)
          {
            new_vector.set( new_SCALAR_indices[i], old_vector(old_SCALAR_indices[i]) );
          }
        }
  }

  new_vector.close();

  // If the old vector was serial, we probably need to send our values
  // to other processors
  //
  // FIXME: I'm not sure how to make a NumericVector do that without
  // creating a temporary parallel vector to use localize! - RHS
  if (old_v.type() == SERIAL)
    {
      AutoPtr<NumericVector<Number> > dist_v = NumericVector<Number>::build(this->comm());
      dist_v->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      dist_v->close();

      for (dof_id_type i=0; i!=dist_v->size(); i++)
        if (new_vector(i) != 0.0)
          dist_v->set(i, new_vector(i));

      dist_v->close();

      dist_v->localize (new_v, this->get_dof_map().get_send_list());
      new_v.close();
    }
  // If the old vector was parallel, we need to update it
  // and free the localized copies
  else if (old_v.type() == PARALLEL)
    {
      // We may have to set dof values that this processor doesn't
      // own in certain special cases, like LAGRANGE FIRST or
      // HERMITE THIRD elements on second-order meshes
      for (dof_id_type i=0; i!=new_v.size(); i++)
        if (new_vector(i) != 0.0)
          new_v.set(i, new_vector(i));
      new_v.close();
    }

  this->get_dof_map().enforce_constraints_exactly(*this, &new_v);

#else

  // AMR is disabled: simply copy the vector
  new_v = old_v;

#endif // #ifdef LIBMESH_ENABLE_AMR

  STOP_LOG("project_vector()", "System");
}



/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (Number fptr(const Point& p,
                                           const Parameters& parameters,
                                           const std::string& sys_name,
                                           const std::string& unknown_name),
                               Gradient gptr(const Point& p,
                                             const Parameters& parameters,
                                             const std::string& sys_name,
                                             const std::string& unknown_name),
                               const Parameters& parameters) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->project_solution(&f, &g);
}


/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (FunctionBase<Number> *f,
                               FunctionBase<Gradient> *g) const
{
  this->project_vector(*solution, f, g);

  solution->localize(*current_local_solution, _dof_map->get_send_list());
}


/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (FEMFunctionBase<Number> *f,
                               FEMFunctionBase<Gradient> *g) const
{
  this->project_vector(*solution, f, g);

  solution->localize(*current_local_solution, _dof_map->get_send_list());
}


/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (Number fptr(const Point& p,
                                         const Parameters& parameters,
                                         const std::string& sys_name,
                                         const std::string& unknown_name),
                             Gradient gptr(const Point& p,
                                           const Parameters& parameters,
                                           const std::string& sys_name,
                                           const std::string& unknown_name),
                             const Parameters& parameters,
                             NumericVector<Number>& new_vector) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->project_vector(new_vector, &f, &g);
}

/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (NumericVector<Number>& new_vector,
                             FunctionBase<Number> *f,
                             FunctionBase<Gradient> *g) const
{
  START_LOG ("project_vector()", "System");

  Threads::parallel_for
    (ConstElemRange (this->get_mesh().active_local_elements_begin(),
		     this->get_mesh().active_local_elements_end() ),
     ProjectSolution(*this, f, g,
                     this->get_equation_systems().parameters,
		     new_vector)
    );

  // Also, load values into the SCALAR dofs
  // Note: We assume that all SCALAR dofs are on the
  // processor with highest ID
  if(this->processor_id() == (this->n_processors()-1))
  {
    // We get different scalars as different
    // components from a new-style f functor.
    DenseVector<Number> fout(this->n_components());
    bool filled_fout = false;

    const DofMap& dof_map = this->get_dof_map();
    for (unsigned int var=0; var<this->n_vars(); var++)
      if(this->variable(var).type().family == SCALAR)
        {
          if (!filled_fout)
            {
              (*f) (Point(), this->time, fout);
	      filled_fout = true;
	    }

          std::vector<dof_id_type> SCALAR_indices;
          dof_map.SCALAR_dof_indices (SCALAR_indices, var);
          const unsigned int n_SCALAR_dofs =
	    libmesh_cast_int<unsigned int>(SCALAR_indices.size());

          for (unsigned int i=0; i<n_SCALAR_dofs; i++)
          {
            const dof_id_type global_index = SCALAR_indices[i];
            const unsigned int component_index =
              this->variable_scalar_number(var,i);
            new_vector.set(global_index, fout(component_index));
          }
        }
  }

  new_vector.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  this->get_dof_map().enforce_constraints_exactly(*this, &new_vector);
#endif

  STOP_LOG("project_vector()", "System");
}


/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (NumericVector<Number>& new_vector,
                             FEMFunctionBase<Number> *f,
                             FEMFunctionBase<Gradient> *g) const
{
  START_LOG ("project_fem_vector()", "System");

  Threads::parallel_for
    (ConstElemRange (this->get_mesh().active_local_elements_begin(),
		     this->get_mesh().active_local_elements_end() ),
     ProjectFEMSolution(*this, f, g, new_vector)
    );

  // Also, load values into the SCALAR dofs
  // Note: We assume that all SCALAR dofs are on the
  // processor with highest ID
  if(this->processor_id() == (this->n_processors()-1))
  {
    // FIXME: Do we want to first check for SCALAR vars before building this? [PB]
    FEMContext context( *this );

    const DofMap& dof_map = this->get_dof_map();
    for (unsigned int var=0; var<this->n_vars(); var++)
      if(this->variable(var).type().family == SCALAR)
        {
	  // FIXME: We reinit with an arbitrary element in case the user
	  //        doesn't override FEMFunctionBase::component. Is there
	  //        any use case we're missing? [PB]
	  Elem *el = const_cast<Elem *>(*(this->get_mesh().active_local_elements_begin()));
	  context.pre_fe_reinit( *this, el );

          std::vector<dof_id_type> SCALAR_indices;
          dof_map.SCALAR_dof_indices (SCALAR_indices, var);
          const unsigned int n_SCALAR_dofs =
	    libmesh_cast_int<unsigned int>(SCALAR_indices.size());

          for (unsigned int i=0; i<n_SCALAR_dofs; i++)
          {
            const dof_id_type global_index = SCALAR_indices[i];
            const unsigned int component_index =
              this->variable_scalar_number(var,i);

            new_vector.set(global_index, f->component(context, component_index, Point(), this->time));
          }
        }
  }

  new_vector.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  this->get_dof_map().enforce_constraints_exactly(*this, &new_vector);
#endif

  STOP_LOG("project_fem_vector()", "System");
}


/**
 * This method projects components of an arbitrary boundary function
 * onto the solution via L2 projections and nodal interpolations on
 * each element.
 */
void System::boundary_project_solution
  (const std::set<boundary_id_type> &b,
   const std::vector<unsigned int> &variables,
   Number fptr(const Point& p,
               const Parameters& parameters,
               const std::string& sys_name,
               const std::string& unknown_name),
   Gradient gptr(const Point& p,
                 const Parameters& parameters,
                 const std::string& sys_name,
                 const std::string& unknown_name),
   const Parameters& parameters)
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->boundary_project_solution(b, variables, &f, &g);
}


/**
 * This method projects an arbitary boundary function onto the
 * solution via L2 projections and nodal interpolations on each
 * element.
 */
void System::boundary_project_solution
  (const std::set<boundary_id_type> &b,
   const std::vector<unsigned int> &variables,
   FunctionBase<Number> *f,
   FunctionBase<Gradient> *g)
{
  this->boundary_project_vector(b, variables, *solution, f, g);

  solution->localize(*current_local_solution);
}





/**
 * This method projects an arbitrary boundary function via L2
 * projections and nodal interpolations on each element.
 */
void System::boundary_project_vector
  (const std::set<boundary_id_type> &b,
   const std::vector<unsigned int> &variables,
   Number fptr(const Point& p,
               const Parameters& parameters,
               const std::string& sys_name,
               const std::string& unknown_name),
   Gradient gptr(const Point& p,
                 const Parameters& parameters,
                 const std::string& sys_name,
                 const std::string& unknown_name),
   const Parameters& parameters,
   NumericVector<Number>& new_vector) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->boundary_project_vector(b, variables, new_vector, &f, &g);
}

/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::boundary_project_vector
  (const std::set<boundary_id_type> &b,
   const std::vector<unsigned int> &variables,
   NumericVector<Number>& new_vector,
   FunctionBase<Number> *f,
   FunctionBase<Gradient> *g) const
{
  START_LOG ("boundary_project_vector()", "System");

  Threads::parallel_for
    (ConstElemRange (this->get_mesh().active_local_elements_begin(),
		     this->get_mesh().active_local_elements_end() ),
     BoundaryProjectSolution(b, variables, *this, f, g,
                             this->get_equation_systems().parameters,
		             new_vector)
    );

  // We don't do SCALAR dofs when just projecting the boundary, so
  // we're done here.

  new_vector.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  this->get_dof_map().enforce_constraints_exactly(*this, &new_vector);
#endif

  STOP_LOG("boundary_project_vector()", "System");
}



#ifndef LIBMESH_ENABLE_AMR
void ProjectVector::operator()(const ConstElemRange &) const
{
  libmesh_error();
}
#else
void ProjectVector::operator()(const ConstElemRange &range) const
{
  START_LOG ("operator()","ProjectVector");

  // A vector for Lagrange element interpolation, indicating if we
  // have visited a DOF yet.  Note that this is thread-local storage,
  // hence shared DOFS that live on thread boundaries may be doubly
  // computed.  It is expected that this will be more efficient
  // than locking a thread-global version of already_done, though.
  //
  // FIXME: we should use this for non-Lagrange coarsening, too
  std::vector<bool> already_done (system.n_dofs(), false);

  // The number of variables in this system
  const unsigned int n_variables = system.n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = system.get_mesh().mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      const Variable& variable = dof_map.variable(var);

      const FEType& base_fe_type = variable.type();

      if (base_fe_type.family == SCALAR)
        continue;

      // Get FE objects of the appropriate type
      AutoPtr<FEBase> fe (FEBase::build(dim, base_fe_type));
      AutoPtr<FEBase> fe_coarse (FEBase::build(dim, base_fe_type));

      // Create FE objects with potentially different p_level
      FEType fe_type, temp_fe_type;

      // Prepare variables for non-Lagrange projection
      AutoPtr<QBase> qrule     (base_fe_type.default_quadrature_rule(dim));
      AutoPtr<QBase> qedgerule (base_fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (base_fe_type.default_quadrature_rule(dim-1));
      std::vector<Point> coarse_qpoints;

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> >& phi_values =
	fe->get_phi();
      const std::vector<std::vector<Real> >& phi_coarse =
	fe_coarse->get_phi();

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();

      // The XYZ locations of the quadrature points on the
      // child element
      const std::vector<Point>& xyz_values =
	fe->get_xyz();


      // The global DOF indices
      std::vector<dof_id_type> new_dof_indices, old_dof_indices;
      // Side/edge local DOF indices
      std::vector<unsigned int> new_side_dofs, old_side_dofs;

      // Iterate over the elements in the range
      for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
	{
	  const Elem* elem = *elem_it;
          // If this element doesn't have an old_dof_object with dofs for the
	  // current system, then it must be newly added, so the user
	  // is responsible for setting the new dofs.

          // ... but we need a better way to test for that; the code
          // below breaks on any FE type for which the elem stores no
          // dofs.
	  // if (!elem->old_dof_object || !elem->old_dof_object->has_dofs(system.number()))
	  //  continue;
	  const Elem* parent = elem->parent();

          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          // Adjust the FE type for p-refined elements
          fe_type = base_fe_type;
          fe_type.order = static_cast<Order>(fe_type.order +
                                             elem->p_level());

          // We may need to remember the parent's p_level
          unsigned int old_parent_level = 0;

	  // Update the DOF indices for this element based on
          // the new mesh
	  dof_map.dof_indices (elem, new_dof_indices, var);

	  // The number of DOFs on the new element
	  const unsigned int new_n_dofs =
	    libmesh_cast_int<unsigned int>(new_dof_indices.size());

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(new_n_dofs, false); // bools
          std::vector<int> free_dof(new_n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (new_n_dofs); Ue.zero();

	  // Update the DOF indices based on the old mesh.
	  // This is done in one of three ways:
	  // 1.) If the child was just refined then it was not
	  //     active in the previous mesh & hence has no solution
	  //     values on it.  In this case simply project or
	  //     interpolate the solution from the parent, who was
          //     active in the previous mesh
	  // 2.) If the child was just coarsened, obtaining a
	  //     well-defined solution may require doing independent
	  //     projections on nodes, edges, faces, and interiors
	  // 3.) If the child was active in the previous
	  //     mesh, we can just copy coefficients directly
	  if (elem->refinement_flag() == Elem::JUST_REFINED)
	    {
	      libmesh_assert(parent);
              old_parent_level = parent->p_level();

              // We may have done p refinement or coarsening as well;
              // if so then we need to reset the parent's p level
              // so we can get the right DoFs from it
              if (elem->p_refinement_flag() == Elem::JUST_REFINED)
                {
                  libmesh_assert_greater (elem->p_level(), 0);
                  (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() - 1);
                }
              else if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
                {
                  (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() + 1);
                }

	      dof_map.old_dof_indices (parent, old_dof_indices, var);
            }
	  else
	    {
	      dof_map.old_dof_indices (elem, old_dof_indices, var);

              if (elem->p_refinement_flag() == Elem::DO_NOTHING)
	        libmesh_assert_equal_to (old_dof_indices.size(), new_n_dofs);
              if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
	        libmesh_assert (elem->has_children());
	    }

	  unsigned int old_n_dofs =
	    libmesh_cast_int<unsigned int>(old_dof_indices.size());

          if (fe_type.family != LAGRANGE) {

	    // For refined non-Lagrange elements, we do an L2
            // projection
            // FIXME: this will be a suboptimal and ill-defined
            // result if we're using non-nested finite element
            // spaces or if we're on a p-coarsened element!
	    if (elem->refinement_flag() == Elem::JUST_REFINED)
	      {
	        // Update the fe object based on the current child
                fe->attach_quadrature_rule (qrule.get());
	        fe->reinit (elem);

	        // The number of quadrature points on the child
	        const unsigned int n_qp = qrule->n_points();

                FEInterface::inverse_map (dim, fe_type, parent,
					  xyz_values, coarse_qpoints);

                fe_coarse->reinit(parent, &coarse_qpoints);

	        // Reinitialize the element matrix and vector for
	        // the current element.  Note that this will zero them
	        // before they are summed.
	        Ke.resize (new_n_dofs, new_n_dofs); Ke.zero();
	        Fe.resize (new_n_dofs); Fe.zero();


	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // The solution value at the quadrature point
	            Number val = libMesh::zero;

		    // Sum the function values * the DOF values
		    // at the point of interest to get the function value
		    // (Note that the # of DOFs on the parent need not be the
		    //  same as on the child!)
		    for (unsigned int i=0; i<old_n_dofs; i++)
		      {
		        val += (old_vector(old_dof_indices[i])*
			        phi_coarse[i][qp]);
		      }

	            // Now \p val contains the solution value of variable
	            // \p var at the qp'th quadrature point on the child
	            // element.  It has been interpolated from the parent
	            // in case the child was just refined.  Next we will
	            // construct the L2-projection matrix for the element.

	            // Construct the Mass Matrix
	            for (unsigned int i=0; i<new_n_dofs; i++)
		      for (unsigned int j=0; j<new_n_dofs; j++)
		        Ke(i,j) += JxW[qp]*phi_values[i][qp]*phi_values[j][qp];

	            // Construct the RHS
	            for (unsigned int i=0; i<new_n_dofs; i++)
		      Fe(i) += JxW[qp]*phi_values[i][qp]*val;

	          } // end qp loop

                Ke.cholesky_solve(Fe, Ue);

                // Fix up the parent's p level in case we changed it
                (const_cast<Elem *>(parent))->hack_p_level(old_parent_level);
	      }
            else if (elem->refinement_flag() == Elem::JUST_COARSENED)
	      {
                FEBase::coarsened_dof_values(old_vector, dof_map,
					     elem, Ue, var, true);
              }
	    // For unrefined uncoarsened elements, we just copy DoFs
	    else
	      {
                // FIXME - I'm sure this function would be about half
                // the size if anyone ever figures out how to improve
                // the DofMap interface... - RHS
                if (elem->p_refinement_flag() == Elem::JUST_REFINED)
                  {
                    libmesh_assert_greater (elem->p_level(), 0);
                    // P refinement of non-hierarchic bases will
                    // require a whole separate code path
                    libmesh_assert (fe->is_hierarchic());
                    temp_fe_type = fe_type;
                    temp_fe_type.order =
                      static_cast<Order>(temp_fe_type.order - 1);
                    unsigned int old_index = 0, new_index = 0;
                    for (unsigned int n=0; n != n_nodes; ++n)
                      {
		        const unsigned int nc =
		          FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                                       elem_type, n);
                        for (unsigned int i=0; i != nc; ++i)
                          {
                            Ue(new_index + i) =
                              old_vector(old_dof_indices[old_index++]);
                          }
                        new_index +=
		          FEInterface::n_dofs_at_node (dim, fe_type,
                                                       elem_type, n);
                      }
		    const unsigned int nc =
		      FEInterface::n_dofs_per_elem (dim, temp_fe_type,
                                                    elem_type);
                    for (unsigned int i=0; i != nc; ++i)
                      {
                        Ue(new_index++) =
                          old_vector(old_dof_indices[old_index+i]);
                      }
                  }
                else if (elem->p_refinement_flag() ==
                         Elem::JUST_COARSENED)
                  {
                    // P coarsening of non-hierarchic bases will
                    // require a whole separate code path
                    libmesh_assert (fe->is_hierarchic());
                    temp_fe_type = fe_type;
                    temp_fe_type.order =
                      static_cast<Order>(temp_fe_type.order + 1);
                    unsigned int old_index = 0, new_index = 0;
                    for (unsigned int n=0; n != n_nodes; ++n)
                      {
		        const unsigned int nc =
		          FEInterface::n_dofs_at_node (dim, fe_type,
                                                       elem_type, n);
                        for (unsigned int i=0; i != nc; ++i)
                          {
                            Ue(new_index++) =
                              old_vector(old_dof_indices[old_index+i]);
                          }
                        old_index +=
		          FEInterface::n_dofs_at_node (dim, temp_fe_type,
                                                       elem_type, n);
                      }
		    const unsigned int nc =
		      FEInterface::n_dofs_per_elem (dim, fe_type,
                                                    elem_type);
                    for (unsigned int i=0; i != nc; ++i)
                      {
                        Ue(new_index++) =
                          old_vector(old_dof_indices[old_index+i]);
                      }
                  }
                else
                  // If there's no p refinement, we can copy every DoF
		  for (unsigned int i=0; i<new_n_dofs; i++)
		    Ue(i) = old_vector(old_dof_indices[i]);
	      }
          }
	  else { // fe type is Lagrange
	    // Loop over the DOFs on the element
	    for (unsigned int new_local_dof=0;
	         new_local_dof<new_n_dofs; new_local_dof++)
	      {
	        // The global DOF number for this local DOF
	        const dof_id_type new_global_dof =
		  new_dof_indices[new_local_dof];

	        // The global DOF might lie outside of the bounds of a
	        // distributed vector.  Check for that and possibly continue
	        // the loop
	        if ((new_global_dof <  new_vector.first_local_index()) ||
		    (new_global_dof >= new_vector.last_local_index()))
		  continue;

	        // We might have already computed the solution for this DOF.
	        // This is likely in the case of a shared node, particularly
	        // at the corners of an element.  Check to see if that is the
	        // case
	        if (already_done[new_global_dof] == true)
		  continue;

		already_done[new_global_dof] = true;

	        if (elem->refinement_flag() == Elem::JUST_REFINED)
		  {
		    // The location of the child's node on the parent element
		    const Point point =
		      FEInterface::inverse_map (dim, fe_type, parent,
					        elem->point(new_local_dof));

		    // Sum the function values * the DOF values
		    // at the point of interest to get the function value
		    // (Note that the # of DOFs on the parent need not be the
		    //  same as on the child!)
		    for (unsigned int old_local_dof=0;
		         old_local_dof<old_n_dofs; old_local_dof++)
		      {
		        const dof_id_type old_global_dof =
			  old_dof_indices[old_local_dof];

		        Ue(new_local_dof) +=
                          (old_vector(old_global_dof)*
			  FEInterface::shape(dim, fe_type, parent,
					     old_local_dof, point));
		      }
		  }
	        else
		  {
		    // Get the old global DOF index
		    const dof_id_type old_global_dof =
		      old_dof_indices[new_local_dof];

		    Ue(new_local_dof) = old_vector(old_global_dof);
		  }
              } // end local DOF loop

            // We may have to clean up a parent's p_level
	    if (elem->refinement_flag() == Elem::JUST_REFINED)
              (const_cast<Elem *>(parent))->hack_p_level(old_parent_level);
          }  // end fe_type if()

	  // Lock the new_vector since it is shared among threads.
	  {
	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

	    for (unsigned int i = 0; i < new_n_dofs; i++)
	      if (Ue(i) != 0.)
		new_vector.set(new_dof_indices[i], Ue(i));
	  }
        }  // end elem loop
    } // end variables loop

  STOP_LOG ("operator()","ProjectVector");
}
#endif // LIBMESH_ENABLE_AMR



void BuildProjectionList::unique()
{
  // Sort the send list.  After this duplicated
  // elements will be adjacent in the vector
  std::sort(this->send_list.begin(),
	    this->send_list.end());

  // Now use std::unique to remove duplicate entries
  std::vector<dof_id_type>::iterator new_end =
    std::unique (this->send_list.begin(),
		 this->send_list.end());

  // Remove the end of the send_list.  Use the "swap trick"
  // from Effective STL
  std::vector<dof_id_type>
    (this->send_list.begin(), new_end).swap (this->send_list);
}



#ifndef LIBMESH_ENABLE_AMR
void BuildProjectionList::operator()(const ConstElemRange &)
{
  libmesh_error();
}
#else
void BuildProjectionList::operator()(const ConstElemRange &range)
{
  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  const dof_id_type first_old_dof = dof_map.first_old_dof();
  const dof_id_type end_old_dof   = dof_map.end_old_dof();

  // We can handle all the variables at once.
  // The old global DOF indices
  std::vector<dof_id_type> di;

  // Iterate over the elements in the range
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
    {
      const Elem* elem = *elem_it;
      // If this element doesn't have an old_dof_object with dofs for the
      // current system, then it must be newly added, so the user
      // is responsible for setting the new dofs.

      // ... but we need a better way to test for that; the code
      // below breaks on any FE type for which the elem stores no
      // dofs.
      // if (!elem->old_dof_object || !elem->old_dof_object->has_dofs(system.number()))
      //  continue;
      const Elem* parent = elem->parent();

      if (elem->refinement_flag() == Elem::JUST_REFINED)
	{
	  libmesh_assert(parent);
	  unsigned int old_parent_level = parent->p_level();

	  if (elem->p_refinement_flag() == Elem::JUST_REFINED)
	    {
	      // We may have done p refinement or coarsening as well;
	      // if so then we need to reset the parent's p level
	      // so we can get the right DoFs from it
	      libmesh_assert_greater (elem->p_level(), 0);
	      (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() - 1);
	    }
	  else if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
	    {
	      (const_cast<Elem *>(parent))->hack_p_level(elem->p_level() + 1);
	    }

	  dof_map.old_dof_indices (parent, di);

	  // Fix up the parent's p level in case we changed it
	  (const_cast<Elem *>(parent))->hack_p_level(old_parent_level);
	}
      else if (elem->refinement_flag() == Elem::JUST_COARSENED)
	{
          std::vector<dof_id_type> di_child;
          di.clear();
          for (unsigned int c=0; c != elem->n_children(); ++c)
            {
	      dof_map.old_dof_indices (elem->child(c), di_child);
              di.insert(di.end(), di_child.begin(), di_child.end());
            }
        }
      else
	dof_map.old_dof_indices (elem, di);

      for (unsigned int i=0; i != di.size(); ++i)
        if (di[i] < first_old_dof || di[i] >= end_old_dof)
          this->send_list.push_back(di[i]);
    }  // end elem loop
}
#endif // LIBMESH_ENABLE_AMR



void BuildProjectionList::join(const BuildProjectionList &other)
{
  // Joining simply requires I add the dof indices from the other object
  this->send_list.insert(this->send_list.end(),
			 other.send_list.begin(),
			 other.send_list.end());
}


void ProjectSolution::operator()(const ConstElemRange &range) const
{
  // We need data to project
  libmesh_assert(f.get());

  /**
   * This method projects an arbitrary solution to the current
   * mesh.  The input function \p f gives the arbitrary solution,
   * while the \p new_vector (which should already be correctly sized)
   * gives the solution (to be computed) on the current mesh.
   */

  // The number of variables in this system
  const unsigned int n_variables = system.n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = system.get_mesh().mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables in the system
  for (unsigned int var=0; var<n_variables; var++)
    {
      const Variable& variable = dof_map.variable(var);

      const FEType& fe_type = variable.type();

      if (fe_type.family == SCALAR)
        continue;

      const unsigned int var_component =
        system.variable_scalar_number(var, 0);

      // Get FE objects of the appropriate type
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

      // Prepare variables for projection
      AutoPtr<QBase> qrule     (fe_type.default_quadrature_rule(dim));
      AutoPtr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> >& phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient> > *dphi = NULL;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          // We'll need gradient data for a C1 projection
          libmesh_assert(g.get());

          const std::vector<std::vector<RealGradient> >&
            ref_dphi = fe->get_dphi();
          dphi = &ref_dphi;
        }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();

      // The XYZ locations of the quadrature points
      const std::vector<Point>& xyz_values =
	fe->get_xyz();

      // The global DOF indices
      std::vector<dof_id_type> dof_indices;
      // Side/edge DOF indices
      std::vector<unsigned int> side_dofs;

      // Iterate over all the elements in the range
      for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
	{
	  const Elem* elem = *elem_it;

          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

	  // Update the DOF indices for this element based on
          // the current mesh
	  dof_map.dof_indices (elem, dof_indices, var);

	  // The number of DOFs on the element
	  const unsigned int n_dofs =
	    libmesh_cast_int<unsigned int>(dof_indices.size());

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(n_dofs, false); // bools
          std::vector<int> free_dof(n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating nodes, then
          // hold those fixed and project edges, then
          // hold those fixed and project faces, then
          // hold those fixed and project interiors

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
	      const unsigned int nc =
		FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
                                             n);
              if (!elem->is_vertex(n))
                {
                  current_dof += nc;
                  continue;
                }
              if (cont == DISCONTINUOUS)
                {
                  libmesh_assert_equal_to (nc, 0);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function
              else if (cont == C_ZERO)
                {
                  libmesh_assert_equal_to (nc, 1);
		  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              // The hermite element vertex shape functions are weird
              else if (fe_type.family == HERMITE)
                {
                  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  // x derivative
                  Ue(current_dof) = grad(0);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Point nxminus = elem->point(n),
                            nxplus = elem->point(n);
                      nxminus(0) -= TOLERANCE;
                      nxplus(0) += TOLERANCE;
		      Gradient gxminus = g->component(var_component,
						      nxminus,
                                                      system.time);
                      Gradient gxplus = g->component(var_component,
                                                     nxplus,
                                                     system.time);
                      // y derivative
                      Ue(current_dof) = grad(1);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      // xy derivative
                      Ue(current_dof) = (gxplus(1) - gxminus(1))
                                        / 2. / TOLERANCE;
                      dof_is_fixed[current_dof] = true;
                      current_dof++;

                      if (dim > 2)
                        {
                          // z derivative
                          Ue(current_dof) = grad(2);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // xz derivative
                          Ue(current_dof) = (gxplus(2) - gxminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // We need new points for yz
                          Point nyminus = elem->point(n),
                                nyplus = elem->point(n);
                          nyminus(1) -= TOLERANCE;
                          nyplus(1) += TOLERANCE;
                          Gradient gyminus = g->component(var_component,
                                                          nyminus,
                                                          system.time);
                          Gradient gyplus = g->component(var_component,
                                                         nyplus,
                                                         system.time);
                          // xz derivative
                          Ue(current_dof) = (gyplus(2) - gyminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem->point(n),
                                nxmyp = elem->point(n),
                                nxpym = elem->point(n),
                                nxpyp = elem->point(n);
                          nxmym(0) -= TOLERANCE;
                          nxmym(1) -= TOLERANCE;
                          nxmyp(0) -= TOLERANCE;
                          nxmyp(1) += TOLERANCE;
                          nxpym(0) += TOLERANCE;
                          nxpym(1) -= TOLERANCE;
                          nxpyp(0) += TOLERANCE;
                          nxpyp(1) += TOLERANCE;
                          Gradient gxmym = g->component(var_component,
                                                        nxmym,
                                                        system.time);
                          Gradient gxmyp = g->component(var_component,
                                                        nxmyp,
                                                        system.time);
                          Gradient gxpym = g->component(var_component,
                                                        nxpym,
                                                        system.time);
                          Gradient gxpyp = g->component(var_component,
                                                        nxpyp,
                                                        system.time);
                          Number gxzplus = (gxpyp(2) - gxmyp(2))
                                         / 2. / TOLERANCE;
                          Number gxzminus = (gxpym(2) - gxmym(2))
                                          / 2. / TOLERANCE;
                          // xyz derivative
                          Ue(current_dof) = (gxzplus - gxzminus)
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              // Assume that other C_ONE elements have a single nodal
              // value shape function and nodal gradient component
              // shape functions
              else if (cont == C_ONE)
                {
                  libmesh_assert_equal_to (nc, 1 + dim);
		  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  for (unsigned int i=0; i!= dim; ++i)
                    {
		      Ue(current_dof) = grad(i);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                    }
                }
              else
                libmesh_error();
            }

          // In 3D, project any edge values next
          if (dim > 2 && cont != DISCONTINUOUS)
            for (unsigned int e=0; e != elem->n_edges(); ++e)
              {
		FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                          side_dofs);

                // Some edge dofs are on nodes and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;

                // There may be nothing to project
                if (!free_dofs)
                  continue;

	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new edge coefficients
                DenseVector<Number> Uedge(free_dofs);

                // Initialize FE data on the edge
                fe->attach_quadrature_rule (qedgerule.get());
	        fe->edge_reinit (elem, e);
	        const unsigned int n_qp = qedgerule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point
	            Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
	            // solution grad at the quadrature point
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
                                              xyz_values[qp],
                                              system.time);

                    // Form edge projection matrix
                    for (unsigned int sidei=0, freei=0;
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                                JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
                                                      (*dphi)[j][qp])
                                                      * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uedge);

                // Transfer new edge solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uedge(i)) < TOLERANCE);
                    ui = Uedge(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

	  // Project any side values (edges in 2D, faces in 3D)
          if (dim > 1 && cont != DISCONTINUOUS)
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              {
		FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                          side_dofs);

		// Some side dofs are on nodes/edges and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;

                // There may be nothing to project
                if (!free_dofs)
                  continue;

	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new side coefficients
                DenseVector<Number> Uside(free_dofs);

                // Initialize FE data on the side
                fe->attach_quadrature_rule (qsiderule.get());
	        fe->reinit (elem, s);
	        const unsigned int n_qp = qsiderule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point
	            Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
	            // solution grad at the quadrature point
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
                                              xyz_values[qp],
                                              system.time);

                    // Form side projection matrix
                    for (unsigned int sidei=0, freei=0;
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
			      Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                               JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
						      (*dphi)[j][qp])
                                                     * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uside);

                // Transfer new side solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uside(i)) < TOLERANCE);
                    ui = Uside(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

	  // Project the interior values, finally

	  // Some interior dofs are on nodes/edges/sides and
          // already fixed, others are free to calculate
          unsigned int free_dofs = 0;
          for (unsigned int i=0; i != n_dofs; ++i)
            if (!dof_is_fixed[i])
              free_dof[free_dofs++] = i;

          // There may be nothing to project
          if (free_dofs)
            {

	  Ke.resize (free_dofs, free_dofs); Ke.zero();
	  Fe.resize (free_dofs); Fe.zero();
          // The new interior coefficients
          DenseVector<Number> Uint(free_dofs);

          // Initialize FE data
          fe->attach_quadrature_rule (qrule.get());
	  fe->reinit (elem);
	  const unsigned int n_qp = qrule->n_points();

	  // Loop over the quadrature points
	  for (unsigned int qp=0; qp<n_qp; qp++)
	    {
	      // solution at the quadrature point
	      Number fineval = f->component(var_component,
                                            xyz_values[qp],
                                            system.time);
	      // solution grad at the quadrature point
	      Gradient finegrad;
              if (cont == C_ONE)
                finegrad = g->component(var_component,
                                        xyz_values[qp],
                                        system.time);

              // Form interior projection matrix
              for (unsigned int i=0, freei=0; i != n_dofs; ++i)
                {
                  // fixed DoFs aren't test functions
                  if (dof_is_fixed[i])
                    continue;
		  for (unsigned int j=0, freej=0; j != n_dofs; ++j)
                    {
                      if (dof_is_fixed[j])
			Fe(freei) -= phi[i][qp] * phi[j][qp] * JxW[qp]
                                     * Ue(j);
                      else
			Ke(freei,freej) += phi[i][qp] * phi[j][qp] *
                                           JxW[qp];
                      if (cont == C_ONE)
                        {
                          if (dof_is_fixed[j])
			    Fe(freei) -= ((*dphi)[i][qp] *
					 (*dphi)[j][qp]) * JxW[qp] *
                                         Ue(j);
                          else
			    Ke(freei,freej) += ((*dphi)[i][qp] *
						(*dphi)[j][qp]) *
                                               JxW[qp];
                        }
                      if (!dof_is_fixed[j])
                        freej++;
                    }
		  Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                  if (cont == C_ONE)
                    Fe(freei) += (finegrad * (*dphi)[i][qp]) * JxW[qp];
                  freei++;
                }
	    }
          Ke.cholesky_solve(Fe, Uint);

          // Transfer new interior solutions to element
	  for (unsigned int i=0; i != free_dofs; ++i)
            {
              Number &ui = Ue(free_dof[i]);
              libmesh_assert(std::abs(ui) < TOLERANCE ||
                     std::abs(ui - Uint(i)) < TOLERANCE);
              ui = Uint(i);
              dof_is_fixed[free_dof[i]] = true;
            }

            } // if there are free interior dofs

          // Make sure every DoF got reached!
	  for (unsigned int i=0; i != n_dofs; ++i)
            libmesh_assert(dof_is_fixed[i]);

	  const dof_id_type
	    first = new_vector.first_local_index(),
	    last  = new_vector.last_local_index();

	  // Lock the new_vector since it is shared among threads.
	  {
	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

	    for (unsigned int i = 0; i < n_dofs; i++)
	      // We may be projecting a new zero value onto
	      // an old nonzero approximation - RHS
	      // if (Ue(i) != 0.)
	      if ((dof_indices[i] >= first) &&
		  (dof_indices[i] <  last))
                {
                  new_vector.set(dof_indices[i], Ue(i));
                }
	  }
        }  // end elem loop
    } // end variables loop
}


void ProjectFEMSolution::operator()(const ConstElemRange &range) const
{
  // We need data to project
  libmesh_assert(f.get());

  /**
   * This method projects an arbitrary solution to the current
   * mesh.  The input function \p f gives the arbitrary solution,
   * while the \p new_vector (which should already be correctly sized)
   * gives the solution (to be computed) on the current mesh.
   */

  FEMContext context( system );

  // The number of variables in this system
  const unsigned int n_variables = context.n_vars();

  // The dimensionality of the current mesh
  const unsigned int dim = context.get_dim();

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;

  // FIXME: Need to generalize this to vector-valued elements. [PB]
  FEBase* fe = NULL;
  FEBase* side_fe = NULL;
  FEBase* edge_fe = NULL;

  // First, loop over all variables and make sure the shape functions phi will
  // get built as well as the shape function gradients if the gradient functor
  // is supplied.
  for (unsigned int var=0; var<n_variables; var++)
    {
      context.get_element_fe( var, fe );
      if (fe->get_fe_type().family == SCALAR)
	continue;
      if( dim > 1 )
	context.get_side_fe( var, side_fe );
      if( dim > 2 )
	context.get_edge_fe( var, edge_fe );

      fe->get_phi();
      if( dim > 1 )
	side_fe->get_phi();
      if( dim > 2 )
	edge_fe->get_phi();

      const FEContinuity cont = fe->get_continuity();
      if(cont == C_ONE)
	{
	  libmesh_assert(g.get());
	  fe->get_dphi();
	  side_fe->get_dphi();
	  if( dim > 2 )
	    edge_fe->get_dphi();
	}
    }

  // Now initialize any user requested shape functions
  f->init_context(context);
  if(g.get())
    g->init_context(context);

  std::vector<unsigned int> side_dofs;

  // Iterate over all the elements in the range
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
    {
      const Elem* elem = *elem_it;

      context.pre_fe_reinit(system, elem);

      // Loop over all the variables in the system
      for (unsigned int var=0; var<n_variables; var++)
	{
	  const Variable& variable = dof_map.variable(var);

	  const FEType& fe_type = variable.type();

	  if (fe_type.family == SCALAR)
	    continue;

	  // Per-subdomain variables don't need to be projected on
	  // elements where they're not active
	  if (!variable.active_on_subdomain(elem->subdomain_id()))
	    continue;

	  const FEContinuity cont = fe->get_continuity();

	  const unsigned int var_component =
	    system.variable_scalar_number(var, 0);

	  const std::vector<dof_id_type>& dof_indices =
	    context.get_dof_indices(var);

	  // The number of DOFs on the element
	  const unsigned int n_dofs =
	    libmesh_cast_int<unsigned int>(dof_indices.size());

	  // Fixed vs. free DoFs on edge/face projections
	  std::vector<char> dof_is_fixed(n_dofs, false); // bools
	  std::vector<int> free_dof(n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating nodes, then
          // hold those fixed and project edges, then
          // hold those fixed and project faces, then
          // hold those fixed and project interiors

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
	      const unsigned int nc =
		FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
                                             n);
              if (!elem->is_vertex(n))
                {
                  current_dof += nc;
                  continue;
                }
              if (cont == DISCONTINUOUS)
                {
                  libmesh_assert_equal_to (nc, 0);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function
              else if (cont == C_ZERO)
                {
                  libmesh_assert_equal_to (nc, 1);
		  Ue(current_dof) = f->component(context,
						 var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              // The hermite element vertex shape functions are weird
              else if (fe_type.family == HERMITE)
                {
                  Ue(current_dof) = f->component(context,
						 var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(context,
					       var_component,
                                               elem->point(n),
                                               system.time);
                  // x derivative
                  Ue(current_dof) = grad(0);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Point nxminus = elem->point(n),
                            nxplus = elem->point(n);
                      nxminus(0) -= TOLERANCE;
                      nxplus(0) += TOLERANCE;
		      Gradient gxminus = g->component(context,
						      var_component,
						      nxminus,
                                                      system.time);
                      Gradient gxplus = g->component(context,
						     var_component,
                                                     nxplus,
                                                     system.time);
                      // y derivative
                      Ue(current_dof) = grad(1);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      // xy derivative
                      Ue(current_dof) = (gxplus(1) - gxminus(1))
                                        / 2. / TOLERANCE;
                      dof_is_fixed[current_dof] = true;
                      current_dof++;

                      if (dim > 2)
                        {
                          // z derivative
                          Ue(current_dof) = grad(2);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // xz derivative
                          Ue(current_dof) = (gxplus(2) - gxminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // We need new points for yz
                          Point nyminus = elem->point(n),
                                nyplus = elem->point(n);
                          nyminus(1) -= TOLERANCE;
                          nyplus(1) += TOLERANCE;
                          Gradient gyminus = g->component(context,
							  var_component,
                                                          nyminus,
                                                          system.time);
                          Gradient gyplus = g->component(context,
							 var_component,
                                                         nyplus,
                                                         system.time);
                          // xz derivative
                          Ue(current_dof) = (gyplus(2) - gyminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem->point(n),
                                nxmyp = elem->point(n),
                                nxpym = elem->point(n),
                                nxpyp = elem->point(n);
                          nxmym(0) -= TOLERANCE;
                          nxmym(1) -= TOLERANCE;
                          nxmyp(0) -= TOLERANCE;
                          nxmyp(1) += TOLERANCE;
                          nxpym(0) += TOLERANCE;
                          nxpym(1) -= TOLERANCE;
                          nxpyp(0) += TOLERANCE;
                          nxpyp(1) += TOLERANCE;
                          Gradient gxmym = g->component(context,
							var_component,
                                                        nxmym,
                                                        system.time);
                          Gradient gxmyp = g->component(context,
							var_component,
                                                        nxmyp,
                                                        system.time);
                          Gradient gxpym = g->component(context,
							var_component,
                                                        nxpym,
                                                        system.time);
                          Gradient gxpyp = g->component(context,
							var_component,
                                                        nxpyp,
                                                        system.time);
                          Number gxzplus = (gxpyp(2) - gxmyp(2))
                                         / 2. / TOLERANCE;
                          Number gxzminus = (gxpym(2) - gxmym(2))
                                          / 2. / TOLERANCE;
                          // xyz derivative
                          Ue(current_dof) = (gxzplus - gxzminus)
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              // Assume that other C_ONE elements have a single nodal
              // value shape function and nodal gradient component
              // shape functions
              else if (cont == C_ONE)
                {
                  libmesh_assert_equal_to (nc, 1 + dim);
		  Ue(current_dof) = f->component(context,
						 var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(context,
					       var_component,
                                               elem->point(n),
                                               system.time);
                  for (unsigned int i=0; i!= dim; ++i)
                    {
		      Ue(current_dof) = grad(i);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                    }
                }
              else
                libmesh_error();
            }

          // In 3D, project any edge values next
          if (dim > 2 && cont != DISCONTINUOUS)
	    {
	      const std::vector<Point>& xyz_values = edge_fe->get_xyz();
	      const std::vector<Real>& JxW = edge_fe->get_JxW();

	      const std::vector<std::vector<Real> >& phi = edge_fe->get_phi();
	      const std::vector<std::vector<RealGradient> >* dphi = NULL;
	      if (cont == C_ONE)
		dphi = &(edge_fe->get_dphi());

	      for (unsigned int e=0; e != elem->n_edges(); ++e)
		{
		  context.edge = e;
		  context.edge_fe_reinit();

                  const QBase& qedgerule = context.get_edge_qrule();
                  const unsigned int n_qp = qedgerule.n_points();

		  FEInterface::dofs_on_edge(elem, dim, fe_type, e,
					    side_dofs);

		  // Some edge dofs are on nodes and already
		  // fixed, others are free to calculate
		  unsigned int free_dofs = 0;
		  for (unsigned int i=0; i != side_dofs.size(); ++i)
		    if (!dof_is_fixed[side_dofs[i]])
		      free_dof[free_dofs++] = i;

		  // There may be nothing to project
		  if (!free_dofs)
		    continue;

		  Ke.resize (free_dofs, free_dofs); Ke.zero();
		  Fe.resize (free_dofs); Fe.zero();
		  // The new edge coefficients
		  DenseVector<Number> Uedge(free_dofs);

		  // Loop over the quadrature points
		  for (unsigned int qp=0; qp<n_qp; qp++)
		    {
		      // solution at the quadrature point
		      Number fineval = f->component(context,
						    var_component,
						    xyz_values[qp],
						    system.time);
		      // solution grad at the quadrature point
		      Gradient finegrad;
		      if (cont == C_ONE)
			finegrad = g->component(context,
						var_component,
						xyz_values[qp],
						system.time);

		      // Form edge projection matrix
		      for (unsigned int sidei=0, freei=0;
			   sidei != side_dofs.size(); ++sidei)
			{
			  unsigned int i = side_dofs[sidei];
			  // fixed DoFs aren't test functions
			  if (dof_is_fixed[i])
			    continue;
			  for (unsigned int sidej=0, freej=0;
			       sidej != side_dofs.size(); ++sidej)
			    {
			      unsigned int j = side_dofs[sidej];
			      if (dof_is_fixed[j])
				Fe(freei) -= phi[i][qp] * phi[j][qp] *
				  JxW[qp] * Ue(j);
			      else
				Ke(freei,freej) += phi[i][qp] *
				  phi[j][qp] * JxW[qp];
			      if (cont == C_ONE)
				{
				  if (dof_is_fixed[j])
				    Fe(freei) -= ( (*dphi)[i][qp] *
						   (*dphi)[j][qp] ) *
				      JxW[qp] * Ue(j);
				  else
				    Ke(freei,freej) += ( (*dphi)[i][qp] *
							 (*dphi)[j][qp] )
				      * JxW[qp];
				}
			      if (!dof_is_fixed[j])
				freej++;
			    }
			  Fe(freei) += phi[i][qp] * fineval * JxW[qp];
			  if (cont == C_ONE)
			    Fe(freei) += (finegrad * (*dphi)[i][qp] ) *
			      JxW[qp];
			  freei++;
			}
		    }

		  Ke.cholesky_solve(Fe, Uedge);

		  // Transfer new edge solutions to element
		  for (unsigned int i=0; i != free_dofs; ++i)
		    {
		      Number &ui = Ue(side_dofs[free_dof[i]]);
		      libmesh_assert(std::abs(ui) < TOLERANCE ||
				     std::abs(ui - Uedge(i)) < TOLERANCE);
		      ui = Uedge(i);
		      dof_is_fixed[side_dofs[free_dof[i]]] = true;
		    }
		}
	    } // end if (dim > 2 && cont != DISCONTINUOUS)

	  // Project any side values (edges in 2D, faces in 3D)
          if (dim > 1 && cont != DISCONTINUOUS)
	    {
	      const std::vector<Point>& xyz_values = side_fe->get_xyz();
	      const std::vector<Real>& JxW = side_fe->get_JxW();

	      const std::vector<std::vector<Real> >& phi = side_fe->get_phi();
	      const std::vector<std::vector<RealGradient> >* dphi = NULL;
	      if (cont == C_ONE)
		dphi = &(side_fe->get_dphi());

	      for (unsigned int s=0; s != elem->n_sides(); ++s)
		{
		  FEInterface::dofs_on_side(elem, dim, fe_type, s,
					    side_dofs);

		  // Some side dofs are on nodes/edges and already
		  // fixed, others are free to calculate
		  unsigned int free_dofs = 0;
		  for (unsigned int i=0; i != side_dofs.size(); ++i)
		    if (!dof_is_fixed[side_dofs[i]])
		      free_dof[free_dofs++] = i;

		  // There may be nothing to project
		  if (!free_dofs)
		    continue;

		  Ke.resize (free_dofs, free_dofs); Ke.zero();
		  Fe.resize (free_dofs); Fe.zero();
		  // The new side coefficients
		  DenseVector<Number> Uside(free_dofs);

		  context.side = s;
		  context.side_fe_reinit();

                  const QBase& qsiderule = context.get_side_qrule();
                  const unsigned int n_qp = qsiderule.n_points();

		  // Loop over the quadrature points
		  for (unsigned int qp=0; qp<n_qp; qp++)
		    {
		      // solution at the quadrature point
		      Number fineval = f->component(context,
						    var_component,
						    xyz_values[qp],
						    system.time);
		      // solution grad at the quadrature point
		      Gradient finegrad;
		      if (cont == C_ONE)
			finegrad = g->component(context,
						var_component,
						xyz_values[qp],
						system.time);

		      // Form side projection matrix
		      for (unsigned int sidei=0, freei=0;
			   sidei != side_dofs.size(); ++sidei)
			{
			  unsigned int i = side_dofs[sidei];
			  // fixed DoFs aren't test functions
			  if (dof_is_fixed[i])
			    continue;
			  for (unsigned int sidej=0, freej=0;
			       sidej != side_dofs.size(); ++sidej)
			    {
			      unsigned int j = side_dofs[sidej];
			      if (dof_is_fixed[j])
				Fe(freei) -= phi[i][qp] * phi[j][qp] *
				  JxW[qp] * Ue(j);
			      else
				Ke(freei,freej) += phi[i][qp] *
				  phi[j][qp] * JxW[qp];
			      if (cont == C_ONE)
				{
				  if (dof_is_fixed[j])
				    Fe(freei) -= ( (*dphi)[i][qp] *
						   (*dphi)[j][qp] ) *
				      JxW[qp] * Ue(j);
				  else
				    Ke(freei,freej) += ( (*dphi)[i][qp] *
							 (*dphi)[j][qp] )
				      * JxW[qp];
				}
			      if (!dof_is_fixed[j])
				freej++;
			    }
			  Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
			  if (cont == C_ONE)
			    Fe(freei) += (finegrad * (*dphi)[i][qp] ) *
			      JxW[qp];
			  freei++;
			}
		    }

		  Ke.cholesky_solve(Fe, Uside);

		  // Transfer new side solutions to element
		  for (unsigned int i=0; i != free_dofs; ++i)
		    {
		      Number &ui = Ue(side_dofs[free_dof[i]]);
		      libmesh_assert(std::abs(ui) < TOLERANCE ||
				     std::abs(ui - Uside(i)) < TOLERANCE);
		      ui = Uside(i);
		      dof_is_fixed[side_dofs[free_dof[i]]] = true;
		    }
		}
	    }// end if (dim > 1 && cont != DISCONTINUOUS)

	  // Project the interior values, finally

	  // Some interior dofs are on nodes/edges/sides and
          // already fixed, others are free to calculate
          unsigned int free_dofs = 0;
          for (unsigned int i=0; i != n_dofs; ++i)
            if (!dof_is_fixed[i])
              free_dof[free_dofs++] = i;

          // There may be nothing to project
          if (free_dofs)
            {
	      context.elem_fe_reinit();

	      const QBase& qrule = context.get_element_qrule();
	      const unsigned int n_qp = qrule.n_points();
	      const std::vector<Point>& xyz_values = fe->get_xyz();
	      const std::vector<Real>& JxW = fe->get_JxW();

	      const std::vector<std::vector<Real> >& phi = fe->get_phi();
	      const std::vector<std::vector<RealGradient> >* dphi = NULL;
	      if (cont == C_ONE)
		dphi = &(side_fe->get_dphi());

	      Ke.resize (free_dofs, free_dofs); Ke.zero();
	      Fe.resize (free_dofs); Fe.zero();
	      // The new interior coefficients
	      DenseVector<Number> Uint(free_dofs);

	      // Loop over the quadrature points
	      for (unsigned int qp=0; qp<n_qp; qp++)
		{
		  // solution at the quadrature point
		  Number fineval = f->component(context,
						var_component,
						xyz_values[qp],
						system.time);
		  // solution grad at the quadrature point
		  Gradient finegrad;
		  if (cont == C_ONE)
		    finegrad = g->component(context,
					    var_component,
					    xyz_values[qp],
					    system.time);

		  // Form interior projection matrix
		  for (unsigned int i=0, freei=0; i != n_dofs; ++i)
		    {
		      // fixed DoFs aren't test functions
		      if (dof_is_fixed[i])
			continue;
		      for (unsigned int j=0, freej=0; j != n_dofs; ++j)
			{
			  if (dof_is_fixed[j])
			    Fe(freei) -= phi[i][qp] * phi[j][qp] * JxW[qp]
			      * Ue(j);
			  else
			    Ke(freei,freej) += phi[i][qp] * phi[j][qp] *
			      JxW[qp];
			  if (cont == C_ONE)
			    {
			      if (dof_is_fixed[j])
				Fe(freei) -= ( (*dphi)[i][qp] *
					       (*dphi)[j][qp] ) * JxW[qp] *
				  Ue(j);
			      else
				Ke(freei,freej) += ( (*dphi)[i][qp] *
						     (*dphi)[j][qp] ) *
				  JxW[qp];
			    }
			  if (!dof_is_fixed[j])
			    freej++;
			}
		      Fe(freei) += phi[i][qp] * fineval * JxW[qp];
		      if (cont == C_ONE)
			Fe(freei) += (finegrad * (*dphi)[i][qp] ) * JxW[qp];
		      freei++;
		    }
		}
	      Ke.cholesky_solve(Fe, Uint);

	      // Transfer new interior solutions to element
	      for (unsigned int i=0; i != free_dofs; ++i)
		{
		  Number &ui = Ue(free_dof[i]);
		  libmesh_assert(std::abs(ui) < TOLERANCE ||
				 std::abs(ui - Uint(i)) < TOLERANCE);
		  ui = Uint(i);
		  dof_is_fixed[free_dof[i]] = true;
		}

            } // if there are free interior dofs

          // Make sure every DoF got reached!
	  for (unsigned int i=0; i != n_dofs; ++i)
            libmesh_assert(dof_is_fixed[i]);

	  const numeric_index_type
	    first = new_vector.first_local_index(),
	    last  = new_vector.last_local_index();

	  // Lock the new_vector since it is shared among threads.
	  {
	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

	    for (unsigned int i = 0; i < n_dofs; i++)
	      // We may be projecting a new zero value onto
	      // an old nonzero approximation - RHS
	      // if (Ue(i) != 0.)
	      if ((dof_indices[i] >= first) &&
		  (dof_indices[i] <  last))
                {
                  new_vector.set(dof_indices[i], Ue(i));
                }
	  }
        }  // end variables loop
    } // end elem loop
}



void BoundaryProjectSolution::operator()(const ConstElemRange &range) const
{
  // We need data to project
  libmesh_assert(f.get());

  /**
   * This method projects an arbitrary boundary solution to the current
   * mesh.  The input function \p f gives the arbitrary solution,
   * while the \p new_vector (which should already be correctly sized)
   * gives the solution (to be computed) on the current mesh.
   */

  // The dimensionality of the current mesh
  const unsigned int dim = system.get_mesh().mesh_dimension();

  // The DofMap for this system
  const DofMap& dof_map = system.get_dof_map();

  // Boundary info for the current mesh
  const BoundaryInfo& boundary_info = *system.get_mesh().boundary_info;

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables we've been requested to project
  for (unsigned int v=0; v!=variables.size(); v++)
    {
      const unsigned int var = variables[v];

      const Variable& variable = dof_map.variable(var);

      const FEType& fe_type = variable.type();

      if (fe_type.family == SCALAR)
        continue;

      const unsigned int var_component =
        system.variable_scalar_number(var, 0);

      // Get FE objects of the appropriate type
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

      // Prepare variables for projection
      AutoPtr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> >& phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient> > *dphi = NULL;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          // We'll need gradient data for a C1 projection
          libmesh_assert(g.get());

          const std::vector<std::vector<RealGradient> >&
            ref_dphi = fe->get_dphi();
          dphi = &ref_dphi;
        }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();

      // The XYZ locations of the quadrature points
      const std::vector<Point>& xyz_values =
	fe->get_xyz();

      // The global DOF indices
      std::vector<dof_id_type> dof_indices;
      // Side/edge DOF indices
      std::vector<unsigned int> side_dofs;

      // Iterate over all the elements in the range
      for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
	{
	  const Elem* elem = *elem_it;

          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

	  // Find out which nodes, edges and sides are on a requested
	  // boundary:
          std::vector<bool> is_boundary_node(elem->n_nodes(), false),
                            is_boundary_edge(elem->n_edges(), false),
                            is_boundary_side(elem->n_sides(), false);
          for (unsigned char s=0; s != elem->n_sides(); ++s)
            {
              // First see if this side has been requested
              const std::vector<boundary_id_type>& bc_ids =
                boundary_info.boundary_ids (elem, s);
              bool do_this_side = false;
              for (unsigned int i=0; i != bc_ids.size(); ++i)
                if (b.count(bc_ids[i]))
                  {
                    do_this_side = true;
                    break;
                  }
              if (!do_this_side)
                continue;

              is_boundary_side[s] = true;

              // Then see what nodes and what edges are on it
              for (unsigned int n=0; n != elem->n_nodes(); ++n)
                if (elem->is_node_on_side(n,s))
                  is_boundary_node[n] = true;
              for (unsigned int e=0; e != elem->n_edges(); ++e)
                if (elem->is_edge_on_side(e,s))
                  is_boundary_edge[e] = true;
            }

	  // Update the DOF indices for this element based on
          // the current mesh
	  dof_map.dof_indices (elem, dof_indices, var);

	  // The number of DOFs on the element
	  const unsigned int n_dofs =
	    libmesh_cast_int<unsigned int>(dof_indices.size());

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(n_dofs, false); // bools
          std::vector<int> free_dof(n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
	  // solution.  We start by interpolating boundary nodes, then
	  // hold those fixed and project boundary edges, then hold
	  // those fixed and project boundary faces,

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
	      const unsigned int nc =
		FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
                                             n);
              if (!elem->is_vertex(n) || !is_boundary_node[n])
                {
                  current_dof += nc;
                  continue;
                }
              if (cont == DISCONTINUOUS)
                {
                  libmesh_assert_equal_to (nc, 0);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function
              else if (cont == C_ZERO)
                {
                  libmesh_assert_equal_to (nc, 1);
		  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              // The hermite element vertex shape functions are weird
              else if (fe_type.family == HERMITE)
                {
                  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  // x derivative
                  Ue(current_dof) = grad(0);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Point nxminus = elem->point(n),
                            nxplus = elem->point(n);
                      nxminus(0) -= TOLERANCE;
                      nxplus(0) += TOLERANCE;
                      Gradient gxminus = g->component(var_component,
                                                      nxminus,
                                                      system.time);
                      Gradient gxplus = g->component(var_component,
                                                     nxplus,
                                                     system.time);
                      // y derivative
                      Ue(current_dof) = grad(1);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      // xy derivative
                      Ue(current_dof) = (gxplus(1) - gxminus(1))
                                        / 2. / TOLERANCE;
                      dof_is_fixed[current_dof] = true;
                      current_dof++;

                      if (dim > 2)
                        {
                          // z derivative
                          Ue(current_dof) = grad(2);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // xz derivative
                          Ue(current_dof) = (gxplus(2) - gxminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // We need new points for yz
                          Point nyminus = elem->point(n),
                                nyplus = elem->point(n);
                          nyminus(1) -= TOLERANCE;
                          nyplus(1) += TOLERANCE;
                          Gradient gyminus = g->component(var_component,
                                                          nyminus,
                                                          system.time);
                          Gradient gyplus = g->component(var_component,
                                                         nyplus,
                                                         system.time);
                          // xz derivative
                          Ue(current_dof) = (gyplus(2) - gyminus(2))
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem->point(n),
                                nxmyp = elem->point(n),
                                nxpym = elem->point(n),
                                nxpyp = elem->point(n);
                          nxmym(0) -= TOLERANCE;
                          nxmym(1) -= TOLERANCE;
                          nxmyp(0) -= TOLERANCE;
                          nxmyp(1) += TOLERANCE;
                          nxpym(0) += TOLERANCE;
                          nxpym(1) -= TOLERANCE;
                          nxpyp(0) += TOLERANCE;
                          nxpyp(1) += TOLERANCE;
                          Gradient gxmym = g->component(var_component,
                                                        nxmym,
                                                        system.time);
                          Gradient gxmyp = g->component(var_component,
                                                        nxmyp,
                                                        system.time);
                          Gradient gxpym = g->component(var_component,
                                                        nxpym,
                                                        system.time);
                          Gradient gxpyp = g->component(var_component,
                                                        nxpyp,
                                                        system.time);
                          Number gxzplus = (gxpyp(2) - gxmyp(2))
                                         / 2. / TOLERANCE;
                          Number gxzminus = (gxpym(2) - gxmym(2))
                                          / 2. / TOLERANCE;
                          // xyz derivative
                          Ue(current_dof) = (gxzplus - gxzminus)
                                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              // Assume that other C_ONE elements have a single nodal
              // value shape function and nodal gradient component
              // shape functions
              else if (cont == C_ONE)
                {
                  libmesh_assert_equal_to (nc, 1 + dim);
		  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  for (unsigned int i=0; i!= dim; ++i)
                    {
		      Ue(current_dof) = grad(i);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                    }
                }
              else
                libmesh_error();
            }

          // In 3D, project any edge values next
          if (dim > 2 && cont != DISCONTINUOUS)
            for (unsigned int e=0; e != elem->n_edges(); ++e)
              {
                if (!is_boundary_edge[e])
                  continue;

		FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                          side_dofs);

                // Some edge dofs are on nodes and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;

                // There may be nothing to project
                if (!free_dofs)
                  continue;

	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new edge coefficients
                DenseVector<Number> Uedge(free_dofs);

                // Initialize FE data on the edge
                fe->attach_quadrature_rule (qedgerule.get());
	        fe->edge_reinit (elem, e);
	        const unsigned int n_qp = qedgerule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point
	            Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
	            // solution grad at the quadrature point
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
                                              xyz_values[qp],
                                              system.time);

                    // Form edge projection matrix
                    for (unsigned int sidei=0, freei=0;
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                                JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
                                                      (*dphi)[j][qp])
                                                      * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uedge);

                // Transfer new edge solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uedge(i)) < TOLERANCE);
                    ui = Uedge(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

	  // Project any side values (edges in 2D, faces in 3D)
          if (dim > 1 && cont != DISCONTINUOUS)
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              {
                if (!is_boundary_side[s])
                  continue;

		FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                          side_dofs);

		// Some side dofs are on nodes/edges and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (unsigned int i=0; i != side_dofs.size(); ++i)
                  if (!dof_is_fixed[side_dofs[i]])
                    free_dof[free_dofs++] = i;

                // There may be nothing to project
                if (!free_dofs)
                  continue;

	        Ke.resize (free_dofs, free_dofs); Ke.zero();
	        Fe.resize (free_dofs); Fe.zero();
                // The new side coefficients
                DenseVector<Number> Uside(free_dofs);

                // Initialize FE data on the side
                fe->attach_quadrature_rule (qsiderule.get());
	        fe->reinit (elem, s);
	        const unsigned int n_qp = qsiderule->n_points();

	        // Loop over the quadrature points
	        for (unsigned int qp=0; qp<n_qp; qp++)
	          {
	            // solution at the quadrature point
	            Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
	            // solution grad at the quadrature point
	            Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
                                              xyz_values[qp],
                                              system.time);

                    // Form side projection matrix
                    for (unsigned int sidei=0, freei=0;
                         sidei != side_dofs.size(); ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
			for (unsigned int sidej=0, freej=0;
                             sidej != side_dofs.size(); ++sidej)
                          {
                            unsigned int j = side_dofs[sidej];
                            if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                           JxW[qp] * Ue(j);
                            else
			      Ke(freei,freej) += phi[i][qp] *
						 phi[j][qp] * JxW[qp];
                            if (cont == C_ONE)
                              {
                                if (dof_is_fixed[j])
                                  Fe(freei) -= ((*dphi)[i][qp] *
					        (*dphi)[j][qp]) *
                                               JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
						      (*dphi)[j][qp])
                                                     * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                        if (cont == C_ONE)
			  Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                                       JxW[qp];
                        freei++;
                      }
	          }

                Ke.cholesky_solve(Fe, Uside);

                // Transfer new side solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number &ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uside(i)) < TOLERANCE);
                    ui = Uside(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

	  const dof_id_type
	    first = new_vector.first_local_index(),
	    last  = new_vector.last_local_index();

	  // Lock the new_vector since it is shared among threads.
	  {
	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

	    for (unsigned int i = 0; i < n_dofs; i++)
	      if (dof_is_fixed[i] &&
                  (dof_indices[i] >= first) &&
		  (dof_indices[i] <  last))
                {
                  new_vector.set(dof_indices[i], Ue(i));
                }
	  }
        }  // end elem loop
    } // end variables loop
}


} // namespace libMesh
