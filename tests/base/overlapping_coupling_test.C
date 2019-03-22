// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/ghosting_functor.h>
#include <libmesh/coupling_matrix.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_modification.h>


#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;


// This coupling functor is attempting to mimic in a basic way
// a particular use case. The idea is that the mesh has overlapping domains
// (we assume only two in this test) and that the elements are not
// topologically connected between the overlapping domains, but when the user
// solves their problem, some dofs on are coupled between the overlapping
// elements. So we want to exercise the GhostingFunctor for this scenario,
// in particular both algebraic ghosting (dof ghosting) and coupling
// ghosting (sparsity pattern). In the use case, the coupling happens at
// quadrature points, so we look for the overlap based on the quadrature points.
//
// Additionally, the use case also has an issue that the coupling functor
// cannot be added until after the System is initialized because the overlap
// depends on the solution. So we will exercise that case as well.
class OverlappingCouplingFunctor : public GhostingFunctor
{
  public:

  OverlappingCouplingFunctor( System & system )
    : GhostingFunctor(),
      _system(system),
      _mesh(system.get_mesh()),
      _point_locator(system.get_mesh().sub_point_locator()),
      _coupling_matrix(nullptr)
  {
    // This is a highly specalized coupling functor where we are assuming
    // only two subdomains that are overlapping so make sure there's only two.
    CPPUNIT_ASSERT_EQUAL(2, (int)_mesh.n_subdomains());

    // We're assuming a specific set of subdomain ids
    std::set<subdomain_id_type> ids;
    _mesh.subdomain_ids(ids);
    CPPUNIT_ASSERT( ids.find(1) != ids.end() );
    CPPUNIT_ASSERT( ids.find(2) != ids.end() );

    _point_locator->enable_out_of_mesh_mode();
  }

  virtual ~OverlappingCouplingFunctor(){};

  void set_coupling_matrix (std::unique_ptr<CouplingMatrix> & coupling_matrix)
  { _coupling_matrix = std::move(coupling_matrix); }

  virtual void operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   std::unordered_map<const Elem *,const CouplingMatrix*> & coupled_elements) override
  {

    for( const auto & elem : as_range(range_begin,range_end) )
      {
        std::set<subdomain_id_type> allowed_subdomains;
        unsigned int var = 0;
        if( elem->subdomain_id() == 1 )
          {
            var = _system.variable_number("V");
            allowed_subdomains.insert(2);
          }
        else if( elem->subdomain_id() == 2 )
          {
            var = _system.variable_number("U");
            allowed_subdomains.insert(1);
          }
        else
          CPPUNIT_FAIL("Something bad happpend");

        unsigned int dim = 2;
        FEType fe_type = _system.variable_type(var);
        std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
        QGauss qrule (dim, fe_type.default_quadrature_order());
        fe->attach_quadrature_rule (&qrule);

        const std::vector<libMesh::Point> & qpoints = fe->get_xyz();

        fe->reinit(elem);

        for ( const auto & qp : qpoints )
          {
            const Elem * overlapping_elem = (*_point_locator)( qp, &allowed_subdomains );

            if( overlapping_elem && overlapping_elem->processor_id() != p )
              coupled_elements.insert( std::make_pair(overlapping_elem,_coupling_matrix.get()) );
          }
      }
  }

private:

  System & _system;

  const MeshBase & _mesh;

  std::unique_ptr<PointLocatorBase> _point_locator;

  std::unique_ptr<CouplingMatrix> _coupling_matrix;
};

