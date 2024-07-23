#ifndef NUMERIC_VECTOR_TEST_H
#define NUMERIC_VECTOR_TEST_H

// test includes
#include "test_comm.h"

// libMesh includes
#include <libmesh/parallel.h>
#include <libmesh/fuzzy_equals.h>

#include "libmesh_cppunit.h"

#include <memory>

#define NUMERICVECTORTEST                       \
  CPPUNIT_TEST( testLocalize );                 \
  CPPUNIT_TEST( testLocalizeBase );             \
  CPPUNIT_TEST( testLocalizeIndices );          \
  CPPUNIT_TEST( testLocalizeIndicesBase );      \
  CPPUNIT_TEST( testLocalizeToOne );            \
  CPPUNIT_TEST( testLocalizeToOneBase );        \
  CPPUNIT_TEST( testNorms );                    \
  CPPUNIT_TEST( testNormsBase );                \
  CPPUNIT_TEST( testOperations );               \
  CPPUNIT_TEST( testOperationsBase );


template <class DerivedClass>
class NumericVectorTest : public CppUnit::TestCase {

protected:
  libMesh::Parallel::Communicator *my_comm = nullptr;

  std::string libmesh_suite_name;

  unsigned int block_size, local_size, global_size;

public:
  void setUp()
  {
    // By default we'll use the whole communicator in parallel;
    // Serial-only NumericVector subclasses will need to override
    // this and set something else first.
    if (!my_comm)
      my_comm = TestCommWorld;

    block_size  = 10;

    // a different size on each processor.
    local_size  = block_size +
      static_cast<unsigned int>(my_comm->rank());

    global_size = 0;
    for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
      global_size += (block_size + static_cast<unsigned int>(p));
  }

  void tearDown()
  {}

  template <class Base, class Derived>
  void Operations()
  {
    auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
    Base & v = *v_ptr;

    const libMesh::dof_id_type
      first = v.first_local_index(),
      last  = v.last_local_index();

    for (libMesh::dof_id_type n=first; n != last; n++)
      v.set (n, static_cast<libMesh::Number>(n+1));
    v.close();

    auto v_clone = v.clone();
    auto & vorig = *v_clone;

    CPPUNIT_ASSERT(relative_fuzzy_equals(v, vorig));

    v += v;

    CPPUNIT_ASSERT(!relative_fuzzy_equals(v, vorig));

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(2*n+2),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v *= v;

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(4*n*n+8*n+4),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v -= vorig;

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(4*n*n+7*n+3),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v /= vorig;

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(4*n+3),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v.add(1);
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(4*n+4),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v.add(2, vorig);
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(6*n+6),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v.scale(1/libMesh::Real(3));
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(2*n+2),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v = 4;
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(4),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v = vorig;
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              libMesh::Real(n+1),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v.reciprocal();
    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v(n)),
                              1/libMesh::Real(n+1),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v.dot(vorig)),
                            libMesh::Real(global_size),
                            libMesh::TOLERANCE*libMesh::TOLERANCE);
  }

  template <class Base, class Derived>
  void Norms()
  {
    auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
    Base & v = *v_ptr;

    const libMesh::dof_id_type
      first = v.first_local_index(),
      last  = v.last_local_index();

    for (libMesh::dof_id_type n=first; n != last; n++)
      v.set (n, -static_cast<libMesh::Number>(n));
    v.close();
    for (libMesh::dof_id_type n=first; n != last; n++)
      v.add (n, -static_cast<libMesh::Number>(n));
    v.close();

    const libMesh::Real exact_l1 =
      global_size * (global_size-libMesh::Real(1));
    const libMesh::Real exact_l2 =
      std::sqrt(libMesh::Real(global_size-1) *
                (2*global_size) *
                (2*global_size-1) / 3);
    LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v.sum()), -exact_l1,
                            libMesh::TOLERANCE*libMesh::TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(v.l1_norm(), exact_l1,
                            libMesh::TOLERANCE*libMesh::TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(v.l2_norm(), exact_l2,
                            libMesh::TOLERANCE*libMesh::TOLERANCE);
    LIBMESH_ASSERT_FP_EQUAL(v.linfty_norm(), 2*(global_size-libMesh::Real(1)),
                            libMesh::TOLERANCE*libMesh::TOLERANCE);

    auto u_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
    Base & u = *u_ptr;
    for (libMesh::dof_id_type n=first; n != last; n++)
      u.set (n, -static_cast<libMesh::Number>(n * n));
    u.close();

    Derived diff_derived(*my_comm, global_size, local_size);
    Base & diff = diff_derived;
    diff = u;
    diff -= v;
    diff.close();

    // Use a relative tolerance here; the norms are O(1e5) when our
    // processor count gets big enough, at which point O(1e-12) is
    // testing for exact equality...
    const auto diff_norm = diff.l2_norm();
    const auto norm_diff = u.l2_norm_diff(v);
    LIBMESH_ASSERT_FP_EQUAL(diff_norm, norm_diff,
                            (std::abs(diff_norm)+std::abs(norm_diff)) *
                            libMesh::TOLERANCE*libMesh::TOLERANCE);
  }

  template <class Base, class Derived>
  void Localize(bool to_one=false)
  {
    const libMesh::processor_id_type root_pid = 0;

    {
      auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
      Base & v = *v_ptr;
      std::vector<libMesh::Number> l(global_size);

      const libMesh::dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (libMesh::dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<libMesh::Number>(n));
      v.close();
      for (libMesh::dof_id_type n=first; n != last; n++)
        v.add (n, static_cast<libMesh::Number>(n));
      v.close();

      if (!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if (!to_one || my_comm->rank() == root_pid)
        // Yes I really mean v.size()
        for (libMesh::dof_id_type i=0; i<v.size(); i++)
          LIBMESH_ASSERT_FP_EQUAL(2.*libMesh::libmesh_real(i),
                                  libMesh::libmesh_real(l[i]),
                                  libMesh::TOLERANCE*libMesh::TOLERANCE);

      for (libMesh::dof_id_type n=first; n != last; n++)
      {
        const auto value = static_cast<libMesh::Number>(n);
        v.insert (&value, std::vector<libMesh::numeric_index_type>({n}));
      }
      v.close();

      if (!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if (!to_one || my_comm->rank() == root_pid)
        // Yes I really mean v.size()
        for (libMesh::dof_id_type i=0; i<v.size(); i++)
          LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(i),
                                  libMesh::libmesh_real(l[i]),
                                  libMesh::TOLERANCE*libMesh::TOLERANCE);
    }
  }


  template <class Base, class Derived>
  void LocalizeIndices()
  {
    {
      auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
      Base & v = *v_ptr;

      // Let's try pulling the same number of entries from each processor
      std::vector<libMesh::Number> values(block_size * my_comm->size());
      std::vector<libMesh::dof_id_type> indices;
      indices.reserve(block_size * my_comm->size());

      const libMesh::dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (libMesh::dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<libMesh::Number>(n));
      v.close();

      libMesh::dof_id_type end_index = 0;
      for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
        {
          end_index += block_size + p;
          for (unsigned int j = 0; j != block_size; ++j)
            indices.push_back(end_index-j-1);
        }

      v.localize(values, indices);

      end_index = 0;
      for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
        {
          end_index += block_size + p;
          for (unsigned int j = 0; j != block_size; ++j)
            LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(values[p*block_size+j]),
                                    libMesh::libmesh_real(end_index-j-1),
                                    libMesh::TOLERANCE*libMesh::TOLERANCE);
        }
    }
  }

  void testLocalize()
  {
    LOG_UNIT_TEST;

    Localize<DerivedClass,DerivedClass>();
  }

  void testLocalizeBase()
  {
    LOG_UNIT_TEST;

    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }

  void testLocalizeToOne()
  {
    LOG_UNIT_TEST;

    Localize<DerivedClass,DerivedClass >(true);
  }

  void testLocalizeToOneBase()
  {
    LOG_UNIT_TEST;

    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>(true);
  }

  void testLocalizeIndices()
  {
    LOG_UNIT_TEST;

    LocalizeIndices<DerivedClass,DerivedClass >();
  }

  void testLocalizeIndicesBase()
  {
    LOG_UNIT_TEST;

    LocalizeIndices<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }

  void testNorms()
  {
    LOG_UNIT_TEST;

    Norms<DerivedClass,DerivedClass >();
  }

  void testNormsBase()
  {
    LOG_UNIT_TEST;

    Norms<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }

  void testOperations()
  {
    LOG_UNIT_TEST;

    Operations<DerivedClass,DerivedClass >();
  }

  void testOperationsBase()
  {
    LOG_UNIT_TEST;

    Operations<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }
};

#endif
