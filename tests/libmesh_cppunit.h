
#ifndef LIBMESH_LIBMESH_CPPUNIT_H
#define LIBMESH_LIBMESH_CPPUNIT_H

// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/libmesh_logging.h>

#if defined(LIBMESH_DEFAULT_QUADRUPLE_PRECISION) || \
    defined(LIBMESH_DEFAULT_TRIPLE_PRECISION)
# define LIBMESH_ASSERT_FP_EQUAL(expected,actual,tolerance) \
         CPPUNIT_ASSERT_DOUBLES_EQUAL(double(expected-actual),0,double(tolerance))
#else
# define LIBMESH_ASSERT_FP_EQUAL(expected,actual,tolerance) \
         CPPUNIT_ASSERT_DOUBLES_EQUAL(expected,actual,tolerance)
#endif

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  This should thus be the *last* header included in our test
// suite files.
#include <libmesh/ignore_warnings.h>
#ifdef __clang__
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif // GCC > 4.5
#endif // __GNUC__ && !__INTEL_COMPILER


extern libMesh::PerfLog * unitlog;

// Our regular libMesh::PerfItem uses fast_push/pop, so doesn't work
// well with strings like TestCase::getName() that get destructed
// before we're done printing logs.  Per-unit-test timings are
// coarse-grained so we can enable them regardless of --enable-perflog
// setting.
class UnitPerfItem
{
public:
  UnitPerfItem(std::string label,
               std::string header,
               libMesh::PerfLog & my_perflog) :
    _label(std::move(label)),
    _header(std::move(header)),
    _perflog(my_perflog)
  {
    _perflog.push(_label, _header);
  }

  ~UnitPerfItem()
  {
    _perflog.pop(_label, _header);
  }

private:
  const std::string _label;
  const std::string _header;
  libMesh::PerfLog & _perflog;
};

#define LIBMESH_CPPUNIT_TEST_SUITE(ATestFixtureType) \
  protected: \
  std::string libmesh_suite_name = #ATestFixtureType; \
  CPPUNIT_TEST_SUITE(ATestFixtureType)


#define LOG_UNIT_TEST UnitPerfItem TOKENPASTE2(perf_item_, __LINE__)(__func__,this->libmesh_suite_name,*unitlog)


#endif // LIBMESH_LIBMESH_CPPUNIT_H
