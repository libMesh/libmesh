// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/getpot.h>

using namespace libMesh;

class GetPotTest : public CppUnit::TestCase {

public:
  CPPUNIT_TEST_SUITE( GetPotTest );

  CPPUNIT_TEST( testVariables );
  CPPUNIT_TEST( testSections );

  CPPUNIT_TEST_SUITE_END();

protected:

  std::stringstream inputfile;

  GetPot input;

public:

  void setup_inputfile()
  {

    const char* text =
"[Section1]\n"
"\n"
"  var1 = '5.0'\n"
"\n"
"  [./SubSection1]\n"
"\n"
"      var2 = 'blah'\n"
"\n"
"[]\n"
"\n"
"[Section2]\n"
"\n"
"   #var3 = '5'\n"
"\n"
"   [./Subsection2]\n"
"\n"
"       var4 = 'false'\n"
"\n"
"       [./Subsection3]\n"
"\n"
"           var5 = '6.02e23'\n"
"\n"
"   [../../Subsection4]\n"
"\n"
"           var6 = '42'\n"
"\n"
"[]\n"
"\n"
"[Section3]\n"
"\n"
"    unused_var = 'not_used'\n";

      inputfile << text;
  }

  void setUp()
  {
    this->setup_inputfile();
    input.parse_input_stream(inputfile);
  }

  void testVariables()
  {
    CPPUNIT_ASSERT( input.have_variable("Section1/var1") );
    Real var1 = input("Section1/var1", 1.0);
    CPPUNIT_ASSERT_EQUAL( var1, Real(5));

    CPPUNIT_ASSERT( input.have_variable("Section1/SubSection1/var2") );
    std::string var2 = input("Section1/SubSection1/var2", "DIE!");
    CPPUNIT_ASSERT_EQUAL( var2, std::string("blah") );

    // This variable is commented out in the input file so we shouldn't find it.
    CPPUNIT_ASSERT( !input.have_variable("Section2/var3") );
    int var3 = input("Section2/var3", 314);
    CPPUNIT_ASSERT_EQUAL( var3, 314 );

    CPPUNIT_ASSERT( input.have_variable("Section2/Subsection2/var4") );
    bool var4 = input("Section2/Subsection2/var4", true);
    CPPUNIT_ASSERT_EQUAL( var4, false );

    CPPUNIT_ASSERT( input.have_variable("Section2/Subsection2/Subsection3/var5") );
    Real var5 = input("Section2/Subsection2/Subsection3/var5", 3.1415);
    CPPUNIT_ASSERT_EQUAL( var5, Real(6.02e23));

    CPPUNIT_ASSERT( input.have_variable("Section2/Subsection4/var6") );
    unsigned int var6 = input("Section2/Subsection4/var6", 21);
    CPPUNIT_ASSERT_EQUAL( var6, (unsigned int)42 );

    // We didn't use Section3/unused_var so it should be a UFO
    std::vector<std::string> ufos = input.unidentified_variables();
    CPPUNIT_ASSERT_EQUAL( ufos.size(),  (std::vector<std::string>::size_type)1 );
    CPPUNIT_ASSERT_EQUAL( ufos[0], std::string("Section3/unused_var") );
  }

  void testSections()
  {
    // GetPot stores the '/' at the end of each section name
    CPPUNIT_ASSERT(input.have_section("Section1/"));
    CPPUNIT_ASSERT(input.have_section("Section1/SubSection1/"));
    CPPUNIT_ASSERT(input.have_section("Section2/Subsection2/"));
    CPPUNIT_ASSERT(input.have_section("Section2/Subsection2/Subsection3/"));
    CPPUNIT_ASSERT(input.have_section("Section2/Subsection4/"));
    CPPUNIT_ASSERT(input.have_section("Section3/"));

    // But we don't need to supply the trailing '/'
    CPPUNIT_ASSERT(input.have_section("Section1"));
    CPPUNIT_ASSERT(input.have_section("Section1/SubSection1"));
    CPPUNIT_ASSERT(input.have_section("Section2/Subsection2/Subsection3"));

    // No such thing as this section
    CPPUNIT_ASSERT(!input.have_section("ImNotASection/"));
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( GetPotTest );
