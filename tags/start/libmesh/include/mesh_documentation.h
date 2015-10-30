//---------------------------------------------------
// Main page documentation
/**
 * \mainpage libMesh Documentation
 *
 * \section intro Introduction
 *
 * The \p libMesh library (to be renamed at a later date...) is a C++
 * framework for the numerical simulation of partial differential
 * equations on serial and parallel platforms.  Development began in
 * March 2002 with the intent of providing a friendly interface to a
 * number of high-quality software packages that are currently
 * available.  Currently the library supports 2D and 3D steady and
 * transient finite element simulations.  \p PETSc
 * (http://www-fp.mcs.anl.gov/petsc) is currently used for the solution of
 * linear systems on both serial and parallel platforms, however the
 * extensibility of the library allows for other solvers to be added
 * with ease.
 *
 * A major goal of the library is to provide support for adaptive mesh
 * refinement (AMR) computations in parallel while allowing a research
 * scientist to focus on the phyiscs they are modelling.  The library
 * currently offers:
 *
 *   - Partitioning Algorithms
 *     - Metis K-Way weighted graph partitioning
 *     - Hilbert and Morton-ordered space filling curves
 *
 *  - Generic 2D Finite Elements
 *     - 3 and 6 noded triangles (\p Tri3, \p Tri6)
 *     - 4, 8, and 9 noded quadrilaterals (\p Quad4, \p Quad8, \p Quad9)
 *     - 4 and 6 noded infinite quadrilaterals (\p InfQuad4, \p InfQuad6)
 *
 *  - Generic 3D Finite Elements
 *     - 4 and 10 noded tetrahedrals (\p Tet4, \p Tet10)
 *     - 8, 20, and 27 noded hexahedrals (\p Hex8, \p Hex20, \p Hex27)
 *     - 6 noded prisms (\p Prism6)
 *     - 5 noded pyramids (\p Pyramid5)
 *     - 8, 16, and 18 noded infinite hexahedrals (\p InfHex8,
 *         \p InfHex16, \p InfHex18)
 *     - 6 and 12 noded infinite prisms (\p InfPrism6, \p InfPrism12)
 *
 *  - Generic Finite Element Families
 *     - Lagrange
 *     - Hierarchic
 *     - Discontinuous Monomials
 *
 *  - Mesh IO & Format Translation Utilities
 *     - Ideas Universal (UNV) format (.unv)
 *     - Sandia National Labs ExodusII format (.exd)
 *     - Amtec Engineering's Tecplot binary format (.plt)
 *     - Amtec Engineering's Tecplot ascii format (.dat)
 *     - Los Alamos National Labs GMV format (.gmv)
 *     - AVS Unstructured UCD format (.ucd)
 *
 *
 *
 * \section install Installation
 *
 * The installation of the library is straightforward. The GNU
 * autoconf package is used to determine site-specific configuration
 * parameters. A standard build will occur after typing
 * \verbatim
     ./configure
     make
   \endverbatim
 * in the top-level project directory.  To see all the configuration
 * options type  
 * \verbatim
     ./configure --help
   \endverbatim     
 *
 * \subsection conf Configuration
 *
 * \subsection build Building the Library
 *
 * To build the library you need GNU Make and a reasonably current
 * C++ compiler. Currently, the library is known to work with the
 * following compilers: 
 *
 * - GNU GCC
 *   - gcc-3.2
 *   - gcc-2.95.3
 *   - gcc-2.96 (RedHat's compiler in the 7.x series)
 *
 * - Intel ICC/ECC
 *   - icc/ifc 7.0
 *   - icc/ifc 6.0
 *   - Earlier versions (<= 5.0) not supported.
 *
 * - SGI MIPSPro Compilers
 *   - Version 7.30
 *   - Not tested (but will likely work) with others
 *
 * - HP aCC
 * - IBM xlC version 5.0
 * - Compaq CXX 6.3.9.6
 *
 *
 * \subsection link Linking
 */
