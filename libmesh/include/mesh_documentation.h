//---------------------------------------------------
// Main page documentation
/**
 * \mainpage libMesh Documentation
 *
 * \section intro Introduction
 *
 * The \p libMesh library is a C++ framework for the numerical
 * simulation of partial differential equations on serial and parallel
 * platforms.  Development began in March 2002 with the intent of
 * providing a friendly interface to a number of high-quality software
 * packages that are currently available.  Currently the library
 * supports 2D and 3D steady and transient finite element simulations.
 * \p PETSc (http://www-fp.mcs.anl.gov/petsc) is currently used for
 * the solution of linear systems on both serial and parallel
 * platforms, however the extensibility of the library allows for
 * other solvers to be added with ease.
 *
 * A major goal of the library is to provide support for adaptive mesh
 * refinement (AMR) computations in parallel while allowing a research
 * scientist to focus on the physics they are modelling.  The library
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
 *  - Mesh Creation & Modification Utilities
 *     - refine or coarsen a mesh: prescribed, level-one-compatible, or uniform
 *     - build equispaced n-cubes out of \p Edge2, \p Tri3, \p Tri6, 
 *          \p Quad4, \p Quad8, \p Quad9, \p Hex8, \p Hex20, \p Hex27
 *     - build circles/spheres out of \p Tri3, \p Tri6, \p Quad4,
 *          \p Quad8, \p Quad9, \p Hex8
 *     - add infinite elements to a volume-based mesh, handle symmetry planes
 *     - convert \p Quad4, \p Quad8, \p Quad9 to \p Tri3, \p Tri6
 *     - distort/translate/rotate/scale a mesh
 *     - determine bounding boxes/spheres
 *     - extract the mesh boundary for BC handling or as a separate mesh
 *
 * \p libMesh is actively developed at The University ot Texas at
 * Austin in the CFDLab and at Technische Universit&auml;t Hamburg-Harburg,
 * Mechanics and Ocean Engineering in Germany.  Many thanks to SourceForge 
 * for hosting the project at http://sourceforge.net/projects/libmesh
 * You can find out what is currently happening in the development branch
 * by checking out the online CVS logs at
 * http://libmesh.sourceforge.net/cvshtml
 *
 *
 * \section install Installation
 *
 *
 *
 * \subsection getsoftware Getting the Software
 *
 * The \p libMesh source can be downloaded from the project's SourceForge
 * homepage: http://sourceforge.net/projects/libmesh  Stable releases are
 * located there as compressed tar archives. You may also access the CVS
 * source tree for the latest code.  You can get read-only access to the
 * CVS repository via:
 * \verbatim
     cvs -d:pserver:anonymous@cvs.sourceforge.net:/cvsroot/libmesh co libmesh \endverbatim
 * If you would like to contribute to the project you will need a
 * SourceForge developer account.
 *
 *
 *
 * \subsection compilers Compilers
 *
 * The library is known to work with the following compilers:
 *
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
 *
 * \subsection conf Configuration
 *
 * Configuring the library is straightforward. The GNU
 * autoconf package is used to determine site-specific configuration
 * parameters. A standard build will occur after typing
 * \verbatim
     ./configure
     make \endverbatim
 * in the top-level project directory.  To see all the configuration
 * options type  
 * \verbatim
     ./configure --help \endverbatim     
 * The configure script will find your compilers and create the \p
 * Make.common file with the configuration for your site. If you want
 * to use different compilers than those found by configure you can
 * specify them in environment variables.  For example, the following
 * will build with the \p MIPS compilers on an SGI:
 * \verbatim
     CXX=CC CC=cc F77=f77 ./configure \endverbatim
 * Note that the FORTRAN compiler is not actually used to compile any
 * part of the library, but \p configure uses it to find out how to
 * link FORTRAN libraries with C++ code.
 *
 *
 *
 * \subsection build Building the Library
 *
 * To build the library you need GNU Make and a supported compiler, as
 * listed in section \ref conf.  After the library is configured
 * simply type \p make to build the library.  Typing \p make \p
 * meshtool will build a mesh translation tool using the library.
 *
 * The Makefiles distributed with the library look at the shell
 * environment variable \p METHOD to determine what mode the library
 * should be built in.  Valid values for \p METHOD are \p opt
 * (optimized mode, the default if \p METHOD is empty), \p debug
 * (build with debug symbols), and \p pro (build with profiling
 * support for use with \p gprof).
 *
 * \subsection link Linking
 */
