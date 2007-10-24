///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator         //
//                                                                           //
// Version 1.3                                                               //
// June 13, 2004                                                             //
//                                                                           //
// Copyright 2002, 2004                                                      //
// Hang Si                                                                   //
// Rathausstr. 9, 10178 Berlin, Germany                                      //
// si@wias-berlin.de                                                         //
//                                                                           //
// You can obtain TetGen via internet: http://tetgen.berlios.de.  It may be  //
//   freely copied, modified, and redistributed under the copyright notices  //
//   given in the file LICENSE.                                              //
//                                                                           //
// TetGen is a program for generating quality tetrahedral meshes and three-  //
//   dimensional Delaunay triangulations.  It currently computes exact       //
//   Delaunay tetrahedralizations, constrained Delaunay tetrahedralizations, //
//   and quality tetrahedral meshes. The latter are nicely graded and whose  //
//   tetrahedra have radius-edge ratio bounded, and are conforming Delaunay  //
//   if there are no input angles smaller than 60 degree.                    //
//                                                                           //
// TetGen incorporates a suit of geometrical and mesh generation algorithms. //
//   A brief description of these algorithms used in TetGen can be found in  //
//   the first section of the user's manual.  References are given for users //
//   who are interesting in these approaches. Below are the main references: //
//                                                                           //
//   The efficient Delaunay tetrahedralization algorithm is: H. Edelsbrunner //
//   and N. R. Shah, "Incremental Topological Flipping Works for Regular     //
//   Triangulations". Algorithmica 15: 223-241, 1996.                        //
//                                                                           //
//   The constrained Delaunay tetrahedralization algorithm is described in:  //
//   H. Si and K. Gaertner, "An Algorithm for Three-Dimensional Constrained  //
//   Delaunay Triangulations". Proc. of the 4th International Conference on  //
//   Engineering Computational Technology, Lisbon, September 2004.           //
//                                                                           //
//   The Delaunay refinement algorithm is from: J. R. Shewchuk, "Tetrahedral //
//   Mesh Generation by Delaunay Refinement". Proceedings of the 14th Annual //
//   Symposium on Computational Geometry, pages 86-95, 1998.                 //
//                                                                           //
// The mesh data structure of TetGen is a combination of two types of mesh   //
//   data structures.  The tetrahedron-based mesh data structure introduced  //
//   by Shewchuk is eligible to implement algorithms of generating Delaunay  //
//   tetrahedralizations. However, constrained Delaunay tetrahedralization   //
//   and quality mesh generation algorithms require other mesh elements      //
//   (subfaces, subsegments) be handled at the same time.  The triangle-edge //
//   data structure from Muecke is adopted for this purpose. Handling        //
//   these data types together is through a set of fast mesh manipulation    //
//   primitives.  References of these two data structures are found below:   //
//                                                                           //
//   J. R. Shewchuk, "Delaunay Refinement Mesh Generation". PhD thesis,      //
//   Carnegie Mellon University, 1997.                                       //
//                                                                           //
//   E. P. Muecke, "Shapes and Implementations in Three-Dimensional          //
//   Geometry". PhD thesis, Univ. of Illinois, Urbana, Illinois, 1993.       //
//                                                                           //
// The research of mesh generation is definitly on the move. A lot of state- //
//   of-the-art algorithms need to be implemented and evaluated.  I heartily //
//   welcome new algorithms especially for quality conforming Delaunay mesh  //
//   generation and anisotropic conforming Delaunay mesh generation. If you  //
//   have any idea on new approaches, please please kindly let me know.      //
//                                                                           //
// TetGen is supported by the "pdelib" project of Weierstrass Institute for  //
//   Applied Analysis and Stochastics (WIAS) in Berlin.  It is a collection  //
//   of software components for solving non-linear partial differential      //
//   equations including 2D and 3D mesh generators, sparse matrix solvers,   //
//   and scientific visualization tools, etc.  For more information please   //
//   see:   http://www.wias-berlin.de/software/pdelib.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgen.h                                                                  //
//                                                                           //
// Header file of the TetGen library. Also is the user-level header file.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen Library Overview                                                   //
//                                                                           //
// TetGen library is comprised by several data types and global functions.   //
//                                                                           //
// If you quickly go through this file, you will find there are only three   //
// main data types defined,  which are "tetgenio", "tetgenbehavior", and     //
// "tetgenmesh". Tetgenio is used to pass data into and out of mesh routines //
// of the library;  tetgenbehavior sets the command line options selected by //
// user and thus controls the behaviors of TetGen;  tetgenmesh, the biggest  //
// data type I've ever defined,  contains everything for creating Delaunay   //
// tetrahedralizations and tetrahedral meshes. These data types are defined  //
// as C++ classes.                                                           //
//                                                                           //
// There are few global functions as well.  "tetrahedralize()" is the (only) //
// user interface for calling TetGen from other programs.  Two functions     //
// "orient3d()" and "insphere()" are incorporated from a public C code for   //
// performing exact geometrical tests.                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef tetgenH
#define tetgenH

// To compile TetGen as a library (e.g. libtet.a) but not as an executable
//   program, define the TETLIBRARY symbol.  The library of TetGen can be
//   linked with programs which want to call TetGen as a function.

// #define TETLIBRARY

// Uncomment the following line to disable assert macros. These macros are
//   inserted in places where I hope to catch bugs. Somewhat, they slow down
//   the speed of TetGen.  They can be ignored by adding the -DNDEBUG
//   compiler switch or uncomment the following line 

// #define NDEBUG

// To insert lots of self-checks for internal errors, define the SELF_CHECK
//   symbol.  This will slow down the program significantly. 

// #define SELF_CHECK

// For single precision ( which will save some memory and reduce paging ),
//   define the symbol SINGLE by using the -DSINGLE compiler switch or by
//   writing "#define SINGLE" below.
//
// For double precision ( which will allow you to refine meshes to a smaller
//   edge length), leave SINGLE undefined.

// #define SINGLE

#ifdef SINGLE
  #define REAL float
#else
  #define REAL double
#endif 	// not defined SINGLE

// Here is the most general used head files for all C/C++ codes

#include <stdio.h>            // Standard IO: FILE, NULL, EOF, printf(), ...
#include <stdlib.h>        // Standard lib: abort(), system(), getenv(), ...
#include <string.h>         // String lib: strcpy(), strcat(), strcmp(), ...
#include <math.h>                     // Math lib: sin(), sqrt(), pow(), ...
#include <assert.h>
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenio data type                                                    //
//                                                                           //
// Used to pass data into and out of the library of TetGen.                  //
//                                                                           //
// If you want to program with the library of TetGen, it's necessary for you //
// to understand the tetgenio data type, while the other two data types can  //
// be hidden through calling the global function "tetrahedralize()".  As you //
// will see below, that basically tetgenio is nothing more than a collection //
// of arrays. These arrays are used to store points, tetrahedra, (triangular)//
// faces, boundary markers, and so forth.  They are used to describe data in //
// input & output files of TetGen.  If you understand TetGen's file formats, //
// then it is straighforward for you to understand these arrays.  The file   //
// formats of TetGen are described in the third section of the user's manual.//
//                                                                           //
// Once you create an object of tetgenio, all arrays are initialized to NULL.//
// This is done by routine "initialize()", it is automatically called by the //
// constructor. Before you set data into these arrays, you need to allocate  //
// enough memory for them.  After you use the object, you need to clear the  //
// memory occupied by these arrays.  Routine "deinitialize()" will be auto-  //
// matically called on deletion of the object. It will clear the memory in   //
// each array if it is not a NULL. However, it assumes that the memory is    //
// allocated by C++ operator 'new'. If you use malloc() to allocate memory,  //
// you should free them yourself, after they're freed, call "initialize()"   //
// once to disable "deinitialize()".                                         //
//                                                                           //
// In all cases, the first item in any array is stored starting at index [0].//
// However, that item is item number `firstnumber' which may be '0' or '1'.  //
// Be sure to set the 'firstnumber' be '1' if your indices pointing into the //
// pointlist is starting from '1'. Default, it is initialized be '0'.        //
//                                                                           //
// Tetgenio also contains routines for reading and writing TetGen's files as //
// well.  Both the library of TetGen and TetView use these routines to parse //
// input files, i.e., .node, .poly, .smesh, .ele, .face, and .edge files.    //
// Other routines are provided mainly for debugging purpose.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenio {

  public:

    // Maximum number of characters in a file name (including the null).
    enum {FILENAMESIZE = 1024};

    // Maxi. numbers of chars in a line read from a file (incl. the null).
    enum {INPUTLINESIZE = 1024};

    // The polygon data structure.  A "polygon" is a planar polygon. It can
    //   be arbitrary shaped (convex or non-convex) and bounded by non-
    //   crossing segments, i.e., the number of vertices it has indictes the
    //   same number of edges.
    // 'vertexlist' is a list of vertex indices (integers), its length is
    //   indicated by 'numberofvertices'.  The vertex indices are odered in
    //   either counterclockwise or clockwise way.
    typedef struct {
      int *vertexlist;
      int numberofvertices;
    } polygon;

    static void init(polygon* p) {
      p->vertexlist = (int *) NULL;
      p->numberofvertices = 0;
    }

    // The facet data structure.  A "facet" is a planar facet.  It is used
    //   to represent a planar straight line graph (PSLG) in two dimension.
    //   A PSLG contains a list of polygons. It also may conatin holes in it,
    //   indicated by a list of hole points (their coordinates).
    typedef struct {
      polygon *polygonlist;
      int numberofpolygons;
      REAL *holelist;
      int numberofholes;
    } facet;

    static void init(facet* f) {
      f->polygonlist = (polygon *) NULL;
      f->numberofpolygons = 0;
      f->holelist = (REAL *) NULL;
      f->numberofholes = 0;
    }

  public:

    // Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
    int firstnumber; 

    // Dimension of the mesh (2 or 3), default is 3.
    int mesh_dim;

    // `pointlist':  An array of point coordinates.  The first point's x
    //   coordinate is at index [0] and its y coordinate at index [1], its
    //   z coordinate is at index [2], followed by the coordinates of the
    //   remaining points.  Each point occupies three REALs. 
    // `pointattributelist':  An array of point attributes.  Each point's
    //   attributes occupy `numberofpointattributes' REALs. 
    // 'addpointlist':  An array of additional point coordinates.
    // `pointmarkerlist':  An array of point markers; one int per point.
    REAL *pointlist;
    REAL *pointattributelist;
    REAL *addpointlist;
    int *pointmarkerlist;
    int numberofpoints;
    int numberofpointattributes;
    int numberofaddpoints;
 
    // `elementlist':  An array of element (triangle or tetrahedron) corners. 
    //   The first element's first corner is at index [0], followed by its
    //   other corners in counterclockwise order, followed by any other
    //   nodes if the element represents a nonlinear element.  Each element
    //   occupies `numberofcorners' ints.
    // `elementattributelist':  An array of element attributes.  Each
    //   element's attributes occupy `numberofelementattributes' REALs.
    // `elementconstraintlist':  An array of constraints, i.e. triangle's
    //   area or tetrahedron's volume; one REAL per element.  Input only.
    // `neighborlist':  An array of element neighbors; 3 or 4 ints per
    //   element.  Output only.
    int *tetrahedronlist;
    REAL *tetrahedronattributelist;
    REAL *tetrahedronvolumelist;
    int *neighborlist;
    int numberoftetrahedra;
    int numberofcorners;
    int numberoftetrahedronattributes;

    // `facetlist':  An array of facets.  Each entry is a structure of facet.
    // `facetmarkerlist':  An array of facet markers; one int per facet.
    facet *facetlist;
    int *facetmarkerlist;
    int numberoffacets;

    // `holelist':  An array of holes.  The first hole's x, y and z
    //   coordinates  are at indices [0], [1] and [2], followed by the
    //   remaining holes. Three REALs per hole. 
    REAL *holelist;
    int numberofholes;

    // `regionlist':  An array of regional attributes and area or volume
    //   constraints. The first constraint's x, y and z coordinates are at
    //   indices [0], [1] and [2], followed by the regional attribute and
    //   index [3], followed by the maximum area or volume at index [4],
    //   followed by the remaining area or volume constraints. Five REALs
    //   per constraint. 
    // Note that each regional attribute is used only if you select the `A'
    //   switch, and each constraint is used only if you select the `a'
    //   switch (with no number following), but omitting one of these
    //   switches does not change the memory layout.
    REAL *regionlist;
    int numberofregions;

    // `trifacelist':  An array of triangular face endpoints.  The first
    //   face's endpoints are at indices [0], [1] and [2], followed by the
    //   remaining faces.  Three ints per face.
    // `trifacemarkerlist':  An array of face markers; one int per face.
    int *trifacelist;
    int *trifacemarkerlist;
    int numberoftrifaces;

    // `edgelist':  An array of edge endpoints.  The first edge's endpoints
    //   are at indices [0] and [1], followed by the remaining edges.  Two
    //   ints per edge.
    // `edgemarkerlist':  An array of edge markers; one int per edge.
    int *edgelist;
    int *edgemarkerlist;
    int numberofedges;

  public:

    // Initialize routine.
    void initialize();
    void deinitialize();

    // Input & output routines.
    bool load_node_call(FILE* infile, int markers, char* nodefilename);
    bool load_node(char* filename);
    bool load_addnodes(char* filename);
    bool load_poly(char* filename);
    bool load_off(char* filename);
    bool load_ply(char* filename);
    bool load_stl(char* filename);
    bool load_medit(char* filename);
    bool load_plc(char* filename, int object);
    bool load_tetmesh(char* filename);
    void save_nodes(char* filename);
    void save_elements(char* filename);
    void save_faces(char* filename);
    void save_edges(char* filename);
    void save_neighbors(char* filename);
    void save_poly(char* filename);

    // Read line and parse string functions.
    char *readline(char* string, FILE* infile, int *linenumber);
    char *findnextfield(char* string);
    char *readnumberline(char* string, FILE* infile, char* infilename);
    char *findnextnumber(char* string);

    // Constructor and destructor.
    tetgenio() {initialize();}
    ~tetgenio() {deinitialize();}
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenbehavior data type                                              //
//                                                                           //
// Used to parse command line switches and file names.                       //
//                                                                           //
// It includes a list of variables corresponding to the commandline switches //
// for control the behavior of TetGen.  These varibales are all initialized  //
// to their default values.                                                  //
//                                                                           //
// Routine "parse_commandline()" defined in this data type is used to change //
// the vaules of the variables. This routine accepts the standard parameters //
// ('argc' and 'argv') that pass to C/C++ main() function. It also accepts a //
// string which contains the command line options.                           //
//                                                                           //
// You don't need to understand this data type. It can be implicitly called  //
// by the global function "tetrahedralize()" defined below.  The necessary   //
// thing you need to know is the meaning of command line switches of TetGen. //
// They are described in the third section of the user's manual.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenbehavior {

  public:

    // Labels define the objects which are acceptable by TetGen. They are 
    //   recognized by TetGen from the file extensions.
    //   - NODES, a list of nodes (.node); 
    //   - POLY, a piecewise linear complex (.poly or .smesh); 
    //   - OFF, a polyhedron (.off, Geomview's file format); 
    //   - PLY, a polyhedron (.ply, file format from gatech);
    //   - STL, a surface mesh (.stl, stereolithography format);
    //   - MEDIT, a surface mesh (.mesh, Medit's file format); 
    //   - MESH, a tetrahedral mesh (.ele).
    //   If no extension is available, the imposed commandline switch
    //   (-p or -r) implies the object. 

    enum objecttype {NONE, NODES, POLY, OFF, PLY, STL, MEDIT, MESH};

    // Variables of command line switches.  After each variable are the 
    //   corresponding switch and its default value.  Read the user's manul
    //   or online instructions to find out the meaning of these switches.

    int plc;                                              // '-p' switch, 0.
    int refine;                                           // '-r' switch, 0.
    int quality;                                          // '-q' switch, 0.
    REAL minratio;                         // number after '-q' switch, 2.0.
    REAL goodratio;               // number calculated from 'minratio', 0.0. 
    REAL minangle;                             // minimum angle bound, 20.0.
    REAL goodangle;                      // cosine squared of minangle, 0.0.
    int varvolume;                         // '-a' switch without number, 0.
    int fixedvolume;                          // '-a' switch with number, 0.
    REAL maxvolume;                       // number after '-a' switch, -1.0.
    int removesliver;                                     // '-s' switch, 0.
    REAL maxdihedral;                      // number after '-s' switch, 0.0.
    int insertaddpoints;                                  // '-i' switch, 0.
    int regionattrib;                                     // '-A' switch, 0.
    REAL epsilon;                       // number after '-T' switch, 1.0e-8.
    int nomerge;           // not merge two coplanar facets, '-M' switch, 0.
    int detectinter;                                      // '-d' switch, 0.
    int checkclosure;                                     // '-c' switch, 0.
    int zeroindex;                                        // '-z' switch, 0.
    int order;             // element order, specified after '-o' switch, 1.
    int facesout;                                         // '-f' switch, 0.
    int edgesout;                                         // '-e' switch, 0.
    int neighout;                                         // '-n' switch, 0.
    int meditview;                                        // '-g' switch, 0.
    int gidview;                                          // '-G' switch, 0.
    int geomview;                                         // '-O' switch, 0.
    int nobound;                                          // '-B' switch, 0.
    int nonodewritten;                                    // '-N' switch, 0.
    int noelewritten;                                     // '-E' switch, 0.
    int nofacewritten;                                    // '-F' switch, 0.
    int noiterationnum;                                   // '-I' switch, 0.
    int nobisect;          // count of how often '-Y' switch is selected, 0.
    int noflip;                     // do not perform flips. '-Y' switch. 0.
    int steiner;                             // number after '-S' switch. 0.
    int dopermute;                        // do permutation. '-P' switch, 0.
    int srandseed;          // number of a random seed after '-P' switch, 1.
    int docheck;                                          // '-C' switch, 0.
    int quiet;                                            // '-Q' switch, 0.
    int verbose;           // count of how often '-V' switch is selected, 0.
    int useshelles;              // '-p', '-r', '-q', 'd', or 'c' switch, 0.
    enum objecttype object;         // determined by -p, or -r switch. NONE.

    // Variables used to save command line switches and in/out file names.
    char commandline[1024];
    char infilename[1024];
    char outfilename[1024];

    // Default initialize and de-initialize functions.
    tetgenbehavior();
    ~tetgenbehavior() {}

    void versioninfo();
    void syntax();
    void usage();

    // Command line parse routine.
    bool parse_commandline(int argc, char **argv);
    bool parse_commandline(char *switches) {
      return parse_commandline(0, &switches);
    }
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Geometric predicates                                                      //
//                                                                           //
// Return one of the values +1, 0, and -1 on basic geometric questions such  //
// as the orientation of point sets, in-circle, and in-sphere tests.  They   //
// are basic units for the composition of geometric algorithms.  TetGen uses //
// two 3D geometric predicates, which are the orientation test and the in-   //
// sphere test (e.g. the locally Deklaunay test).                            //
//                                                                           //
// Orientation test:  let a, b, c be a sequence of 3 non-collinear points in //
// R^3.  They defines a unique hypeplane H.  Let H+ and H- be the two spaces //
// separated by H, which are defined as follows (using the left-hand rule):  //
// make a fist using your left hand in such a way that your fingers follow   //
// the order of a, b and c, then your thumb is pointing to H+.  Given any    //
// point d in R^2, the orientation test returns +1 if d lies in H+, -1 if d  //
// lies in H-, or 0 if d lies on H.                                          //
//                                                                           //
// In-sphere test:  let a, b, c, d be 4 non-coplanar points in R^3.  They    //
// defines a unique circumsphere S.  Given any point e in R^3, the in-sphere //
// test returns +1 if e lies inside S, or -1 if e lies outside S, or 0 if e  //
// lies on S.                                                                //
//                                                                           //
// The correctness of geometric predicates is crucial for the control flow   //
// and hence for the correctness and robustness of an implementation of a    //
// geometric algorithm.  The following routines use arbitrary precision      //
// floating-point arithmetic. They are fast and robust. It is provided by J. //
// Schewchuk in public domain (http://www.cs.cmu.edu/~quake/robust.html).    //
// The source code are found in a separate file "predicates.cxx".            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL exactinit();
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenmesh data type                                                  //
//                                                                           //
// Includes data types and mesh routines for creating tetrahedral meshes and //
// Delaunay tetrahedralizations, mesh input & output, and so on.             //
//                                                                           //
// An object of tetgenmesh can be used to store a triangular or tetrahedral  //
// mesh and its settings. TetGen's functions operates on one mesh each time. //
// This type allows reusing of the same function for different meshes.       //
//                                                                           //
// The mesh data structure (tetrahedron-based and triangle-edge data struct- //
// ures) are declared. There are other accessary data type defined as well,  //
// for efficient memory management and link list operations, etc.            //
//                                                                           //
// All algorithms TetGen used are implemented in this data type as member    //
// functions. References of these algorithms can be found in user's manual.  //
//                                                                           //
// It's not necessary to understand this type. There is a global function    //
// "tetrahedralize()" (defined at the end of this file) implicitly creates   //
// the object and calls its member functions according to the command line   //
// switches you specified.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenmesh {

  public:

    // Maximum number of characters in a file name (including the null).
    enum {FILENAMESIZE = 1024};

    // For efficiency, a variety of data structures are allocated in bulk.
    //   The following constants determine how many of each structure is
    //   allocated at once.
    enum {VERPERBLOCK = 4092, SUBPERBLOCK = 4092, ELEPERBLOCK = 8188};

    // Used for the point location scheme of Mucke, Saias, and Zhu, to
    //   decide how large a random sample of tetrahedra to inspect.
    enum {SAMPLEFACTOR = 11};

    // Labels that signify two edge rings of a triangle defined in Muecke's
    //   triangle-edge data structure, one (CCW) traversing edges in count-
    //   erclockwise direction and one (CW) in clockwise direction.
    enum {CCW = 0, CW = 1};

    // Labels that signify whether a record consists primarily of pointers
    //   or of floating-point words.  Used to make decisions about data
    //   alignment.
    enum wordtype {POINTER, FLOATINGPOINT};

    // Labels that signify the type of a vertex. An UNUSEDVERTEX is a vertex
    //   read from input (.node file or tetgenio structure) or an isolated
    //   vertex (outside the mesh).  It is the default type for a newpoint.
    enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, NONACUTEVERTEX,
                   ACUTEVERTEX, FREESEGVERTEX, FACETVERTEX, FREESUBVERTEX,
                   FREEVOLVERTEX, DEADVERTEX = -32768};
 
    // Labels that signify the type of a subsegment or a subface.  An input
    //   (sub)segment may have type NONSHARPSEGMENT, SHARPSEGMENT. A subface
    //   is SKINNYSUB is it has two segments form a small angle (less than
    //   10 degree), such faces are unrefineable.
    enum shestype {NONSHARPSEGMENT, SHARPSEGMENT, NONSKINNYSUB, SKINNYSUB};

    // Labels that signify the type of a degenerate face or tetrahedron.
    //   A degenerate face has its three vertices collinear, DEGENEEDLE;
    //   a degenerate tetrahedron may have it three vertices collinear,
    //   DEGETRI, or four vertices collinear, DEGENEEDLE, or four vertices
    //   coplanar, DEGEKITE, or DEGECAP.
    enum degetype {NONDEGENERATE, DEGENEEDLE, DEGETRI, DEGEKITE, DEGECAP};

    // Labels that signify the type of flips can be applied on a face.
    //   A flipable face has the one of the types T23, T32, T22, and T44.
    //   Types UNFLIPABLE, NONCONVEX are unflipable.
    enum fliptype {T23, T32, T22, T44, UNFLIPABLE, FORBIDDENFACE,
                   FORBIDDENEDGE, NONCONVEX};

    // Labels that signify the result of triangle-triangle intersection.
    //   The result indicates that two triangles t1 and t2 are completely
    //   DISJOINT, or adjoint only at a vertex SHAREVERTEX, or adjoint at
    //   an edge SHAREEDGE, or coincident SHAREFACE, or INTERSECT.
    enum intersectresult {DISJOINT, SHAREVERTEX, SHAREEDGE, SHAREFACE,
                          INTERSECT};

    // Labels that signify the result of point location.  The result of a
    //   search indicates that the point falls inside a tetrahedron, inside
    //   a triangle, on an edge, on a vertex, or outside the mesh. 
    enum locateresult {INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX, OUTSIDE};

    // Labels that signify the result of vertex insertion.  The result
    //   indicates that the vertex was inserted with complete success, was
    //   inserted but encroaches upon a subsegment, was not inserted because
    //   it lies on a segment, or was not inserted because another vertex
    //   occupies the same location.
    enum insertsiteresult {SUCCESSINTET, SUCCESSONFACE, SUCCESSONEDGE,
                           ENCROACHINGPOINT, DUPLICATEPOINT, OUTSIDEPOINT};

    // Labels that signify the result of direction finding.  The result
    //   indicates that a segment connecting the two query points accross
    //   an edge of the direction triangle/tetrahedron, across a face of
    //   the direction tetrahedron, along the left edge of the direction
    //   triangle/tetrahedron, along the right edge of the direction
    //   triangle/tetrahedron, or along the top edge of the tetrahedron.
    enum finddirectionresult {ACROSSEDGE, ACROSSFACE, LEFTCOLLINEAR,
                              RIGHTCOLLINEAR, TOPCOLLINEAR, BELOWHULL};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The basic mesh element data structures                                    //
//                                                                           //
// There are four types of mesh elements: tetrahedra, subfaces, subsegments, //
// and points,  where subfaces and subsegments are triangles and edges which //
// appear on boundaries.  A tetrahedralization of a 3D point set comprises   //
// tetrahedra and points;  a surface mesh of a 3D domain comprises subfaces  //
// (triangles), subsegments and points;  and it is the elements of all the   //
// four types consist of a tetrahedral mesh of a 3D domain.  However, TetGen //
// uses three data types: 'tetrahedron', 'shellface', and 'point' to repres- //
// ent the basic mesh elements.  A 'tetrahedron' is a tetrahedron; while a   //
// 'shellface' can represent either a subface or a subsegment; and a 'point' //
// represent a point.  Theese three data types, linked by pointers comprise  //
// a tetrahedralization or a mesh.                                           //
//                                                                           //
// The data type 'tetrahedron' primarily consists of a list of four pointers //
// to its corners, a list of four pointers to its adjoining tetrahedra, a    //
// list of four pointers to its adjoining subfaces(when subfaces are needed).//
// Optinoally, (depending on the selected switches), it may contain an arbi- //
// trary number of user-defined floating-point attributes,  an optional max- //
// imum volume constraint (-a switch), and a pointer to a list of high-order //
// nodes (-o2 switch). Because the size of a tetrahedron is not determined   //
// until running time, it is not simply declared as a structure.             //
//                                                                           //
// For purpose of storing geometric information, it is important to know the //
// ordering of the vertices of a tetrahedron.  Let v0, v1, v2, and v3 be the //
// four nodes corresponding to the order of their storage in a tetrahedron.  //
// v3 always has negative orientation with respect to v0, v1, v2 (in other   //
// words, v3 lies above the oriented plane passes through v0, v1, v2).  Let  //
// the four faces of the tetrahedron be f0, f1, f2, and f3. Vertices of each //
// face are stipulated as follows: f0 (v0, v1, v2), f1 (v0, v3, v1), f2 (v1, //
// v3, v2), f3 (v2, v3, v0).  Adjoining tetrahedra as well as subfaces are   //
// stored in the order of its faces, e.g., the first adjoining tetrahedra is //
// the neighbor at f0, and so on.                                            //
//                                                                           //
// A subface is represented by the data type 'shellface'.  It has three      //
// pointers to its vertices, three pointers to its adjoining subfaces, three //
// pointers to subsegments, two pointers to its adjoining tetrahedra, and a  //
// boundary marker (an integer). Furthermore, the pointers to vertices,      //
// adjoining subfaces, and subsegments are ordered in a way that indicates   //
// their geometric relation.  Let the three vertices according to the order  //
// of their storage be v0, v1 and v2, respectively, and e0, e1 and e2 be the //
// three edges, then we have: e0 (v0, v1), e1 (v1, v2), e2 (v2, v0). Adjoin- //
// ing subfaces and subsegments are stored in the order of its edges.        //
//                                                                           //
// A subsegment is also represented by a 'shellface'. It has exactly the     //
// same data fields as a subface has, but only uses some of them. It has two //
// pointers to its endpoints, two pointers to its adjoining (and collinear)  //
// subsegments, one pointer to a subface containing it (there may exist any  //
// number of subfaces having it, choose one of them arbitrarily). The geome- //
// tric relation between its endpoints and adjoining (collinear) subsegments //
// is kept with respect to the storing order of its endpoints. The adjoining //
// subsegment at the first endpoint is saved ahead of the other.             //
//                                                                           //
// The data type 'point' is relatively simple. A point is a list of floating //
// -point numbers, starting with the x, y, and z coordinates, followed by an //
// arbitrary number of optional user-defined floating-point attributes, an   //
// integer boundary marker, an integer for the point type, and a pointer to  //
// a tetrahedron. The latter is used for speeding up point location during   //
// the mesh generation.                                                      //
//                                                                           //
// For a tetrahedron on a boundary (or a hull) of the mesh, some or all of   //
// the adjoining tetrahedra may not be present. For an interior tetrahedron, //
// often no neighboring subfaces are present,  Such absent tetrahedra and    //
// subfaces are never represented by the NULL pointers; they are represented //
// by two special records: `dummytet', the tetrahedron fills "outer space",  //
// and `dummysh',  the vacuous subfaces which are omnipresent.               //
//                                                                           //
// Tetrahedra and adjoining subfaces are glued together through the pointers //
// saved in each data fields of them. Subfaces and adjoining subsegments are //
// connected in the same fashion.  However, there are no pointers directly   //
// gluing tetrahedra and adjoining subsegments.  For the purpose of saving   //
// space, the connections between tetrahedra and subsegments are entirely    //
// mediated through subfaces.  The following part is an explaination of how  //
// subfaces are connected in TetGen.                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The subface-subface and subface-subsegment connections                    //
//                                                                           //
// Adjoining subfaces sharing a common edge are connected in such a way that //
// they form a face ring around the edge. It is in deed a single linked list //
// which is cyclic, e.g., one can start from any subface in it and traverse  //
// back. When the edge is not a subsegment, the ring only has two coplanar   //
// subfaces which are pointing to each other. Otherwise, the face ring may   //
// have any number of subfaces (and are not all coplanar).                   //
//                                                                           //
// The face ring around a subsegment is formed as follows.  Let s be a sub-  //
// segment, and f be a subface containing s as an edge.  The direction of s  //
// is stipulated from its first endpoint to its second (the first and second //
// endpoints are according to their storage in s).  When the direction of s  //
// is determined, other two edges of f are oriented following this direction.//
// The "directional normal" of f is a ray starts from any point in f, points //
// to the direction of the cross product of any of two edge vectors of f.    //
//                                                                           //
// The face ring of s is a cyclic ordered set of subfaces containing s, i.e.,//
// F(s) = {f1, f2, ..., fn}, n >= 1.  Where the order is defined as follows: //
// let fi, fj be two faces in F(s), the "normal-angle", nangle(i,j) (range   //
// from 0 to 360 degree) is the angle between the directional normals of fi  //
// and fj;  that fi is in front of fj (or symbolically, fi < fj) if there    //
// exists another fk in F(s), and nangle(k, i) < nangle(k, j). The face ring //
// of s can be represented as: f1 < f2 < ... < fn < f1.                      //
//                                                                           //
// The easiest way to imagine how a face ring is formed is to use the right- //
// hand rule.  Make a fist using your right hand with the thumb pointing to  //
// the direction of the subsegment. The face ring is connected following the //
// direction of your fingers.                                                //
//                                                                           //
// The subface and subsegment are also connected through pointers stored in  //
// their own data fields.  Every subface has a pointer ti its adjoining sub- //
// segment. However, a subsegment only has one pointer to a subface which is //
// containing it. Such subface can be choosn arbitrarily, other subfaces are //
// found through the face ring.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // The tetrahedron data structure.  Fields of a tetrahedron contains:
    //   - a list of four adjoining tetrahedra;
    //   - a list of four vertices;
    //   - a list of four subfaces (optional, used for -p switch);
    //   - a list of user-defined floating-point attributes (optional);
    //   - a volume constraint (optional, used for -a switch);
    //   - a pointer to a list of high-ordered nodes (optional, -o2 switch);

    typedef REAL **tetrahedron;

    // The shellface data structure.  Fields of a shellface contains:
    //   - a list of three adjoining subfaces;
    //   - a list of three vertices;
    //   - a list of two adjoining tetrahedra;
    //   - a list of three adjoining subsegments;
    //   - a pointer to a badface containing it (optional, used for -q);
    //   - an area constraint (optional, used for -q);
    //   - an integer for boundary marker;
    //   - an integer for type: SHARPSEGMENT, NONSHARPSEGMENT, ...;

    typedef REAL **shellface;

    // The point data structure.  It is actually an array of REALs:
    //   - x, y and z coordinates;
    //   - a list of user-defined point attributes (optional);
    //   - a pointer to a simplex (tet, tri, edge, or vertex);
    //   - a pointer to a parent point (optional, used for -q);
    //   - an integer for boundary marker;
    //   - an integer for verttype: INPUTVERTEX, FREEVERTEX, ...;

    typedef REAL *point;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The mesh handle (triface, face) data types                                //
//                                                                           //
// Two special data types, 'triface' and 'face' are defined for maintaining  //
// and updating meshes. They are like pointers (or handles), which allow you //
// to hold one particular part of the mesh, i.e., a tetrahedron, a triangle, //
// an edge and a vertex.  However, these data types do not themselves store  //
// any part of the mesh. The mesh is made of the data types defined above.   //
//                                                                           //
// Muecke's "triangle-edge" data structure is the prototype for these data   //
// types.  It allows a universal representation for every tetrahedron,       //
// triangle, edge and vertex.  For understanding the following descriptions  //
// of these handle data structures,  readers are required to read both the   //
// introduction and implementation detail of "triangle-edge" data structure  //
// in Muecke's thesis.                                                       //
//                                                                           //
// A 'triface' represents a face of a tetrahedron and an oriented edge of    //
// the face simultaneously.  It has a pointer 'tet' to a tetrahedron, an     //
// integer 'loc' (range from 0 to 3) as the face index, and an integer 'ver' //
// (range from 0 to 5) as the edge version. A face of the tetrahedron can be //
// uniquly determined by the pair (tet, loc), and an oriented edge of this   //
// face can be uniquly determined by the triple (tet, loc, ver).  Therefore, //
// different usages of one triface are possible.  If we only use the pair    //
// (tet, loc), it refers to a face, and if we add the 'ver' additionally to  //
// the pair, it is an oriented edge of this face.                            //
//                                                                           //
// A 'face' represents a subface and an oriented edge of it simultaneously.  //
// It has a pointer 'sh' to a subface, an integer 'shver'(range from 0 to 5) //
// as the edge version.  The pair (sh, shver) determines a unique oriented   //
// edge of this subface.  A 'face' is also used to represent a subsegment,   //
// in this case, 'sh' points to the subsegment, and 'shver' indicates the    //
// one of two orientations of this subsegment, hence, it only can be 0 or 1. //
//                                                                           //
// Mesh navigation and updating are accomplished through a set of mesh       //
// manipulation primitives which operate on trifaces and faces.  They are    //
// introduced below.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class triface {

      public:

        tetrahedron* tet;
        int loc, ver;

        // Constructors;
        triface() : tet(0), loc(0), ver(0) {}
        // Operators;
        triface& operator=(const triface& t) {
          tet = t.tet; loc = t.loc; ver = t.ver;
          return *this;
        }
        bool operator==(triface& t) {
          return tet == t.tet && loc == t.loc && ver == t.ver;
        }
        bool operator!=(triface& t) {
          return tet != t.tet || loc != t.loc || ver != t.ver;
        }
    };

    class face {

      public:

        shellface *sh;
        int shver;

        // Constructors;
        face() : sh(0), shver(0) {}
        // Operators;
        face& operator=(const face& s) {
          sh = s.sh; shver = s.shver;
          return *this;
        }
        bool operator==(face& s) {return (sh == s.sh) && (shver == s.shver);}
        bool operator!=(face& s) {return (sh != s.sh) || (shver != s.shver);}
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The badface structure                                                     //
//                                                                           //
// A multiple usages structure. Despite of its name, a 'badface' can be used //
// to represent the following objects:                                       //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;             //
//   - an encroached subsegment or subface;                                  //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;      //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;    //
//   - a degenerate tetrahedron (see routine checkdegetet()).                //
//   - a recently flipped face (saved for undoing the flip later).           //
//                                                                           //
// It has the following fields:  'tt' holds a tetrahedron; 'ss' holds a sub- //
// segment or subface; 'cent' is the circumcent of 'tt' or 'ss', 'key' is a  //
// special value depending on the use, it can be either the square of the    //
// radius-edge ratio of 'tt' or the flipped type of 'tt';  'forg', 'fdest',  //
// 'fapex', and 'foppo' are vertices saved for checking the object in 'tt'   //
// or 'ss' is still the same when it was stored; 'noppo' is the fifth vertex //
// of a degenerate point set.  'previtem' and 'nextitem' implement a double  //
// link for managing many basfaces.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    struct badface {
      triface tt; 
      face ss; 
      REAL key;
      REAL cent[3];
      point forg, fdest, fapex, foppo;
      point noppo;
      struct badface *previtem, *nextitem; 
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The list, link and queue data structures                                  //
//                                                                           //
// These data types are used to manipulate a set of (same-typed) data items. //
// For a given set S = {a, b, c, ...}, a list stores the elements of S in a  //
// piece of continuous memory. It allows quickly accessing each element of S,//
// thus is suitable for storing a fix-sized set.  While a link stores its    //
// elements incontinuously. It allows quickly inserting or deleting an item, //
// thus is suitable for storing a size-variable set.  A queue is basically a //
// special case of a link where one data element joins the link at the end   //
// and leaves in an ordered fashion at the other end.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // The compfunc data type.  "compfunc" is a pointer to a linear-order
    //   function, which takes two 'void*' arguments and returning an 'int'. 
    //   
    // A function: int cmp(const T &, const T &),  is said to realize a
    //   linear order on the type T if there is a linear order <= on T such
    //   that for all x and y in T satisfy the following relation:
    //                 -1  if x < y.
    //   comp(x, y) =   0  if x is equivalent to y.
    //                 +1  if x > y.
    typedef int (*compfunc) (const void *, const void *);

    // The predefined compare functions for primitive data types.  They
    //   take two pointers of the corresponding date type, perform the
    //   comparation, and return -1, 0 or 1 indicating the default linear
    //   order of them.

    // Compare two 'integers'.
    static int compare_2_ints(const void* x, const void* y);
    // Compare two 'longs'. 
    static int compare_2_longs(const void* x, const void* y);
    // Compare two 'unsigned longs'. 
    static int compare_2_unsignedlongs(const void* x, const void* y);

    // The function used to determine the size of primitive data types and
    //   set the corresponding predefined linear order functions for them.
    static void set_compfunc(char* str, int* itembytes, compfunc* pcomp);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// List data structure.                                                      //
//                                                                           //
// A 'list' is an array of items with automatically reallocation of memory.  //
// It behaves like an array.                                                 //
//                                                                           //
// 'base' is the starting address of the array;  The memory unit in list is  //
//   byte, i.e., sizeof(char). 'itembytes' is the size of each item in byte, //
//   so that the next item in list will be found at the next 'itembytes'     //
//   counted from the current position.                                      //
//                                                                           //
// 'items' is the number of items stored in list.  'maxitems' indicates how  //
//   many items can be stored in this list. 'expandsize' is the increasing   //
//   size (items) when the list is full.                                     //
//                                                                           //
// 'comp' is a pointer pointing to a linear order function for the list.     //
//   default it is set to 'NULL'.                                            //
//                                                                           //
// The index of list always starts from zero, i.e., for a list L contains    //
//   n elements, the first element is L[0], and the last element is L[n-1].  //
//   This feature lets lists likes C/C++ arrays.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class list {

      public:

        char *base;
        int  itembytes;
        int  items, maxitems, expandsize;
        compfunc comp;

      public:

        list(int itbytes, compfunc pcomp, int mitems = 256, int exsize = 128) {
          listinit(itbytes, pcomp, mitems, exsize);
        }
        list(char* str, int mitems = 256, int exsize = 128) {
          set_compfunc(str, &itembytes, &comp);
          listinit(itembytes, comp, mitems, exsize);
        }
        ~list() { free(base); }

        void *operator[](int i) { return (void *) (base + i * itembytes); }

        void listinit(int itbytes, compfunc pcomp, int mitems, int exsize);
        void setcomp(compfunc compf) { comp = compf; }    
        void clear() { items = 0; }
        int  len() { return items; }
        void *append(void* appitem);
        void *insert(int pos, void* insitem);
        void del(int pos);
        int  hasitem(void* checkitem);
        int  remove(void* remitem);
        void sort();
    }; 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Memorypool data structure.                                                //
//                                                                           //
// A type used to allocate memory.  (It is incorporated from Shewchuk's      //
// Triangle program)                                                         //
//                                                                           //
// firstblock is the first block of items. nowblock is the block from which  //
//   items are currently being allocated. nextitem points to the next slab   //
//   of free memory for an item. deaditemstack is the head of a linked list  //
//   (stack) of deallocated items that can be recycled.  unallocateditems is //
//   the number of items that remain to be allocated from nowblock.          //
//                                                                           //
// Traversal is the process of walking through the entire list of items, and //
//   is separate from allocation.  Note that a traversal will visit items on //
//   the "deaditemstack" stack as well as live items.  pathblock points to   //
//   the block currently being traversed.  pathitem points to the next item  //
//   to be traversed.  pathitemsleft is the number of items that remain to   //
//   be traversed in pathblock.                                              //
//                                                                           //
// itemwordtype is set to POINTER or FLOATINGPOINT, and is used to suggest   //
//   what sort of word the record is primarily made up of.  alignbytes       //
//   determines how new records should be aligned in memory.  itembytes and  //
//   itemwords are the length of a record in bytes (after rounding up) and   //
//   words.  itemsperblock is the number of items allocated at once in a     //
//   single block.  items is the number of currently allocated items.        //
//   maxitems is the maximum number of items that have been allocated at     //
//   once; it is the current number of items plus the number of records kept //
//   on deaditemstack.                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class memorypool {

      public:

        void **firstblock, **nowblock;
        void *nextitem;
        void *deaditemstack;
        void **pathblock;
        void *pathitem;
        wordtype itemwordtype;
        int  alignbytes;
        int  itembytes, itemwords;
        int  itemsperblock;
        long items, maxitems;
        int  unallocateditems;
        int  pathitemsleft;

      public:

        memorypool();
        memorypool(int, int, enum wordtype, int);
        ~memorypool();
    
        void poolinit(int, int, enum wordtype, int);
        void restart();
        void *alloc();
        void dealloc(void*);
        void traversalinit();
        void *traverse();
    };  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Link data structure.                                                      //
//                                                                           //
// A 'link' is a double linked nodes. It uses the memorypool data structure  //
// for memory management.  Following is an image of a link.                  //
//                                                                           //
//   head-> ____0____      ____1____      ____2____      _________<-tail     //
//         |__next___|--> |__next___|--> |__next___|--> |__NULL___|          //
//         |__NULL___|<-- |__prev___|<-- |__prev___|<-- |__prev___|          //
//         |         |    |_       _|    |_       _|    |         |          //
//         |         |    |_ Data1 _|    |_ Data2 _|    |         |          //
//         |_________|    |_________|    |_________|    |_________|          //
//                                                                           //
// The unit size for storage is size of pointer, which may be 4-byte (in 32- //
//   bit machine) or 8-byte (in 64-bit machine). The real size of an item is //
//   stored in 'linkitembytes'.                                              //
//                                                                           //
// 'head' and 'tail' are pointers pointing to the first and last nodes. They //
//   do not conatin data (See above).                                        //
//                                                                           //
// 'nextlinkitem' is a pointer pointing to a node which is the next one will //
//   be traversed. 'curpos' remembers the position (1-based) of the current  //
//   traversing node.                                                        //
//                                                                           //
// 'linkitems' indicates how many items in link. Note it is different with   //
//   'items' of memorypool.                                                  //
//                                                                           //
// The index of link starts from 1, i.e., for a link K contains n elements,  //
//   the first element of the link is K[1], and the last element is K[n].    //
//   See the above figure.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class link : public memorypool {

      public:

        void **head, **tail;
        void *nextlinkitem;
        int  linkitembytes;
        int  linkitems;
        int  curpos;
        compfunc comp;

      public:

        link(int _itembytes, compfunc _comp, int itemcount) {
          linkinit(_itembytes, _comp, itemcount);
        }
        link(char* str, int itemcount) {
          set_compfunc(str, &linkitembytes, &comp);
          linkinit(linkitembytes, comp, itemcount);
        }

        void linkinit(int _itembytes, compfunc _comp, int itemcount);
        void setcomp(compfunc compf) { comp = compf; }
        void rewind() { nextlinkitem = *head; curpos = 1; }
        void goend() { nextlinkitem = *(tail + 1); curpos = linkitems; }    
        long len() { return linkitems; }
        void clear();
        bool move(int numberofnodes);
        bool locate(int pos);
        void *add(void* newitem);
        void *insert(int pos, void* insitem);
        void *del(void* delitem);
        void *del(int pos);
        void *getitem();
        void *getnitem(int pos);
        int  hasitem(void* checkitem);
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Queue data structure.                                                     //
//                                                                           //
// A 'queue' is a basically a link.  Following is an image of a queue.       //
//              ___________     ___________     ___________                  //
//   Pop() <-- |_         _|<--|_         _|<--|_         _| <-- Push()      //
//             |_  Data0  _|   |_  Data1  _|   |_  Data2  _|                 //
//             |___________|   |___________|   |___________|                 //
//              queue head                       queue tail                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class queue : public link {

      public:

        queue(int bytes, int count = 256) : link(bytes, NULL, count) {}
        queue(char* str, int count = 256) : link(str, count) {}

        int  empty() { return linkitems == 0; }
        void *push(void* newitem) { return link::add(newitem); }
        void *bot() { return link::getnitem(1); }
        void *pop() { return link::del(1); }
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Following are variables used in 'tetgenmesh' for miscellaneous purposes.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // Pointer to an object of 'tetgenio', which contains input data.
    tetgenio *in;

    // Pointer to an object of 'tetgenbehavor', which contains the user-
    //   defined command line swithes and filenames.
    tetgenbehavior *b;

    // Variables used to allocate and access memory for tetrahedra, subfaces
    //   subsegments, points, encroached subfaces, encroached subsegments,
    //   bad-quality tetrahedra, and so on.
    memorypool *tetrahedrons;
    memorypool *subfaces;
    memorypool *subsegs;
    memorypool *points;   
    memorypool *badsubsegs;
    memorypool *badsubfaces;
    memorypool *badtetrahedrons;
    memorypool *flipstackers;

    // Pointer to a recently visited tetrahedron. Improves point location
    //   if proximate points are inserted sequentially.
    triface recenttet;

    // Pointer to the 'tetrahedron' that occupies all of "outer space".
    tetrahedron *dummytet;
    tetrahedron *dummytetbase; // Keep base address so we can free it later.

    // Pointer to the omnipresent subface.  Referenced by any tetrahedron,
    //   or subface that isn't connected to a subface at that location.
    shellface *dummysh;
    shellface *dummyshbase;    // Keep base address so we can free it later.

    // List of lifting points of facets used for surface triangulation.
    REAL *liftpointarray;

    // List used for Delaunay refinement algorithm. 
    list *qualchecktetlist;

    // Queues that maintain the bad (badly-shaped or too large) tetrahedra.
    //   The tails are pointers to the pointers that have to be filled in to
    //   enqueue an item.  The queues are ordered from 63 (highest priority)
    //   to 0 (lowest priority).
    badface *subquefront[2], **subquetail[2];
    badface *tetquefront[64], **tetquetail[64];

    // Array (size = numberoftetrahedra * 6) for storing high-order nodes of
    //   tetrahedra (only used when -o2 switch is selected).
    point *highordertable;

    REAL xmax, xmin, ymax, ymin, zmax, zmin;      // Bounding box of points.
    REAL longest;                       // The longest possible edge length.
    long hullsize;                        // Number of faces of convex hull.
    long insegment;                             // Number of input segments.
    int steinerleft;               // Number of Steiner points not yet used.
    int pointmarkindex;         // Index to find boundary marker of a point.
    int point2simindex;      // Index to find a simplex adjacent to a point.
    int highorderindex; // Index to find extra nodes for highorder elements.
    int elemattribindex;       // Index to find attributes of a tetrahedron.
    int volumeboundindex;    // Index to find volume bound of a tetrahedron.
    int shmarkindex;          // Index to find boundary marker of a subface.
    int areaboundindex;            // Index to find area bound of a subface.
    int checksubfaces;                // Are there subfaces in the mesh yet?
    int checkquality;                // Has quality triangulation begun yet?
    int nonconvex;                            // Is current mesh non-convex?
    int dupverts;                          // Are there duplicated vertices?
    int unuverts;                              // Are there unused vertices?
    long samples;            // Number of random samples for point location.
    unsigned long randomseed;                 // Current random number seed.
    REAL macheps;                                    // The machine epsilon.
    long flip23s, flip32s, flip22s, flip44s;   // Number of flips performed.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Fast lookup tables for mesh manipulation primitives.                      //
//                                                                           //
// Mesh manipulation primitives (given below) are basic operations on mesh   //
// data structures. They answer basic queries on mesh handles, such as "what //
// is the origin (or destination, or apex) of the face?", "what is the next  //
// (or previous) edge in the edge ring?", and "what is the next face in the  //
// face ring?", and so on.                                                   //
//                                                                           //
// The implementation of basic queries can take advangtage of the fact that  //
// the mesh data structures additionally store geometric informations.  For  //
// example, we have ordered the four vertices (from 0 to 3) and four faces   //
// (from 0 to 3) of a tetrahedron,  and for each face of the tetrahedron, a  //
// sequence of vertices has stipulated,  therefore the origin of any face of //
// the tetrahedron can be quickly determined by a table 'locver2org', which  //
// takes the index of the face and the edge version as inputs.  A list of    //
// fast lookup tables are defined below. They're just like global variables. //
// All tables are initialized once at the runtime and used by all objects of //
// tetgenmesh.                                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // For enext() primitive, uses 'ver' as the index. 
    static int ve[6];

    // For org(), dest() and apex() primitives, uses 'ver' as the index.
    static int vo[6], vd[6], va[6];

    // For org(), dest() and apex() primitives, uses 'loc' as the first
    //   index and 'ver' as the second index.
    static int locver2org[4][6];
    static int locver2dest[4][6];
    static int locver2apex[4][6];

    // For oppo() primitives, uses 'loc' as the index.
    static int loc2oppo[4];

    // For fnext() primitives, uses 'loc' as the first index and 'ver' as
    //   the second index,  returns an array containing a new 'loc' and a
    //   new 'ver'. Note: Only valid for 'ver' equals one of {0, 2, 4}.
    static int locver2nextf[4][6][2];

    // For enumerating three edges of a triangle.
    static int plus1mod3[3];
    static int minus1mod3[3];

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh manipulation primitives                                              //
//                                                                           //
// A serial of mesh operations such as topological maintenance,  navigation, //
// local modification, etc.,  is accomplished through a set of mesh manipul- //
// ation primitives. These primitives are indeed very simple functions which //
// take one or two handles ('triface's and 'face's) as parameters,  perform  //
// basic operations such as "glue two tetrahedra at a face",  "return the    //
// origin of a tetrahedron", "return the subface adjoining at the face of a  //
// tetrahedron", and so on.                                                  //
//                                                                           //
// In the following, symbols t, t1, and t2 denote handles of type 'triface', //
// i.e., t is a face of a tetrahedron. Likewise, handles of type 'face' are  //
// denoted by s, s1, s2; e denotes an oriented edge, and v denotes a vertex. //
//                                                                           //
// The basic primitives for tetrahedra are:                                  //
//                                                                           //
//   sym(t1, t2)      t1 and t2 refer to the same face but point to two      //
//                    different tetrahedra respectively.                     //
//   bond(t1, t2)     Bonds t1 and t2 together.  t1 and t2 should refer to   //
//                    the same face.                                         //
//   dissolve(t)      Detaches the adjoining tetrahedron from t. t bonds to  //
//                    'dummytet' after this operation.                       //
//                                                                           //
//   v = org(t)       v is the origin of t.                                  //
//   v = dest(t)      v is the destination of t.                             //
//   v = apex(t)      v is the apex of t.                                    //
//   v = oppo(t)      v is the opposite of t.                                //
//                                                                           //
//   esym(t1, t2)     t2 is the inversed edge of t1, i.e., t1 and t2 are two //
//                    directed edges of the same undirected edge.            //
//   enext(t1, t2)    t2 is the successor of t1 in the edge ring.            //
//   enext2(t1, t2)   t2 is the precessor of t1 in the edge ring.            //
//                                                                           //
//   fnext(t1, t2)    t2 is the successor of t1 in the face ring.            //
//                                                                           //
// The basic primitives for subfaces (as well as subsegments) are:           //
//                                                                           //
//   spivot(s1, s2)   s1 and s2 refer to the same edge but point to two      //
//                    different subfaces respectively.                       //
//   sbond(s1, s2)    Bonds s1 and s2 together (at an edge).                 //
//   sbond1(s1, s2)   Only bonds s2 to s1 (but not s1 to s2).  It is used    //
//                    for creating the face ring.                            //
//   sdissolve(s)     Detaches the adjoining subface from s. s bonds to      //
//                    'dummysh' after this operation.                        //
//                                                                           //
//   v = sorg(s)      v is the origin of s.                                  //
//   v = sdest(s)     v is the destination of s.                             //
//   v = sapex(s)     v is the apex of s.                                    //
//                                                                           //
//   sesym(s1, s2)    s2 is the inversed edge of s1.                         //
//   senext(s1, s2)   s2 is the successor of s1 in the edge ring.            //
//   senext2(s1, s2)  s2 is the precessor of s1 in the edge ring..           //
//                                                                           //
// For interacting tetrahedra and subfaces:                                  //
//                                                                           //
//   tspivot(t, s)    Returns the adjoining subface of t in s. s may hold    //
//                    'dummysh' when t is an internal face.                  //
//   stpivot(s, t)    Returns the adjoining tetrahedron of s in t. t may be  //
//                    'dummytet'.                                            //
//   tsbond(t, s)     Bond t and s together.  t and s must represent the     //
//                    same face.                                             //
//   tsdissolve(t)    Detaches the adjoining subface from t.                 //
//   stdissolve(s)    Detaches the adjoining tetrahedron from s.             //
//                                                                           //
// For interacting subfaces and subsegments:                                 //
//                                                                           //
//   sspivot(s, e)    Returns the adjoining subsegment of s in e.            //
//   ssbond(s, e)     Bond s and e together.  s and e must represent the     //
//                    same edge.                                             //
//   ssdissolve(s)    Detaches the adjoining subsegment from s.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // Primitives for tetrahedra.
    inline void decode(tetrahedron ptr, triface& t);
    inline tetrahedron encode(triface& t);
    inline void sym(triface& t1, triface& t2);
    inline void symself(triface& t);
    inline void bond(triface& t1, triface& t2);
    inline void dissolve(triface& t);
    inline point org(triface& t);
    inline point dest(triface& t);
    inline point apex(triface& t);
    inline point oppo(triface& t);
    inline void setorg(triface& t, point pointptr);
    inline void setdest(triface& t, point pointptr);
    inline void setapex(triface& t, point pointptr);
    inline void setoppo(triface& t, point pointptr);
    inline void esym(triface& t1, triface& t2);
    inline void esymself(triface& t);
    inline void enext(triface& t1, triface& t2);
    inline void enextself(triface& t);
    inline void enext2(triface& t1, triface& t2);
    inline void enext2self(triface& t);
    inline bool fnext(triface& t1, triface& t2);
    inline bool fnextself(triface& t);
    inline void enextfnext(triface& t1, triface& t2);
    inline void enextfnextself(triface& t);
    inline void enext2fnext(triface& t1, triface& t2);
    inline void enext2fnextself(triface& t);
    inline void infect(triface& t);
    inline void uninfect(triface& t);
    inline bool infected(triface& t);
    inline REAL elemattribute(tetrahedron* ptr, int attnum);
    inline void setelemattribute(tetrahedron* ptr, int attnum, REAL value);
    inline REAL volumebound(tetrahedron* ptr);
    inline void setvolumebound(tetrahedron* ptr, REAL value);
 
    // Primitives for subfaces and subsegments.
    inline void sdecode(shellface sptr, face& s);
    inline shellface sencode(face& s);
    inline void spivot(face& s1, face& s2);
    inline void spivotself(face& s);
    inline void sbond(face& s1, face& s2);
    inline void sbond1(face& s1, face& s2);
    inline void sdissolve(face& s);
    inline point sorg(face& s);
    inline point sdest(face& s);
    inline point sapex(face& s);
    inline void setsorg(face& s, point pointptr);
    inline void setsdest(face& s, point pointptr);
    inline void setsapex(face& s, point pointptr);
    inline void sesym(face& s1, face& s2);
    inline void sesymself(face& s);
    inline void senext(face& s1, face& s2);
    inline void senextself(face& s);
    inline void senext2(face& s1, face& s2);
    inline void senext2self(face& s);
    inline void sfnext(face&, face&);
    inline void sfnextself(face&);
    inline badface* shell2badface(face& s);
    inline void setshell2badface(face& s, badface* value);
    inline REAL areabound(face& s);
    inline void setareabound(face& s, REAL value);
    inline int shellmark(face& s);
    inline void setshellmark(face& s, int value);
    inline enum shestype shelltype(face& s);
    inline void setshelltype(face& s, enum shestype value); 
    inline void sinfect(face& s);
    inline void suninfect(face& s);
    inline bool sinfected(face& s);

    // Primitives for interacting tetrahedra and subfaces.
    inline void tspivot(triface& t, face& s);
    inline void stpivot(face& s, triface& t);
    inline void tsbond(triface& t, face& s);
    inline void tsdissolve(triface& t);    
    inline void stdissolve(face& s);

    // Primitives for interacting subfaces and subsegs.
    inline void sspivot(face& s, face& edge);
    inline void ssbond(face& s, face& edge);
    inline void ssdissolve(face& s);

    // Primitives for points.
    inline int  pointmark(point pt);
    inline void setpointmark(point pt, int value);
    inline enum verttype pointtype(point pt);
    inline void setpointtype(point pt, enum verttype value);
    inline tetrahedron point2tet(point pt);
    inline void setpoint2tet(point pt, tetrahedron value);
    inline shellface point2sh(point pt);
    inline void setpoint2sh(point pt, shellface value);
    inline point point2pt(point pt);
    inline void setpoint2pt(point pt, point value);
    inline point point2ppt(point pt);
    inline void setpoint2ppt(point pt, point value);
    inline point getliftpoint(int facetmark);

    // Advanced primitives.
    inline void adjustedgering(triface& t, int direction);
    inline void adjustedgering(face& s, int direction);
    inline bool isdead(triface* t);
    inline bool isdead(face* s);
    inline bool isfacehaspoint(face* t, point testpoint);
    inline bool isfacehasedge(face* s, point tend1, point tend2);
    inline bool issymexist(triface* t);
    bool getnextface(triface*, triface*);
    void getnextsface(face*, face*);
    void tsspivot(triface*, face*);
    void sstpivot(face*, triface*);   
    bool findorg(triface* t, point dorg);
    bool findorg(face* s, point dorg);
    void findedge(triface* t, point eorg, point edest);
    void findedge(face* s, point eorg, point edest);
    void findface(triface *fface, point forg, point fdest, point fapex);
    void getonextseg(face* s, face* lseg);
    void getseghasorg(face* sseg, point dorg);
    point getsubsegfarorg(face* sseg);
    point getsubsegfardest(face* sseg);
    void printtet(triface*);
    void printsh(face*);    

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Triangle-triangle intersection test                                       //
//                                                                           //
// The triangle-triangle intersection test is implemented with exact arithm- //
// etic. It exactly tells whether or not two triangles in three dimensions   //
// intersect.  Before implementing this test myself,  I tried two C codes    //
// (implemented by Thomas Moeller and Philippe Guigue, respectively), which  //
// are all public available. However both of them failed frequently. Another //
// unconvenience is that both codes only tell whether or not two triangles   //
// are intersecting and not distinguish the cases whether they are exactly   //
// intersecting in interior or they share a vertex, or share an edge.  All   //
// theese cases are acceptable and should return not intersection in TetGen. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    enum intersectresult edge_vertex_collinear_inter(REAL*, REAL*, REAL*);
    enum intersectresult edge_edge_coplanar_inter(REAL*, REAL*, REAL*,
                                                  REAL*, REAL*);
    enum intersectresult triangle_vertex_coplanar_inter(REAL*, REAL*, REAL*,
                                                        REAL*, REAL*);
    enum intersectresult triangle_edge_coplanar_inter(REAL*, REAL*, REAL*,
                                                      REAL*, REAL*, REAL*);
    enum intersectresult triangle_edge_inter_tail(REAL*, REAL*, REAL*, REAL*,
                                                  REAL*, REAL, REAL);
    enum intersectresult triangle_edge_inter(REAL*, REAL*, REAL*, REAL*,
                                             REAL*);
    enum intersectresult triangle_triangle_inter(REAL*, REAL*, REAL*, REAL*,
                                                 REAL*, REAL*);

    // Degenerate cases tests
    bool iscollinear(REAL*, REAL*, REAL*, REAL epspp);
    bool iscoplanar(REAL*, REAL*, REAL*, REAL*, REAL vol6, REAL epspp);
    bool iscospheric(REAL*, REAL*, REAL*, REAL*, REAL*, REAL epspp);

    // Linear algebra functions
    inline REAL dot(REAL* v1, REAL* v2);
    inline void cross(REAL* v1, REAL* v2, REAL* n);
    void initm44(REAL a00, REAL a01, REAL a02, REAL a03,
                 REAL a10, REAL a11, REAL a12, REAL a13,
                 REAL a20, REAL a21, REAL a22, REAL a23,
                 REAL a30, REAL a31, REAL a32, REAL a33, REAL M[4][4]);
    void m4xm4(REAL m1[4][4], REAL m2[4][4]);
    void m4xv4(REAL v2[4], REAL m[4][4], REAL v1[4]);
    bool lu_decmp(REAL lu[3][3], int n, int* ps, REAL* d, int N);
    void lu_solve(REAL lu[3][3], int n, int* ps, REAL* b, int N);

    // Geometric quantities calculators.
    inline REAL distance(REAL* p1, REAL* p2);
    REAL shortdistance(REAL* p, REAL* e1, REAL* e2);
    REAL interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n);
    void projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj);
    void projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj);
    void facenormal(REAL* pa, REAL* pb, REAL* pc, REAL* n, REAL* nlen);
    void edgeorthonormal(REAL* e1, REAL* e2, REAL* op, REAL* n);
    REAL facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2);
    void tetalldihedral(point, point, point, point, REAL dihed[6]);
    bool circumsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
    void inscribedsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
    void rotatepoint(REAL* p, REAL rotangle, REAL* p1, REAL* p2);
    void spherelineint(REAL* p1, REAL* p2, REAL* C, REAL R, REAL p[7]);
    void linelineint(REAL *p1,REAL *p2, REAL *p3, REAL *p4, REAL p[7]);

    // Memory managment routines.
    void dummyinit(int, int);
    void initializepointpool();
    void initializetetshpools();
    void tetrahedrondealloc(tetrahedron*);
    tetrahedron *tetrahedrontraverse();
    void shellfacedealloc(memorypool*, shellface*);
    shellface *shellfacetraverse(memorypool*);
    void badfacedealloc(memorypool*, badface*);
    badface *badfacetraverse(memorypool*);
    void pointdealloc(point);
    point pointtraverse();
    void maketetrahedron(triface*);
    void makeshellface(memorypool*, face*);
    void makepoint(point*);

    // Mesh items searching routines.
    void makepoint2tetmap();
    void makeindex2pointmap(point*& idx2verlist);
    void makesegmentmap(int*& idx2seglist, shellface**& segsperverlist);
    void makesubfacemap(int*& idx2facelist, shellface**& facesperverlist);
    void maketetrahedronmap(int*& idx2tetlist, tetrahedron**& tetsperverlist);

    // Point location routines.
    unsigned long randomnation(unsigned int choices);
    REAL distance2(tetrahedron* tetptr, point p);
    enum locateresult preciselocate(point searchpoint, triface* searchtet);
    enum locateresult locate(point searchpoint, triface* searchtet);
    enum locateresult adjustlocate(point searchpoint, triface* searchtet,
                                   enum locateresult precise, REAL epspp);    

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh Local Transformation Operators                                       //
//                                                                           //
// These operators (including flips, insert & remove vertices and so on) are //
// used to transform (or replace) a set of mesh elements into another set of //
// mesh elements.                                                            //
//                                                                           //
// Flip Types                                                                //
//                                                                           //
// If abc is a hull face, it is unflipable, and is locally Delaunay.  In the //
// following, we assume abc is an interior face, and the other tetrahedron   //
// adjoining at abc is bace.                                                 //
//                                                                           //
// If the convex hull CH of the set {a, b, c, d, e} only has four vertices,  //
// i.e., one vertex lies inside CH, then abc is unflipable, and is locally   //
// Delaunay. If CH is the vertex set itself, we have the following cases to  //
// determine whether abc is flipable or not.                                 //
//                                                                           //
// If no four points of {a, b, c, d, e} are coplanar, a 2-to-3 flip can be   //
// applied to abc if the edge de crosses the triangle abc; a 3-to-2 flip can //
// be applied to abc if ab crosses cde, and abde exists, otherwise, face abc //
// is unflipable, i.e., the tetrahedron abde is not present.                 //
//                                                                           //
// If four points of {a, b, c, d, e} are coplanar (two faces are coplanar).  //
// Assume faces abd and abe are coplanar (it is impossible be abc). If a, b, //
// d, e form a non-convex quadrilateral, then abc is unflipable, furthermore,//
// it is locally Delaunay.  Assume they are convex quadrilateral, if abd and //
// abe are hull faces, a 2-to-2 flip can be applied to abc;  if abd and abe  //
// are interior faces,  assume two tetrahedra adjoining abd and abe at the   //
// opposite sides are abdg and abef, respectively.  If g = f, a 4-to-4 flip  //
// can be applied to abc, otherwise, abc is unflipable.                      //
//                                                                           //
// There are other cases which can cause abc unflipable. If abc is a subface,//
// a 2-to-3 flip is forbidden;  if ab is a subsegment, flips 3-to-2, 2-to-2, //
// and 4-to-4 are forbidden.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // Mesh transformation routines.
    enum fliptype categorizeface(triface& horiz);
    void enqueueflipface(triface& checkface, queue* flipqueue);
    void enqueueflipedge(face& checkedge, queue* flipqueue);
    void flip23(triface* flipface, queue* flipqueue);
    void flip32(triface* flipface, queue* flipqueue);
    void flip22(triface* flipface, queue* flipqueue);
    void flip22sub(face* flipedge, queue* flipqueue);
    long flip(queue* flipqueue, badface **plastflip);
    void undoflip(badface *lastflip);

    void splittetrahedron(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetrahedron(triface* splittet);
    void splittetface(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetface(triface* splittet);
    void splitsubface(point newpoint, face* splitface, queue* flipqueue);
    void unsplitsubface(face* splitsh);
    void splittetedge(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetedge(triface* splittet);
    void splitsubedge(point newpoint, face* splitsh, queue* flipqueue);
    void unsplitsubedge(face* splitsh);
    enum insertsiteresult insertsite(point newpoint, triface* searchtet,
                                     bool approx, queue* flipqueue);
    void undosite(enum insertsiteresult insresult, triface* splittet, 
                  point torg, point tdest, point tapex, point toppo);
    void inserthullsite(point inspoint, triface* horiz, queue* flipqueue,
                        link* hulllink, int* worklist);
    void collectcavtets(point newpoint, list* cavtetlist);

    // Incremental flip Delaunay triangulation routines.
    void incrflipinit(queue* insertqueue);
    long incrflipdelaunay();

    // Surface triangulation routines.
    enum locateresult locatesub(point searchpt, face* searchsh, int stopatseg);
    long flipsub(queue* flipqueue);
    bool incrflipinitsub(int facetidx, list* ptlist, point* idx2verlist);
    void collectvisiblesubs(int facetidx, point inspoint, face* horiz,
                            queue* flipqueue);
    void incrflipdelaunaysub(int facetidx, list* ptlist, point* idx2verlist,
                             queue* flipqueue);
    enum finddirectionresult finddirectionsub(face* searchsh, point tend);
    void insertsubseg(face* tri);
    bool scoutsegmentsub(face* searchsh, point tend);
    void delaunayfixup(face* fixupsh, int leftside);
    void constrainededge(face* startsh, point tend);
    void insertsegmentsub(point tstart, point tend);
    void infecthullsub(memorypool* viri);
    void plaguesub(memorypool* viri);
    void carveholessub(int holes, REAL* holelist);
    void triangulatefacet(int facetidx, list* ptlist, list* conlist,
                          point* idx2verlist, queue* flipqueue);
    void unifysegments();
    void mergefacets(queue* flipqueue);
    long meshsurface();

    // Detect intersecting facets of PLC.
    void interecursive(shellface** subfacearray, int arraysize, int axis,
                       REAL bxmin, REAL bxmax, REAL bymin, REAL bymax,
                       REAL bzmin, REAL bzmax, int* internum);
    void detectinterfaces(); 

    // Segments recovery routines.
    void markacutevertices(REAL acuteangle);
    enum finddirectionresult finddirection(triface* searchtet, point tend);
    void getsearchtet(point p1, point p2, triface* searchtet, point* tend);
    bool isedgeencroached(point p1, point p2, point testpt, bool degflag);
    point scoutrefpoint(triface* searchtet, point tend);
    point getsegmentorigin(face* splitseg);
    point getsplitpoint(face* splitseg, point refpoint);
    void delaunizesegments();

    // Perturbation routines.
    REAL randgenerator(REAL range);
    void tallcocirsubs(REAL eps, queue* cocirque);
    void tallcosphtets(REAL eps, queue* cosphque);
    void perturbcocirsub(face* cocirsub, queue* flipqueue);
    bool coplanartest(triface* cosphtet1, triface* cosphtet2, REAL eps);
    void perturbcosphtet(triface* cosphtet, queue* flipqueue);
    void incrperturbvertices(REAL eps);

    // Constrained Delaunay triangulation routines.
    bool insertsubface(face* insertsh, triface* searchtet);
    bool tritritest(triface* checktet, point p1, point p2, point p3);
    void initializecavity(list* floorlist, list* ceillist, list* floorptlist,
                          list* ceilptlist, link* frontlink, link* ptlink);
    bool reducecavity(link* frontlink, link* ptlink, queue* flipqueue);
    bool reducecavity1(link* frontlink, queue* flipqueue);
    void triangulatecavity(list* floorlist, list* ceillist, list* floorptlist,
                           list* ceilptlist);
    void formmissingregion(face* missingsh, list* missingshlist,
                           list* equatptlist, int* worklist);
    bool scoutcrossingedge(list* missingshlist, list* boundedgelist,
                           list* crossedgelist, int* worklist);
    void rearrangesubfaces(list* missingshlist, list* boundedgelist,
                           list* equatptlist, int* worklist);
    void recoversubfaces(list* missingshlist, list* crossedgelist,
                         list* equatptlist, int* worklist);
    void constrainedfacets();

    // Carving out holes and concavities routines.
    void indenthull();
    void infecthull(memorypool *viri);
    void plague(memorypool *viri);
    void regionplague(memorypool *viri, REAL attribute, REAL volume);
    void carveholes();

    // Mesh update rotuines.
    long reconstructmesh();
    void insertaddpoints();

    // Mesh repair routines.
    enum degetype checkdegetet(triface* testtet, REAL eps, REAL dihed,
                               list* degetetlist);
    void removetetbypeeloff(triface *badtet, queue* flipqueue);
    void removetetbyflip32(triface *badtet, queue* flipqueue);
    bool removekite(triface* akite, queue* flipqueue);
    void talldegetets(REAL eps, REAL dihed, list* degetetlist);
    void repairdegetets(REAL eps, REAL dihed);

    // Delaunay refinement routines.
    void initializerpsarray(REAL* rpsarray);
    void marksharpfacets(int* idx2facetlist, REAL dihedbound);
    void enqueuebadtet(triface *instet, REAL ratio, point insorg,
                       point insdest, point insapex, point insoppo,
                       point inscent);
    badface* dequeuebadtet();
    bool checkseg4encroach(face* testseg, point testpt, bool enqueueflag);
    bool checksub4encroach(face* testsub, point testpt, bool enqueueflag);
    bool checksub4badqual(face* testsub);
    bool checktet4badqual(triface* testtet);
    bool checkseg4splitting(face* testseg, REAL* rpsarray, bool bqual);
    bool checksub4splitting(face* testsub);
    void doqualchecktetlist();
    bool tallencsegs(point testpt, list *cavtetlist);
    bool tallencsubs(point testpt, list *cavtetlist);
    void tallbadtetrahedrons();
    void repairencsegs(REAL* rpsarray, bool bqual, queue* flipqueue);
    void repairencsubs(REAL* rpsarray, int* idx2facetlist, list* cavtetlist,
                       queue* flipqueue);
    void repairbadtets(REAL* rpsarray, int* idx2facetlist, list* cavtetlist,
                       queue* flipqueue);
    void enforcequality();

    // I/O routines
    void transfernodes();
    void jettisonnodes();
    void highorder();
    void outnodes(tetgenio* out);
    void outelements(tetgenio* out);
    void outfaces(tetgenio* out);
    void outhullfaces(tetgenio* out);
    void outsubfaces(tetgenio* out);
    void outsubsegments(tetgenio* out);
    void outneighbors(tetgenio* out);
    void outsmesh(char* smfilename);
    void outmesh2medit(char* mfilename);
    void outmesh2gid(char* gfilename);
    void outmesh2off(char* ofilename);

    // User interaction routines.
    void internalerror();
    void checkmesh();
    void checkshells();
    void checkdelaunay(queue* flipqueue);
    void checkconforming();
    void qualitystatistics();
    void statistics();

  public:

    // Constructor and destructor.
    tetgenmesh();
    ~tetgenmesh();

};                                               // End of class tetgenmesh.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    Interface for using TetGen's library to generate      //
//                     Delaunay tetrahedralizations, constrained Delaunay    //
//                     tetrahedralizations, quality tetrahedral meshes.      //
//                                                                           //
// Two functions (interfaces) are available. The difference is only the way  //
// of passing switches.  One directly accepts an object of 'tetgenbehavior', //
// the other accepts a string which is the same as one can used in command   //
// line.  The latter may be more convenient for users who don't know the     //
// 'tetgenbehavir' structure.                                                //
//                                                                           //
// 'in' is the input object containing a PLC or a list of points. It should  //
// not be a NULL.  'out' is for outputting the mesh or tetrahedralization    //
// created by TetGen. If it is NULL, the output will be redirect to file(s). //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(tetgenbehavior *b, tetgenio *in, tetgenio *out);
void tetrahedralize(char *switches, tetgenio *in, tetgenio *out);

#endif // #ifndef tetgenH
