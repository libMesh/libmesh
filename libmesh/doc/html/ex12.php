<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex12",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
Example 12 -- The \p MeshData class

<br><br>The previous examples covered the certainly involved
aspects of simulating multiple equation systems, and prior 
to this even using adaptive refinement.  I.e.: cool stuff.

<br><br>This example now reduces brain effort significantly,
and presents some supplements concerning data I/O,
especially how to handle the \p MeshData object.
The \p MeshData class may be used for handling input
data for a simulation, like actual material properties, 
(not indicators, like in the \p BoundaryInfo class),
or mode shapes that serve as input for acoustic radiation
simulations.  The use of the \p MeshData during simulation
is straightforward:  the 

<br><br>- \p Number               \p MeshData::operator()(const Node*, int),
get the i-th floating-point value associated with the node
- \p bool                 \p MeshData::has_data(const Node*),
verify whether a certain node has data associated
- \p std::vector&lt;Number&gt;&amp; \p MeshData::get_data (const Node* node)
to get read-access to the data associated with this node (better
make sure first that this node @e has data, see \p has_data() )
- iterator for nodal data \p MeshData::const_node_data_iterator
to directly iterate over the set of nodes that hold data, 
instead of asking through \p has_data() each time again.

<br><br>(and corresponding methods for const Elem*) provide access to
the floating-point data associated with nodes/elements.
This example does @e not show how to use these aforementioned
methods, this is straightforward.  Simply check if the current 
Elem* has a node with data. If so, add this contribution to the 
RHS, or whatever. Or ask the \p MeshData for each Elem* for its 
material properties...

<br><br>Here, only the actual file I/O necessary to handle such 
nodal/element data is presented.  Lean back, relax...


<br><br>C++ include files that we need.
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Functions to initialize the library
and provide some further features (e.g. 
our own \p pi)
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh and 
\p MeshData functionality.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        
</pre>
</div>
<div class = "comment">
Function prototype for creating artificial nodal data
that can be inserted into a \p MeshData object.
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh&amp; mesh,
                                     std::map&lt;const Node*, std::vector&lt;Number&gt; &gt;&amp; art_data);
        
</pre>
</div>
<div class = "comment">
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize the library.
</div>

<div class ="fragment">
<pre>
          libMesh::init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Force all our objects to have local scope.  
</div>

<div class ="fragment">
<pre>
          {    
</pre>
</div>
<div class = "comment">
Check for proper usage. The program is designed to be run
as follows:

<br><br>./ex12 -d 3 in_mesh.unv

<br><br>where in_mesh.unv should better be a Universal file.
</div>

<div class ="fragment">
<pre>
            if (argc &lt; 4)
              {
                std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -d &lt;dim&gt; in_mesh.unv"
                          &lt;&lt; std::endl;
                
                error();
              }
        
</pre>
</div>
<div class = "comment">
Get the dimensionality of the mesh from argv[2]
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = atoi(argv[2]);
            
</pre>
</div>
<div class = "comment">
The filename of the mesh
</div>

<div class ="fragment">
<pre>
            const std::string mesh_file = argv[3];
        
</pre>
</div>
<div class = "comment">
The following example makes currently sense
only with a Universal file
</div>

<div class ="fragment">
<pre>
            if (mesh_file.rfind(".unv") &gt;= mesh_file.size())
              {
        
                std::cerr &lt;&lt; "ERROR:  This example works only properly with a Universal mesh file!"
                          &lt;&lt; std::endl;
                
                error();
              }
          
</pre>
</div>
<div class = "comment">
Some notes on \p MeshData:

<br><br>- The \p MeshData is @e not a mesh!  Consult the \p Mesh,
\p MeshBase, and \p BoundaryMesh classes for details on
handling meshes!

<br><br>- The \p MeshData is an object handling arbitrary floating-point 
data associated with mesh entities (nodes, elements), currently
most features only available for nodal data.  However,
extending to element-data (does @e anybody need this?)
is straightforward.

<br><br>- Currently std::vector&lt;Number&gt; is the data (associated
with nodes/elements) that can be handled by \p MeshData,

<br><br>- In order to provide @e full functionality, the \p MeshData
@e has to be activated prior to reading in a mesh.  However,
there is a so-called compatibility mode available when you 
(intentionally) forgot to "turn on" the \p MeshData.

<br><br>- It is possible to provide user-defined nodal data that
may subsequently be used or written to file.

<br><br>- Translate the nodal-/element-associated data to formats
suitable for visual inspection, e.g. GMV files.

<br><br>- Two I/O file formats are currently supported: the Universal 
file format with dataset 2414, and an extension of 
the libMesh-own xda/xdr format, named xtr, xta.
Some details on this:

<br><br>- The xtr/xta format is simply an extension of the
xdr/xda format to read/write also nodal/element-associated
floating-point data.  The xtr interface of the \p MeshData
creates stand-alone files that may be binary or ASCII.
You cannot append files created by the \p MeshData I/O methods
to a mesh file (xdr/xda).
The xtr files may be seen as a "backdrop", especially when
binary files are preferred.  Note that unv files are @e always
ASCII and may become pretty big!

<br><br>- The unv format is an extremely versatile text-based file format
for arbitrary data.  Its functionality is @e large, and \p libMesh
supports only a small share: namely the I/O for nodes, elements and
arbitrary node-/element-associated data, stored in the 
so-called datasets "2411", "2412", and "2414", respectively.
Here, only the last dataset is covered.  The two first datasets are
implemented in the Universal support for I/O of meshes.  A single
unv file may hold @e multiple datasets of type "2414".  To
distinguish data, each dataset has a header.  The \p libMesh pendant
to this header is the \p MeshDataUnvHeader class.  When you
read/write unv files using the \p MeshData, you @always
automatically also read/write such headers.

<br><br>Enough babble for now.  Examples. 
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
Create a mesh with the requested dimension.
Below we require a 3-dim mesh, therefore assert
it.
</div>

<div class ="fragment">
<pre>
              assert (dim == 3);
              Mesh mesh(dim);
        
</pre>
</div>
<div class = "comment">
Activate the \p MeshData of the mesh, so that
we can remember the node and element ids used
in the file.  When we do not activate the \p MeshData,
then there is no chance to import node- or element-
associated data.
</div>

<div class ="fragment">
<pre>
              mesh.data.activate();
        
</pre>
</div>
<div class = "comment">
Now we can safely read the input mesh.  Note that
this should better be a .unv or .xdr/xda file, otherwise
we cannot load/save associated data.
</div>

<div class ="fragment">
<pre>
              mesh.read (mesh_file);
            
</pre>
</div>
<div class = "comment">
Print information about the mesh and the data
to the screen.  Obviously, there is no data
(apart from node &amp; element ids) in it, yet.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Finished reading the mesh.  MeshData is active but empty:" &lt;&lt; std::endl
                        &lt;&lt; "---------------------------------------------------------" &lt;&lt; std::endl;
              
              mesh.print_info();
              mesh.data.print_info();
            
</pre>
</div>
<div class = "comment">
Create some artificial node-associated data and store
it in the mesh's \p MeshData.  Use a \p std::map for this. 
Let the function create_artificial_data do the work.
</div>

<div class ="fragment">
<pre>
              {
                std::map&lt;const Node*, std::vector&lt;Number&gt; &gt; artificial_data;
                
                create_artificial_data (mesh, artificial_data);
                
</pre>
</div>
<div class = "comment">
Before we let the map go out of scope, insert it into
the \p MeshData.  Note that by default the element-associated
data containers are closed, so that the \p MeshData is
ready for use.
</div>

<div class ="fragment">
<pre>
                mesh.data.insert_node_data(artificial_data);
        
</pre>
</div>
<div class = "comment">
Let \p artificial_data go out of scope
</div>

<div class ="fragment">
<pre>
              }
              
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "After inserting artificial data into the MeshData:" &lt;&lt; std::endl
                        &lt;&lt; "--------------------------------------------------" &lt;&lt; std::endl;
              
              mesh.data.print_info();
        
</pre>
</div>
<div class = "comment">
Write the artificial data into a universal
file.  Use a default header for this.
</div>

<div class ="fragment">
<pre>
              std::string first_out_data="data_first_out_with_default_header.unv";
              std::cout &lt;&lt; "Writing MeshData to: " &lt;&lt; first_out_data &lt;&lt; std::endl;
              mesh.data.write(first_out_data);
        
</pre>
</div>
<div class = "comment">
Alternatively, create your own header.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Attach our own MeshDataUnvHeader to the MeshData:" &lt;&lt; std::endl
                        &lt;&lt; "-------------------------------------------------" &lt;&lt; std::endl
                        &lt;&lt; " (note the warning: the number of values per node in" &lt;&lt; std::endl
                        &lt;&lt; "  my_header is not correct)" &lt;&lt; std::endl &lt;&lt; std::endl;
              MeshDataUnvHeader my_header;
        
</pre>
</div>
<div class = "comment">
Specify an int for this dataset.
This is particularly helpful
when there are multiple datasets in
a file and you want to find a specific one.
For this, use MeshDataUnvHeader::which_dataset(int),
which is not covered here.
</div>

<div class ="fragment">
<pre>
              my_header.dataset_label = 3;
        
</pre>
</div>
<div class = "comment">
Specify some text that helps the @e user 
identify the data.  This text is @e not used for
finding a specific dataset.  Leave default for
the remaining 2 lines.
</div>

<div class ="fragment">
<pre>
              my_header.id_lines_1_to_5[0] = "Artificial data";
              my_header.id_lines_1_to_5[1] = "sin curve in z-direction";
              my_header.id_lines_1_to_5[2] = "line in x-direction";
              
</pre>
</div>
<div class = "comment">
Also some float data, not associated with nodes,
can be stored in the header.
</div>

<div class ="fragment">
<pre>
              my_header.record_12[0] = libMesh::pi;
              
</pre>
</div>
<div class = "comment">
Now attach this header to the \p MeshData, and write
the same file again, but with the personalized header.
</div>

<div class ="fragment">
<pre>
              mesh.data.set_unv_header(&amp;my_header);
        
</pre>
</div>
<div class = "comment">
Write again to file.
</div>

<div class ="fragment">
<pre>
              std::string second_out_data="data_second_with_header_out.unv";
              std::cout &lt;&lt; "Writing MeshData to: " &lt;&lt; second_out_data &lt;&lt; std::endl;
              mesh.data.write(second_out_data);
              
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Before clearing the MeshData:" &lt;&lt; std::endl
                        &lt;&lt; "-----------------------------" &lt;&lt; std::endl;
              mesh.data.print_info();
        
</pre>
</div>
<div class = "comment">
Now clear only the data associated with nodes/elements,
but keep the node/element ids used in the mesh.  To
clear also these ids, use MeshData::slim() (not used here).
</div>

<div class ="fragment">
<pre>
              mesh.data.clear();
        
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "After clearing the MeshData:" &lt;&lt; std::endl
                        &lt;&lt; "----------------------------" &lt;&lt; std::endl;
              mesh.data.print_info();
        
</pre>
</div>
<div class = "comment">
Now the \p MeshData is open again to read data from
file.  Read the file that we created first.
</div>

<div class ="fragment">
<pre>
              mesh.data.read(first_out_data);
              
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "After re-reading the first file:" &lt;&lt; std::endl
                        &lt;&lt; "--------------------------------" &lt;&lt; std::endl;
              mesh.data.print_info();
        
</pre>
</div>
<div class = "comment">
Let the mesh, the unv_header etc go out of scope, and
do another example.
</div>

<div class ="fragment">
<pre>
            }
        
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "----------------------------------------------" &lt;&lt; std::endl
                      &lt;&lt; "---------- next example with MeshData --------" &lt;&lt; std::endl
                      &lt;&lt; "----------------------------------------------" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Create a new mesh, read it again -- but this time
with de-activated \p MeshData.  Then we are
able to use the compatibility mode not only for
handling \p MeshData dat, but also to @e write 
a mesh in unv format.
The libMesh-internal node and element ids are used.
</div>

<div class ="fragment">
<pre>
            {
              Mesh mesh(dim);
        
</pre>
</div>
<div class = "comment">
Read the input mesh, but with deactivated \p MeshData.
</div>

<div class ="fragment">
<pre>
              mesh.read (mesh_file);
            
</pre>
</div>
<div class = "comment">
Print information about the mesh and the data
to the screen.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "De-activated MeshData:" &lt;&lt; std::endl
                        &lt;&lt; "----------------------" &lt;&lt; std::endl;
              mesh.print_info();
              mesh.data.print_info();
         
</pre>
</div>
<div class = "comment">
Write the @e mesh (not the MeshData!) as .unv file.
In general, the \p MeshBase interface for .unv I/O
needs an active \p MeshData.  However, use compatibility
mode to at least write a .unv file, but with the 
ids from libMesh.
</div>

<div class ="fragment">
<pre>
              const std::string out_mesh = "mesh_with_libmesh_ids.unv";
              std::cout &lt;&lt; "Writing _Mesh_ to: " &lt;&lt; out_mesh &lt;&lt; std::endl
                        &lt;&lt; "Try 'diff " &lt;&lt; out_mesh &lt;&lt; " " &lt;&lt; mesh_file &lt;&lt; "'" &lt;&lt; std::endl
                        &lt;&lt; "to see the differences in node numbers." &lt;&lt; std::endl
                        &lt;&lt; "---------------------------------------" &lt;&lt; std::endl
                        &lt;&lt; std::endl;
              mesh.write(out_mesh);
        
</pre>
</div>
<div class = "comment">
Again create some artificial node-associated data,
as before.
</div>

<div class ="fragment">
<pre>
              {
                std::map&lt;const Node*, std::vector&lt;Number&gt; &gt; artificial_data;
                create_artificial_data (mesh, artificial_data);
                mesh.data.insert_node_data(artificial_data);
              }
        
</pre>
</div>
<div class = "comment">
Note that even with (only) compatibility mode MeshData
data can be written.  But again, the ids from libMesh
are used.  Consult the warning messages issued in
DEBUG mode.  And the user @e has to specify that the
\p MeshData should change to compatibility mode.
</div>

<div class ="fragment">
<pre>
              mesh.data.enable_compatibility_mode();
        
</pre>
</div>
<div class = "comment">
Now that compatibility mode is used, data can be written.
_Without_ explicitly enabling compatibility mode, we
would get an error message and (with the following write()
statement).
</div>

<div class ="fragment">
<pre>
              std::string mesh_data_file = "data_third_with_libmesh_ids_out.unv";
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Writing MeshData to: " &lt;&lt; mesh_data_file &lt;&lt; std::endl
                        &lt;&lt; "----------------------------------------------------------" 
                        &lt;&lt; std::endl &lt;&lt; std::endl;
              mesh.data.write (mesh_data_file);
        
        #ifdef HAVE_ZLIB_H
        
</pre>
</div>
<div class = "comment">
As may already seen, UNV files are text-based, so the may
become really big.  When \p ./configure found \p zlib.h,
then we may also @e read or @e write \p .unv files in gzip'ed 
format! -- Pretty cool, and also pretty fast, due to zlib.h.

<br><br>Note that this works also for mesh files, not only for 
meshdata files.

<br><br>In order to write a ".unv.gz" file instead of a ".unv" file,
simply provide the full name with ".gz" appended to the
write method; it will then figure out whether this file should
be gzip'ed or not.
</div>

<div class ="fragment">
<pre>
              std::string packed_mesh_data_file =  "packed_" + mesh_data_file + ".gz";
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Writing gzip'ed MeshData to: " &lt;&lt; packed_mesh_data_file &lt;&lt; std::endl
                        &lt;&lt; "---------------------------------------------------------------------------" &lt;&lt; std::endl
                        &lt;&lt; " To verify the integrity of the packed version, type:" &lt;&lt; std::endl &lt;&lt; std::endl
                        &lt;&lt; "   gunzip " &lt;&lt; packed_mesh_data_file &lt;&lt; "; " &lt;&lt; std::endl
                        &lt;&lt; "   diff packed_" &lt;&lt; mesh_data_file &lt;&lt; " " 
                        &lt;&lt; mesh_data_file &lt;&lt; std::endl &lt;&lt; std::endl;
              
              mesh.data.write (packed_mesh_data_file);
              
        #endif
        
</pre>
</div>
<div class = "comment">
And now a last gimmick: The \p MeshData::translate()
conveniently converts the nodal- or element-associated
data (currently only nodal) to vectors that may be used 
for writing a mesh with "solution" vectors, where the solution 
vector contains the data from the \p MeshData.  Particularly
useful for @e inspecting the data contained in \p MeshData.

<br><br>And even better: the user can choose the mesh for which
to export the data.  E.g. not only use the \p mesh
itself, but alternatively use the \p BoundaryMesh 
(any mesh that uses the @e same nodes, i.e. the \p Node* 
have to be the same.  Only exception that will not work:
A mesh created using Mesh::create_submesh()
actually will @e not work with \p MeshData::translate() ).

<br><br>All in all not bad, hm?
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
have a vector for the actual values and a vector
for the names of the data available.
</div>

<div class ="fragment">
<pre>
                std::vector&lt;Number&gt; translated_data;
                std::vector&lt;std::string&gt; data_names;
                
</pre>
</div>
<div class = "comment">
Use the \p mesh itself.  Alternatively, use the 
\p BoundaryMesh of \p mesh.
</div>

<div class ="fragment">
<pre>
                mesh.data.translate (mesh,
                                     translated_data,
                                     data_names);
        
        
</pre>
</div>
<div class = "comment">
And write the data to a GMV file
</div>

<div class ="fragment">
<pre>
                const std::string gmv_file = "data_and_mesh_out.gmv";
                std::cout &lt;&lt; std::endl 
                          &lt;&lt; "Writing the data from the MeshData to the GMV file " 
                          &lt;&lt; gmv_file &lt;&lt; std::endl
                          &lt;&lt; "------------------------------------------------------------------------" 
                          &lt;&lt; std::endl;
                
                mesh.write_gmv_binary (gmv_file,
                                       &amp;translated_data,
                                       &amp;data_names);
                
</pre>
</div>
<div class = "comment">
Let the vectors with translated data
go out of scope.
</div>

<div class ="fragment">
<pre>
              }
              
</pre>
</div>
<div class = "comment">
Let the second mesh go out of scope.
</div>

<div class ="fragment">
<pre>
            }
          }
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return libMesh::close();
        }
        
</pre>
</div>
<div class = "comment">
This function creates the data to populate the \p MeshData object
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh&amp; mesh,
                                     std::map&lt;const Node*, std::vector&lt;Number&gt; &gt;&amp; art_data)
        {
</pre>
</div>
<div class = "comment">
get the bounding box to have some sensible data
</div>

<div class ="fragment">
<pre>
          std::pair&lt;Point, Point&gt; b_box = mesh.bounding_box();
        
          const Real z_min = b_box.first (2);
          const Real z_max = b_box.second(2);
          assert (fabs(z_max-z_min) &gt; TOLERANCE);
        
          const Real x_min = b_box.first (0);
          const Real x_max = b_box.second(0);
          assert (fabs(x_max-x_min) &gt; TOLERANCE);
        
        
          const_node_iterator node_it = mesh.nodes_begin();
          const const_node_iterator node_end = mesh.nodes_end();
        
          for (; node_it != node_end; ++node_it)
            {
</pre>
</div>
<div class = "comment">
All the vectors in \p artificial_data @e have to have the
same size.  Here we use only two entries per node,
but theoretically arbitrary size is possible.
</div>

<div class ="fragment">
<pre>
              std::vector&lt;Number&gt; node_data;
              node_data.resize(2);
        
</pre>
</div>
<div class = "comment">
Use a sin curve in z-direction and a linear
increase in x-direction
</div>

<div class ="fragment">
<pre>
              const Point&amp; p = **node_it;
                
              const Real z_normalized = (p(2)-z_min)/(z_max-z_min);
              const Real x_normalized = (p(0)-x_min)/(x_max-x_min);
        
              node_data[0] = sin(2*libMesh::pi*z_normalized);
              node_data[1] = x_normalized;
              
</pre>
</div>
<div class = "comment">
Insert the current data together with a pointer to
the current node in the map.
</div>

<div class ="fragment">
<pre>
              art_data.insert (std::make_pair(*node_it,node_data));
            }
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;math.h&gt;
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> create_artificial_data (<FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh,
  			     std::map&lt;<FONT COLOR="#228B22"><B>const</FONT></B> Node*, std::vector&lt;Number&gt; &gt;&amp; art_data);
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
  
    {    
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 4)
        {
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -d &lt;dim&gt; in_mesh.unv&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	error();
        }
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = atoi(argv[2]);
      
      <FONT COLOR="#228B22"><B>const</FONT></B> std::string mesh_file = argv[3];
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (mesh_file.rfind(<FONT COLOR="#BC8F8F"><B>&quot;.unv&quot;</FONT></B>) &gt;= mesh_file.size())
        {
  
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR:  This example works only properly with a Universal mesh file!&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	error();
        }
    
      {
        assert (dim == 3);
        Mesh mesh(dim);
  
        mesh.data.activate();
  
        mesh.read (mesh_file);
      
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Finished reading the mesh.  MeshData is active but empty:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;---------------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
        
        mesh.print_info();
        mesh.data.print_info();
      
        {
  	std::map&lt;<FONT COLOR="#228B22"><B>const</FONT></B> Node*, std::vector&lt;Number&gt; &gt; artificial_data;
  	
  	create_artificial_data (mesh, artificial_data);
  	
  	mesh.data.insert_node_data(artificial_data);
  
        }
        
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;After inserting artificial data into the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;--------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
        
        mesh.data.print_info();
  
        std::string first_out_data=<FONT COLOR="#BC8F8F"><B>&quot;data_first_out_with_default_header.unv&quot;</FONT></B>;
        std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; first_out_data &lt;&lt; std::endl;
        mesh.data.write(first_out_data);
  
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Attach our own MeshDataUnvHeader to the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;-------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; (note the warning: the number of values per node in&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  my_header is not correct)&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl;
        MeshDataUnvHeader my_header;
  
        my_header.dataset_label = 3;
  
        my_header.id_lines_1_to_5[0] = <FONT COLOR="#BC8F8F"><B>&quot;Artificial data&quot;</FONT></B>;
        my_header.id_lines_1_to_5[1] = <FONT COLOR="#BC8F8F"><B>&quot;sin curve in z-direction&quot;</FONT></B>;
        my_header.id_lines_1_to_5[2] = <FONT COLOR="#BC8F8F"><B>&quot;line in x-direction&quot;</FONT></B>;
        
        my_header.record_12[0] = libMesh::pi;
        
        mesh.data.set_unv_header(&amp;my_header);
  
        std::string second_out_data=<FONT COLOR="#BC8F8F"><B>&quot;data_second_with_header_out.unv&quot;</FONT></B>;
        std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; second_out_data &lt;&lt; std::endl;
        mesh.data.write(second_out_data);
        
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Before clearing the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;-----------------------------&quot;</FONT></B> &lt;&lt; std::endl;
        mesh.data.print_info();
  
        mesh.data.clear();
  
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;After clearing the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;----------------------------&quot;</FONT></B> &lt;&lt; std::endl;
        mesh.data.print_info();
  
        mesh.data.read(first_out_data);
        
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;After re-reading the first file:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;--------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
        mesh.data.print_info();
  
      }
  
      std::cout &lt;&lt; std::endl 
  	      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;----------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
  	      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;---------- next example with MeshData --------&quot;</FONT></B> &lt;&lt; std::endl
  	      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;----------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
  
      {
        Mesh mesh(dim);
  
        mesh.read (mesh_file);
      
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;De-activated MeshData:&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;----------------------&quot;</FONT></B> &lt;&lt; std::endl;
        mesh.print_info();
        mesh.data.print_info();
   
        <FONT COLOR="#228B22"><B>const</FONT></B> std::string out_mesh = <FONT COLOR="#BC8F8F"><B>&quot;mesh_with_libmesh_ids.unv&quot;</FONT></B>;
        std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing _Mesh_ to: &quot;</FONT></B> &lt;&lt; out_mesh &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Try 'diff &quot;</FONT></B> &lt;&lt; out_mesh &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> &lt;&lt; mesh_file &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;'&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;to see the differences in node numbers.&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;---------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; std::endl;
        mesh.write(out_mesh);
  
        {
  	std::map&lt;<FONT COLOR="#228B22"><B>const</FONT></B> Node*, std::vector&lt;Number&gt; &gt; artificial_data;
  	create_artificial_data (mesh, artificial_data);
  	mesh.data.insert_node_data(artificial_data);
        }
  
        mesh.data.enable_compatibility_mode();
  
        std::string mesh_data_file = <FONT COLOR="#BC8F8F"><B>&quot;data_third_with_libmesh_ids_out.unv&quot;</FONT></B>;
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; mesh_data_file &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;----------------------------------------------------------&quot;</FONT></B> 
  		&lt;&lt; std::endl &lt;&lt; std::endl;
        mesh.data.write (mesh_data_file);
  
  #ifdef HAVE_ZLIB_H
  
        std::string packed_mesh_data_file =  <FONT COLOR="#BC8F8F"><B>&quot;packed_&quot;</FONT></B> + mesh_data_file + <FONT COLOR="#BC8F8F"><B>&quot;.gz&quot;</FONT></B>;
        std::cout &lt;&lt; std::endl 
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing gzip'ed MeshData to: &quot;</FONT></B> &lt;&lt; packed_mesh_data_file &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;---------------------------------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; To verify the integrity of the packed version, type:&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;   gunzip &quot;</FONT></B> &lt;&lt; packed_mesh_data_file &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;; &quot;</FONT></B> &lt;&lt; std::endl
  		&lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;   diff packed_&quot;</FONT></B> &lt;&lt; mesh_data_file &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> 
  		&lt;&lt; mesh_data_file &lt;&lt; std::endl &lt;&lt; std::endl;
        
        mesh.data.write (packed_mesh_data_file);
        
  #endif
  
        {
  	std::vector&lt;Number&gt; translated_data;
  	std::vector&lt;std::string&gt; data_names;
  	
  	mesh.data.translate (mesh,
  			     translated_data,
  			     data_names);
  
  
  	<FONT COLOR="#228B22"><B>const</FONT></B> std::string gmv_file = <FONT COLOR="#BC8F8F"><B>&quot;data_and_mesh_out.gmv&quot;</FONT></B>;
  	std::cout &lt;&lt; std::endl 
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Writing the data from the MeshData to the GMV file &quot;</FONT></B> 
  		  &lt;&lt; gmv_file &lt;&lt; std::endl
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;------------------------------------------------------------------------&quot;</FONT></B> 
  		  &lt;&lt; std::endl;
  	
  	mesh.write_gmv_binary (gmv_file,
  			       &amp;translated_data,
  			       &amp;data_names);
  	
        }
        
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close();
  }
  
  <FONT COLOR="#228B22"><B>void</FONT></B> create_artificial_data (<FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh,
  			     std::map&lt;<FONT COLOR="#228B22"><B>const</FONT></B> Node*, std::vector&lt;Number&gt; &gt;&amp; art_data)
  {
    std::pair&lt;Point, Point&gt; b_box = mesh.bounding_box();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real z_min = b_box.first (2);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real z_max = b_box.second(2);
    assert (fabs(z_max-z_min) &gt; TOLERANCE);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x_min = b_box.first (0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x_max = b_box.second(0);
    assert (fabs(x_max-x_min) &gt; TOLERANCE);
  
  
    const_node_iterator node_it = mesh.nodes_begin();
    <FONT COLOR="#228B22"><B>const</FONT></B> const_node_iterator node_end = mesh.nodes_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (; node_it != node_end; ++node_it)
      {
        std::vector&lt;Number&gt; node_data;
        node_data.resize(2);
  
        <FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p = **node_it;
  	
        <FONT COLOR="#228B22"><B>const</FONT></B> Real z_normalized = (p(2)-z_min)/(z_max-z_min);
        <FONT COLOR="#228B22"><B>const</FONT></B> Real x_normalized = (p(0)-x_min)/(x_max-x_min);
  
        node_data[0] = sin(2*libMesh::pi*z_normalized);
        node_data[1] = x_normalized;
        
        art_data.insert (std::make_pair(*node_it,node_data));
      }
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex12
***************************************************************
 

 Universal file Interface:
  Counting nodes and elements
  Convert from "D" to "e"
  Nodes   : 3977
  Elements: 3520
  Reading nodes
  Reading elements
  Finished.


Finished reading the mesh.  MeshData is active but empty:
---------------------------------------------------------
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=3977
  n_elem()=3520
   n_local_elem()=3520
   n_active_elem()=3520
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 MeshData Information:
  object activated.


After inserting artificial data into the MeshData:
--------------------------------------------------
 MeshData Information:
  object activated.
  Element associated data initialized.
   n_val_per_elem()=0
   n_elem_data()=0
  Node associated data initialized.
   n_val_per_node()=2
   n_node_data()=3977

Writing MeshData to: data_first_out_with_default_header.unv

Attach our own MeshDataUnvHeader to the MeshData:
-------------------------------------------------
 (note the warning: the number of values per node in
  my_header is not correct)

Writing MeshData to: data_second_with_header_out.unv
WARNING: nvaldc=3 of attached MeshDataUnvHeader object not valid!
         re-set nvaldc to 2

Before clearing the MeshData:
-----------------------------
 MeshData Information:
  object activated.
  Element associated data initialized.
   n_val_per_elem()=0
   n_elem_data()=0
  Node associated data initialized.
   n_val_per_node()=2
   n_node_data()=3977


After clearing the MeshData:
----------------------------
 MeshData Information:
  object activated.


After re-reading the first file:
--------------------------------
 MeshData Information:
  object activated.
  Element associated data initialized.
   n_val_per_elem()=0
   n_elem_data()=0
  Node associated data initialized.
   n_val_per_node()=2
   n_node_data()=3977


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 4Elem reference count information:
| Creations:    24984
| Destructions: 24984
| 4Node reference count information:
| Creations:    3977
| Destructions: 3977
 ---------------------------------------------------------------------------- 

----------------------------------------------
---------- next example with MeshData --------
----------------------------------------------

 Universal file Interface:
  Counting nodes and elements
  Convert from "D" to "e"
  Nodes   : 3977
  Elements: 3520
  Reading nodes
  Reading elements
  Finished.


De-activated MeshData:
----------------------
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=3977
  n_elem()=3520
   n_local_elem()=3520
   n_active_elem()=3520
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 MeshData neither active nor in compatibility mode.

Writing _Mesh_ to: mesh_with_libmesh_ids.unv
Try 'diff mesh_with_libmesh_ids.unv ../../examples/ex8/pipe-mesh.unv'
to see the differences in node numbers.
---------------------------------------


 Universal file Interface:

*************************************************************************
* WARNING: MeshData neither active nor in compatibility mode.           *
*          Enable compatibility mode for MeshData.  Use this Universal  *
*          file with caution: libMesh node and element ids are used.    *
*************************************************************************

  Writing 3977 nodes
  Writing elements
  Finished writing 3520 elements

Writing MeshData to: data_third_with_libmesh_ids_out.unv
----------------------------------------------------------

WARNING: MeshData in compatibility mode.  Node and element ids
         written to file may differ from libMesh numbering
         next time this file is read!

Writing gzip'ed MeshData to: packed_data_third_with_libmesh_ids_out.unv.gz
---------------------------------------------------------------------------
 To verify the integrity of the packed version, type:

   gunzip packed_data_third_with_libmesh_ids_out.unv.gz; 
   diff packed_data_third_with_libmesh_ids_out.unv data_third_with_libmesh_ids_out.unv

WARNING: MeshData in compatibility mode.  Node and element ids
         written to file may differ from libMesh numbering
         next time this file is read!

Writing the data from the MeshData to the GMV file data_and_mesh_out.gmv
------------------------------------------------------------------------
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
Option left: name:-d value: 3

 ----------------------------------------------------------------------------
| Time:           Mon Nov 10 15:42:16 2003
| OS:             Linux
| HostName:       hactar
| OS Release      2.4.20-19.9smp
| OS Version:     #1 SMP Tue Jul 15 17:04:18 EDT 2003
| Machine:        i686
| Username:       benkirk
 ----------------------------------------------------------------------------
 ----------------------------------------------------------------------------
| libMesh Performance: Alive time=1.69033, Active time=1.52033
 ----------------------------------------------------------------------------
| Event                         nCalls  Total       Avg         Percent of   |
|                                       Time        Time        Active Time  |
|----------------------------------------------------------------------------|
|                                                                            |
|                                                                            |
| Mesh                                                                       |
|   read()                      2       0.9003      0.450145    59.22        |
|   write()                     1       0.1339      0.133948    8.81         |
|                                                                            |
| MeshBase                                                                   |
|   find_neighbors()            2       0.1889      0.094473    12.43        |
|   renumber_nodes_and_elem()   2       0.0036      0.001777    0.23         |
|                                                                            |
| MeshData                                                                   |
|   read()                      1       0.1333      0.133305    8.77         |
|   translate()                 1       0.0017      0.001714    0.11         |
|   write()                     4       0.1586      0.039644    10.43        |
 ----------------------------------------------------------------------------
| Totals:                       13      1.5203                  100.00       |
 ----------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex12
***************************************************************
</pre>
</div>
<?php make_footer() ?>
</body>
</html>
<?php if (0) { ?>
\#Local Variables:
\#mode: html
\#End:
<?php } ?>
