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
<h1>Example 12 - The <code>MeshData</code> class</h1>

<br><br>The previous examples covered the certainly involved
aspects of simulating multiple equation systems, and prior 
to this even using adaptive refinement.  I.e.: cool stuff.

<br><br>This example now reduces brain effort significantly,
and presents some supplements concerning data I/O,
especially how to handle the <code> MeshData </code> object.
The <code> MeshData </code> class may be used for handling input
data for a simulation, like actual material properties, 
(not indicators, like in the <code> BoundaryInfo </code> class),
or mode shapes that serve as input for acoustic radiation
simulations.  The use of the <code> MeshData </code> during simulation
is straightforward:  the 
<ul>
<li>
<code>Number MeshData::operator()(const Node*, int)</code>,
get the i-th floating-point value associated with the node
</li>
<li>
<code> bool MeshData::has_data(const Node*)</code>,
verify whether a certain node has data associated
</li>
<li>
<code>std::vector<Number>& MeshData::get_data (const Node* node)</code>
to get read-access to the data associated with this node (better
make sure first that this node <i>has</i> data, see <code> has_data() )
</li>
<li>
iterator for nodal data <code>MeshData::const_node_data_iterator</code>
to directly iterate over the set of nodes that hold data, 
instead of asking through <code>has_data()</code> each time again.
</li>
</ul>
(and corresponding methods for <code> const Elem*</code>) provide access to
the floating-point data associated with nodes/elements.
This example does <i>not</i> show how to use these aforementioned
methods, this is straightforward.  Simply check if the current 
<code>Elem*</code> has a node with data. If so, add this contribution to the 
RHS, or whatever. Or ask the <code> MeshData </code> for each <code>Elem*</code> for its 
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
our own pi)
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh and 
<code> MeshData </code> functionality.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        #include "gmv_io.h"
        
        
        
</pre>
</div>
<div class = "comment">
Function prototype for creating artificial nodal data
that can be inserted into a <code>MeshData</code> object.
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh&amp; mesh,
        			     std::map<const Node*, std::vector<Number> >& art_data);
        
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
<pre>
$ ./ex12 -d 3 in_mesh.unv
</pre>
where in_mesh.unv should be a Universal file.
</div>

<div class ="fragment">
<pre>
            if (argc &lt; 4)
              {
        	std::cerr << "Usage: " << argv[0] << " -d <dim> in_mesh.unv"
        		  << std::endl;
        	
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
        	std::cerr << "ERROR:  This example works only properly with a Universal mesh file!"
        		  << std::endl;
        	
        	error();
              }
          
</pre>
</div>
<div class = "comment">
Some notes on <code>MeshData</code>:
<ul>
<li>
The <code>MeshData</code> is <i>not</i> a mesh!  Consult the <code>Mesh</code>,
<code>MeshBase</code>, and <code>BoundaryMesh</code> classes for details on
handling meshes!
</li>
<li>
The <code>MeshData</code> is an object handling arbitrary floating-point 
data associated with mesh entities (nodes, elements), currently
most features only available for nodal data.  However,
extending to element-data (does <i>anybody</i> need this?)
is straightforward.
</li>
<li>
Currently std::vector<Number> is the data (associated
with nodes/elements) that can be handled by <code>MeshData</code>,
</li>
<li>
In order to provide <i>full</i> functionality, the <code>MeshData</code>
<i>has</i> to be activated prior to reading in a mesh.  However,
there is a so-called compatibility mode available when you 
(intentionally) forgot to "turn on" the <code>MeshData</code>.
</li>
<li>
It is possible to provide user-defined nodal data that
may subsequently be used or written to file.
</li>
<li>
Translate the nodal-/element-associated data to formats
suitable for visual inspection, e.g. GMV files.
</li>
<li>
Two I/O file formats are currently supported: the Universal 
file format with dataset 2414, and an extension of 
the libMesh-own xda/xdr format, named xtr, xta.
Some details on this:
</li>
<li>
The xtr/xta format is simply an extension of the
xdr/xda format to read/write also nodal/element-associated
floating-point data.  The xtr interface of the <code>MeshData</code>
creates stand-alone files that may be binary or ASCII.
You cannot append files created by the <code>MeshData</code> I/O methods
to a mesh file (xdr/xda).
The xtr files may be seen as a "backdrop", especially when
binary files are preferred.  Note that unv files are <i>always</i>
ASCII and may become pretty big!
</li>
<li>
The unv format is an extremely versatile text-based file format
for arbitrary data.  Its functionality is <i>large</i>, and <code>libMesh</code>
supports only a small share: namely the I/O for nodes, elements and
arbitrary node-/element-associated data, stored in the 
so-called datasets "2411", "2412", and "2414", respectively.
Here, only the last dataset is covered.  The two first datasets are
implemented in the Universal support for I/O of meshes.  A single
unv file may hold <i>multiple</i> datasets of type "2414".  To
distinguish data, each dataset has a header.  The <code>libMesh</code> pendant
to this header is the <code>MeshDataUnvHeader</code> class.  When you
read/write unv files using the <code>MeshData</code>, you <i>always</i>
automatically also read/write such headers.
</li>
</ul>
Enough babble for now.  Examples. 
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
Activate the <code>MeshData</code> of the mesh, so that
we can remember the node and element ids used
in the file.  When we do not activate the <code>MeshData</code>,
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
this should be a .unv or .xdr/xda file, otherwise
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
(apart from node & element ids) in it, yet.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
        		<< "Finished reading the mesh.  MeshData is active but empty:" << std::endl
        		<< "---------------------------------------------------------" << std::endl;
              
              mesh.print_info();
              mesh.data.print_info();
            
</pre>
</div>
<div class = "comment">
Create some artificial node-associated data and store
it in the mesh's <code>MeshData</code>.  Use a <code>std::map</code> for this. 
Let the function <code>create_artificial_data()</code> do the work.
</div>

<div class ="fragment">
<pre>
              {
        	std::map<const Node*, std::vector<Number> > artificial_data;
        	
        	create_artificial_data (mesh, artificial_data);
        	
</pre>
</div>
<div class = "comment">
Before we let the map go out of scope, insert it into
the <code>MeshData</code>.  Note that by default the element-associated
data containers are closed, so that the <code>MeshData</code> is
ready for use.
</div>

<div class ="fragment">
<pre>
                mesh.data.insert_node_data(artificial_data);
        
</pre>
</div>
<div class = "comment">
Let <code>artificial_data()</code> go out of scope
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
        		<< "After inserting artificial data into the MeshData:" << std::endl
        		<< "--------------------------------------------------" << std::endl;
              
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
              std::cout << "Writing MeshData to: " << first_out_data << std::endl;
              mesh.data.write(first_out_data);
        
</pre>
</div>
<div class = "comment">
Alternatively, create your own header.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
        		<< "Attach our own MeshDataUnvHeader to the MeshData:" << std::endl
        		<< "-------------------------------------------------" << std::endl
        		<< " (note the warning: the number of values per node in" << std::endl
        		<< "  my_header is not correct)" << std::endl << std::endl;
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
Specify some text that helps the <i>user</i> 
identify the data.  This text is <i>not</i> used for
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
Now attach this header to the <code>MeshData</code>, and write
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
              std::cout << "Writing MeshData to: " << second_out_data << std::endl;
              mesh.data.write(second_out_data);
              
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; std::endl 
        		<< "Before clearing the MeshData:" << std::endl
        		<< "-----------------------------" << std::endl;
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
        		<< "After clearing the MeshData:" << std::endl
        		<< "----------------------------" << std::endl;
              mesh.data.print_info();
        
</pre>
</div>
<div class = "comment">
Now the <code>MeshData</code> is open again to read data from
file.  Read the file that we created first.
</div>

<div class ="fragment">
<pre>
              mesh.data.read(first_out_data);
              
              std::cout << std::endl 
        		<< "After re-reading the first file:" << std::endl
        		<< "--------------------------------" << std::endl;
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
        
            std::cout << std::endl 
        	      << "----------------------------------------------" << std::endl
        	      << "---------- next example with MeshData --------" << std::endl
        	      << "----------------------------------------------" << std::endl;
        
</pre>
</div>
<div class = "comment">
Create a new mesh, read it again -- but this time
with de-activated <code>MeshData</code>.  Then we are
able to use the compatibility mode not only for
handling <code>MeshData</code> dat, but also to <i>write</i> 
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
Read the input mesh, but with deactivated <code>MeshData</code>.
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
        		<< "De-activated MeshData:" << std::endl
        		<< "----------------------" << std::endl;
              mesh.print_info();
              mesh.data.print_info();
         
</pre>
</div>
<div class = "comment">
Write the <i>mesh</i> (not the MeshData!) as .unv file.
In general, the <code>MeshBase</code> interface for .unv I/O
needs an active <code>MeshData</code>.  However, use compatibility
mode to at least write a .unv file, but with the 
ids from libMesh.
</div>

<div class ="fragment">
<pre>
              const std::string out_mesh = "mesh_with_libmesh_ids.unv";
              std::cout << "Writing _Mesh_ to: " << out_mesh << std::endl
        		<< "Try 'diff " << out_mesh << " " << mesh_file << "'" << std::endl
        		<< "to see the differences in node numbers." << std::endl
        		<< "---------------------------------------" << std::endl
        		<< std::endl;
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
        	std::map<const Node*, std::vector<Number> > artificial_data;
        	create_artificial_data (mesh, artificial_data);
        	mesh.data.insert_node_data(artificial_data);
              }
        
</pre>
</div>
<div class = "comment">
Note that even with (only) compatibility mode MeshData
data can be written.  But again, the ids from libMesh
are used.  Consult the warning messages issued in
DEBUG mode.  And the user <i>has</i> to specify that the
<code>MeshData</code> should change to compatibility mode.
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
              std::cout << std::endl 
        		<< "Writing MeshData to: " << mesh_data_file << std::endl
        		<< "----------------------------------------------------------" 
        		<< std::endl << std::endl;
              mesh.data.write (mesh_data_file);
        
        #ifdef HAVE_ZLIB_H
        
</pre>
</div>
<div class = "comment">
As may already seen, UNV files are text-based, so the may
become really big.  When <code>./configure</code> found <code>zlib.h</code>,
then we may also <i>read</i> or <i>write</i> <code>.unv</code> files in gzip'ed 
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
              std::cout << std::endl 
        		<< "Writing gzip'ed MeshData to: " << packed_mesh_data_file << std::endl
        		<< "---------------------------------------------------------------------------" << std::endl
        		<< " To verify the integrity of the packed version, type:" << std::endl << std::endl
        		<< "   gunzip " << packed_mesh_data_file << "; " << std::endl
        		<< "   diff packed_" << mesh_data_file << " " 
        		<< mesh_data_file << std::endl << std::endl;
              
              mesh.data.write (packed_mesh_data_file);
              
        #endif
        
</pre>
</div>
<div class = "comment">
And now a last gimmick: The <code>MeshData::translate()</code>
conveniently converts the nodal- or element-associated
data (currently only nodal) to vectors that may be used 
for writing a mesh with "solution" vectors, where the solution 
vector contains the data from the <code>MeshData</code>.  Particularly
useful for <i>inspecting</i> the data contained in <code> MeshData</code>.

<br><br>And even better: the user can choose the mesh for which
to export the data.  E.g. not only use the <code> mesh</code>
itself, but alternatively use the <code> BoundaryMesh </code>
(any mesh that uses the <i>same nodes</i>, i.e. the <code> Node*</code> 
have to be the same.  Only exception that will not work:
A mesh created using <code> Mesh::create_submesh() </code>
actually will <i>not</i> work with <code> MeshData::translate() </code>).

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
        	std::vector<std::string> data_names;
        	
</pre>
</div>
<div class = "comment">
Use the <code> mesh</code> itself.  Alternatively, use the 
<code> BoundaryMesh</code> of <code> mesh</code>.
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
        	std::cout << std::endl 
        		  << "Writing the data from the MeshData to the GMV file " 
        		  << gmv_file << std::endl
        		  << "------------------------------------------------------------------------" 
        		  << std::endl;
        	
        	GMVIO(mesh).write_nodal_data (gmv_file,
        				      translated_data,
        				      data_names);
        	
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
This function creates the data to populate the <code> MeshData</code> object
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh&amp; mesh,
        			     std::map<const Node*, std::vector<Number> >& art_data)
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
          assert (fabs(z_max-z_min) > TOLERANCE);
        
          const Real x_min = b_box.first (0);
          const Real x_max = b_box.second(0);
          assert (fabs(x_max-x_min) > TOLERANCE);
        
        
          const_node_iterator node_it = mesh.nodes_begin();
          const const_node_iterator node_end = mesh.nodes_end();
        
          for (; node_it != node_end; ++node_it)
            {
</pre>
</div>
<div class = "comment">
All the vectors in <code> artificial</code>_data <i>have</i> to have the
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
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  
  
  
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
  	
  	GMVIO(mesh).write_nodal_data (gmv_file,
  				      translated_data,
  				      data_names);
  	
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
Compiling C++ (in debug mode) ex12.C...
Linking ex12...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
***************************************************************
* Running Example  ./ex12
***************************************************************
 
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
| Creations:    24951
| Destructions: 24951
| 4Node reference count information:
| Creations:    3977
| Destructions: 3977
 ---------------------------------------------------------------------------- 

----------------------------------------------
---------- next example with MeshData --------
----------------------------------------------
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


*************************************************************************
* WARNING: MeshData neither active nor in compatibility mode.           *
*          Enable compatibility mode for MeshData.  Use this Universal  *
*          file with caution: libMesh node and element ids are used.    *
*************************************************************************


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
