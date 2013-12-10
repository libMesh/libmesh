<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>Download and Installation</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("download",$root)?>

<div class="content">
<h1>Installation Instructions</h1>


<a name="getsoftware"></a><h2>Getting the Software</h2>

The <code>libMesh</code> source can be <a href="https://github.com/libMesh/libmesh/releases">downloaded from our GitHub release page</a>.
Stable releases are located there as compressed tar archives.

You may also access the Git source tree for the latest code. You can get read-only access
to the Git repository via:
<br>

<div class="fragment">
  <pre>git clone git://github.com/libMesh/libmesh.git </pre>
</div>

<br>
If you would like to contribute to the project you will need a GitHub
developer account, or you can contribute patches. To create a patch from a
modified Git tree simply do:
<br>
<div class="fragment">
  <pre>git diff &gt; patch </pre>
</div>

<br>
in the top-level directory. You can then submit the file <code>patch</code>.

<a name="compilers"></a><h2>Compilers</h2>

<code>libMesh</code> makes extensive use of the standard C++ library,
so you will need a decent, standards-compliant compiler. We have tried
very hard to make the code completely compiler-agnostic by avoiding
questionable (but legal) constructs. If you have a compiler that won't
build the code please let us know. You will also need a decent C compiler
if you want to build some of the contributed packages that add functionality
to the library.

<p>
The library is continuously tested with the following compilers:

<br>
<ul>
  <li>GNU GCC</li>
    <ul>
      <li><code>gcc</code> 4.2</li>
      <li><code>gcc</code> 4.4</li>
      <li><code>gcc</code> 4.5</li>
      <li><code>gcc</code> 4.6</li>
      <li><code>gcc</code> 4.7</li>
    </ul>
  <li>Intel ICC</li>
    <ul>
      <li><code>icc</code> 11.1</li>
      <li><code>icc</code> 12.0</li>
      <li><code>icc</code> 12.1</li>
      <li><code>icc</code> 13.0</li>
    </ul>
  <li>Clang</li>
    <ul>
      <li><code>clang++</code> 2.9</li>
      <li><code>clang++</code> 3.0</li>
      <li><code>clang++</code> 3.1</li>
      <li><code>clang++</code> 3.2</li>
    </ul>
  <li>Sun Studio/Oracle</li>
    <ul>
      <li><code>CC</code> 12.3</li>
      <tt>  $ ./configure --disable-fparser CXXFLAGS=-library=stlport4 --disable-unordered-containers</tt>
    </ul>
  <li> Portland Group</li>
    <ul>
      <li><code>pgi</code> 11.7</li>
      <li><code>pgi</code> 12.9</li>
      <li><code>pgi</code> 13.4</li>
      <tt>  $ ./configure --disable-unordered-containers --disable-fparser --enable-static --disable-shared</tt>
    </ul>
</ul>


<a name="conf"></a><h2>Configuration</h2>

Configuring the library is straightforward. The GNU autoconf package is used
to determine site-specific configuration parameters. A standard build will
occur after typing
<br>

<div class="fragment">
<pre>./configure
make
</pre>
</div>

<br>
in the top-level project directory. To see all the configuration options type
<p>
<div class="fragment">
<pre>./configure --help</pre>
</div>

<br>
The configure script will find your compilers and
create <code>Makefile</code>s with the configuration specific for your
site. If you want to use different compilers than those found by
configure you can specify them in environment variables. For example,
the following will build with the macports <code>Clang</code>
compilers, and also specifies nonstandard search paths for a number of
optional libraries:

<br>
<div class="fragment">
<pre>./configure --prefix=/tmp/foo \
       --with-glpk-include=/opt/local/include \
       --with-glpk-lib=/opt/local/lib \
       --with-vtk-include=/opt/local/include/vtk-5.10 \
       --with-vtk-lib=/opt/local/lib/vtk-5.10 \
       --with-eigen-include=/opt/local/include/eigen3 \
       --with-cxx=clang++-mp-3.2 --with-cc=clang-mp-3.2 --disable-fortran
</div>

<br>
Note that the Fortran compiler is not actually used to compile any part of the library,
but <code>configure</code> uses it to find out how to link Fortran libraries with C++ code, and it is possible to compile <code>libMesh</code> without a Fortran compiler.


<a name="build"></a><h2>Building the Library</h2>

To build the library you need <code>GNU</code> <code>Make</code> and a supported compiler,
as listed in the <a href="installation.php#compilers">Compiler</a> section. After the library
is configured simply type <code>make</code> to build the library.

<br>
The <code>./configure</code> script distributed with the library looks at the shell
environment variable <code>METHODS</code> to determine what modes the library should be built in.
Valid values for <code>METHOD</code> are <code>opt</code> (optimized mode), <code>dbg</code> (build with debug symbols),
and <code>pro</code> (build with profiling support for use with <code>gprof)</code>.
Once the library is configured you can build it simply by typing

<div class="fragment">
<pre>make</pre>
</div>


<a name="test"></a><h2>Testing the Library</h2>
<h3>Running the Examples</h3>
<code>libMesh</code> includes a number of examples in the <code>examples</code>
directory. From the top-level directory you can build and run the example programs
by typing
<div class="fragment">
<pre>make check</pre>
</div>

<br>
<!-- Note that the example programs all create output in the <code>GMV</code> format, -->
<!-- since you can <a href="http://laws.lanl.gov/XCM/gmv/GMVHome.html">download GMV</a> -->
<!-- for free from Los Alamos National Lab. It is a simple matter to change the source -->
<!-- in the example to write a different format, just replace the <code>write_gmv</code> -->
<!-- function call with whatever you like. -->

Note that many of the the example programs create output in the <code>ExodusII</code> format,
since you can <a href="http://www.paraview.org">download Paraview</a>
for free, and it is a highly capable postprocessing tool. It is a simple matter to change the source
in the example to write a different formats, however.


<h3>Unit Tests</h3>
The source tree contains a <code>tests</code> entry in the main trunk
that contains a series of unit tests which can be used to validate a <code>libMesh</code>
installation.  These unit tests require <a href="https://sourceforge.net/apps/mediawiki/cppunit/index.php?title=Main_Page">CPPUnit</a>
to run properly.  To run the unit test suite, simply do

<pre>make -C test check </pre>


<a name="link"></a><h2>Linking With Your Application</h2>

Since <code>libMesh</code> can be configured with many additional packages we recommend
including the <code>Make.common</code> file created in the top-level directory in the
<code>Makefile</code> of any application you want to use with the library. This will
properly set the <code>libmesh_INCLUDE</code> and <code>libmesh_LIBS</code> variables, which you can
append to with your own stuff.

<br>
For testing simple programs you may want to use the <code>libmesh-config</code> script
included in the <code>contrib/bin</code> directory instead of creating a <code>Makefile</code>.
This script may be used to determine the relevant compilation and linking flags
used by <code>libMesh</code>. For example, you could build the application <code>foo</code> from
<code>foo.C</code> like this:
<div class="fragment">
<pre> `libmesh-config --cxx` -o foo foo.C `libmesh-config --cxxflags --include --ldflags`</pre>
</div>


<br>


</div>

<?php make_footer() ?>

</body>
</html>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
