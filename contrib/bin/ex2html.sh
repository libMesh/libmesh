#!/bin/bash

orig_pwd=`pwd`
base=$1
output_file=$orig_pwd/$base.html
cd $2

# First we put our special headings in
cat <<EOF > $output_file
<!doctype html>
<html lang="en-US">
<head>
  <meta http-equiv="Content-Type" content="text/html;charset=utf-8">
  <title>libMesh - A C++ Finite Element Library</title>
  <meta name="author" content="libMesh development team">
  <link rel="stylesheet" type="text/css" media="all" href="../styles.css">
  <link rel="stylesheet" type="text/css" media="all" href="../doxygen_stylesheet.css">
</head>

<body>
  <nav id="fixedbar">
    <ul id="fixednav">
      <li><a href="../index.html">Home</a></li>
      <li><a href="../support.html">About Us</a></li>
      <li><a href="../publications.html">Publications</a></li>
      <li><a href="https://github.com/libMesh/libmesh/graphs/contributors">Developers</a></li>
      <li><a href="../installation.html">Installation</a></li>
      <li><a href="../examples.html">Examples</a></li>
     <li><a href="https://mooseframework.inl.gov/docs/doxygen/libmesh/index.html">Documentation</a></li>
    </ul>
  </nav>

  <div id="w">
    <header id="logo"><a href="../index.html"><span id="logobg">SomeWebsiteLogo</span></a></header>

    <nav id="navigation">
      <ul>
        <li><a href="../index.html">Home</a></li>
        <li><a href="../support.html">About Us</a></li>
        <li><a href="../publications.html">Publications</a></li>
        <li><a href="https://github.com/libMesh/libmesh/graphs/contributors">Developers</a></li>
        <li><a href="../installation.html">Installation</a></li>
        <li><a href="../examples.html">Examples</a></li>
        <li><a href="https://mooseframework.inl.gov/docs/doxygen/libmesh/index.html">Documentation</a></li>
      </ul>
    </nav>

<div id="content">

EOF

# Put a link to the lastest version of this example on GitHub.  This
# is much easier to maintain than what we were doing previously with
# perl, enscript, etc.

# Directory name minus the foo_ex1 bit.
one_path_stripped=$(dirname "$orig_pwd")

# Now basename returns just the last directory name, which for
# $one_path_stripped is just the name of ".."
last_two_paths=$(basename "$one_path_stripped")/$(basename "$orig_pwd")

echo "<br> <h1> Link to the source code for this example: </h1>" >> $output_file
echo "<a href=\"https://github.com/libMesh/libmesh/tree/master/examples/$last_two_paths\" target=\"_blank\">Open $base in new tab.</a>" >> $output_file

# Now add the stdout.log if it exists
cd $orig_pwd
if (test -f stdout.log); then
    echo "<a name=\"output\"></a> "                           >> $output_file;
    echo "<br><br><br> <h1> The console output of the program: </h1> " >> $output_file;
    echo "<pre>" >> $output_file;
    cat stdout.log >> $output_file;
    echo "</pre>" >> $output_file;
fi;

# Now put our special php footer in
cat <<EOF  >> $output_file
</div>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>
<script type="text/javascript">
\$(document).ready(function(){
  \$(window).on('scroll',function() {
    var scrolltop = \$(this).scrollTop();

    if(scrolltop >= 215) {
      \$('#fixedbar').fadeIn(250);
    }

    else if(scrolltop <= 210) {
      \$('#fixedbar').fadeOut(250);
    }
  });
});
</script>

<!-- Google Analytics stuff -->
<script type="text/javascript">
  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-24978333-1']);
  _gaq.push(['_trackPageview']);
  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();
</script>

</body>
</html>
EOF

# Local Variables:
# mode: shell-script
# End:
