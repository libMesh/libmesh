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
      <li><a href="../developers.html">Developers</a></li>
      <li><a href="../installation.html">Installation</a></li>
      <li><a href="../examples.html">Examples</a></li>
     <li><a href="../doxygen/index.html">Documentation</a></li>
    </ul>
  </nav>

  <div id="w">
    <header id="logo"><a href="../index.html"><span id="logobg">SomeWebsiteLogo</span></a></header>

    <nav id="navigation">
      <ul>
        <li><a href="../index.html">Home</a></li>
        <li><a href="../support.html">About Us</a></li>
        <li><a href="../publications.html">Publications</a></li>
        <li><a href="../developers.html">Developers</a></li>
        <li><a href="../installation.html">Installation</a></li>
        <li><a href="../examples.html">Examples</a></li>
        <li><a href="../doxygen/index.html">Documentation</a></li>
      </ul>
    </nav>

<div id="content">

EOF

# First generate the html with the nice comments in it.
for input_file in *.h *.C ; do
    if (test -f $input_file); then
	echo "<a name=\"comments\"></a> "                             >> $output_file
	echo "<br><br><br> <h1> The source file $input_file with comments: </h1> " >> $output_file
	program2html.pl $input_file >> $output_file
    fi
done

# Now generate the code with no comments
for input_file in *.h *.C ; do
    if (test -f $input_file); then
        # Now put some kind of separating message
	echo "<a name=\"nocomments\"></a> "                           >> $output_file
	echo "<br><br><br> <h1> The source file $input_file without comments: </h1> " >> $output_file
        # Now bust out the magic.  We are going to use
        # two perl scripts and enscript to get sweet looking
        # C code into our html.
	stripcomments.pl $input_file | enscript -q --header="" --pretty-print=cpp --color --title="" --language=html --output=- | stripenscript.pl >> $output_file;
    fi
done


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

# # Put in a few emacs comments to force syntax highlighting
# echo "<?php if (0) { ?>" >> $output_file;
# echo "\#Local Variables:" >> $output_file;
# echo "\#mode: html" >> $output_file;
# echo "\#End:" >> $output_file;
# echo "<?php } ?>" >> $output_file;

# Local Variables:
# mode: shell-script
# End:
