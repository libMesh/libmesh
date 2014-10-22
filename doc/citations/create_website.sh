#!/bin/bash

# This script creates a single html file of all the .bib files
# passed in the year_names array below.  Make sure there is a
# corresponding entry in the year_numbers array as well.  Why
# not just name the bib files 2009.bib etc?  For some reason
# multibib can't handle filenames with numbers in them.

# This script requires bibtex2html (see associated Makefile and
# libmesh_html rule)
if [ "x`which bibtex2html`" == "x" ]; then
  echo "bibtex2html is required to create the website!";
  exit 1;
fi

# Declare array of year names and numbers.  Zero-based indexing!
year_names=( libmesh theses preprints fifteen fourteen thirteen twelve eleven ten nine eight seven six five four )
year_numbers=( 'Please use the following citation to reference libMesh' 'Dissertations & Theses' Preprints 2015 2014 2013 2012 2011 2010 2009 2008 2007 2006 2005 2004 )
link_names=( 'skip me' 'Dissertations' 'Preprints' 'Articles' )

# Length of the arrays
N=${#year_names[@]}
# echo "N=$N"

# Remove previous file, if it exists
if [ -f master.html ]; then
  rm master.html
fi

for ((i=0; i < $N ; i++))
do
  # echo $i
  make libmesh_html SOURCE=${year_names[$i]}.bib

  # Add header
  # echo "<h2>${year_numbers[$i]}</h2>" >> master.html

  # Ben added links to some of these headers so they can
  # be accessed from the left navigation pane on the actual
  # site...
  if [ $i -gt 0 ]; then
    if [ $i -lt 4 ]; then
      echo "<a name=\"${link_names[$i]}\"></a>" >> master.html
    fi
  fi

  echo "<h2>${year_numbers[$i]}</h2>" >> master.html

  # Combine
  cat ${year_names[$i]}.html >> master.html
done

# The name of the final file that we will produce
final_file=publications.html

# Create publications.html by adding header/footer
if [ -f $final_file ]; then
  rm $final_file
fi

# This is the new header we use, but it probably shouldn't be
# hard-coded?  Not sure how Ben is auto-generating these now...
echo "<!doctype html>" >> $final_file
echo "<html lang=\"en-US\">" >> $final_file
echo "<head>" >> $final_file
echo "  <meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\">" >> $final_file
echo "  <title>libMesh - A C++ Finite Element Library</title>" >> $final_file
echo "  <meta name=\"author\" content=\"Benjamin S. Kirk\">" >> $final_file
echo "  <link rel=\"stylesheet\" type=\"text/css\" media=\"all\" href=\"styles.css\">" >> $final_file
echo "</head>" >> $final_file
echo "" >> $final_file
echo "<body>" >> $final_file
echo "  <nav id=\"fixedbar\">" >> $final_file
echo "    <ul id=\"fixednav\">" >> $final_file
echo "      <li><a href=\"index.html\">Home</a></li>" >> $final_file
echo "      <li><a href=\"support.html\">About Us</a></li>" >> $final_file
echo "      <li><a href=\"publications.html\">Publications</a></li>" >> $final_file
echo "      <li><a href=\"developers.html\">Developers</a></li>" >> $final_file
echo "      <li><a href=\"installation.html\">Installation</a></li>" >> $final_file
echo "      <li><a href=\"examples.html\">Examples</a></li>" >> $final_file
echo "      <li><a href=\"doxygen/index.html\">Documentation</a></li>" >> $final_file
echo "    </ul>" >> $final_file
echo "  </nav>" >> $final_file
echo "" >> $final_file
echo "  <div id=\"w\">" >> $final_file
echo "    <header id=\"logo\"><a href=\"index.html\"><span id=\"logobg\">SomeWebsiteLogo</span></a></header>" >> $final_file
echo "" >> $final_file
echo "    <nav id=\"navigation\">" >> $final_file
echo "      <ul>" >> $final_file
echo "        <li><a href=\"index.html\">Home</a></li>" >> $final_file
echo "        <li><a href=\"support.html\">About Us</a></li>" >> $final_file
echo "        <li><a href=\"publications.html\">Publications</a></li>" >> $final_file
echo "        <li><a href=\"developers.html\">Developers</a></li>" >> $final_file
echo "        <li><a href=\"installation.html\">Installation</a></li>" >> $final_file
echo "        <li><a href=\"examples.html\">Examples</a></li>" >> $final_file
echo "        <li><a href=\"doxygen/index.html\">Documentation</a></li>" >> $final_file
echo "      </ul>" >> $final_file
echo "    </nav>" >> $final_file
echo "" >> $final_file
echo "<div id=\"content\">" >> $final_file

# Inject all the references
cat master.html >> $final_file

# Inject Ben's footer
echo "<br>" >> $final_file
echo "</div>" >> $final_file
echo "<script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\"></script>" >> $final_file
echo "<script type=\"text/javascript\">" >> $final_file
echo "\$(document).ready(function(){" >> $final_file
echo "  \$(window).on('scroll',function() {" >> $final_file
echo "    var scrolltop = \$(this).scrollTop();" >> $final_file
echo "" >> $final_file
echo "    if(scrolltop >= 215) {" >> $final_file
echo "      \$('#fixedbar').fadeIn(250);" >> $final_file
echo "    }" >> $final_file
echo "" >> $final_file
echo "    else if(scrolltop <= 210) {" >> $final_file
echo "      \$('#fixedbar').fadeOut(250);" >> $final_file
echo "    }" >> $final_file
echo "  });" >> $final_file
echo "});" >> $final_file
echo "</script>" >> $final_file
echo "<!-- Google Analytics stuff -->" >> $final_file
echo "<script type=\"text/javascript\">" >> $final_file
echo "  var _gaq = _gaq || [];" >> $final_file
echo "  _gaq.push(['_setAccount', 'UA-24978333-1']);" >> $final_file
echo "  _gaq.push(['_trackPageview']);" >> $final_file
echo "  (function() {" >> $final_file
echo "    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;" >> $final_file
echo "    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';" >> $final_file
echo "    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);" >> $final_file
echo "  })();" >> $final_file
echo "</script>" >> $final_file
echo "</body>" >> $final_file
echo "</html>" >> $final_file
