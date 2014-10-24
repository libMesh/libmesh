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
final_file=../html/src/publications.html

# Create publications.html
if [ -f $final_file ]; then
  rm $final_file
fi

# Inject all the references
cat master.html >> $final_file

# bibtex2html unhelpfully inserts trailing whitespace, so kill that off
perl -pli -e 's/\s+$//' $final_file
