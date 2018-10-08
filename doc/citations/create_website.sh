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

# Process command line arguments.
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
count_only=false
while [[ $# -ge 1 ]]
do
    key="$1"
    case $key in
        -c|--count-only)
            count_only=true
            shift # past argument
            ;;
        *)
            # unknown option
            shift
            ;;
    esac
done

# Declare array of year names and numbers.  Zero-based indexing!
year_names=( libmesh preprints theses nineteen eighteen seventeen sixteen fifteen fourteen thirteen twelve eleven ten nine eight seven six five four )
year_numbers=( 'Please use the following citation to reference libMesh' 'Preprints' 'Dissertations & Theses' 2019 2018 2017 2016 2015 2014 2013 2012 2011 2010 2009 2008 2007 2006 2005 2004 )
link_names=( 'skip me' 'Preprints' 'Dissertations' 'Articles' )
# The short_numbers are used for making Python strings of the years.
short_numbers=( 'na' 'na' \'T\' \'\\\'19\' \'\\\'18\' \'\\\'17\' \'\\\'16\' \'\\\'15\' \'\\\'14\' \'\\\'13\' \'\\\'12\' \'\\\'11\' \'\\\'10\' \'\\\'09\' \'\\\'08\' \'\\\'07\' \'\\\'06\' \'\\\'05\' \'\\\'04\')

# Stores strings to be output in reverse
output_strings=()

# Length of the arrays
N=${#year_names[@]}

# Remove previous file, if it exists
if [ -f master.html ]; then
  rm master.html
fi

for ((i=0; i < $N ; i++))
do
  # Store the bibtex2html output in $output_filename so we can comb through it later.
  input_filename=${year_names[$i]}.bib
  output_filename=${year_names[$i]}.b2h
  make libmesh_html SOURCE=$input_filename &> $output_filename

  # This strips the number of entries from the bibtex2thml
  # output. There is probably a better, cleaner way to do this with
  # perl, but I have not bothered to look it up yet.
  n_entries=`cat $output_filename | grep -m 1 "entries)" | cut -d"(" -f 2 | cut -d" " -f 1`
  output_strings+=("${short_numbers[$i]}, $n_entries,")

  # Clean up the temporary file.
  rm $output_filename

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

# Create publications.html
if [ "$count_only" = false ] ; then
    # The name of the final file that we will produce
    final_file=../html/src/publications.html

    if [ -f $final_file ]; then
        rm $final_file
    fi

    # Inject all the references
    cat master.html >> $final_file

    # bibtex2html unhelpfully inserts trailing whitespace, so kill that off
    perl -pli -e 's/\s+$//' $final_file
fi

# Print the citation count results by year, in reverse, and we don't
# care about the first two categories.
for ((i=$[$N-1]; i > 1 ; i--))
do
    echo ${output_strings[$i]}
done
