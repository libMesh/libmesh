#!/bin/bash

# The first argument is the C file name (ex1.C) which
# is going to be converted to html.
# The output file name is implied by it.
input_file=$1;
output_file=`basename $input_file | cut -d"." -f1`.html;

# First generate the html with the nice comments in it.
perl program2html.pl $input_file > $output_file;

# Now put some kind of separating message
echo "<br> <h1> The program without comments: </h1> </br>" >> $output_file;

# Now bust out the magic.  We are going to use
# two perl scripts and enscript to get sweet looking
# C code into our html.
perl stripcomments.pl $input_file | enscript -q --header="" --pretty-print=cpp --color --title="" -W html --output=- | perl stripenscript.pl >> $output_file;