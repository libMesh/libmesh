#!/bin/bash

orig_pwd=`pwd`
base=$1
output_file=$orig_pwd/$base.php;
cd $2

# First we put our special php headings in
echo "<?php \$root=\"\"; ?>" > $output_file;
echo "<?php require(\$root.\"navigation.php\"); ?>" >> $output_file;
echo "<html>" >> $output_file;
echo "<head>" >> $output_file;
echo "  <?php load_style(\$root); ?>" >> $output_file;
echo "</head>" >> $output_file;
echo " " >> $output_file;
echo "<body>" >> $output_file;
echo " " >> $output_file;
echo "<?php make_navigation(\"$base\",\$root)?>" >> $output_file;
echo " " >> $output_file;
echo "<div class=\"content\">" >> $output_file;

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
echo "</div>" >> $output_file;
echo "<?php make_footer() ?>" >> $output_file;
echo "</body>" >> $output_file;
echo "</html>" >> $output_file;

# Put in a few emacs comments to force syntax highlighting
echo "<?php if (0) { ?>" >> $output_file;
echo "\#Local Variables:" >> $output_file;
echo "\#mode: html" >> $output_file;
echo "\#End:" >> $output_file;
echo "<?php } ?>" >> $output_file;

# Local Variables:
# mode: shell-script
# End:
