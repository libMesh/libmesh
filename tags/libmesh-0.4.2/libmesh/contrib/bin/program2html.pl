# Copyright (C) 1999, 2000, 2001, 2002 by Wolfgang Bangerth, University of Heidelberg
# With modifications for libmesh by John W. Peterson, 2003.

# print "<a name=\"CommProg\"></a>\n";
# print "<h1> The commented program</h1>\n";

# ignore comments at the start of the program. this includes CVS
# tags and copyright notices.
$_ = <>;
while ( m!^/[/\*]!  ||  m/^$/ ) {
    $_ = <>;
}

# have two states, in which the program can be:
# comment-mode and program-mode
$comment_mode = 0;
$program_mode = 1;
$state =  $comment_mode;

print "<div class = \"comment\">\n";

do {
    # substitute special characters
    # s/&/&amp;/g;
#     s/</&lt;/g;
#     s/>/&gt;/g;
#     s/\t/        /g;

    if (($state == $program_mode) && m!^\s*//!)
    {     
	$state = $comment_mode;
	# End the preceeding code div
	print "</pre>\n";
	print "</div>\n";
	
	# Begin a new comment division
	print "<div class = \"comment\">\n";
    }
    
    # if in comment mode and no comment line: toggle state.
    # don't do so, if only a blank line
    elsif (($state == $comment_mode) && !m!^\s*//! && !m!^\s*$!)
    {     
	$state = $program_mode;
	# End the comment division
	print "</div>\n\n";

	# Begin a code fragment div
	print "<div class =\"fragment\">\n";

	# Also start a <pre> block
	print "<pre>\n";

	# substitute special characters
	s/&/&amp;/g;
	s/</&lt;/g;
	s/>/&gt;/g;
	s/\t/        /g;
    }
    
    if ($state == $comment_mode) 
    {
	# in comment mode: first skip leading whitespace and
	# comment // signs
	s!\s*//\s*(.*)\n!$1!;

	# second, replace section headers, and generate addressable
	# anchor...don't need this
	# if ( /\@sect/ ) {
# 	   s!\@sect(\d)\{(.*)\}\s*$!<h$1>$2</h$1>!g;
# 	   $sect_name = $2;
# 	   $sect_name =~ s/\s/_/g;
# 	   $_ = "\n<a name=\"$sect_name\"></a>" . $_;
# 	}

	# finally print this line
	print $_, "\n";

	# if empty line, introduce paragraph break
	print "<br><br>" if  $1 =~ m!^\s*$!;
    }
    else
    {
	print "        $_";
    }
} while (<>);


# Close off the last code section
if ($state == $program_mode) {
    print "</pre>\n";
    print "</div>\n\n";
}

