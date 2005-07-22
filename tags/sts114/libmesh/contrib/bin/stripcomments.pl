# Modifications to program2plain originally
# Copyright (C) 1999, 2000, 2001, 2002 by Wolfgang Bangerth, University of Heidelberg
# by John W. Peterson.
# Now it doesn't put any html into the output.

# ignore comments at the start of the program. this includes CVS
# tags and copyright notices.
$_ = <>;
while ( m!^/[/\*]!  ||  m/^$/ ) {
    $_ = <>;
}



do {
    # ignore comment lines
    if ( ! m!^\s*//! ) {
	print "  $_";
    }
} while (<>);


