# This perl script strips out the junk that enscript puts before and
# after the colorized code that it generates.  It cuts out everything
# before <PRE> and after </PRE>.  This way the output can be stuck
# into an existing html file without worrying about having too many
# </html> commands!
#
# Author: John W. Peterson

# Some nice variables to tell us where to start and stop
$found_begin = 0;
$found_end   = 0;

# You can put in your own header if you want.  This should
# probably be commented out!
# print "<html> \n";

do {

    if ( m!\</PRE\>! ) {
	$found_end = 1;
	print "</pre> \n";
    }

    if ( ($found_begin == 1) && ($found_end == 0)) {
	print "$_"
	}
    
    if ( m!\<PRE\>! ) {
	$found_begin = 1;
	print "<pre> \n";
    }
    
} while (<>);

# Now print your own footer if you want.  You probably
# want to leave this since we are now at the end of the file.
# print "</html> \n";

