#!/usr/bin/perl

# This script currently takes exactly 1 argument which is the name of
# the file to open and read. This script does not actually change any
# lines, it just reports "suspicious" lines so that an actual human
# can fix them, this is to prevent spurious changes from bad regexes
# accidentally triggering massive recompiles. To run this script on
# every .C file in the library, for example, you can use a bash loop:
#
# for i in src/*/*.C; do ./contrib/bin/operator_spacing.pl $i; done
use warnings;
use strict;

# Check for expected number of command line args.
die "Usage: $0 filename\n" if @ARGV != 1;

my $filename = $ARGV[0];

# Open file for reading
open(my $fh, "<", $filename) or die "Cannot open < $filename: $!";

# Code to skip multiline C style comments from:
# http://www.perlmonks.org/?node_id=1132036
my $comment = 0;

# Keep track of what line number we are on so we can report it to the user
my $line_number = 0;

while (my $line = <$fh>) {
    ++$line_number;

    # Skip C-style comments.
    if ($line =~ m#/\*#) {
        $comment = 1;
        next;
    }

    if ($line =~ m#\*/#) {
        $comment = 0;
        next;
    }

    next if $comment;

    # If we made it here, we aren't in the middle of a C-style
    # comment, so check for a C++ style comment.
    if ($line =~ m!^\s*//!) {
        next;
    }

    # Print lines that match "int &x;" for several user-defined types.
    if ($line =~ m!(void|int|double|Real|float) &[a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }

    # Print lines that match "int& x;" for several user-defined types.
    if ($line =~ m!(void|int|double|Real|float)& [a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }
    # Print lines that match "CamelCase& camel;"
    if ($line =~ m![A-Z][A-Z|a-z]+& [a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }
    # Print lines that match "CamelCase &camel;"
    if ($line =~ m![A-Z][A-Z|a-z]+ &[a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }

    # Print lines that match "int *x;" for several user-defined types.
    if ($line =~ m!(void|int|double|Real|float) \*[a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }

    # Print lines that match "int* x;" for several user-defined types.
    if ($line =~ m!(void|int|double|Real|float)\* [a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }

    # Print lines that match "CamelCase *camel;"
    if ($line =~ m![A-Z][A-Z|a-z]+ \*[a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }

    # Print lines that match "CamelCase* camel;"
    if ($line =~ m![A-Z][A-Z|a-z]+\* [a-z]!) {
        print "$filename:$line_number: $line";
        next;
    }
}
