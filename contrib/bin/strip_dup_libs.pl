#!/usr/bin/perl -w

#
# Courtesy of Trilinos v9
#
# This perl script removes duplicate libraries from the right to the left and
# removes duplicate -L library paths from the left to the right
use strict;

my @all_libs = @ARGV;
#
# Move from left to right and remove duplicate -l libraries
#
my @cleaned_up_libs_first;
foreach( reverse @all_libs ) {
        $_ = remove_rel_paths($_);
        if( $_=~/-L/ ) {
                unshift @cleaned_up_libs_first, $_;
        }
        else {
                if( !entry_exists($_,\@cleaned_up_libs_first) ) {
                        unshift @cleaned_up_libs_first, $_;
                }
        }
}

#
# Move from right to left and remove duplicate -Wl linker run paths
#
my @cleaned_up_libs_second;
foreach( @cleaned_up_libs_first ) {
        $_ = remove_rel_paths($_);
        if( !($_=~/-Wl/) ) {
                push @cleaned_up_libs_second, $_;
        }
        elsif( !entry_exists($_,\@cleaned_up_libs_second) ) {
                push @cleaned_up_libs_second, $_;
        }
}

#
# Move from right to left and remove duplicate -L library paths
#
my @cleaned_up_libs_third;
foreach( @cleaned_up_libs_second ) {
        $_ = remove_rel_paths($_);
        if( !($_=~/-L/) ) {
                push @cleaned_up_libs_third, $_;
        }
        elsif( !entry_exists($_,\@cleaned_up_libs_third) ) {
                push @cleaned_up_libs_third, $_;
        }
}

#
# Move system library paths to the end of the link line.  This way,
# there should be less risk of libmesh linking against something in
# /usr/lib by accident.
#
my @cleaned_up_libs_fourth;
my @system_lib_paths;
foreach( @cleaned_up_libs_third ) {
    $_ = remove_rel_paths($_);
    # Discard libary paths starting with /usr.  This includes e.g.
    # /usr/lib and /usr/local/lib, which are assumed to be system
    # library paths.
    if( ($_=~/-L\/usr/) )
    {
        push @system_lib_paths, $_;
    }
    # Discard -Wl,* flags starting with /usr
    elsif( ($_=~/-Wl,.*\/usr/) )
    {
        push @system_lib_paths, $_;
    }
    # Otherwise, keep the path.
    else
    {
        push @cleaned_up_libs_fourth, $_;
    }
}

#
# Push system lib paths onto the end of the list
#
foreach( @system_lib_paths )
{
    push @cleaned_up_libs_fourth, $_;
}

#
# Print the new list of libraries and paths
#
print join( " ", @cleaned_up_libs_fourth );


#
# Subroutines
#
sub entry_exists {
        my $entry = shift; # String
        my $list  = shift; # Reference to an array
        foreach( @$list ) {
                if( $entry eq $_ ) { return 1; }
        }
        return 0;
}
#
sub remove_rel_paths {
        my $entry_in = shift;
        if ($entry_in=~/-L\.\./) {
                return $entry_in;
        }
        my @paths = split("/",$entry_in);
        my @new_paths;
        foreach( @paths ) {
                if( !($_=~/\.\./) ) {
                        push @new_paths, $_;
                }
                else {
                        pop @new_paths
                }
        }
        return join("/",@new_paths);
}
