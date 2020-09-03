#!/usr/bin/perl -w

#
# Courtesy of Trilinos v9
#
# This perl script removes duplicate libraries from the right to the left and
# removes duplicate -L library paths from the left to the right
use strict;
use warnings;

my @all_libs = @ARGV;
#
# Move from left to right and remove duplicate -l libraries
#
my @cleaned_up_libs_first;
foreach( reverse @all_libs ) {
        $_ = remove_rel_paths($_);
        if( ($_=~/^-l/) ) {
                if( !entry_exists($_,\@cleaned_up_libs_first) ) {
                        unshift @cleaned_up_libs_first, $_;
                }
        }
        else {
                unshift @cleaned_up_libs_first, $_;
        }
}

#
# Move from right to left and remove duplicate -Wl linker run paths
#
my @cleaned_up_libs_second;
foreach( @cleaned_up_libs_first ) {
        $_ = remove_rel_paths($_);
        if( ($_=~/-Wl/) ) {
                if ( !entry_exists($_,\@cleaned_up_libs_second) ) {
                        push @cleaned_up_libs_second, $_;
                }
        }
        else {
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
# separate system paths (assumed here to be things in /usr/ or /lib/) from
# non-system paths:
#
my @system_lib_paths;
my @other_lib_paths;
my @link_libraries;
foreach( @cleaned_up_libs_third ) {
    $_ = remove_rel_paths($_);
    if( ($_=~/^-L\/usr/) )
    {
        push @system_lib_paths, $_;
    }
    elsif( ($_=~/^-Wl,[^\/]*\/usr/) )
    {
        push @system_lib_paths, $_;
    }
    elsif( ($_=~/^-L\/lib/) )
    {
        push @system_lib_paths, $_;
    }
    elsif( ($_=~/^-Wl,[^\/]*\/lib/) )
    {
        push @system_lib_paths, $_;
    }
    elsif( ( $_=~/^-L\// ) )
    {
        push @other_lib_paths, $_;
    }
    elsif( ( $_=~/^-Wl,[^\/]*\// ) )
    {
        push @other_lib_paths, $_;
    }
    else
    {
        push @link_libraries, $_;
    }
}

#
# Push non-system linker paths, system linker paths, then library link commands:
#
my @cleaned_up_libs_fourth;
foreach( @other_lib_paths )
{
    push @cleaned_up_libs_fourth, $_;
}

foreach( @system_lib_paths )
{
    push @cleaned_up_libs_fourth, $_;
}

foreach( @link_libraries )
{
    push @cleaned_up_libs_fourth, $_;
}

#
# Print the new list of libraries and paths
#
print join( " ", @cleaned_up_libs_fourth );
print "\n";

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
