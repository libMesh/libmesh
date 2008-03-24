#!/usr/bin/perl -w

use strict;

my @paramnames;
my @paramtypes;  # = "range", "list", whatever
my @paramlists;  # = [start, interval, stop] for ranges
my @paramcomments;
my @paramvalcount; # i.e. how many values to instantiate?

# Scan through each config file first to find param* to expand
foreach my $configfile (@ARGV) {
	open(CONFFILE, $configfile) or 
		die "Can't open input file $configfile : $!\n";

	while (<CONFFILE>) {
		if (m/[ \t]*([^# \t]*)[ \t]*= *paramrange *\((.*?),(.*?),(.*?)\)(.*)/) {
			push @paramnames, $1;
			push @paramtypes, "range";
			my @varlist = ($2, $3, $4);
			push @paramlists, \@varlist;
			push @paramcomments, $5;
			my $valcount = 0;
			for (my $value = $varlist[0]; $value <= $varlist[2] 
			        + $varlist[1]/2; $value += $varlist[1]) {
				$valcount++;
			}
			push @paramvalcount, $valcount;
		} elsif (m/[ \t]*([^# \t]*)[ \t]*= *paramlist *\((.*?)\)(.*)/) {
			push @paramnames, $1;
			push @paramtypes, "list";
			my @varlist = split(/ *\, */, $2);
			push @paramlists, \@varlist;
			push @paramcomments, $3;
			my $valcount = 0;
			foreach my $value (@varlist) {
				$valcount++;
			}
			push @paramvalcount, $valcount;
		}
	}

	close CONFFILE;
}

# How many directories are we about to make?
my $totalcount = 1;
foreach my $valcount (@paramvalcount) {
	$totalcount *= $valcount;
}

print "About to create $totalcount parameterized directories.\n";

# Create directories and write config files

sub parameterizeonevariable {
	my $fixednamesref = shift;
	my $fixedvaluesref = shift;
	my $fixedcommentsref = shift;
	my $varnamesref = shift;
	my $vartypesref = shift;
	my $varlistsref = shift;
	my $varcommentsref = shift;

	my $fixedsubdir = "parameterized/";
	mkdir $fixedsubdir;
	opendir SUBDIR, $fixedsubdir or die "Can't mkdir $fixedsubdir\n";
	closedir SUBDIR;
	my $index = 0;

	foreach my $fixedname (@$fixednamesref) {
		my $fixedvalue = $$fixedvaluesref[$index++];
		$fixedvalue =~ s/ /_/g;
		$fixedvalue =~ s/\'/_/g;
		$fixedsubdir .= $fixedname . "_" . $fixedvalue . "/";
	}

	if (scalar(@$varnamesref) > 0) {
		my $nextname = shift @$varnamesref;
		my $nexttype = shift @$vartypesref;
		my $nextlistref = shift @$varlistsref;
		my @nextlist = @$nextlistref;
		my $nextcomment = shift @$varcommentsref;
		push @$fixednamesref, $nextname;
		push @$fixedcommentsref, $nextcomment;
		if ($nexttype eq "range") {
			for (my $value = $nextlist[0]; $value <= $nextlist[2] 
			        + $nextlist[1]/2; $value += $nextlist[1]) {
		        	my $finalvalue = eval('$value ' . $nextcomment);
				push @$fixedvaluesref, $finalvalue;
				$finalvalue =~ s/ /_/g;
				$finalvalue =~ s/\'/_/g;
				my $newsubdir = $fixedsubdir . $nextname . "_"
			  	. $finalvalue . "/";
				mkdir $newsubdir;
			       	opendir NEWSUBDIR, $newsubdir
				  or die "Can't mkdir $newsubdir\n";
				closedir NEWSUBDIR;
				parameterizeonevariable($fixednamesref,
					$fixedvaluesref, $fixedcommentsref,
					$varnamesref, $vartypesref,
					$varlistsref, $varcommentsref);
				pop @$fixedvaluesref;
			}
		} elsif ($nexttype eq "list") {
			foreach my $value (@nextlist) {
		        	my $finalvalue = eval('$value ' . $nextcomment);
				push @$fixedvaluesref, $finalvalue;
				$finalvalue =~ s/ /_/g;
				$finalvalue =~ s/\'/_/g;
				my $newsubdir = $fixedsubdir . $nextname . "_"
			  	. $finalvalue . "/";
				mkdir $newsubdir;
			       	opendir NEWSUBDIR, $newsubdir
				  or die "Can't mkdir $newsubdir\n";
				closedir NEWSUBDIR;
				parameterizeonevariable($fixednamesref,
					$fixedvaluesref, $fixedcommentsref,
					$varnamesref, $vartypesref,
					$varlistsref, $varcommentsref);
				pop @$fixedvaluesref;
			}
		} else {
			die "Found unexpected type";
		}
		pop @$fixednamesref;
		pop @$fixedcommentsref;
		unshift @$varnamesref, $nextname;
		unshift @$vartypesref, $nexttype;
		unshift @$varlistsref, $nextlistref;
		unshift @$varcommentsref, $nextcomment;
	} else {
		my $index = 0;
		foreach my $configfile (@ARGV) {
			open(INPUTFILE, $configfile) or 
				die "Can't open input file $configfile : $!\n";

			my $configbasename = $configfile;
			$configbasename =~ s#.*/(.*?)#$1#;

			open(OUTPUTFILE, ">$fixedsubdir$configbasename") or 
				die "Can't open output file $fixedsubdir$configbasename : $!\n";

			while (<INPUTFILE>) {
				if (m/[ \t]*([^# \t]*)[ \t]*= *paramrange *\((.*),(.*),(.*)\)(.*)/ or m/[ \t]*([^# \t]*)[ \t]*= *paramlist *\((.*?)\)(.*)/) {
					my $varname = $1;
					if ($varname ne $$fixednamesref[$index]) {
						print "Expected $varname eq $$fixednamesref[$index] at index $index\n";
						die;
					}

					print OUTPUTFILE $varname, " = ",
					$$fixedvaluesref[$index], ";";
				        if ($$fixedcommentsref[$index] ne ";") {
						print OUTPUTFILE " #", $$fixedcommentsref[$index];
					}
					print OUTPUTFILE "\n";
					$index++;
				} else {
					print OUTPUTFILE $_;
				}
			}

			close OUTPUTFILE;
			close INPUTFILE;
		}
	}
}

parameterizeonevariable([], [], [], \@paramnames, \@paramtypes,
	\@paramlists, \@paramcomments);
