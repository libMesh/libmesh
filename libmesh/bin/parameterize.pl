#!/usr/bin/perl -w

use strict;

my @paramnames;
my @paramtypes;  # = "range", "list", whatever
my @paramlists;  # = [start, interval, stop] for ranges
my @paramprefixes;
my @paramcomments;
my @paramvalcount; # i.e. how many values to instantiate?

# Scan through each config file first to find param* to expand
foreach my $configfile (@ARGV) {
	open(CONFFILE, $configfile) or 
		die "Can't open input file $configfile : $!\n";

	while (<CONFFILE>) {
		if (m/^[ \t]*([^# \t]*)[ \t]*=(.*?)paramrange *\((.*?),(.*?),(.*?)\)(.*)/) {
			push @paramnames, $1;
			push @paramprefixes, $2;
			push @paramtypes, "range";
			my @varlist = ($3, $4, $5);
			push @paramlists, \@varlist;
			push @paramcomments, $6;
			my $valcount = 0;
			for (my $value = $varlist[0]; $value <= $varlist[2] 
			        + $varlist[1]/2; $value += $varlist[1]) {
				$valcount++;
			}
			push @paramvalcount, $valcount;
		} elsif (m/^[ \t]*([^# \t]*)[ \t]*=(.*?)paramlist *\((.*?)\)(.*)/) {
			push @paramnames, $1;
			push @paramprefixes, $2;
			push @paramtypes, "list";
			my @varlist = split(/ *\, */, $3);
			push @paramlists, \@varlist;
			push @paramcomments, $4;
			my $valcount = 0;
			foreach my $value (@varlist) {
				$valcount++;
			}
			push @paramvalcount, $valcount;
		} elsif (m/^[ \t]*([^# \t]*)[ \t]*=(.*?)paramvalue *\((.*?)\)(.*)/) {
			push @paramnames, $1;
			push @paramprefixes, $2;
			push @paramtypes, "value";
			my @varlist = ($3,);
			push @paramlists, \@varlist;
			push @paramcomments, $4;
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
my $totaldepth = 0;
foreach my $valcount (@paramvalcount) {
	$totalcount *= $valcount;
	$totaldepth++;
}

# print "About to create $totalcount parameterized directories, $totaldepth deep.\n";

# Create directories and write config files

sub parameterizeonevariable {
	my $fixednamesref = shift;
	my $fixedvaluesref = shift;
	my $fixedprefixesref = shift;
	my $fixedcommentsref = shift;
	my $varnamesref = shift;
	my $vartypesref = shift;
	my $varlistsref = shift;
	my $varprefixesref = shift;
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
		my $nextprefix = shift @$varprefixesref;
		my $nextcomment = shift @$varcommentsref;
		push @$fixednamesref, $nextname;
		push @$fixedprefixesref, $nextprefix;
		push @$fixedcommentsref, $nextcomment;
		if ($nexttype eq "range") {
			for (my $value = $nextlist[0]; $value <= $nextlist[2] 
			        + $nextlist[1]/2; $value += $nextlist[1]) {
		        	my $finalvalue = eval($nextprefix .
					              '$value ' .
						      $nextcomment);
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
					$fixedvaluesref,
					$fixedprefixesref,
					$fixedcommentsref,
					$varnamesref, $vartypesref,
					$varlistsref, $varprefixesref,
					$varcommentsref);
				pop @$fixedvaluesref;
			}
		} elsif ($nexttype eq "list") {
			foreach my $valuename (@nextlist) {
		        	my $finalvalue = eval($nextprefix .
					              '$value ' .
						      $nextcomment);
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
					$fixedvaluesref,
					$fixedprefixesref,
					$fixedcommentsref,
					$varnamesref, $vartypesref,
					$varlistsref, $varprefixesref,
					$varcommentsref);
				pop @$fixedvaluesref;
			}
		} elsif ($nexttype eq "value") {
			foreach my $valuename (@nextlist) {
				my $value = 0;
				for (my $i = 0;
				     $i != $#$fixednamesref + 1; ++$i) {
					my $fixedname = @$fixednamesref[$i];
					if ($fixedname eq $valuename) {
						$value = @$fixedvaluesref[$i];
					}
				}
		        	my $finalvalue = eval($nextprefix .
						      '$value ' .
						      $nextcomment);
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
					$fixedvaluesref,
					$fixedprefixesref,
					$fixedcommentsref,
					$varnamesref, $vartypesref,
					$varlistsref, $varprefixesref,
					$varcommentsref);
				pop @$fixedvaluesref;
			}
		} else {
			die "Found unexpected type";
		}
		pop @$fixednamesref;
		pop @$fixedprefixesref;
		pop @$fixedcommentsref;
		unshift @$varnamesref, $nextname;
		unshift @$vartypesref, $nexttype;
		unshift @$varlistsref, $nextlistref;
		unshift @$varprefixesref, $nextprefix;
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
				if (m/[ \t]*([^# \t]*)[ \t]*= *paramrange *\((.*),(.*),(.*)\)(.*)/ or 
				    m/[ \t]*([^# \t]*)[ \t]*= *paramlist *\((.*?)\)(.*)/ or
				    m/[ \t]*([^# \t]*)[ \t]*= *paramvalue *\((.*?)\)(.*)/) {
					my $varname = $1;
					if ($varname ne $$fixednamesref[$index]) {
						print "Expected $varname eq $$fixednamesref[$index] at index $index\n";
						die;
					}

					print OUTPUTFILE $varname, " = ",
					$$fixedvaluesref[$index], ";";
				        if ($$fixedcommentsref[$index] ne ";" or
					    $$fixedprefixesref[$index] ne " ") {
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

parameterizeonevariable([], [], [], [], \@paramnames, \@paramtypes,
	\@paramlists, \@paramprefixes, \@paramcomments);
