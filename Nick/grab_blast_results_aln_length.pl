#! /usr/bin/perl

# Written by Tom de Man

use strict;

#my $minlen = 140000;
my $minlen = 300;
my %querycovs;
my $blast_path = shift;

my @blasts = &get_files("blast");
#sleep(15);
foreach (@blasts) {
	#print "$_\n";
	my $blastf = $_;
	my $totallen = 0;
	my $bp = "$blast_path$_";
	my $biggest_contig=0;
	#print "$bp";
	open FILE, "$bp";
	while (<FILE>) {
		chomp;
		my @split = split('\t', $_);
		if ($split[3] > $biggest_contig) {
			$biggest_contig=$split[3];
		}
		if ($split[3] > $minlen) {
			#print "Adding $split[3]";
			#sleep(10);
			$totallen += $split[3];
		}
		else {
			#print "Not Adding $split[3]";
		}
	}
	close FILE;
	my $totallen_rel = ($totallen/178154)*100;
	$totallen_rel = "$totallen_rel	$biggest_contig";
	$querycovs{$blastf} = $totallen_rel;
}

my @filesblast = keys(%querycovs);
foreach my $element (@filesblast) {
	#print "$element\n";
}

foreach my $el (sort hashSort (keys(%querycovs))) {
	print "$querycovs{$el}\t$el\n";
}

sub hashSort {
	$querycovs{$b} <=> $querycovs{$a};
}

sub get_files {
	my $ext = qr/$_[0]/;
	my @blastfiles;
	opendir(DIR, $blast_path) or die "cannot open $blast_path \n";
	my @files = readdir(DIR);
	close DIR;
	
	foreach my $file (@files){
		#print "gf-$file";
		next if (!($file =~ /\.$ext$/));
		push @blastfiles, $file;
	}
	return @blastfiles;
}
