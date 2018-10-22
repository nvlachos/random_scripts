#!/usr/bin/perl

# Written by Tom de Man

use strict;
use warnings;

my $sam_path = $ARGV[0];

my @sams = &getSam;
&main(@sams);

sub main {
        foreach my $file (@_) {
        	my @read_dbs;
        	open FILE, "$file";
        	while (<FILE>) {
        		chomp;
        		my @split_list = split(' ', $_);
				if ($split_list[2] eq '*') {
        			push @read_dbs, $split_list[0];
				}
			}
			close FILE;

			my $amount_of_reads = scalar @read_dbs;
			print "amount of unmapped reads in the SAM $file: $amount_of_reads\n";

		}
}

sub getSam {
        my @sam_docs;
        opendir(DIR, $sam_path) or die "Cannot open $sam_path \n";
        my @sam_files = readdir(DIR);
        close DIR;

        foreach my $sam_file (@sam_files) {
                next if (!($sam_file =~ /\.sam$/));
                push @sam_docs, $sam_file;
        }
        return @sam_docs;
}
         
