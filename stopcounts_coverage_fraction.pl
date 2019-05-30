#!/usr/bin/perl -w
#
# Perl script to calaculate the fraction of stop counts and coverage for each transcript
#

use strict;
use List::Util qw(sum);

if ($#ARGV < 0)
{
	print "\tInput count file missing\n";
	print "\tUsage: perl $0 <.cnt file> > output_file\n";
	exit(0);
}


open (FH1, "$ARGV[0]") or die ("Can't open input file\n");
my $line_count=1;
my $id;
my $stop_counts;
my $coverage;

print ("ID\tStop_count_sum\tCoverage_sum\tLength\tStop_count_sum/Length\tCoverage_sum/Length\n");

while (<FH1>)
{
	chomp $_;
	if ($line_count ==1)
	{
		$id = $_;
		$line_count++;
		print "$id\t";
	}
	elsif ($line_count == 2)
	{
		$stop_counts=$_;
		$line_count++;
		print "Stpc $stop_counts\n";
	}
	elsif ($line_count == 3)
	{
		$coverage = $_;
		$line_count++;
		print "Coverage $coverage\n";
	}
	elsif ($line_count==4)
	{
		my @stpc = split(/\t/,$stop_counts);
		my @cov = split(/\t/,$coverage);
		
		my $length = scalar(@stpc);
		
		my $stpcounts = join(',', @stpc);
		#print "$length\t$stpcounts\n";
		
		my $stpcount_sum = sum(0,@stpc);
		my $coverage_sum = sum(0,@cov);

		my $stp_frac = sprintf("%.3f", $stpcount_sum/$length);
		my $cov_frac = sprintf("%.3f", $coverage_sum/$length);
		print("$id\t$stpcount_sum\t$coverage_sum\t$length\t$stp_frac\t$cov_frac\n");
		$line_count =1;
	}
}
		
