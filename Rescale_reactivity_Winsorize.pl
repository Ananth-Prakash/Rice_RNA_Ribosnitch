#! /usr/bin/perl -w
#
# Perl script to rescale reactivity values between 0 to 2
# by Winsorization
# Normalise all other raw reactivity values by the 95th percentile value
# and mulitply by 2, so that all rescaled values are in the range 0 to 2.
#

use strict();
use List::Util qw(min max);

if ($#ARGV < 1)
{
	print "Input list file or output directory path missing\n";
	print "Usage: perl $0 <Raw_reactivity_files.list> <Directory/path/for/out/file>\n";
	exit();
}

my @file_list;

open (FH1, "<$ARGV[0]") or die ("Can't open input file\n");
while (<FH1>)
{
	chomp $_;
	push (@file_list, $_);
}
close FH1;

foreach my $file(@file_list)
{
	open (FH2, "<$file") or die ("Can't open reactivity file: $file\n");
	my @nucl;
	my @reac;

	while(<FH2>)
	{
		chomp $_;
		
		my @inp = split(/\s+/,$_);
		
		push(@nucl, $inp[0]);	
		push(@reac, $inp[1]); 
	}
	close FH2;
	
	my @sorted_reactivity = sort {$a <=> $b} @reac;

	my $NinetyFifth_percentile_value = &percentile(95, \@sorted_reactivity);

	my $winsorized_reactivity = &winsorize($NinetyFifth_percentile_value, \@reac);

	#print "Raw array\n";
	#print "@reac\n";
	#print "Sorted raw array\n";
	#print "@sorted_reactivity\n";
	#print "95th percentile value: $NinetyFifth_percentile_value\n";
	#print "Winsorized array\n";
	#print "@$winsorized_reactivity\n";

	my $rescaled_reactivity = &rescale($winsorized_reactivity);

	#print "Rescaled reactivity\n";
	#printf ("@$rescaled_reactivity\n\n");

	$file =~ s/.*\///g;

	open (FW1, ">$ARGV[1]/$file.rescaled");

	for (my $i=0; $i<= $#nucl; $i++)
	{
		print FW1 ("$nucl[$i]  ");
		printf FW1 ("%.4f",@$rescaled_reactivity[$i]);
		print FW1 "\n";
	}
	close FW1;
}

sub percentile
{
	my($percentile, $sorted_reactivity) = @_;
	
	my $index = int($percentile * $#{$sorted_reactivity}/100);

	my $NinetyFifth_percentile_value = @$sorted_reactivity[$index];

	return($NinetyFifth_percentile_value);	

}

sub winsorize
{
	my ($NinetyFifth_percentile_value, $reac) = @_;

	my @winsorized_reactivity;
	foreach my $val(@$reac)
	{
		my $raw_reactivity = $val;
		if ($raw_reactivity > $NinetyFifth_percentile_value)
		{
			$raw_reactivity = $NinetyFifth_percentile_value;
		}

		push (@winsorized_reactivity, $raw_reactivity);
	}

	return (\@winsorized_reactivity);
}

sub rescale
{
	my ($winsorized_reactivity) = @_;

	my @rescaled_reactivity;

	my $max_winsorized_value = max(@$winsorized_reactivity);

	#print "Maximum winsorized value: $max_winsorized_value\n";

	foreach my $val(@$winsorized_reactivity)
	{
		my $winsorized_reactivity = $val;
		my $rescaled_val = ($winsorized_reactivity/$max_winsorized_value)*2;
		
		push(@rescaled_reactivity, $rescaled_val);
	}

	return(\@rescaled_reactivity);
}
