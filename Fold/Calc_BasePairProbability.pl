#!/usr/bin/perl -w
#
# Perl script to extract base pair accessibility values for all base pairs in *_rss.ps file
# then convert to Base Pairing Probability (BPP) = 1-accesibility
#

use strict();

if($#ARGV < 0)
{
	print "\tList of input *_rss.ps files missing\n";
	print "\tUsage: perl $0 path_to_fold_rss.ps_files.list > base_pair_probability.out\n";
	exit();
}


print "Transcript_id\tBase_pair_probability\n";

open(FH1, "<$ARGV[0]") or die ("Can't open input file $ARGV[0]\n");
while (<FH1>)
{
	chomp $_;

	my $input_file_path = $_;

	#./Relplot/UNCONSTR/UNCONSTR_LOC_Os02g30630.1_SNP10_Nipponbare_LOC_Os02g30630.1.Nipponbare_rss.ps	
	#./Relplot/UNCONSTR/UNCONSTR_LOC_Os02g30630.1_SNP10_9311_LOC_Os02g30630.1.9311_rss.ps	

	$input_file_path =~ m/(LOC_\S+SNP\d+)/g;
	my $id = $1;

	my $start_block=0;
	my $end_block=0;
	my @access=();

	open(FN, "<$input_file_path") or die ("Can't open input file $input_file_path \n");
	LINE: while (<FN>)
	{
		chomp $_;

		# Beginning of base pair accessibility field
		if($_ =~ m/^\/S \[/g)
		{
			$start_block=1;
			next LINE;
		}

		# End of base pair accessibility field
		if( ($start_block) && ($_ =~ /^\] def/g) )
		{
			$end_block=1;
		}

		if( ($start_block) && ($end_block == 0) )
		{
			$_ =~ s/\s+//g;
			push(@access, $_);
		}
	}

	# Calculate base pairing probability
	my @bpp=();
	foreach my $acc(@access)
	{
		my $bpp_value = 1-$acc;
		push(@bpp, sprintf("%.5f", $bpp_value));
	}
	
	my $bpp_string = join(',', @bpp); 
	
	print "$id\t$bpp_string\n";

}
