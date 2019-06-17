#! /usr/bin/perl -w
#
# Perl script to map geneomic and transcriptomic positions
# such that the nucleotide conservation values (in genomic coords) can be assigned to their correct transcript position
#

if ($#ARGV < 1)
{
	print "One of more input files missing\n";
	print "Usage perl $0 /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/Resources/all_gff3_MSU7.gff3 SNPConservation_All_3K_Accessions.out > SNPConservation_All_3K_Accessions_transcript_positions_mapped.out\n";
	exit();
}

my @snp_freq;

# Read the input conservation score file with genomic cordinates
open(FH2, "<$ARGV[1]") or die ("Can't open input file $ARGV[1]\n");
while (<FH2>)
{
	#   Chr	Genomic_SNP_pos	A	C	G	T	A/G	C/T	G/C	A/T	G/T	A/C	C/G/T	A/G/T	A/C/T	A/C/G	Anybase	Gap
	#   Chr1	1026	0.53	0	0.00	0	0	0	0	0	0	0	0	0	0	0	0	0.47
	#   Chr1	1033	0.56	0	0.00	0	0	0	0	0	0	0	0	0	0	0	0	0.44
	chomp $_;
	my $inp = $_;
	$inp =~ s/^\s+//g;
	$inp =~ s/\s+/\t/g;
	push(@snp_freq, $inp);
}
close FH2;


print "Chr\tStrand\tTranscript\tFeature\tGenomic_SNP_position\tTranscript_SNP_position\tA\tC\tG\tT\tA/G\tC/T\tG/C\tA/T\tG/T\tA/C\tC/G/T\tA/G/T\tA/C/T\tA/C/G\tAnybase\tGap\n";

my $transcript_length = 0;
my $prev_id = "Start";
my @coords;

# Read gff3 file
open(FH1, "<$ARGV[0]") or die ("Can't open input file $ARGV[0]\n");
LINE: while(<FH1>)
{
	#Chr1	MSU_osa1r7	gene	309459	313170	.	-	.	ID=LOC_Os01g01620;Name=LOC_Os01g01620;Note=kinase%2C%20pfkB%20family%2C%20putative%2C%20expressed
	#Chr1	MSU_osa1r7	gene	319713	322248	.	+	.	ID=LOC_Os01g01640;Name=LOC_Os01g01640;Note=Rer1%20protein%2C%20putative%2C%20expressed

	chomp $_;

	my @inp = split(/\t/, $_);
	my $candidate_chr = $inp[0];
	my $feature = $inp[2];
	
	goto LINE if( ($feature eq "gene") || ($feature eq "mRNA") || ($feature eq "exon") || ($feature =~ m/^$/g) );

	my $feature_start = $inp[3];
	my $feature_end = $inp[4];
	my $feature_length = $feature_end - $feature_start+1;
	my $strand = $inp[6]; 
	my $id = $inp[8];
	$id =~ s/^ID=//g;
	$id =~ s/[:;].*//g;


	foreach my $val(@snp_freq)
	{
		my @line = split(/\t/, $val);
		my $snpfreq_chr = $line[0];
		my $snp = $line[1];
		
		#   Chr	Genomic_SNP_pos	A	C	G	T	A/G	C/T	G/C	A/T	G/T	A/C	C/G/T	A/G/T	A/C/T	A/C/G	Anybase	Gap
		my $fraction_A = $line[2];
		my $fraction_C = $line[3];
		my $fraction_G = $line[4];
		my $fraction_T = $line[5];
		my $fraction_AorG = $line[6];
		my $fraction_CorT = $line[7];
		my $fraction_GorC = $line[8];
		my $fraction_AorT = $line[9];
		my $fraction_GorT = $line[10];
		my $fraction_AorC = $line[11];
		my $fraction_CorGorT = $line[12];
		my $fraction_AorGorT = $line[13];
		my $fraction_AorCorT = $line[14];
		my $fraction_AorCorG = $line[15];
		my $fraction_Anybase = $line[16];
		my $fraction_Gap = $line[17];


		if ($candidate_chr eq $snpfreq_chr)
		{
			if($snp >= $feature_start && $snp <= $feature_end)
			{
				my $nucleotide_transcript_pos;

				if($id ne $prev_id)
				{
					$transcript_length=0;
				}
				
				# Minu strand
				if($strand eq "-")
				{
					$nucleotide_transcript_pos = $feature_end-$snp+1;	
					$nucleotide_transcript_pos = $nucleotide_transcript_pos+$transcript_length;	
				}
				# Plus strand
				else
				{
					$nucleotide_transcript_pos = $snp-$feature_start+1;
					$nucleotide_transcript_pos = $nucleotide_transcript_pos+$transcript_length;
				}

				print "$candidate_chr\t$strand\t$id\t$feature\t$snp\t$nucleotide_transcript_pos\t$fraction_A\t$fraction_C\t$fraction_G\t$fraction_T\t$fraction_AorG\t$fraction_CorT\t$fraction_GorC\t$fraction_AorT\t$fraction_GorT\t$fraction_AorC\t$fraction_CorGorT\t$fraction_AorGorT\t$fraction_AorCorT\t$fraction_AorCorG\t$fraction_Anybase\t$fraction_Gap\n";
			}
		}
	}
	
	if( ($prev_id eq "Start") || ($id eq $prev_id) )
	{
		$transcript_length = $transcript_length + $feature_length;
		$prev_id = $id;
	}
	else
	{
		$transcript_length = $feature_length;
		$prev_id = $id;
	}
	
} 
