#! /usr/bin/perl -w
#
# Perl script to extract a slice of the genome form two reference genomes: Nipponbare and the newly
# generated Alternate reference for 93-11. The input is a VCF file that has Chromosome number (ex. Chr1.vcf),
# and the SNP position and the paths to two reference genome files. The script slices 15 nucleotides 
# upstream and downstream of the SNP in both genomes and outputs the sequences for comparison.
#

use strict;
use Bio::SeqIO;


if ($#ARGV < 2)
{
	print "Input file missing\n";
	print "Usage: perl $0 <input_vcf_coordinate_file> <Genome1.fasta> <Genome2.fasta> > outfile\n";
	print "       perl $0 FILTERED-SNPs-forref.vcf 93-11_AltReferenceGenome.fasta Nipponbare_ref_genome.fasta > genome_slices.out\n";
	exit(0);
}

open (FH1, "$ARGV[0]") or die("Can't open input file $ARGV[0]\n");
while (<FH1>)
{
	chomp $_;
	my $input = $_;

	# Skip comment lines
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  9mNAI_B1        9mNAI_B2
	#Chr1    149929  .       T       C       70.69   PASS    AC=4;AF=1.00;AN=4;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=60.00;QD=10.10;SOR=4.174 GT:AD:DP:GQ:PL  1/1:0,4:4:10:38,10,0    1/1:0,3:3:9:59,9,0	
	
	if ($input !~ m/^#/g)
	{
		my @values = split(/\t/, $input);
		my $chr = $values[0];
		my $snp_position = $values[1];

		my @genome_slice_1 = &getGenomeSlice($chr, $snp_position, $ARGV[1]);
		my @genome_slice_2 = &getGenomeSlice($chr, $snp_position, $ARGV[2]);

		my $subseq_1 = $genome_slice_1[0];
		my $subseq_2 = $genome_slice_2[0];

		my $upstream_pos_1 = $genome_slice_1[1];
		my $downstream_pos_1 = $upstream_pos_1+30;

		my $upstream_pos_2 = $genome_slice_2[1];
		my $downstream_pos_2 = $upstream_pos_2+30;

		my $file_name_1 = $ARGV[1];
		my $file_name_2 = $ARGV[2];

		$file_name_1 =~ s/.*\///g;	
		$file_name_2 =~ s/.*\///g;	

		$subseq_1 =~ m/\[(\w)\]/g;
		my $snp1 = $1;
		$subseq_2 =~ m/\[(\w)\]/g;
		my $snp2 = $1;


		my $filter;
		if ($snp1 ne $snp2)
		{
			$filter = "PASSED";
		}
		else
		{
			$filter = "FAILED";
		}

		print "$file_name_1\t$chr\t$upstream_pos_1\t$subseq_1\t$downstream_pos_1\n";
		print "$file_name_2\t\t\t$chr\t$upstream_pos_2\t$subseq_2\t$downstream_pos_2\t$filter\n\n";
	}
}

sub getGenomeSlice
{
	my ($chr, $snp_position, $file_path) = @_;
	
	my $in  = Bio::SeqIO->new(-file => "$file_path",-format => 'Fasta');
	my $seq_object;
	my $seq_slice;
	my $index_pos;

	while ( $seq_object = $in->next_seq() )
	{
		my $Chr_id = $seq_object->display_id(); 
	
		if ($Chr_id eq $chr)
		{
			#$seq->seq();
			my $seq = $seq_object->seq();
			
			# perl substr() function offset for getting upstream sequence
			# from desired postion.
			# Here the offset is set at -16 positions upstream from the snp
			# position and from there 15 downstream nucleotides are selected.
			#
			my $minus_offset = $snp_position-15-1;
					
			my $minus_15 = substr($seq, $minus_offset, 15);
			my $plus_15 = substr($seq, $snp_position, 15);
			my $snp_position_nucleotide = substr($seq, $snp_position-1, 1);

			my $seq_string = $minus_15.$snp_position_nucleotide.$plus_15;
			# Index position is the start position of the substring within the chromosome sequence.
			$index_pos = index($seq, $seq_string);
			$index_pos = $index_pos+1;

			$seq_slice = $minus_15." [".$snp_position_nucleotide."] ".$plus_15;
		}
	}
	my @return_string;
	push (@return_string, $seq_slice);
	push (@return_string, $index_pos);

	return (@return_string);
}
