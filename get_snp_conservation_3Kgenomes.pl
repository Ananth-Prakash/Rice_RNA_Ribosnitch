#!/usr/bin/perl

# Script to compute the fraction (or frequency) of each nucleotide at SNP position in the multiallelic SNP dataset
# File downloaded from :http://snp-seek.irri.org/_download.zul
# Dataset:	Full 3K RG SNPs Dataset (includes multi-allelic SNPs)
#		https://s3.amazonaws.com/3kricegenome/reduced/3k_RG_32mio_All_multiallelic_biallelic_SNP_dataset.zip

use strict();
use warnings();
use IO::Handle;

if($#ARGV < 2)
{
	print "One or more input files missing\n";
	print "Usage: perl $0 Mapped_Rice3K_accessions_varieties_types.out 3kall_snpposition_map.tsv Universe_matrix_geno_NB\n";
	exit();
}

# Read RiceAcession variety index and group
# variety index is useful later, since this maps to genomes in Universe_matrix_geno_NB file.
#
my %group_type;

my $time = localtime;
print STDERR "Started processing at $time\n";

open(FH1, "<$ARGV[0]") or die("Can't open input file $ARGV[0]\n");
my $variety_header = <FH1>;
while(<FH1>)
{
	chomp $_;
#	Accession	Group	VARIETY_INDEX	NAME	BOX_CODE	GS_ACCESSION
#	IRIS 313-15896	Indica	1	FEDEARROZ 50G1	AB01	126957
#	IRIS 313-15897	Indica	2	SANHUANGZHAN NO 2G1	AB02	126968
#	IRIS 313-15898	Indica	3	IR 4630-22-2-5-1-3G1	AB03	126962
	my @inp_variety = split(/\t/, $_);
	my $accession = $inp_variety[0];
	my $group = $inp_variety[1];
	my $variety_index = $inp_variety[2];
	my $variety_name = $inp_variety[3];

	$group_type{$group}{$variety_index} = ();
}
close FH1;

$time = localtime;
print STDERR "Finished reading file1 $ARGV[0] at $time\n";



# Read genomic SNP positions
my @SNP_genomic_positions;
open(FH2, "<$ARGV[1]") or die("Can't open input file $ARGV[1]\n");
my $snppos_header = <FH2>;
while(<FH2>)
{
	chomp $_;

#	SNP_INDEX	CHROMOSOME	POSITION	REFCALL
#	0	1	1026	A
#	1	1	1033	A
	my @snp_inp = split(/\t/, $_);
	my $annotation = join("\t", $snp_inp[1], $snp_inp[2]);
	push(@SNP_genomic_positions, $annotation);
}
close FH2;

$time = localtime;
print STDERR "Finished reading file2 $ARGV[1] at $time\n";



# Read Universal SNPs for 3K genomes
#
my %SNP_column_3K;
my %SNP_column_Indica_subgroup;
my %SNP_column_Japonica_subgroup;
my %SNP_column_TropicalJaponica_subgroup;
my %SNP_column_TemperateJaponica_subgroup;

my $All_accession_counts = 0;
my $Indica_subgroup_accession_counts=0;
my $Japonica_subgroup_accession_counts=0;
my $TemperateJaponica_subgroup_accession_counts=0;
my $TropicalJaponica_subgroup_accession_counts=0;

# Read Universe SNP file
# Read snps at all positions within one accession
# each line is an accession, therefore reading one accession at a time.
#
#
open(FH3, "<$ARGV[2]") or die ("Can't open input file $ARGV[2]\n");
while(<FH3>)
{
	chomp $_;
	my @snp_universe_inp = split(/\t/, $_);
	my $snp_universe_accession = $snp_universe_inp[0];

	my $i=0;
	my $k=0;
	$All_accession_counts++;
	#Coloumn-wise counting each variation
	# (I) Within ALL 3K genomes

	map {$SNP_column_3K{$i++}{$_}, $SNP_column_3K{$k++}{$_}++} split(//, $snp_universe_inp[1]);

	# (II) Within Indica subgroup accessions only
	if(exists($group_type{"Indica"}{$snp_universe_accession}))
	{
		$i=0;
		$k=0;
		$Indica_subgroup_accession_counts++;
		map {$SNP_column_Indica_subgroup{$i++}{$_}, $SNP_column_Indica_subgroup{$k++}{$_}++} split(//, $snp_universe_inp[1]);
	}

	# (III) Within Japonica subgroup accessions only
	if(exists($group_type{"Japonica"}{$snp_universe_accession}))
	{
		$i=0;
		$k=0;
		$Japonica_subgroup_accession_counts++;
		map {$SNP_column_Japonica_subgroup{$i++}{$_}, $SNP_column_Japonica_subgroup{$k++}{$_}++} split(//, $snp_universe_inp[1]);
	}

	# (IV) Within Tropical Japonica subgroup accessions only
	if(exists($group_type{"Tropical japonica"}{$snp_universe_accession}))
	{
		$i=0;
		$k=0;
		$TropicalJaponica_subgroup_accession_counts++;
		map {$SNP_column_TropicalJaponica_subgroup{$i++}{$_}, $SNP_column_TropicalJaponica_subgroup{$k++}{$_}++} split(//, $snp_universe_inp[1]);
	}
	
	# (V) Within Temperate Japonica subgroup accessions only
	if(exists($group_type{"Temperate japonica"}{$snp_universe_accession}))
	{
		$i=0;
		$k=0;
		$TemperateJaponica_subgroup_accession_counts++;
		map {$SNP_column_TemperateJaponica_subgroup{$i++}{$_}, $SNP_column_TemperateJaponica_subgroup{$k++}{$_}++} split(//, $snp_universe_inp[1]);
	}

$time = localtime;
print STDERR "Done reading accession: $snp_universe_accession at $time\n";

}
close FH3;

$time = localtime;
print STDERR "Finished reading file3 $ARGV[2] at $time\n";
print STDERR "Started calculation of SNP conservation\n";



#print "3K\n";
&print_SNPConservation(\%SNP_column_3K, $All_accession_counts, \@SNP_genomic_positions, "SNPConservation_All_3K_Accessions.out");
$time = localtime;
print STDERR "Finished SNP conservation caclulation for 3K Genomes at $time\n";

#print "Indica subgroup only\n";
&print_SNPConservation(\%SNP_column_Indica_subgroup, $Indica_subgroup_accession_counts, \@SNP_genomic_positions, "SNPConservation_Indica_Subgroup_only.out");
$time = localtime;
print STDERR "Finished SNP conservation caclulation for Indica Subgroup at $time\n";

#print "Japonica subgroup only\n";
&print_SNPConservation(\%SNP_column_Japonica_subgroup, $Japonica_subgroup_accession_counts, \@SNP_genomic_positions, "SNPConservation_Japonica_Subgroup_only.out");
$time = localtime;
print STDERR "Finished SNP conservation caclulation for Japonica Subgroup at $time\n";

#print "Japonica Tropical subgroup only\n";
&print_SNPConservation(\%SNP_column_TropicalJaponica_subgroup, $TropicalJaponica_subgroup_accession_counts, \@SNP_genomic_positions, "SNPConservation_TropicalJaponica_Subgroup_only.out");
$time = localtime;
print STDERR "Finished SNP conservation caclulation for Tropical Japonica Subgroup at $time\n";

#print "Japonica Temperate subgroup only\n";
&print_SNPConservation(\%SNP_column_TemperateJaponica_subgroup, $TemperateJaponica_subgroup_accession_counts, \@SNP_genomic_positions, "SNPConservation_TemperateJaponica_Subgroup_only.out");
$time = localtime;
print STDERR "Finished SNP conservation caclulation for Temperate Japonica Subgroup at $time\n";


sub print_SNPConservation
{
	my ($hash_ref_SNP_column, $number_of_accessions, $array_ref_snp_genomic_positions, $outfile) = @_;

	my %SNP_column = %$hash_ref_SNP_column;
	my @SNP_genomic_positions = @$array_ref_snp_genomic_positions;
	my $snp_pos_count=0;
	
	#IUPAC nomenclature: https://www.bioinformatics.org/sms/iupac.html
	my $fractionA=0;
	my $fractionC=0;
	my $fractionG=0;
	my $fractionT=0;
	my $fractionR=0; #A or G
	my $fractionY=0; #C or T
	my $fractionS=0; #G or C
	my $fractionW=0; #A or T
	my $fractionK=0; #G or T
	my $fractionM=0; #A or C
	my $fractionB=0; #C or G or T
	my $fractionD=0; #A or G or T
	my $fractionH=0; #A or C or T
	my $fractionV=0; #A or C or G
	my $fractionN=0; #any base
	my $fraction_questionmark=0; #gap

	open (FW1, ">./$outfile") or die("Can't open output file $outfile\n");
	{
		print FW1 "Chr\tGenomic_SNP_pos\tA\tC\tG\tT\tA/G\tC/T\tG/C\tA/T\tG/T\tA/C\tC/G/T\tA/G/T\tA/C/T\tA/C/G\tAnybase\tGap\n";

		foreach my $nucl_index(sort {$a<=>$b} keys %SNP_column)
		{
			#print "$SNP_genomic_positions[$snp_pos_count]\n";
	
			foreach my $nucl_snp(keys %{$SNP_column{$nucl_index}})
			{
				#print "$nucl_index\t$nucl_snp\t$SNP_column{$nucl_index}{$nucl_snp}\n";
				if($nucl_snp eq "A"){$fractionA = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "C"){$fractionC = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "G"){$fractionG = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "T"){$fractionT = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "R"){$fractionR = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "Y"){$fractionY = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "S"){$fractionS = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "W"){$fractionW = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "K"){$fractionK = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "M"){$fractionM = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "B"){$fractionB = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "D"){$fractionD = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "H"){$fractionH = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "V"){$fractionV = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "N"){$fractionN = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				if($nucl_snp eq "?"){$fraction_questionmark = sprintf("%.2f", $SNP_column{$nucl_index}{$nucl_snp}/$number_of_accessions)}
				#print "$nucl_snp\t$SNP_column{$nucl_index}{$nucl_snp}\t";
			}
			print FW1 "Chr$SNP_genomic_positions[$snp_pos_count]\t$fractionA\t$fractionC\t$fractionG\t$fractionT\t$fractionR\t$fractionY\t$fractionS\t$fractionW\t$fractionK\t$fractionM\t$fractionB\t$fractionD\t$fractionH\t$fractionV\t$fractionN\t$fraction_questionmark\n";
			$snp_pos_count++;
			$fractionA=0;
			$fractionC=0;
			$fractionG=0;
			$fractionT=0;
			$fractionR=0;
			$fractionY=0;
			$fractionS=0;
			$fractionW=0;
			$fractionK=0;
			$fractionM=0;
			$fractionB=0;
			$fractionD=0;
			$fractionH=0;
			$fractionV=0;
			$fractionN=0;
			$fraction_questionmark=0;
		}
	}
	close FW1;
}
