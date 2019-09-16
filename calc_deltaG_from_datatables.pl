#!/usr/bin/perl -w

use strict();
use warnings();
use POSIX;

# Enthlapy (dH) and Free Energy (dG) calculated at 37degC were printed out using RNAstructure package (David Matthews).
# The Tscale() function in src/rna_library.h was modified to prinout dH and dG values as well as the scaled dG value at the new temperature
#
#
# NOTE: The dH values in rna.coaxial.dh and dG valuess in rna.coaxial.dg


if($#ARGV < 0)
{
	print "Input file or parameter missing\n";
	print "Usage: perl $0 <CT_files_101Nts.list>  > output\n";
	exit();
}

# The Enthalpy (dH) calculated at 37deg originally from UV melt experiments can be accessed at RNAstructure/data_tables/rna.coaxial.dh
# The description of the data_tables are in data_tables/description.txt
# Enthalpy values in rna.coaxial.dh are multiplied by 10 within the RNAstructure program as printed in dG_T at src/rna_library.cpp called in dg.stack[a][b][c][d] = Tscale()
# Later on the dGscaled value is divided by conversion factor 10. 
my %enthalpy_dH=(
	"AU_followedby_AU" => -68.0,
	"AU_followedby_CG" => -114.0,
	"AU_followedby_GC" => -105.0,
	"AU_followedby_UA" => -94.0,
	"AU_followedby_GU" => -32.0,
	"AU_followedby_UG" => -88.0,
	"CG_followedby_AU" => -104.0,
	"CG_followedby_CG" => -134.0,
	"CG_followedby_GC" => -106.0,
	"CG_followedby_UA" => -105.0,
	"CG_followedby_GU" => -56.0,
	"CG_followedby_UG" => -121.0,
	"GC_followedby_AU" => -124.0,
	"GC_followedby_CG" => -149.0,
	"GC_followedby_GC" => -134.0,
	"GC_followedby_UA" => -114.0,
	"GC_followedby_GU" => -83.0,
	"GC_followedby_UG" => -126.0,
	"GU_followedby_AU" => -128.0,
	"GU_followedby_CG" => -126.0,
	"GU_followedby_GC" => -121.0,
	"GU_followedby_UA" => -88.0,
	"GU_followedby_GU" => -135.0,
	"GU_followedby_UG" => -146.0,
	"UA_followedby_AU" => -77.0,
	"UA_followedby_CG" => -124.0,
	"UA_followedby_GC" => -104.0,
	"UA_followedby_UA" => -68.0,
	"UA_followedby_GU" => -70.0,
	"UA_followedby_UG" => -128.0,
	"UG_followedby_AU" => -70.0,
	"UG_followedby_CG" => -83.0,
	"UG_followedby_GC" => -56.0,
	"UG_followedby_UA" => -32.0,
	"UG_followedby_GU" => -93.0,
	"UG_followedby_UG" => -135.0
);

# The Free Energy (dG) calculated at 37deg originally from UV melt experiments can be accessed at RNAstructure/data_tables/rna.coaxial.dg
# The description of the data_tables are in data_tables/description.txt
# Free Energy (dG) values in rna.coaxial.dg are multiplied by 10 within the RNAstructure program as printed in dG_T at src/rna_library.cpp called in dg.stack[a][b][c][d] = Tscale()
# Later on the dGscaled value is divided by conversion factor 10. 
my %data_table_dG=(
	"AU_followedby_AU" => -9.0,
	"AU_followedby_CG" => -22.0,
	"AU_followedby_GC" => -21.0,
	"AU_followedby_UA" => -11.0,
	"AU_followedby_GU" => -6.0,
	"AU_followedby_UG" => -14.0,
	"CG_followedby_AU" => -21.0,
	"CG_followedby_CG" => -33.0,
	"CG_followedby_GC" => -24.0,
	"CG_followedby_UA" => -21.0,
	"CG_followedby_GU" => -14.0,
	"CG_followedby_UG" => -21.0,
	"GC_followedby_AU" => -24.0,
	"GC_followedby_CG" => -34.0,
	"GC_followedby_GC" => -33.0,
	"GC_followedby_UA" => -22.0,
	"GC_followedby_GU" => -15.0,
	"GC_followedby_UG" => -25.0,
	"GU_followedby_AU" => -13.0,
	"GU_followedby_CG" => -25.0,
	"GU_followedby_GC" => -21.0,
	"GU_followedby_UA" => -14.0,
	"GU_followedby_GU" => -5.0,
	"GU_followedby_UG" => +13.0,
	"UA_followedby_AU" => -13.0,
	"UA_followedby_CG" => -24.0,
	"UA_followedby_GC" => -21.0,
	"UA_followedby_UA" => -9.0,
	"UA_followedby_GU" => -10.0,
	"UA_followedby_UG" => -13.0,
	"UG_followedby_AU" => -10.0,
	"UG_followedby_CG" => -15.0,
	"UG_followedby_GC" => -14.0,
	"UG_followedby_UA" => -6.0,
	"UG_followedby_GU" => +3.0,
	"UG_followedby_UG" => -5.0
);	

#Since the input CT files were folded at 28degC, the temperature here is set to 301.15K (28degC)
my $temperature = 301.15;

#foreach my $key(keys %enthalpy_dH)
#{
#	print ("Base step: $key\t enthalpy (dH): $enthalpy_dH{$key}\tfree_energy_data_table (dG): $data_table_dG{$key}\n");
#}

print("File\tBase_Pair_Stacks\tGiven_Temp1_Kelvin\tdH_fromDataTables(base/steps)\tdG_fromDataTables(base/steps)\tdG_TempScaled(base/steps)\tdG_TempScaled_total\tEstimatedTemp_Kelvin_to_openBase_to_make_scaled_deltaG_zero: GivenTemperature/(1-(total_dG_data_tables_unconverted/total_dH_data_tables_unconverted))\tTemperature_difference\n");


# Read CT file
# if 51st position is paired as well as 50th and 52nd.
#Read list of input CT files
open(FH1, "<$ARGV[0]") or die("Can't open input file $ARGV[0]\n");
NEXT_CTFILE: while(<FH1>)
{
	chomp $_;
	my $input_CT_file = $_;
	my @ct_data;
	my %ct_hash;
	my @pairs;

	open(FH2, "<$input_CT_file") or die("Can't open input CT file\n");
	my $ct_header = <FH2>;
	while(<FH2>)
	{
		chomp $_;
	  #101 ENERGY =   -12.6    SHAPE_CONSTR_101Nt/m1.9b-0.5_LOC_Os01g01640.1_SNP3_Nipponbare_LOC_Os01g01640.1.Nipponbare
	  #  1 G       0    2    0    1
	  #  2 A       1    3    0    2
	  #  3 U       2    4    0    3
	  # 50 U      49   51   64   50
	  # 51 U      50   52   63   51
	  # 52 G      51   53   62   52
		$_=~ s/^\s+//g;
		$_=~ s/\s+/\t/g;
		@ct_data = split(/\t/, $_);
		# make a hash of nucleotide number and its sequence
		$ct_hash{$ct_data[0]} = $ct_data[1];

		# consider only the 50th, 51st and 52nd bases
		if( ($ct_data[0] == 50) || ($ct_data[0] == 51) || ($ct_data[0] == 52) )
		{
			# check if the 3 bases are paired
			# if atleast one of them is unpaired then don't consider this file and read the next CT file
			if($ct_data[4] == 0) { goto NEXT_CTFILE }
			else
			{
				# get their respective pairing nucleotides
				push(@pairs, join(":", $ct_data[0], $ct_data[4]));
			}
		}
	}
	my $bp_before_snppos;
	my $bp_snppos;
	my $bp_after_snppos;

	foreach my $val (@pairs)
	{
		my @pairing_partners = split(/:/, $val);
		if($pairing_partners[0] == 50)
		{
			$bp_before_snppos = join('',$ct_hash{50}, $ct_hash{$pairing_partners[1]});
		}
		if($pairing_partners[0] == 51)
		{
			$bp_snppos = join('',$ct_hash{51}, $ct_hash{$pairing_partners[1]});
		}
		if($pairing_partners[0] == 52)
		{
			$bp_after_snppos = join('',$ct_hash{52}, $ct_hash{$pairing_partners[1]});
		}
	}

	my $first_basestep = $bp_before_snppos."_followedby_".$bp_snppos;
	my $second_basestep= $bp_snppos."_followedby_".$bp_after_snppos;


	# Scaled_dG = dH - floor( [ (dH - datatable_dG)*temperature/310.15 ]+0.5 )
	# after this divide by conversion factor
	my $scaled_dG_first_base_step =   $enthalpy_dH{$first_basestep} - floor( ($enthalpy_dH{$first_basestep} - $data_table_dG{$first_basestep})*$temperature/310.15 +0.5 );
	my $scaled_dG_second_base_step =  $enthalpy_dH{$second_basestep} - floor( ($enthalpy_dH{$second_basestep} - $data_table_dG{$second_basestep})*$temperature/310.15 +0.5);

	# conversion factor used in RNAstructure. defined in RNAstructure/src/defines.h
	my $conversion_factor= 10;
	#Adding up the scaled dG and dH and dG from datatables for the three basestacks of the structure at a given temeprature. 
	my $total_dG_scaled_basesteps = ($scaled_dG_first_base_step + $scaled_dG_second_base_step)/$conversion_factor;

	#------------------------------------
	#my $dH_data_tables_converted_first_base_step = $enthalpy_dH{$first_basestep}/$conversion_factor;
	#my $dH_data_tables_converted_second_base_step = $enthalpy_dH{$second_basestep}/$conversion_factor;

	#my $dG_data_tables_converted_first_base_step = $data_table_dG{$first_basestep}/$conversion_factor;
	#my $dG_data_tables_converted_second_base_step = $data_table_dG{$second_basestep}/$conversion_factor;
	#------------------------------------

	my $total_dH_data_tables_unconverted = $enthalpy_dH{$first_basestep} + $enthalpy_dH{$second_basestep};
	my $total_dG_data_tables_unconverted = $data_table_dG{$first_basestep} + $data_table_dG{$second_basestep};

	my $total_dH_data_tables_converted = $total_dH_data_tables_unconverted/$conversion_factor;
	my $total_dG_data_tables_converted = $total_dG_data_tables_unconverted/$conversion_factor;


	my $estimated_temp;
	#
	#
	#Yiliang's idea (derivation see my notes: red book) 
	$estimated_temp = $temperature/(1-($total_dG_scaled_basesteps/$total_dH_data_tables_unconverted));

	my $temp_difference = $estimated_temp - $temperature;	
	
	print("$input_CT_file\t5' $bp_before_snppos/$bp_snppos/$bp_after_snppos 3'\t$temperature\t$enthalpy_dH{$first_basestep}/$enthalpy_dH{$second_basestep}\t$data_table_dG{$first_basestep}/$data_table_dG{$second_basestep}\t$scaled_dG_first_base_step/$scaled_dG_second_base_step\t$total_dG_scaled_basesteps\t$estimated_temp\t$temp_difference\n");
}
