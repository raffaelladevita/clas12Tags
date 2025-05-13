use strict;
use warnings;

our %configuration;


sub materials
{
	# uploading the mat definition
	
	# Scintillator
	my %mat = init_mat();
	$mat{"name"}          = "dcgas";
	$mat{"description"}   = "clas12 dc gas";
	$mat{"density"}       = "0.0018";
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "G4_Ar 0.9 G4_CARBON_DIOXIDE 0.1";
	print_mat(\%configuration, \%mat);

	%mat = init_mat();
	$mat{"name"}          = "dcg10";
	$mat{"description"}   = "clas12 dc region 2 frame g10";
	$mat{"density"}       = "1.7";
	$mat{"ncomponents"}   = "4";
	$mat{"components"}    = "G4_Si 0.2805 G4_O 0.3945 G4_C 0.2480 G4_H 0.0770";
	print_mat(\%configuration, \%mat);

}


