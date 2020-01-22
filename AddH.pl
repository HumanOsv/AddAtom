#!/usr/bin/perl -s


#  Cordero, B.; Gómez, V.; Platero-Prats, A. E.; Revés, M.; Echeverría, J.; Cremades, E.; Barragán, F.; Alvarez, S.
#  Covalent Radii Revisited. J. Chem. Soc. Dalt. Trans. 2008, No. 21, 2832–2838. DOI: 10.1039/B801115J
my $other_element = 0.8;

my %Atomic_radii = ( 'H'  => '0.31', 'He' => '0.28', 'Li' => '1.28', 'Be' => '0.96',
                     'B'  => '0.84', 'C'  => '0.76', 'N'  => '0.71', 'O'  => '0.66',
                     'F'  => '0.57', 'Ne' => '0.58', 'Na' => '1.66', 'Mg' => '1.41',
                     'Al' => '1.21', 'Si' => '1.11', 'P'  => '1.07', 'S'  => '1.05',
                     'Cl' => '1.02', 'Ar' => '1.06', 'K'  => '2.03', 'Ca' => '1.77',
                     'Sc' => '1.70', 'Ti' => '1.60', 'V'  => '1.53', 'Cr' => '1.39',
                     'Mn' => '1.39', 'Fe' => '1.32', 'Co' => '1.26', 'Ni' => '1.24',
                     'Cu' => '1.32', 'Zn' => '1.22', 'Ga' => '1.22', 'Ge' => '1.20',
                     'As' => '1.19', 'Se' => '1.20', 'Br' => '1.20', 'Kr' => '1.16',
                     'Rb' => '2.20', 'Sr' => '1.95', 'Y'  => '1.90', 'Zr' => '1.75',
                     'Nb' => '1.64', 'Mo' => '1.54', 'Tc' => '1.47', 'Ru' => '1.46',
                     'Rh' => '1.42', 'Pd' => '1.39', 'Ag' => '1.45', 'Cd' => '1.44',
                     'In' => '1.42', 'Sn' => '1.39', 'Sb' => '1.39', 'Te' => '1.38',
                     'I'  => '1.39', 'Xe' => '1.40', 'Cs' => '2.44', 'Ba' => '2.16',
                     'La' => '2.07', 'Ce' => '2.04', 'Pr' => '2.03', 'Nd' => '2.01',
                     'Pm' => '1.99', 'Sm' => '1.98', 'Eu' => '1.98', 'Gd' => '1.96',
                     'Tb' => '1.94', 'Dy' => '1.92', 'Ho' => '1.92', 'Er' => '1.89',
                     'Tm' => '1.90', 'Yb' => '1.87', 'Lu' => '1.87', 'Hf' => '1.75',
                     'Ta' => '1.70', 'W'  => '1.62', 'Re' => '1.51', 'Os' => '1.44',
                     'Ir' => '1.41', 'Pt' => '1.36', 'Au' => '1.36', 'Hg' => '1.32',
                     'Tl' => '1.45', 'Pb' => '1.46', 'Bi' => '1.48', 'Po' => '1.40',
                     'At' => '1.50', 'Rn' => '1.50', 'Fr' => '2.60', 'Ra' => '2.21',
                     'Ac' => '2.15', 'Th' => '2.06', 'Pa' => '2.00', 'U'  => '1.96',
                     'Np' => '1.90', 'Pu' => '1.87', 'Am' => '1.80', 'Cm' => '1.69' );

my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',
		      '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',
		      '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',
		      '107' => 'Bh', '5'   => 'B',  '35'  => 'Br', '48'  => 'Cd',
		      '20'  => 'Ca', '98'  => 'Cf', '6'   => 'C',  '58'  => 'Ce',
		      '55'  => 'Cs', '17'  => 'Cl', '27'  => 'Co', '29'  => 'Cu',
		      '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
		      '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',
		      '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',
		      '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',
		      '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',
		      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',
		      '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
		      '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',
		      '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',
		      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',
		      '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',
		      '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',
		      '76'  => 'Os', '8'   => 'O',  '46'  => 'Pd', '47'  => 'Ag',
		      '78'  => 'Pt', '82'  => 'Pb', '94'  => 'Pu', '84'  => 'Po',
		      '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',
		      '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',
		      '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
		      '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
		      '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',
		      '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',
		      '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',
		      '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
		      '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
		      '30'  => 'Zn', '40'  => 'Zr' );

#
sub read_cart_file{
	my ($input_file)=@_;
	my @data=();
	open(FILE, "<", $input_file) or die "Can't open";
	my @lines=<FILE>;
	close(FILE);
	foreach $i (@lines){
		chomp($i);
		$data[++$#data]=$i;
	}
	return @data;
}
#
sub gen_xyz {
	my ($radioAct,$Box_x, $Box_y, $Box_z ) = @_;
	# generate a random number in perl in the range box size
	my $lower_limit_x = ( 0.15 + $radioAct + $Box_x ) * -1;
	my $upper_limit_x = ( 0.15 + $radioAct + $Box_x );
	my $lower_limit_y = ( 0.15 + $radioAct + $Box_y ) * -1;
	my $upper_limit_y = ( 0.15 + $radioAct + $Box_y );
	my $lower_limit_z = ( 0.15 + $radioAct + $Box_z ) * -1;
	my $upper_limit_z = ( 0.15 + $radioAct + $Box_z );
	#
	my $x = (rand($upper_limit_x-$lower_limit_x) + $lower_limit_x);
	my $y = (rand($upper_limit_y-$lower_limit_y) + $lower_limit_y);
	my $z = (rand($upper_limit_z-$lower_limit_z) + $lower_limit_z);
	#
	my $x_coord = sprintf '%.6f', $x;
	my $y_coord = sprintf '%.6f', $y;
	my $z_coord = sprintf '%.6f', $z;
	my $coords  = "$x_coord  $y_coord  $z_coord";
	return $coords;
}
###################################
# Fisher-Yates Algorithm
sub fisher_yates_shuffle {
	my $array = shift;
	my $i = @$array;
	while ( --$i ) {
		my $j = int rand( $i+1 );
		@$array[$i,$j] = @$array[$j,$i];
	}
	return @$array;
}
###################################
# Euclidean distance between points
sub Euclidean_distance {
	# array coords basin 1 and basin 2
	my ($p1,$p2,$p3, $axis_x, $axis_y, $axis_z) = @_;
	# variables
	my $x1 = $axis_x;
	my $y1 = $axis_y;
	my $z1 = $axis_z;
	# measure distance between two point
	my $dist = sqrt(
					($x1-$p1)**2 +
					($y1-$p2)**2 +
					($z1-$p3)**2
					);
	return $dist;
}
###################################
# Verification
sub verification_opt {
	my ($a1, $a2, $dist)=@_;
	# hash values
	my $sum = 0.7;
	my $resultado;
	# steric effects if radio1+radio2 < distance
	if($dist <= $sum){
		# Steric problem
		$resultado = 1;
	}else{
		$resultado = 0;
	}
	return $resultado;
}


my @tmparr = (H,H,H,H,H,H);

my ($file) = @ARGV;
if (not defined $file) {
	die "\nCambio_Formato must be run with:\n\nUsage:\n\tCambio_Formato.pl [file-all.xyz]\n";
	exit(1);
}
#read and parse
my @data  = read_cart_file($file);
my $word  = $data[0];
#
my @array_energy = ();
my @linearray_1  = ();
my @linearray_2  = ();
my $count = 0;
my $lines = 0;
#
foreach my $i (@data) {
	if ( $i == $word ) {
		my @array_tabs = ();
		@array_tabs    = split (' ',$i);
		push (@linearray_1,$count);
	}
	$count++;
}
#
#
for (my $x=0; $x < 4 ; $x++) {
#
#
for (my $i=0; $i < scalar(@linearray_1) ; $i++) {
	my $number_atoms = $linearray_1[1] - ($linearray_1[0] + 2);
	my $tmpsum_1     = $linearray_1[$i] + 2;
	my $tmpsum_2     = $linearray_1[$i] + $number_atoms + 1;
	my $count        = 0;
	my @array_2      = ();
	for my $x ($tmpsum_1 .. $tmpsum_2) {
		push (@array_2,$data[$x]);
	}

	my $boolean      = 0;
	while ( $boolean < 1 ) {
		#
		my @newCoords  = ();
		#
		my @total = ();
		#print "24\n\n";
		#
		for my $lol (@array_2) {
			my ($elem,$axisX,$axisY,$axisZ) = split (' ',$lol);
			#
			my $atomic_val;
			if ( exists $Atomic_number{$elem} ) {
				# exists
				$atomic_val = $Atomic_number{$elem};
			} else {
				# not exists
				$atomic_val = $elem;
			}
			#
			#print "$atomic_val  $axisX  $axisY  $axisZ\n";
			push (@total,"$atomic_val  $axisX  $axisY  $axisZ");
			#
			my $radioAct = $Atomic_radii{$atomic_val} || $other_element;
			#print "$atomic_val  $axisX  $axisY  $axisZ\n";
			my $box = gen_xyz ($radioAct,$axisX,$axisY,$axisZ);
			push (@newCoords, $box);
		}
		# Ver impedimento esterico
		# Creacion coordenadas
		my @newCoords_Fi = fisher_yates_shuffle (\@newCoords);
		for ( my $i = 0 ; $i < scalar (@tmparr) ; $i = $i + 1 ){
			#print "$tmparr[$i]  $newCoords_Fi[$i]\n";
			push (@total,"$tmparr[$i]  $newCoords_Fi[$i]");
		}
		#
		#
		#
		#get size
		my $final_trial = 0;
		my $resultado   = 0;
		#
		for (my $i=0; $i < scalar(@total);$i++){
			for (my $j=0; $j < scalar(@total); $j++){
				if ( $i < $j ){
					my ($atom_1,$axis_x_1,$axis_y_1,$axis_z_1) = split ' ', $total[$i];
					my ($atom_2,$axis_x_2,$axis_y_2,$axis_z_2) = split ' ', $total[$j];
					#
					my $distance    = Euclidean_distance($axis_x_1,$axis_y_1,$axis_z_1,$axis_x_2,$axis_y_2,$axis_z_2);
					$final_trial = $final_trial + verification_opt($atom_1,$atom_2,$distance);
					#
					if( $final_trial ==	 1 ) {
						$resultado = 1;
						last;
					}
					#
				}
			}
		}
		# verify for steric impediment, 1 yes, 0 no;
		#
		#print " $resultado \n";
		if ( $resultado == 0 ) {
			$boolean = 1;
			print "18\n mol \n";
			foreach my $u (@total) {
				print "$u\n";
			}
			#
		}
		}
}
#
#
}
