#! c:\perl\bin\perl

my $seq_file = "file4q.txt";
my $seq_dir = "TriTrypDB\\whole\\slices";
my @allfiles;
my $outrnafold;
my $outcentroidfold;
my $rnafold_file;
my $centroidfold_file;
my $rnafold_struct;
my $centroidfold_struct;
my $combo;
my $rspredict_file;
my $rspredict_struct;
my $rspredict_APSI;
my $rspredict_pairing_threshold;
my $rspredict_energy_density;
my $outrspredict;
my $vienna;
my $ct;

# get aln files

unless ( open(SEQS, "$seq_dir\\$seq_file") ) {

	print "Cannot open file \"$seq_dir\\$seq_file\"\n\n";
	exit;
}

my @allfiles = <SEQS>;

# Close the file

close SEQS;

# print header

print "sample type,sample name,rspredict structure,rspredict base pairs,rspredict number of loops,rspredict largest loop,rspredict largest bulge,rspredict hairpin length,rspredict max consecutive base pairs,rspredict length,rspredict APSI,rspredict pairing threshold,rspredict energy density\n";


# read aln file

foreach my $file (@allfiles) {

	chomp($file);

	$rspredict_file = $file;
	$rspredict_file =~ s/-fin/-rspredict/;

	$combo = $file;
	$combo =~ s/-fin.txt//;

	$vienna = $file;
	$vienna =~ s/.txt/.vienna/;

	$ct = $file;
	$ct =~ s/.txt/.ct/;


# run through rnaalifold

	my $outrspredict = `java -jar RSpredict.jar -c $seq_dir\\$file >$seq_dir\\$rspredict_file`;


# parse rspredict structure

	unless ( open(RNA, "$seq_dir\\$rspredict_file") ) {
		print "Cannot open file \"$seq_dir\\$rspredict_file\"\n\n";
		exit;
	}

	
	my @rspredict_lines = <RNA>;

	close RNA;

	unlink "$seq_dir\\$rspredict_file";

	unlink "$seq_dir\\$vienna";

	unlink "$seq_dir\\$ct";

	foreach $line (@rspredict_lines) {

		if ($line =~ /^[\.\(\)]/) {
			$rspredict_struct = $line;

			chomp($rspredict_struct);
		}
		elsif ($line =~ /^The average pairwise sequence identity is:/) {
			$rspredict_APSI = $line;
			$rspredict_APSI =~ s/^The average pairwise sequence identity is:\s//;
			$rspredict_APSI =~ s/%//;
			$rspredict_APSI = $rspredict_APSI/100;
		}
		elsif ($line =~ /^The pairing threshold is:/) {
			$rspredict_pairing_threshold = $line;
			$rspredict_pairing_threshold =~ s/^The pairing threshold is:\s*//;
			$rspredict_pairing_threshold =~ s/\s*$//;
		}
		elsif ($line =~ /^The energy density/) {
			$rspredict_energy_density = $line;
			$rspredict_energy_density =~ s/^The energy density of the predicted structure is:\s*//;
		}

	}


#	print "rspredict $rspredict_struct\n";
	

	my $i = 0;
	my $s0 = "";
	my $s1 = "";
	my $s2 = "";
	my $s3 = "";
	my $s4 = "";
	my $end;
	my $ss_length = 0;
	my $largest_loop_length = 0;
	my $largest_bulge_length = 0;
	my $loop_number = 0;
	my $bulge_number = 0;
	my $blech_number = 0;
	my $consecutive_base_pairs = 0;
	my $max_consecutive_base_pairs = 0;
	my $rspredict_hairpin = "";
	my $rspredict_hairpin_length = 0;
	my $rspredict_length = 0;

# find structure length

	$rspredict_length = length($rspredict_struct);

#	print "rspredict_length $rspredict_length\n";

# find hairpin length

	$rspredict_hairpin = $rspredict_struct;
	$rspredict_hairpin =~ s/^\.*//g;
	$rspredict_hairpin =~ s/\.*$//g;
	$rspredict_hairpin_length = length($rspredict_hairpin);

#	print "rspredict_hairpin_length $rspredict_hairpin_length\n";

	my @rspredict_array = split('',$rspredict_struct);

	my $current = $rspredict_array[$i];
	my $state = "s0";
	my $previous_state = "start";

	while ($i < $rspredict_length) {

#		print "current $current\n";

		if ($state eq "s0" && $current eq ".") {
			$s0 .= $current;
		}
		elsif ($state eq "s0" && $current eq "(") {
			$s1 .= $current;
			$state = "s1";
			$previous_state = "s0";
			$consecutive_base_pairs++;
		}
		elsif ($state eq "s1" && $current eq "(") {
			$s1 .= $current;
			$consecutive_base_pairs++;
		}
		elsif ($state eq "s1" && $current eq ".") {
			$s2 .= $current;
			$state = "s2";
			$previous_state = "s1";
			$ss_length++;
			if ($consecutive_base_pairs >= $max_consecutive_base_pairs) {
				$max_consecutive_base_pairs = $consecutive_base_pairs;
			}
		$consecutive_base_pairs = 0;		
		}
		elsif ($state eq "s2" && $current eq ".") {
			$s2 .= $current;
			$ss_length++;
		}
		elsif ($state eq "s2" && $current eq "(") {
			$s1 .= $current;
			$state = "s1";
			$previous_state = "s2";
#			print "bulge of length $ss_length ends at $i\n";

		# find largest bulge length

			if ($ss_length >= $largest_bulge_length) {
				$largest_bulge_length = $ss_length;
			}

			$bulge_number++;
			$ss_length = 0;
			$consecutive_base_pairs++;
		}
		elsif ($state eq "s2" && $current eq ")") {
			$s3 .= $current;
			$state = "s3";
			$previous_state = "s3";
#			print "loop of length $ss_length ends at $i\n";

			# find largest loop length

			if ($ss_length >= $largest_loop_length) {
				$largest_loop_length = $ss_length;
			}

			$loop_number++;
			$ss_length = 0;
		}
		elsif ($state eq "s3" && $current eq ".") {	
			$s4 .= $current;
			$state = "s4";
			$previous_state = "s3";
			$ss_length++;
		}
		elsif ($state eq "s3" && $current eq "(") {
			$s1 .= $current;
			$state = "s1";
			$previous_state = "s3";
			$ss_length = 0;
		}
		elsif ($state eq "s3" && $current eq ")") {
			$s3 .= $current;
		}
		elsif ($state eq "s4" && $current eq ".") {
			$s4 .= $current;
			$ss_length++;
		}
		elsif ($state eq "s4" && $current eq "(") {
			$s1 .= $current;
			$state = "s1";
			$previous_state = "s4";
#			print "blech of length $ss_length ends at $i\n";
			$blech_number++;
			$ss_length = 0;	
		}
		elsif ($state eq "s4" && $current eq ")") {
			$s3 .= $current;
			if ($ss_length != 0) {
#				print "bulge of length $ss_length ends at $i\n";

				# find largest bulge length

				if ($ss_length >= $largest_bulge_length) {
					$largest_bulge_length = $ss_length;
				}
			$bulge_number++;
			}
			$ss_length = 0;
		}
		$i = $i + 1;
		$current = $rspredict_array[$i];

	}

my $bp = length($s1);

print "positive,$combo,$rspredict_struct,$bp,$loop_number,$largest_loop_length,$largest_bulge_length,$rspredict_hairpin_length,$max_consecutive_base_pairs,$rspredict_length,$rspredict_APSI,$rspredict_pairing_threshold,$rspredict_energy_density";

}


