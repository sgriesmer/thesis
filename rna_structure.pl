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

# get aln files

unless ( open(SEQS, "$seq_dir\\$seq_file") ) {

	print "Cannot open file \"$seq_dir\\$seq_file\"\n\n";
	exit;
}

my @allfiles = <SEQS>;

# Close the file

close SEQS;

# print header

print "sample type,sample name,rnafold structure,rnafold base pairs,rnafold number of loops,rnafold largest loop,rnafold largest bulge,rnafold hairpin length,rnafold max consecutive base pairs,rnafold length,centroidfold structure,centroidfold base pairs,centroidfold number of loops,centroidfold largest loop,centroidfold largest bulge,centroidfold hairpin length,centroidfold max consecutive base pairs,centroidfold length\n";


# read aln file

foreach my $file (@allfiles) {

	chomp($file);

	$rnafold_file = $file;
	$rnafold_file =~ s/-fin/-rnafold/;

	$centroidfold_file = $file;
	$centroidfold_file =~ s/-fin/-centroidfold/;

	$combo = $file;
	$combo =~ s/-fin.txt//;

# run through rnaalifold

	my $outrnafold = `rnaalifold <$seq_dir\\$file >$seq_dir\\$rnafold_file`;

# run through centroidfold

	my $outcentroidfold = `centroid_fold $seq_dir\\$file >$seq_dir\\$centroidfold_file`;

# parse rnaalifold structure

	unless ( open(RNA, "$seq_dir\\$rnafold_file") ) {
		print "Cannot open file \"$seq_dir\\$rnafold_file\"\n\n";
		exit;
	}

	
	my @rnafold_lines = <RNA>;

	close RNA;

	unlink "$seq_dir\\$rnafold_file";

	foreach $line (@rnafold_lines) {

		if ($line =~ /^[\.\(\)]*\s/) {
			$line =~ s/ \(.*\) //;
			$rnafold_struct = $line;

			chomp($rnafold_struct);
		}

	}

	

# parse centroidfold structure

	unless ( open(CENTROID, "$seq_dir\\$centroidfold_file") ) {
		print "Cannot open file \"$seq_dir\\$centroidfold_file\"\n\n";
		exit;
	}


	my @centroidfold_lines = <CENTROID>;

	close CENTROID;

	unlink "$seq_dir\\$centroidfold_file";

	foreach $line (@centroidfold_lines) {

		if ($line =~ /^[\.\(\)]*\s/) {
			$line =~ s/ \(.*\)//;
			$centroidfold_struct = $line;
			
			chomp($centroidfold_struct);
		}

	}

#	print "rnafold $rnafold_struct\n";
#	print "centroidfold $centroidfold_struct\n";

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
	my $rnafold_hairpin = "";
	my $rnafold_hairpin_length = 0;
	my $rnafold_length = 0;

# find structure length

	$rnafold_length = length($rnafold_struct);

#	print "rnafold_length $rnafold_length\n";

# find hairpin length

	$rnafold_hairpin = $rnafold_struct;
	$rnafold_hairpin =~ s/^\.*//g;
	$rnafold_hairpin =~ s/\.*$//g;
	$rnafold_hairpin_length = length($rnafold_hairpin);

#	print "rnafold_hairpin_length $rnafold_hairpin_length\n";

	my @rnafold_array = split('',$rnafold_struct);

	my $current = $rnafold_array[$i];
	my $state = "s0";
	my $previous_state = "start";

	while ($i < $rnafold_length) {

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
		$current = $rnafold_array[$i];

	}

my $bp = length($s1);

print "positive,$combo,$rnafold_struct,$bp,$loop_number,$largest_loop_length,$largest_bulge_length,$rnafold_hairpin_length,$max_consecutive_base_pairs,$rnafold_length,";

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
	my $centroidfold_hairpin = "";
	my $centroidfold_hairpin_length = 0;
	my $centroidfold_length = 0;

# find structure length

	$centroidfold_length = length($centroidfold_struct);

#	print "centroidfold_length $centroidfold_length\n";

# find hairpin length

	$centroidfold_hairpin = $centroidfold_struct;
	$centroidfold_hairpin =~ s/^\.*//g;
	$centroidfold_hairpin =~ s/\.*$//g;
	$centroidfold_hairpin_length = length($centroidfold_hairpin);

#	print "centroidfold_hairpin_length $centroidfold_hairpin_length\n";

	my @centroidfold_array = split('',$centroidfold_struct);

	my $current = $centroidfold_array[$i];
	my $state = "s0";
	my $previous_state = "start";

	while ($i < $rnafold_length) {

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
		$current = $centroidfold_array[$i];

	}

my $bp = length($s1);

print "$centroidfold_struct,$bp,$loop_number,$largest_loop_length,$largest_bulge_length,$centroidfold_hairpin_length,$max_consecutive_base_pairs,$centroidfold_length\n";

}


