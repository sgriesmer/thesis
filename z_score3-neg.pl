#! c:\perl\bin\perl

use Statistics::Basic qw(:all);

# open directory

my $seq_dir = "srprna_ali-fasta_60-blastclust-sub1-nh\\negative";
my $updated;
my $neg_dir = "z_score";
my $err_file = "error";
my $score_file;
my $j;
my $pos_mfe;
my $pos_rnaali;

my $neg_path = ".\\".$seq_dir."\\".$neg_dir;

# print "$neg_path\n";

mkdir $neg_path;


opendir(THISDIR, $seq_dir) or die "Cannot open $!\n";
my @allfiles = grep {/combo_[0-9]*_3_[0-9]*-ran.txt/} readdir THISDIR;
closedir THISDIR;


foreach my $file (@allfiles) {

	# find mfe and structure

	my $pos_file = $file;
	$pos_file =~ s/ran/pos/;

#	print "$pos_file\n";

	my $out_rnaalign = `rnaalifold $seq_dir\\$file >> $seq_dir\\$pos_file`;

	unless ( open(SEQFILE, "$seq_dir\\$pos_file") ) {

		print "Cannot open file \"$seq_dir\\$pos_file\"\n\n";
		exit;
	}

	# Read RNAfold records and store in array

	my @pos_rnafold_lines = <SEQFILE>;

	# Close the file

	close SEQFILE;

	# delete the file

	unlink "$seq_dir\\$pos_file";

	foreach my $line (@pos_rnafold_lines) {

		chomp($line);

#		print "$line\n";

	# discard sequence

		if ($line =~ /^[ACGU]/) {
			next;

		}

	# discard comment line

		elsif($line =~ /^[\.\(\)]*\s/) {

			$pos_rnaali = $line;
			$pos_rnaali =~ s/\s\(.*$//;


			$pos_mfe = $line;
			$pos_mfe =~ s/^[\.\(\)]*\s\(//;
			$pos_mfe =~ s/ =.*$//;

#			print "pos_rnaali $pos_rnaali\n\n";
		
#			print "pos_mfe $pos_mfe\n\n";

		}
	}

	# generate random files for z score

	$updated = $file;
#	$updated =~ s/fin/ran/;
	$score_file = $file;
	$score_file =~ s/fin/score/;
	$updated =~ s/ran/ran_0/;
	
	for (my $i = 1; $i <= 100; $i++) {

		$j = $i - 1;
	
		$updated =~ s/ran_$j/ran_$i/;

#		print "$updated\n\n";

		my $out_rnazran = `rnazRandomizeAln.pl $seq_dir\\$file > $seq_dir\\$neg_dir\\$updated`;

		my $out_rnaalign_ran = `rnaalifold $seq_dir\\$neg_dir\\$updated >> $seq_dir\\$neg_dir\\$score_file`;

		unlink "$seq_dir\\$neg_dir\\$updated";

	}

#	print "$updated\n\n";


	unless ( open(SEQFILE, "$seq_dir\\$neg_dir\\$score_file") ) {

		print "Cannot open file \"$seq_dir\\$neg_dir\\$score_file\"\n\n";
		exit;
	}

	# Read RNAfold records and store in array

	my @rnafold_lines = <SEQFILE>;
	my @amfe;

	# Close the file

	close SEQFILE;

	# Unlink score files

	unlink "$seq_dir\\$neg_dir\\$score_file";

	foreach my $line (@rnafold_lines) {

		chomp($line);

#		print "$line\n";

		# discard sequence

		if ($line =~ /^[ACGU]/) {
			next;

		}

		# discard comment line

		elsif($line =~ /^[\.\(\)]*\s/) {

			my $mfe = $line;
			$mfe =~ s/^[\.\(\)]*\s\(//;
			$mfe =~ s/ =.*$//;
		
#			print "mfe $mfe\n\n";

			push @amfe, $mfe;
		
		}
	}

	my $mean = mean(@amfe);

#	print "mean $mean\n";

	my $stddev = stddev(@amfe);

#	print "stddev $stddev\n";

	my $z = ($pos_mfe - $mean)/$stddev;

#	print "z $z\n";

	my $combo = $pos_file;
	$combo =~ s/-pos.txt//;

	print "$combo,$pos_rnaali,$pos_mfe,$z\n";

}


