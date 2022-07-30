#! c:\perl\bin\perl

use Statistics::Basic qw(:all);

# open directory

my $seq_dir = "RF00029_60-blastclust-sub1-nh";
my $updated;
my $neg_dir = "z_score";
my $err_file = "error";
my $seq_file = "file4.txt";
my $fasta_file;
my $mfe_file;
my @mfe_array;
my $combo;
my $score_file;
my $j;
my $pos_mfe;
my $mfe;
my $pos_rnaali;
my @name;
my $nm;
my @seq;
my $sequence;

my $neg_path = ".\\".$seq_dir."\\".$neg_dir;

# print "$neg_path\n";

mkdir $neg_path;

unless ( open(SEQS, "$seq_dir\\$seq_file") ) {

	print "Cannot open file \"$seq_dir\\$seq_file\"\n\n";
	exit;
}

# Read RNAfold records and store in array

my @allfiles = <SEQS>;

# Close the file

close SEQS;


foreach my $file (@allfiles) {

	# initialize arrays

	@mfe_array = ();
	@seq = ();
	@name = ();
	$sequence = "";

	# find mfe and structure

	chomp($file);

	unless ( open(SEQFILE, "$seq_dir\\$file") ) {

		print "Cannot open file \"$seq_dir\\$file\"\n\n";
		exit;
	}

	# Read Clustal file

	my @pos_rnafold_lines = <SEQFILE>;

	# Close the file

	close SEQFILE;

	my $i = 0;

	foreach my $line (@pos_rnafold_lines) {

		chomp($line);

	# discard blank lines

		if ($line =~ /^$/) {
			$i = 0;
			next;

		}


	# discard header line

		elsif ($line =~ /^CLUSTAL/) {
			next;

		}

	# discard similarity line

		elsif ($line =~ /^\s\*/) {
			next;

		}

	# sequence line

		elsif($line =~ /^[a-z]/) {
			$nm = $line;
			$nm =~ s/\s[-UACG]*$//;
			$name[$i] = $nm;
			$sequence = $line;
			$sequence =~ s/^[a-z].*\s//;
			$seq[$i] .= $sequence;
			$i++;

		}
	}

	$fasta_file = $file;
	$fasta_file =~ s/fin/fasta/;
	$mfe_file = $file;
	$mfe_file =~ s/fin/mfe/;

	unless ( open(FASTA, "> $seq_dir\\$fasta_file") ) {

		print "Cannot open file \"$seq_dir\\$fasta_file\"\n\n";
		exit;
	}



	while ($i > 0) {
		$i--;
		print FASTA ">$name[$i]\n";
		print FASTA "$seq[$i]\n";
	}

	my $outmfe = `rnafold <$seq_dir\\$fasta_file >$seq_dir\\$mfe_file`;

	# Close and delete fasta file

	close FASTA;

	unlink "$seq_dir\\$fasta_file";

	unless ( open(MFE, "$seq_dir\\$mfe_file") ) {

		print "Cannot open file \"$seq_dir\\$mfe_file\"\n\n";
		exit;
	}

	my @mfe_lines = <MFE>;

	# Close and delete MFE file

	close MFE;

	unlink "$seq_dir\\$mfe_file";

	foreach $line (@mfe_lines) {

		if ($line =~ /\(-[0-9]/){
			$mfe = substr($line, -9);
			$mfe =~ s/\)//;
			$mfe =~ s/\(//;
			
		#	print "mfe $mfe\n";

			push @mfe_array, $mfe;
		}

	}

	$combo = $file;
	$combo =~ s/-fin.txt//;

	my $mean = mean(@mfe_array);

	print "$combo,$mean\n";

}


