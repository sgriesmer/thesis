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
my @pair_line;
my $prob;
my $shannon;
my $dbp;
my $sum_shannon;
my $sum_dbp;
my $length;
my $norm_shannon;
my $norm_dbp;

# get aln files

unless ( open(SEQS, "$seq_dir\\$seq_file") ) {

	print "Cannot open file \"$seq_dir\\$seq_file\"\n\n";
	exit;
}

my @allfiles = <SEQS>;

# Close the file

close SEQS;

# print header

print "sample type,sample name,shannon entropy,base-pair distance\n";



# read aln file

foreach my $file (@allfiles) {

	chomp($file);

	$combo = $file;
	$combo =~ s/-fin.txt//;

# initialize sums

	$sum_shannon = 0;
	$sum_dbp = 0;

# run through rnaalifold

	my $outrnafold = `rnaalifold -p <$seq_dir\\$file >$seq_dir\\temp`;

	unlink "$seq_dir\\temp";


# parse rnaalifold structure

	unless ( open(RNA, "alifold.out") ) {
		print "Cannot open file \"alifold.out\"\n\n";
		exit;
	}

	
	my @rnafold_lines = <RNA>;

	close RNA;

	unlink "$alifold.out";

	foreach $line (@rnafold_lines) {
		if ($line =~ /length of alignment/) {
			$length = $line;

#			print "line $line\n";

			$length =~ s/.*length of alignment\s//;
			
#			print "length $length\n";
		}
		elsif ($line =~ /^\s*[0-9]/) {
			@pair_line = split " ", $line;
			$prob = $pair_line[3];
			$prob =~ s/\%//;
			$prob = $prob/100;
			if ($prob != 0) {
				$shannon = (-1)*$prob*log($prob)/0.6931;
			}
			
#			print "pair_line $shannon\n";

			$sum_shannon += $shannon;

			$dbp = $prob - ($prob)**2;
			$sum_dbp += $dbp;
		}

	}

	$norm_shannon = $sum_shannon/$length;
	$norm_dbp = $sum_dbp/$length;

	print "positive,$combo,$norm_shannon,$norm_dbp\n";
}


