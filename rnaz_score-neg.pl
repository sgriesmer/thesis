#! c:\perl\bin\perl

use Statistics::Basic qw(:all);

# open directory

my $seq_dir = "RF00029_60-blastclust-sub1-nh\\negative";
my $updated;
my $neg_dir = "z_score";
my $err_file = "error";
my $seq_file = "file4.txt";
my $score_file;
my $j;
my $pos_mfe;
my $pos_rnaali;

my $rnaz_file = "rnazout-neg-".$seq_file;

unless ( open(SEQS, "$seq_dir\\$seq_file") ) {

	print "Cannot open file \"$seq_dir\\$seq_file\"\n\n";
	exit;
}

# Read RNAfold records and store in array

my @allfiles = <SEQS>;

# Close the file

close SEQS;


foreach my $file (@allfiles) {

	# find mfe and structure

	print ".";

	chomp($file);

	my $out_rnaalign = `rnaz $seq_dir\\$file >> $seq_dir\\$rnaz_file`;

}


