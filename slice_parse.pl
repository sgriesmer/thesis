#! c:\perl\bin\perl

use Statistics::Basic qw(:all);

# open directory

my $slice_dir = "TriTrypDB\\whole";
my $seq_file = "whole-slices-1q.txt";
my $out_dir = "slices";
my $file_no = 1;


unless ( open(SEQS, "$slice_dir\\$seq_file") ) {

	print "Cannot open file \"$slice_dir\\$seq_file\"\n\n";
	exit;
}

# Read RNAfold records and store in array

my @all_lines = <SEQS>;

# Close the file

close SEQS;


foreach my $line (@all_lines) {

	if ($line =~ /CLUSTAL/) {

		if ($file_no > 1) {
			close SEQFILE;
		}

		$outfile = "combo_try_" . $file_no ."-fin" . ".txt";

		unless (open(SEQFILE, ">$slice_dir\\$out_dir\\$outfile")) {
			print "Cannot open file \"$slice_dir\\$out_dir\\$outfile\"\n\n";
			exit;
		}

		print SEQFILE $line;

		$file_no++;

	}
	else {
		print SEQFILE $line;
	}

}



