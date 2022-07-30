#! c:\perl\bin\perl

use Statistics::Basic qw(:all);

# program parameters

my $slice_dir = "srprna_ali-fasta_60-blastclust-sub1-nh\\negative";
my $seq_file = "rnazout-neg-file2.txt";
my $seq_no = 0;
my $ref_seq = 0;
my $brucei = 1;

# RNAz features

my $sequences;
my $reading_direction;
my $apsi;
my $mfe_avg;
my $mfe_consensus;
my $z_score;
my $sci;
my $probability;
my $prediction;
my $seq_id;
my $length;


unless ( open(SEQS, "$slice_dir\\$seq_file") ) {

	print "Cannot open file \"$slice_dir\\$seq_file\"\n\n";
	exit;
}

# Read RNAz records into an array

my @all_lines = <SEQS>;
$length = scalar(@all_lines);

# Close the file

close SEQS;

print "seq_no,seq_id,brucei_seq,sequences,reading_direction,apsi,mfe_avg,mfe_consensus,z_score,sci,probability,prediction\n";

foreach my $line (@all_lines) {

# delete cr

	chomp ($line);

	$length--;

	if (($line =~ /RNAz 1.0/) || $length == 0) {

		if ($seq_no > 0) {
			
			print "$seq_no,$seq_id,$brucei_seq, $sequences,$reading_direction,$apsi,$mfe_avg,$mfe_consensus,$z_score,$sci,$probability,$prediction\n";
		}

		$seq_no++;
		$ref_seq = 0;
	}
	elsif ($line =~ /Sequences:/) {
		$sequences = $line;
		$sequences =~ s/Sequences: //;
		$sequences =~ s/\s//g;
	}
	elsif ($line =~ /Reading direction:/) {
		$reading_direction = $line;
		$reading_direction =~ s/Reading direction://;
		$reading_direction =~ s/\s//g;
	}
	elsif ($line =~ /Mean pairwise identity:/) {
		$apsi = $line;
		$apsi =~ s/Mean pairwise identity://;
		$apsi =~ s/\s//g;
	}
	elsif ($line =~ /Mean single sequence MFE:/) {
		$mfe_avg = $line;
		$mfe_avg =~ s/Mean single sequence MFE://;
		$mfe_avg =~ s/\s//g;
	}
	elsif ($line =~ /Consensus MFE:/) {
		$mfe_consensus = $line;
		$mfe_consensus =~ s/Consensus MFE://;
		$mfe_consensus =~ s/\s//g;
	}
	elsif ($line =~ /Mean z-score:/) {
		$z_score = $line;
		$z_score =~ s/Mean z-score://;
		$z_score =~ s/\s//g;
	}
	elsif ($line =~ /Structure conservation index:/) {
		$sci = $line;
		$sci =~ s/Structure conservation index://;
		$sci =~ s/\s//g;
	}
	elsif ($line =~ /SVM RNA-class probability:/) {
		$probability = $line;
		$probability =~ s/SVM RNA-class probability://;
		$probability =~ s/\s//g;
	}
	elsif ($line =~ /Prediction:/) {
		$prediction = $line;
		$prediction =~ s/Prediction://;
		$prediction =~ s/\s//g;
	}
	elsif ($line =~ /\>/  && $ref_seq == 0) {
		$seq_id = $line;
		$seq_id =~ s/\s//g;
		$ref_seq = 1;
		$brucei = 0;
	}
	elsif ($line =~ /^[AUCG-]*/ && $brucei == 0) {
		$brucei_seq = $line;
		$brucei = 1;
	}
}