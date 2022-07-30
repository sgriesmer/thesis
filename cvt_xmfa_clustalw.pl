#! c:\perl\bin\perl

use Bio::AlignIO;

$in  = Bio::AlignIO->new(-file => "whole.xmfa" ,
                         -format => 'xmfa');
$out = Bio::AlignIO->new(-file => ">whole.aln",
                         -format => 'clustalw');
while ( my $aln = $in->next_aln() ) {
	 $out->write_aln($aln); 
	}

