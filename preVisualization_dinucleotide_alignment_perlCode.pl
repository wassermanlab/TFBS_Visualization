#!/usr/bin/perl -w

####################################
# Rebecca Worsley Hunt  2013.11.21 #
# Wasserman Lab                    #
####################################
#
#
# preVisualization_dinucleotide_alignment_perlCode.pl
# Purpose: to align sequences on a given feature and calculate the frequency of dinucleotides at each position of the 
#          the alignment. The output file is submitted to the R function: dinucleotide_view() found in the R file:
#          Visualization_methods_Rcode.R
# Alignment of provided sequences is anchored on a feature, such as a transcription factor motif in each sequence.
# If the feature has an orientation, i.e. it is stranded, then provision of the relevant strand results in all
# sequences being oriented to place the motif in a 5' to 3' orientation before aligning and counting the dinucleotides
# at each position relative to the motif. If no strand is provided the alignment retains the orientation of the
# sequences as they were provided.
#
#	NOTE: if the existing offset information is the centre of the motif, then the user must convert the offset to be 
#         the start of the motif
#   NOTE: Any dinucleotide containing 'N' has a count of 0 
#
# Required input: 
#		1) a fasta file of sequences containing a feature (motif) of interest
# 	 	2) a containing the offset of the feature relative to the start of the sequence. The feature is the point
#           of alignment.  The offset is expected to be an integer.
#		3) the column number or name that contains the offset data. The first column in the file is number 1.
#		4) the name of the output file (an extension of *.dinuc.align  will be appended)
# Optional input:  
#		Required if orienting all the motifs/features 5'-3'
# 		1) the column number or name in the meta file that contains the strand information. Strand is either + or -
#           The first column in the file is number 1. 
#		2) the width of the motif i.e. length of the PWM. The width is critical for adjusting alignments. The width
#           is an integer.
# Output:
#       a single file, with extension *.dinuc.align
#


use strict;
use Getopt::Long;
use IO::Handle;
use List::Util qw(max);

use constant FALSE => 0;
use constant TRUE => 1;

my $seqFile;
my $metaFile;  		
my $colposition;  		# the column (in meta file) with the motif offset
my $colstrand; 			# the column (in meta file) with the strand of motif. Needed if orienting alignment on motif 5'.
my $width_pfm = -1; 	# needed for strand alignment
my $dinucl_filename;
my $dinucl_map_file; 	# pre-cursor file
my $dinucl_alignment_file; # the output alignment file, that is used by R code


GetOptions(
	'f=s'   => \$seqFile,
	'm=s'	=> \$metaFile,  	 
	'p=i'	=> \$colposition,  	
	's=s'	=> \$colstrand, 
	'w=i'	=> \$width_pfm, 	
	'o=s'   => \$dinucl_filename,	
);
if( !$seqFile || !$metaFile || ! $colposition || !$dinucl_filename || !-e $seqFile || !-e $metaFile){
	print "-f  = fasta seq file\n",
		  "-m  = meta file: motif offset RELATIVE TO START OF SEQ and strand if using it\n",
		  "-p  = column number containing motif offset (a column in -m metafile). First column is 1. \n",
		  "-o  = output name  (output extension will be: *.dinuc.align)\n", 
		  "-s  = (optional) column number containing strand +/- (a column in -m meta file). First column is 1.\n",
		  "-w  = (optional) motif 'width' i.e. length of PWM. Necessary for alignment adjustment between +/- strands.\n";
	exit 1;
}
$dinucl_map_file = join("", $dinucl_filename, ".dinuc.map");
$dinucl_alignment_file = join("", $dinucl_filename, ".dinuc.align");
if( $metaFile && $colposition && $colstrand && $width_pfm == -1 ){
	print "Is position a motif location? (y/n) ";
	my $answer = <>;
	chomp($answer);
	if($answer eq "yes" || $answer eq "y"){
		print "you forgot to include width of the motif\n-w = width of pfm\n";
		exit 1;
	}
}
print "Begin ".scalar localtime(),"\n";
open(FH, "<$seqFile") || die "Error opening input\n";
open(FH_OUT, ">$dinucl_map_file") || die "Error opening outfile\n";

# verify metadata
if( $colposition > 0 && $colposition =~ /\d+/ ){
	my $testmeta = `sed '/^#/d; /^\$/d' $metaFile | cut -f $colposition | grep -v ^[1-9] | wc -l`;
	if( $testmeta != 0){
		print "there are offsets not starting with [1-9]. Are there -ve numbers? Integers must be positive and cannot be 0.\n";
		exit 1;
	}
}
if( $colstrand > 0 && $colstrand =~ /\d+/ ){
	my $testmeta = `sed '/^#/d; /^\$/d' $metaFile | cut -f $colstrand | grep -v ^[+-] | wc -l`;
	if( $testmeta != 0){
		print "there are invalid strand values. Values must be +/-.\n";
		exit 1;
	}
}

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# dinucleotide number IDs:
#  AA=1, GG=2, CC=3, TT=4, AG=5, AC=6, AT=7, GA=8, GC=9, GT=10, CA=11, CG=12 CT=13, TA=14,TG=15, TC=16

my $firstline = TRUE;
my $seqID = 0;
my $sequence = "";
my $dinuc = "";
my $seqnum = -1;
my $use_strand = FALSE;
my $positionFile;
my $dna_strandFile;
my @strands;
my @pos_strand;
my @offsets;

if( $colstrand =~ /\d+/ && $colstrand >0 ){
	$use_strand = TRUE;
	$dna_strandFile = join("", $seqFile,".strand.tmp");
	@strands = split(" ", `cut -f $colstrand $metaFile` );  
}
if( $colposition =~ /\d+/ && $colstrand =~ /\d+/ ){
	$positionFile = join("", $seqFile,".position.tmp");
	@offsets = split(" ", `cut -f $colposition $metaFile`); 
}

while( my $line = <FH>) {
	chomp($line);
	if( $line =~ m/^>/) {
		if( $firstline == FALSE ) {
			$seqnum++;
			#deal with previous sequence
			if( $use_strand == TRUE ) {
				push(@pos_strand, map_dinucleotides(\*FH_OUT, $seqID, $sequence, 
				     $strands[$seqnum], $offsets[$seqnum], $width_pfm ) ); 
			} else{  
				map_dinucleotides( \*FH_OUT, $seqID, $sequence, FALSE, FALSE, FALSE ) ;
			}
		}
		$seqID = $line;
		$sequence = ""; 
		$firstline = FALSE; 
	}else {
		$sequence = $sequence . $line;
	}
}continue{
	if(eof FH) {
		$seqnum++;
		if( $use_strand == TRUE ) {
			push(@pos_strand, map_dinucleotides( \*FH_OUT, $seqID, $sequence, 
			     $strands[$seqnum], $offsets[$seqnum], $width_pfm ) ); 
		} else{
			map_dinucleotides(\*FH_OUT, $seqID, $sequence, FALSE, FALSE, FALSE ); 
		}
	}
}
close(FH);
close(FH_OUT);
# get count of each dinucleotide occuring relative to a feature position 
#  e.g. relative to peakMax, top scoring motif etc
if( $positionFile ) {
	#if there is a strand provided, use both position and strand. Else only use position
	if( $dna_strandFile ) {
	    print "strand and offset  used\n"; #TTT	
		count_dinucl_relative_to_position( $dinucl_map_file, \@pos_strand );
	}else{
	    print "only offset used, not strand\n"; #TTT	
		count_dinucl_relative_to_position( $dinucl_map_file, \@offsets );
	}
}
unlink $dinucl_map_file; # delete map file
print "Fini ".scalar localtime()."\n\n";

exit;

#::::::::::::::::::::::::::::::::::::::::::::::::::::

sub count_dinucl_relative_to_position
{	
	# strandedness already dealt with in dinucleotide_map
	my ($dinucl_map, $position_ref ) = @_;
	my $max_offset = 0;
	my $max_end = 0;
	my $countseq = -1;
	my $arraylength = 0; 
	my $seqnum = -1;
	my $position = 0;

	$max_offset = max(@$position_ref); 
	# determine from dinucl_map file how big the array needs to be which will hold
	# 	the dinucleotide counts at each aligned position
	open(FH_DINUC, "<$dinucl_map" ) || die "Error, can't open the dinucleotide file\n";
	while( my $line = <FH_DINUC> ) {
		chomp($line);
		if( $line !~ m/^>/) {
			$countseq++;
			my @seq = split(/,/, $line);
			my $end_offset = scalar (@seq) - $position_ref->[ $countseq ];
			if( $end_offset > $max_end ) {
				$max_end = $end_offset;
			}
		}	
	}
	$arraylength = $max_offset + $max_end;

	seek(FH_DINUC, 0, 0); #reset FH 
	# array for each dincleotide. Why didn't I do hash?
	my @AA = ((0) x $arraylength);
	my @GG = ((0) x $arraylength);
	my @CC = ((0) x $arraylength);
	my @TT = ((0) x $arraylength);
	my @AG = ((0) x $arraylength);
	my @AC = ((0) x $arraylength);
	my @AT = ((0) x $arraylength);
	my @GA = ((0) x $arraylength);
	my @GC = ((0) x $arraylength);
	my @GT = ((0) x $arraylength);
	my @CA = ((0) x $arraylength);
	my @CG = ((0) x $arraylength);
	my @CT = ((0) x $arraylength);
	my @TA = ((0) x $arraylength);
	my @TG = ((0) x $arraylength);
	my @TC = ((0) x $arraylength);
	my @count_seqs_pos = ((0) x $arraylength);
	
	while( my $line = <FH_DINUC> ) {
		chomp($line);
		if( $line !~ m/^>/) {
			$seqnum++;
			my @seq = split(/,/, $line);
			my $start = $max_offset - $position_ref->[ $seqnum ];
			my $pos = -1;
			foreach( @seq ) {
				$pos++;
				# dinucleotide IDs: 1-16
				if( $_ == 1 ) { # 1 = AA 
					$AA[$start+$pos] += 1;
				}elsif ($_ == 2 ) { # 2 = GG 
					$GG[$start+$pos] += 1;
				}elsif ( $_ == 3) {
					$CC[$start+$pos] += 1;
				}elsif( $_ == 4){
					$TT[$start+$pos] += 1;
				}elsif ($_ == 5 ) {
					$AG[$start+$pos] += 1;
				}elsif ($_ == 6 ) {
					$AC[$start+$pos] += 1;
				}elsif ($_ == 7 ) {
					$AT[$start+$pos] += 1; 
				}elsif ($_ == 8 ) {
					$GA[$start+$pos] += 1;
				}elsif ($_ == 9) {
					$GC[$start+$pos] += 1; 
				}elsif ( $_ == 10) {
					$GT[$start+$pos] += 1;
				}elsif ($_ == 11) {
					$CA[$start+$pos] += 1;
				}elsif ($_ == 12) {
					$CG[$start+$pos] += 1; 
				}elsif ( $_ == 13) {
					$CT[$start+$pos] += 1;
				}elsif ($_ == 14 ) {
					$TA[$start+$pos] += 1; 
				}elsif ( $_ == 15) {
					$TG[$start+$pos] += 1;
				}elsif ($_ == 16 ) {
					$TC[$start+$pos] += 1;
				}else {
					#do nothing with 'N' containing cases (ID = 0)	
				}
				$count_seqs_pos[$start+$pos] += 1;
			}
		}
	}
	close(FH_DINUC);
	open(FH, ">$dinucl_alignment_file") || die "Error, couldn't open the final outfile\n";
	# print column numbers (used in R), followed by the number of seqs contributing to a position
	# and then the counts for each dinucleotide of the one strand considered
	print FH "position,";
	$position = 0;
	for( $position = 1-$max_offset; $position < $arraylength-$max_offset; $position++){
		print FH $position,",";
	}
	print FH $position,"\n"; #the last column number (already incremented above)
	print FH "countseq.alignPos",$max_offset,",",join(',', @count_seqs_pos),"\n",
			 "AA,",join(',', @AA),"\n",
			 "GG,",join(',', @GG),"\n",
			 "CC,",join(',', @CC),"\n",
			 "TT,",join(',', @TT),"\n",
			 "AG,",join(',', @AG),"\n",
			 "AC,",join(',', @AC),"\n",
			 "AT,",join(',', @AT),"\n",
			 "GA,",join(',', @GA),"\n",
			 "GC,",join(',', @GC),"\n",
			 "GT,",join(',', @GT),"\n",
			 "CA,",join(',', @CA),"\n",
			 "CG,",join(',', @CG),"\n",
			 "CT,",join(',', @CT),"\n",
			 "TA,",join(',', @TA),"\n",
			 "TG,",join(',', @TG),"\n",
			 "TC,",join(',', @TC),"\n";
	close(FH);
}

	
sub map_dinucleotides 
{
	# if using strandedness, reverse compliment -ve strand
	# label each position with its dinucleotide (ID 1-16)
	# store as a string
	# return string of dinuc IDs; comma separated
	my ($fh_out, $seqID_ref, $seq, $dna_strand, $pos_offset, $footprint ) = @_;
	my @dinucl_map;
	my $pos_new = $pos_offset;
	my @seqarray = split(//, $seq); #split up string into individual characters
	my $length = scalar( @seqarray);
	# -ve strands must be reversed (complimented below)
	if( $dna_strand eq "-"){
		@seqarray = reverse( @seqarray );
	}
	#for -ve strand, calculate 5' offset (currently 3')
	if( $dna_strand eq "-" && $pos_offset != FALSE){
		$pos_new = $length - ($pos_new + $footprint -1) +1;
	}
	# -ve strand: compliment the dinucleotide
	# +ve strand leave as is
	for( my $i=0; $i < $length-1; $i++) {
		my $di = $seqarray[$i].$seqarray[$i+1];
		if( $di eq "AA" ) {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 4);
			}else{ 
				push( @dinucl_map, 1);
			}
		}elsif( $di eq "GG") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 3);
			}else{ 
				push( @dinucl_map, 2);
			}	
		}elsif( $di eq "CC") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 2);
			}else{ 
				push( @dinucl_map, 3);
			}
		}elsif( $di eq "TT") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 1);
			}else{ 
				push( @dinucl_map, 4);
			}
		}elsif( $di eq "AG") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 16);
			}else{ 
				push( @dinucl_map, 5);
			}
		}elsif( $di eq "AC") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 15);
			}else{ 
				push( @dinucl_map, 6);
			}
		}elsif( $di eq "AT") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 14);
			}else{ 
				push( @dinucl_map, 7);
			}
		}elsif( $di eq "GA") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 13);
			}else{ 
				push( @dinucl_map, 8);
			}
		}elsif( $di eq "GC") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 12);
			}else{ 
				push( @dinucl_map, 9);
			}
		}elsif( $di eq "GT") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 11);
			}else{ 
				push( @dinucl_map, 10);
			}
		}elsif( $di eq "CA") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 10);
			}else{ 
				push( @dinucl_map, 11);
			}
		}elsif( $di eq "CG") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 9);
			}else{ 
				push( @dinucl_map, 12);
			}
		}elsif( $di eq "CT") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 8);
			}else{ 
				push( @dinucl_map, 13);
			}
		}elsif( $di eq "TA") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 7);
			}else{ 
				push( @dinucl_map, 14);
			}
		}elsif( $di eq "TG") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 6);
			}else{ 
				push( @dinucl_map, 15);
			}
		}elsif( $di eq "TC") {
			if( $dna_strand eq "-" ){
				push( @dinucl_map, 5);
			}else{ 
				push( @dinucl_map, 16);
			}
		}else{
			push( @dinucl_map, 0); # cases containing 'N'
		}
	}
	print $fh_out $seqID_ref," | strand ",$dna_strand,"\n", join( ',', @dinucl_map ), "\n";
	return $pos_new; 
}

__END__
