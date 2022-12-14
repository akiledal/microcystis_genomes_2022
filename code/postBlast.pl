#!/usr/bin/perl

=head1 Description

	extract blast outputs based on any of the following: evalue, bitscore, percent ID, alignment length and % query coverage

=head1 USAGE

	perl postBlast.pl -b <blast output> 

=head2 Options

=head3 Thresholds

	-e	e-value; default=1
	-s	bit score; default=0
	-p	min % id; default=0
	-l	aln len; default=0
	-c	Query Coverage in %;
		requires a number between 0-100; 
		requires '-covCol'
	-covCol	Column number for query coverage in blast output;
	-list	Just get me the queries/subj in the list File;
		if you wish to search in subj column for matches; mention with flag "-subj";
		default= will look in query;
	
=head3 Dependencies
	
	-subj	look for items in list in the subj col; required if using '-list' AND wish to only check subjects for matches to the list.

=head3 Output

	-o	output file name; default= processID.pb.out
	-r	inverse of output file; default= not printed at all.

=head1 Author

	Sunit Jain, January 2010 : sunitj [AT] umich [DOT] edu
	- updated September 2015
	
=cut

use strict;
use Getopt::Long;

my $setEval=1; 
my $setPer=0;
my $setLen=0;
my $setScore=0;
my $covCol;
my $out=$$.".pb.out";
my $version="postBlast.pl\tv0.1.7";
my $wholeSeqPerID=0;

my(
$blast,
$list,
$subj,
$queryFasta,
$setCoverage,
$getCoverage,
$qLenCol,
$remainder,
);

GetOptions(
	'b|blast:s'=>\$blast,
	'e|evalue:f'=>\$setEval,
	's|bit_score:f'=>\$setScore,
	'p|per:i'=>\$setPer,
	'l|len:i'=>\$setLen,
	'o|out:s'=>\$out,
	'r|remainder:s'=>\$remainder,
	'list:s'=>\$list,
	'c|cov:i'=>\$setCoverage,
	'covCol:i'=>\$covCol,
	'qLenCol:i'=>\$qLenCol,
	'f|q|fasta|query:s'=>\$queryFasta,
	'subj'=>\$subj,
	'add_cov'=>\$getCoverage,
	'spid:i'=>\$wholeSeqPerID,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";
die "[ERROR] Blast output file required. See $0 -h for help.\n" if (!$blast);

if (($list) && (!$subj)){ print "# Searching Queries\n";}
elsif (($list) && ($subj)){ print "# Searching Subjects\n";}

## create index file from list;
my %index;
&readListFile if ($list);

## get Lengths and set coverage value
my (%lengths, $printCov);
if ($setCoverage){
	die "Coverage column required. Use -covCol to specify the column number\n" if (! $covCol);
	$printCov=$setCoverage;
#	&getLengths;
}
else{
	$setCoverage=0;
	$printCov=0;
}


open (BLAST, $blast) || die "[ERROR]: $!\n";
open (OUT, ">".$out) || die "[ERROR]: $!\n";
if($remainder){open (REM, ">".$remainder) || die "[ERROR]: $!\n";}
print OUT "#Thresholds: Blast Output=$blast; \%ID=$setPer; aln_len=$setLen; e-value=$setEval; bitScore=$setScore; QueryCoverage=$printCov\%\n";
my $header="#Columns: Query\tSubject\tPercent_ID\tAln_Len\tMismatch\tGap\tQuery_Start\tQuery_Stop\tSubject_Start\tSubject_Stop\tE-Value\tBit_Score";
$header.= ($getCoverage ? "\tCoverage\n" : "\n");
print OUT $header;
while(my $line=<BLAST>){
	next if $line=~ /^#/;
	$line=~ s/\r//;
	chomp($line);
	next unless $line;

	my ($q, $s, $per, $len, $mm, $gap, $q_start, $q_stop, $s_start, $s_stop, $eval, $score, @extra)=split(/\t/, $line);
	my $cov=$setCoverage > 0 ? $extra[$covCol-13] : 0;
	my $qLen=$extra[$qLenCol-13];
	my $id_along_whole_seq = $wholeSeqPerID > 0 ? (($len-$mm-$gap)/$lengths{$q}) : 0;
	my $item=$subj ? $s : $q;
	if (($per >= $setPer) && ($len >= $setLen) && ($eval <= $setEval) && ($score >= $setScore) && ($cov >= $setCoverage)){
		my $printLine="";
		if (($list) && ($index{$item})){
			$printLine = $line
		}
		elsif(($list) && (! $index{$item})){
			next;
		}
		else{
			$printLine = $line;
		}

		if($getCoverage){
			$cov=sprintf("%.3f", ($cov*100));
			$printLine.="\t".$lengths{$q}."\t".$cov."\t".$id_along_whole_seq;
		}

		print OUT $printLine."\n";
	}
	elsif($remainder){
		print REM $line."\n";
	}
}

## Sub-routines ##
sub getLengths{
	open(FASTA, $queryFasta)|| die "[ERROR]: $!\n";
	$/=">";
	while(my $line=<FASTA>){
		chomp $line;
		my ($head, @sequence)= split(/\n/, $line);
		my ($query, @stuff)=split(/\s+/, $head);

		if ($list){
			next unless $index{$query};
		}
		
		my $seqLen=length(join("", @sequence));
		$lengths{$query}=$seqLen;
	}
	$/="\n";
	close FASTA;
	return;
}

sub readListFile{
	open (LIST, $list)|| die "[ERROR]: $!\n";
	while(my $line=<LIST>){
		next if $line=~ /^#/;
		$line=~ s/\r//;
		chomp($line);
		next unless $line;

		$index{$line}++;
	}
	close LIST;
	return;
}
