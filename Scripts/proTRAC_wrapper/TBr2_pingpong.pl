#!/usr/bin/perl
use Getopt::Long;
$|=1;

$script=$0;
$script=~s/^.+[\\\/]//;
$info_text=
"
============================== NGS TOOLBOX ==============================
Program: pingpong
Version: 2.1.1
LAST MODIFIED: 17. November 2015

Usage: perl $script -i input -o output (-xhits -xreads)

Input files are map files as output by SeqMap or sRNAmapper (the small
RNA mapping tool that comes along with proTRAC). Map files that are post
processed with re-al-locate are also accepted. Ideally the FASTA or FASTQ
headers in the file you used for mapping referred to sequence read counts
(use the collapse tool from the NGS-TOOLBOX to convert your FASTA/FASTQ
file into this format). Default output is stdout (usually console). The
option -xhits will ignore the number of genomic hits produced by a
sequence (set to 1). The option -xreads will ignore the number of
sequence read counts for a given sequence (set to 1).

-h or -help will print this information.

Contact:
David Rosenkranz
Institute of Anthropology
Johannes Gutenberg University Mainz, Germany
email: rosenkranz\@uni-mainz.de

";

GetOptions
	(
	"help"=>\$print_info,
	"h"=>\$print_info,
	"i=s"=>\$input_file,
	"o=s"=>\$output_file,
	"xhits"=>\$norm_hits,
	"xreads"=>\$norm_reads,
	);

# print info
if($print_info)
	{
	print$info_text;
	exit;
	}

print"
============================== NGS TOOLBOX ==============================
Program: pingpong
Version: 2.1.1
last modified: 17. November 2015
=========================================================================
";

# process input file
print"\nAnalyze map file $input_file...";
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
if($output_file)
	{
	open(OUT,">$output_file")||die print"\nCould not create output file $output_file.\n$!\n\n";
	}

$max_length=0;
while(<IN>)
	{
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$_=~s/\s*$//;
		@d=split("\t",$_);
		if($d[6]=~/\+/)
			{
			$hits_plus++;
			}
		else
			{
			$hits_minus++;
			}
		if($norm_hits)
			{
			$hits_per_seq{$d[4]}=1;
			}
		else
			{
			$hits_per_seq{$d[4]}++;
			}
		if($norm_reads||$d[3]!~/^\d+$/)
			{
			$reads_per_seq{$d[4]}=1;
			}
		else
			{
			$reads_per_seq{$d[4]}=$d[3];
			}
		if(length$d[4]>$max_length) # save max. sequence length (-> maximum overlap)
			{
			$max_length=length$d[4];
			
			# random test map file format
			if($d[1]!~/^\d+$/||$d[2]!~/^[ATGCUN]+$/||$d[4]!~/^[ATGCUN]+$/||$d[5]!~/^\d+$/||$d[6]!~/^[\+-]$/)
				{
				print"\nBad input file format!\n\n";
				exit;
				}
			}
		}
	}
close IN;
$map_type="This is a standard map file:";
if($d[8])
	{
	$map_type="This is a post processed map file:";
	}
print" done.\n$map_type +hits: $hits_plus / -hits: $hits_minus

Scan map file for 5' overlaps:
0%                                   50%                                   100%
|-------------------------------------|-------------------------------------|
";

# read map with sliding window
@sw=();
$progress=0;
$processed_plus=0;
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
while(<IN>)
	{
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		@d=split("\t",$_);
		if($d[0]ne$prev_loc) # empty sliding window when location changes
			{
			@sw=();
			}
		if($d[6]=~/-/) # add hit to sliding window if current hit is on the minus strand
			{
			push(@sw,[$d[1],$d[4],$d[8]]); # sliding window entries comprise coordinate|sequence|allocated reads (if available)
			while(1) # check if hit(s) must be removed from sliding window
				{
				if(@sw>1&&$d[1]>($sw[0][0])+$max_length)
					{
					shift@sw;
					}
				else
					{
					last;
					}
				}
			}
		else # count overlaps in the sliding window if current hit is on the plus strand
			{
			foreach$i(0..@sw-1)
				{
				if($sw[$i][0]+(length$sw[$i][1])>$d[1]&&$sw[$i][0]+(length$sw[$i][1])<$d[1]+length$d[4]) # check whether hits overlap
					{
					if($d[8]) # post processed map (re-al-locate)
						{
						$overlap{$sw[$i][0]+(length$sw[$i][1])-$d[1]}+=($d[8]*$sw[$i][2]);
						}
					else # standard map
						{
						$overlap{$sw[$i][0]+(length$sw[$i][1])-$d[1]}+=(($reads_per_seq{$d[4]}/$hits_per_seq{$d[4]})*($reads_per_seq{$sw[$i][1]}/$hits_per_seq{$sw[$i][1]}));
						}
					}
				}
			$processed_plus++;
			if($processed_plus>=$hits_plus/77)
				{
				$progress++;
				$processed_plus=0;
				print".";
				}
			}
		$prev_loc=$d[0];
		}
	}
close IN;
foreach($progress..76)
	{
	print".";
	}
print"\n\n";

$results="";
$results.="Total scores for 5' overlap:\noverlap [nt]\tread pairs\n";
foreach(1..($max_length-1))
	{
	if($_>=1&&$_<=9)
		{
		$average_background+=$overlap{$_};
		}
	elsif($_>=11&&$_<=20)
		{
		$average_background+=$overlap{$_};
		}
	$results.="$_\t$overlap{$_}\n";
	}
$average_background=$average_background/19;

$variance=0;
foreach(1..($max_length-1))
	{
	if($_>=1&&$_<=9)
		{
		$variance+=($overlap{$_}-$average_background)**2;
		}
	elsif($_>=11&&$_<=20)
		{
		$variance+=($overlap{$_}-$average_background)**2;
		}
	}
$variance=$variance/19;
$stddeviation=sqrt($variance);
if($stddeviation>0)
	{
	$Zscore=($overlap{10}-$average_background)/$stddeviation;
	}
else
	{
	if($overlap{10}>0)
		{
		$Zscore="infinite";
		}
	else
		{
		$Zscore=0;
		}
	}

$results.="
Z-score calculation:
Average background (overlaps 1-9+11-20): $average_background
Variance background: $variance
Standard deviation: $stddeviation
Ping-Pong Z-Score: $Zscore

Z-score >= 1.6449 (one-tailed hypothesis) -> p < 0.05
Z-score >= 2.3264 (one-tailed hypothesis) -> p < 0.01
";

if($output_file)
	{
	print OUT"Results for $input_file\n$results";
	}
else
	{
	print$results;
	}
exit;