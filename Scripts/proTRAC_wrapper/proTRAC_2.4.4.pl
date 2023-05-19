#!/usr/bin/perl
use Getopt::Long;
use Scalar::Util qw( looks_like_number );
$proTRAC_version="2.4.4";
$proTRAC_date="10. August 2019";
$|=1;

###   SET DEFAULTS   ###
$sliding_window_size=5000;
$sliding_window_increament=1000;
$min_piRNAsize=24;
$max_piRNAsize=32;
$threshold_size=0.75;
$threshold_T1A10=0.75;
$threshold_T1A10_both=0.5;
$threshold_clustersize=1000;
$threshold_hits=0;
$threshold_hits_normalized=0;
$threshold_strandbias=0.75;
$bidirectionality_split=0.05;
$threshold_density_pvalue=0.01;
$threshold_density_absolute=0;
$flank_size=0;
$save_info=0;
$script_name_parsed=$0;
$script_name_parsed=~s/^.*[\\\/]//;
$help_desk="
================================= proTRAC ====================================
VERSION: .......... $proTRAC_version
LAST MODIFIED: .... $proTRAC_date

Usage:
perl $script_name_parsed -map mapfile.sam -genome genome.fas [-option value]

Please cite:
Rosenkranz D, Zischler H. proTRAC - a software for probabilistic piRNA cluster
detection, visualization and analysis. 2012. BMC Bioinformatics 13:5.

Contact:
David Rosenkranz
Institute of Organismic and Molecular Evolutionary Biology, Anthropology
Johannes Gutenberg University Mainz
email: rosenkranz\@uni-mainz.de

You can find the latest proTRAC version at:
http://www.smallRNAgroup.uni-mainz.de
==============================================================================
                              +++ OPTIONS +++
 ftp  = floating point number
 int  = integer
 0..1 = floating point number between 0 and 1. Type 0.0 for 0 and 1.0 for 1.

 -genome OR -g [file]    Name of the file that contains the genomic sequence
                         and that was used for mapping the sequence reads.
 -map OR -m [file]       Name of the file that contains mapped reads in SAM
                         or ELAND3 format. We recommend to use SeqMap with
                         option /output_all_matches or sRNAmapper with default
                         settings to create an appropriate file. When using a
                         more popular mapper you should make sure to allow for
                         multi-mapping.
 -format OR -f [s]       Specify the input format. Allowed values are SAM and
                         ELAND3. This is only required if the input file
			         contains less than 1000 hits.
 -help OR -h             Shows this information.
 -repeatmasker [file]    Name of the file that contains the RepeatMasker
                         annotation. Make sure that the names for the
                         chromosomes/scaffolds are identical in your Repeat-
                         Masker and genome file.
 -geneset [file]         Name of the file that contains gene annotation (GTF-
                         file from Ensembl database). Make sure that the names
                         for the chromosomes/scaffolds are identical in your
                         GTF- and genome file.
 -swsize [int]           Size of the sliding window (default=5000)
 -swincr [int]           Increment of the sliding window (default=1000)
 -nohc                   Do not consider total number of genomic hits for the
                         sequence in question for calculation of hit counts.
 -norc                   Do not consider number of reads for the sequence in
                         question for calculation of hit counts.
 -norpm                  Do not normalize the hit count with the total number
                         of mapped sequence reads (reads per million).
                         Normalization will make the values comparable accross
                         different piRNA libraries and/or organisms.
 -dens [fpt]             Define an absolute minimum number of (normalized) 
                         read counts per kb.
 -pdens [0..1]           Define a p-value for minimum number of (normalized)
                         read counts per kb. A p-value of 0.01 means that the
                         (normalized) read count value in a sliding window must
                         belong to the top 1% of all sliding windows.
 -est [int]              Use that option together with -pdens. Estimate the 
                         required minimum number of (normalized) read counts
                         in a sliding windows based on n random 1kb sliding
                         windows (faster). Without that option proTRAC will
                         scan the map file and calculate the required minimum
                         number of (normalized) read counts in a sliding
                         window based on the observed distribution.
 -pisize [0..1]          Fraction of (normalized) read counts that have
                         the typical piRNA size (default=0.75).
 -pimin [int]            Define the minimum length of a piRNA (default=24).
 -pimax [int]            Define the maximum length of a piRNA (default=32).
 -1Tor10A [0..1]         Fraction of (normalized) read counts that have 1T
                         (1U) or 10A (default=0.75).
 -1Tand10A [0..1]        If the fraction of (normalized) read counts with 1T
                         (1U) OR 10A is below the value defined by -1Tor10A,
                         accept the sliding window if BOTH the 1T (1U) and the
                         10A fraction reach this value (default=0.5).
 -distr [ftp-ftp]        To avoid false positive piRNA cluster annotation of
                         loci with one or few mapped sequences represented by
                         exceptionally many reads. You can e.g. type
                         '-distr 10-75' which means that the TOP 10% of
                         mapped sequences account for max. 75% of the piRNA
                         clusters (normalized) read counts. Otherwise the
                         locus is rejected (default=5-90).
 -spike [ftp-ftp]        Somewhat similar to -distr. With e.g. 50-1000 you
                         allow max. 50% of reads in a cluster to be within
                         any 1000bp sliding window. Otherwise the locus is
                         rejected (default=50-1000).			 
 -clsize [int]           Set the minimum size for a piRNA cluster (default=
                         1000).
 -clhits [int]           Minimum number of sequence hit loci per piRNA cluster
                         (default=0).
 -clhitsn [ftp]          Minimum number of normalized sequence read counts per
                         piRNA cluster (default=0).
 -clstrand [0..1]        Fraction of (normalized) read counts that map to the
                         main strand (default=0.75).
 -clsplit [0..1]         Minimum fraction of (normalized) read counts on the
                         smaller arm of a bi-directional piRNA cluster.
                         Otherwise the cluster will be annotated as
                         mono-directional (default=0.1).
 -nohtml                 Do not output .html files for each piRNA cluster.
 -notable                Do not output a summary table.
 -nofaspi                Do not output a FASTA file comprising mapped piRNAs
                         for each cluster.
 -nofascl                Do not output a FASTA file comprising all piRNA
                         cluster sequences.
 -nogtf                  Do not output a GTF file for predicted piRNA
                         clusters.
 -nomotif                Do not search for transcription factor binding
                         sites.
 -flank [int]            Include n bp of flanking sequence in output files.
 -pti                    Output a file that contains information on mapped
                         sequence reads in a tab-delimited table that
                         comprises sequence, reads, genomic hits e.g:
                         TGGGCACGCAAATTCGAGTATCG   12   4


";
GetOptions
	(
	"help|h"=>\$print_help,# print help/info text
	"genome|g=s"=>\$genome,# genome in fasta format
	"map|m=s"=>\$map,# SAM format or ELAND output from sRNAmapper or SeqMap using the option /output_all_matches
	"format|f=s"=>\$file_format,
	"repeatmasker=s"=>\$RMannotation,# output from RepeatMasker
	"geneset=s"=>\$GeneSet,# Gene Set (GTF-file) from ensembl database
	"swsize=i"=>\$sliding_window_size,
	"swincr=i"=>\$sliding_window_increament,
	"nohc"=>\$normalize_by_hit_number,# normalize hits by number of genomic hits
	"norc"=>\$normalize_by_read_number,# normalize hits by number of reads
	"norpm"=>\$normalize_by_total_number_of_mapped_reads,# normalize hits by number of reads
	"dens=f"=>\$threshold_density_absolute,# set user-defined minimum number of normalized hits per kb
	"pdens=f"=>\$threshold_density_pvalue,# p value for hit density
	"est=i"=>\$accuracy,# number of resampling steps when calculating minimum density based on expected hit density
	"pisize=f"=>\$threshold_size,# fraction of normalized reads with typical piRNA size
	"pimin=i"=>\$min_piRNAsize,# minimum size for a typical piRNA [nt]
	"pimax=i"=>\$max_piRNAsize,# maximum size for a typical piRNA [nt]
	"1Tor10A=f"=>\$threshold_T1A10,# fraction of normalized reads with 1T or 10A (either/or)
	"1Tand10A=f"=>\$threshold_T1A10_both,# fraction of normalized reads with 1T or 10A (both)
	"distr=s"=>\$avoid_peaks,# avoid false positive piRNA cluster prediction due to peaks of few mapped sequences with many reads (x-y which means: x% of top sequences account for max. y% of reads. Default: 5-50);
	"spike=s"=>\$no_spikes, # maximum n% of reads in any x bp sliding window.
	"clsize=i"=>\$threshold_clustersize,# minimum size of a piRNA cluster [bp]
	"clhits=i"=>\$threshold_hits,# minimum number of sequence hits per cluster
	"clhitsn=f"=>\$threshold_hits_normalized,# minimum number of mapped reads (normalized) per cluster
	"clstrand=f"=>\$threshold_strandbias,# fraction of normalized reads on mainstrand
	"clsplit=f"=>\$bidirectionality_split,# minimum fraction of normalized hits per cluster arm (for bi-directional clusters)
	"nohtml"=>\$html_files,# make html image file for each cluster [0=no/1=yes]
	"notable"=>\$results_table,# make a summary table with data for each cluster
	"nofaspi"=>\$fasta_piRNA_files,# make a fasta file for each cluster comprising mapped sequences [0=no/1=yes]
	"nofascl"=>\$fasta_cluster_file,# make one fasta file comprising all cluster sequences [0=no/1=yes]
	"nogtf"=>\$gtf_for_clusters,# make one GTF file for predicted piRNA clusters [0=no/1=yes]
	"nomotif"=>\$search_bindingsites,# search sequence motifs in piRNA clusters (e.g. binding sites for transcription factors like A-MYB)
	"flank=i"=>\$flank_size,# flanking sequence [bp]
	"pti"=>\$save_info,# save information from map file and genome file in ~.pTi for faster future computation
	);
if($print_help||!$map){print$help_desk;exit;}

###   CHECK INPUT FILES   ###
if(!-e$map){print"\nERROR. Map file $map does not exist.\n\n";exit;}
if(!-e$genome){print"\nERROR. Genome file $genome does not exist.\n\n";exit;}
if($RMannotation){if(!-e$RMannotation){print"\nERROR. RepeatMasker file $RMannotation does not exist. Continue nevertheless...\n";}}else{$RMannotation="n.a.";}
if($GeneSet){if(!-e$GeneSet){print"\nERROR. GeneSet file $GeneSet does not exist. Continue nevertheless...\n";}}else{$GeneSet="n.a.";}

###   CHECK OPTIONS   ###
if($threshold_density_absolute>0){$exp_or_obs=0;$options_info="Minimum hit density: $threshold_density_absolute hits/kb";}
if($accuracy){$exp_or_obs=1;$options_info="Significant (p<=$threshold_density_pvalue) hit density will be estimated based\non $accuracy random 1kb sliding windows.";}
elsif($threshold_density_absolute==0){$exp_or_obs=2;$options_info="Significant (p<=$threshold_density_pvalue) hit density will be calculated based\non observed hit distribution.";}
if(!$normalize_by_hit_number){$normalize_by_hit_number=1;$normalization_in_words="normalization: 1*(reads)";}else{$normalize_by_hit_number=0;$normalization_in_words="normalization: 1";}
if(!$normalize_by_read_number){$normalize_by_read_number=1;$normalization_in_words.="/(genomic hits)";}else{$normalize_by_read_number=0;}
if(!$normalize_by_total_number_of_mapped_reads){$normalize_by_total_number_of_mapped_reads=1;$normalization_in_words.="/((total mapped reads)*10^6)";}else{$normalize_by_total_number_of_mapped_reads=0;}$normalization_in_words=~s/1$/no/;
if(!$html_files){$html_files=1;}else{$html_files=0;}
if(!$results_table){$results_table=1;}else{$results_table=0;}
if(!$fasta_piRNA_files){$fasta_piRNA_files=1;}else{$fasta_piRNA_files=0;}
if(!$fasta_cluster_file){$fasta_cluster_file=1;}else{$fasta_cluster_file=0;}
if(!$gtf_for_clusters){$gtf_for_clusters=1;}else{$gtf_for_clusters=0;}
if(!$search_bindingsites){$search_bindingsites=1;}else{$search_bindingsites=0;}
%YES_or_NO=("0"=>"no","1"=>"yes");
if(!$avoid_peaks)
	{
	$avoid_peaks[0]=5;
	$avoid_peaks[1]=90;
	}
else
	{
	if($avoid_peaks=~/^\d+-\d+$/||$avoid_peaks=~/\d+\.\d+-\d+\.\d+/||$avoid_peaks=~/\d+-\d+\.\d+/||$avoid_peaks=~/\d+\.\d+-\d+/)
		{
		@avoid_peaks=split('-',$avoid_peaks);
		}
	else
		{
		print"\nERROR: Option -distr requires format: [>0..<100]-[>0..<100], e.g: 1-90\nValues can be integers or floating point numbers\n\n";
		exit;
		}
	}
if(!$no_spikes)
	{
	$no_spikes[0]=50;
	$no_spikes[1]=1000;
	}
else
	{
	if($no_spikes=~/^\d+-\d+$/||$no_spikes=~/\d+\.\d+-\d+\.\d+/||$no_spikes=~/\d+-\d+\.\d+/||$no_spikes=~/\d+\.\d+-\d+/)
		{
		@no_spikes=split('-',$no_spikes);
		}
	else
		{
		print"\nERROR: Option -spike requires format: [>0..<100]-[>0..inf], e.g: 50-1000\nValues can be integers or floating point numbers\n\n";
		exit;
		}
	}

$proTRAC_runinfo=
"
                                              /\\
                _______________________/\\___ /  \\_______
               I                      /  \\  /    \\      I
               I     pro             /    \\/      \\     I
               I        TRAC        /              \\    I
               I   ________________/________________\\_  I
               I   \\              /                     I
               I    \\            /                      I
               I     \\  /\\      /       V.$proTRAC_version         I
               I      \\/  \\    /                        I
               I___________\\  /_________________________I
                            \\/


================================= proTRAC ====================================
VERSION: .......... $proTRAC_version
LAST MODIFIED: .... $proTRAC_date

Please cite:
Rosenkranz D, Zischler H. proTRAC - a software for probabilistic piRNA cluster
detection, visualization and analysis. 2012. BMC Bioinformatics 13:5.


Contact:
David Rosenkranz
Institute of Organismic and Molecular Evolutionary Biology
Dept. Anthropology, small RNA group
Johannes Gutenberg University Mainz
email: rosenkranz\@uni-mainz.de

You can find the latest proTRAC version at:
http://sourceforge.net/projects/protrac/files
http://www.smallRNAgroup.uni-mainz.de
==============================================================================

PARAMETERS:
Map file: ...............$map
Genome file: ............$genome
RepeatMasker annotation: $RMannotation
GeneSet:.................$GeneSet

$options_info

Sliding window size: ........................................ $sliding_window_size bp
Sliding window increament: .................................. $sliding_window_increament bp
Normalize each hit by number of genomic hits: ............... $YES_or_NO{$normalize_by_hit_number}
Normalize each hit by number of sequence reads: ............. $YES_or_NO{$normalize_by_read_number}
Normalize values (-> per million mapped reads): ............. $YES_or_NO{$normalize_by_total_number_of_mapped_reads}
Min. fraction of hits with 1T(U) or 10A: .................... $threshold_T1A10
Alternatively: Min. fraction of hits with 1T(U) and 10A: .... $threshold_T1A10_both
Min. fraction of hits with typical piRNA length: ............ $threshold_size
Typical piRNA length: ....................................... ${min_piRNAsize}-${max_piRNAsize} nt
Min. size of a piRNA cluster: ............................... $threshold_clustersize bp.
Min. number of hits (absolute): ............................. $threshold_hits
Min. number of hits (normalized): ........................... $threshold_hits_normalized
Min. fraction of hits on the mainstrand: .................... $threshold_strandbias
Top fraction of mapped sequences (in terms of read counts): . $avoid_peaks[0]%
Top fraction accounts for max. n% of sequence reads: ........ $avoid_peaks[1]%
Max. percentage of reads within any x bp sliding window: .... $no_spikes[0]%
Size of sliding window with max x % of reads (see above): ... $no_spikes[1] bp
Min. fraction of hits on each arm of a bidirectional cluster: $bidirectionality_split
Output html file for each cluster: .......................... $YES_or_NO{$html_files}
Output a summary table: ..................................... $YES_or_NO{$results_table}
Output a FASTA file for each cluster (piRNA sequences): ..... $YES_or_NO{$fasta_piRNA_files}
Output a FASTA file comprising cluster sequences: ........... $YES_or_NO{$fasta_cluster_file}
Output a GTF file for predicted piRNA clusters: ..............$YES_or_NO{$gtf_for_clusters}
Search DNA motifs in clusters: .............................. $YES_or_NO{$search_bindingsites}
Output flanking sequences: +/- .............................. $flank_size bp
Output ~.pTi file: .......................................... $YES_or_NO{$save_info}
==============================================================================\n\n\n";
print"\n\n$proTRAC_runinfo";
$proTRAC_runinfo_html=$proTRAC_runinfo;
$proTRAC_runinfo_html=~s/\n/<br>/g;
$proTRAC_runinfo_html=~s/ /&nbsp;/g;

# create folder for results
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$mon++;
$year+=1900;
$map_parsed=$map;
$map_parsed=~s/^.+[\\\/]//;
$out_folder="proTRAC_${map_parsed}_${year}y${mon}m${mday}d${hour}h${min}m${sec}s";
if(-d $out_folder)
	{
	$folder_id=1;
	while(1)
		{
		$folder_id++;
		$out_folder="proTRAC_${map}_${year}y${mon}m${mday}d${hour}h${min}m${sec}s_$folder_id";
		unless(-d$out_folder)
			{
			mkdir("$out_folder")||die print"\nERROR. Could not create results folder.\n$!\n";
			last;
			}
		}
	}
else
	{
	mkdir("$out_folder")||die print"\nERROR. Could not create results folder.\n$!\n";
	}

# check input file format (SAM or ELAND3)
if(!$file_format)
	{
	print"\nCheck input file format...";
	$sam_lines=0;
	$eland3_lines=0;
	$tested_lines=0;
	open(IN,$map)||die print"\nERROR. Unable to read map file $map.\n$!\n\n";
	while(<IN>)
		{
		if($_!~/^trans_id/&&$_!~/^[#\@]/)# skip all the header lines in map files
			{
			$head_is_over=1;
			}
		next if($head_is_over==0);
		$_=~s/\s*$//;
		@d=split("\t",$_);
		$tested_lines++;
		if(looks_like_number($d[1])&&$d[1]>-1&&$d[1]<2**16&&$d[3]=~/^\d+$/&&$d[9]=~/^[ATGCUN]+$/)
			{
			$sam_lines++;
			}
		elsif($d[1]=~/^\d+$/&&$d[2]=~/^[ATGCUN]+$/&&looks_like_number($d[3])&&$d[4]=~/^[ATGCUN]+$/&&$d[5]=~/^\d+$/&&$d[6]=~/^[\+-]$/)
			{
			$eland3_lines++;
			}
		else
			{
			print"
ERROR. Input map file is not valid SAM or ELAND3 format.
Use a popular aligner such as STAR, bowtie or bowtie2 to
create map files in SAM format. Alternatively you can use
sRNAmapper to create valid map files in ELAND3 format.\n\n";
			exit;
			}
		if($sam_lines==1000)
			{
			if($tested_lines>1000)
				{
				print"\nATTENTION: Map file seems to have SAM format but some lines look like ELAND3. Continue nevertheless...";
				}
			$file_format='SAM';
			last;
			}
		elsif($eland3_lines==1000)
			{
			if($tested_lines>1000)
				{
				print"\nATTENTION: Map file seems to have ELAND3 format but some lines look like SAM. Continue nevertheless...";
				}
			$file_format='ELAND3';
			last;
			}
		last if($tested_lines==2000);
		}
	}
if($tested_lines<1000&&!$file_format)
	{
	print"\nERROR. Cannot check file format for $map. Less than 1000 hits).\nUse the option -format to specify input file format (SAM or ELAND3).\n\n";
	exit;
	}
close IN;
print" $file_format format.";

# save flags that refer to '-' or '+' strand and mapped/unmapped reads in SAM files 
if($file_format eq'SAM')
	{foreach$bit1(0,1)
	{foreach$bit2(0,2)
	{foreach$bit3(0,4)
	{foreach$bit4(0,8)
	{foreach$bit5(0,16)
	{foreach$bit6(0,32)
	{foreach$bit7(0,64)
	{foreach$bit8(0,128)
	{foreach$bit9(0,256)
	{foreach$bit10(0,512)
	{foreach$bit11(0,1024)
	{foreach$bit12(0,2048)
		{
		$flag=$bit1+$bit2+$bit3+$bit4+$bit5+$bit6+$bit7+$bit8+$bit9+$bit10+$bit11+$bit12;
		if($bit3==4)
			{
			$sam_mapped{$flag}=0;
			}
		else
			{
			$sam_mapped{$flag}=1;
			}
		if($bit5==16)
			{
			$sam_strand{$flag}='-';
			}
		else
			{
			$sam_strand{$flag}='+';
		}}}}}}}}}}}}}
	
	print"\nCreate sorted ELAND3 map file:";
	
	# check if SAM was made with collapsed data
	open(SAM,$map)||die print"\nERROR. Unable to open input file $input.\n$!\n\n";
	$extracted_sequences_nonid=0;
	$extracted_sequence_reads=0;
	$all_identifiers_are_numeric=1;
	while(<SAM>)
		{
		next if($_=~/^@/); # skip SAM header lines
		@sam=split("\t",$_);
		next if($sam_mapped{$sam[1]}==0||$sam[3]==0); # skip unmapped reads
		if(!looks_like_number($sam[0]))
			{
			$all_identifiers_are_numeric=0;
			}
		# reverse complement if sequence was mapped to the minus strand (SAM gives sequence of the locus, not of the mapped sequence)
		if($sam_strand{$sam[1]}eq'-')
			{
			$sam[9]=reverse$sam[9];
			$sam[9]=~tr/ATGC/TACG/;
			}
		$extracted_sequences{$sam[9]}{$sam[0]}++;
		}
	print"\n  -> Identifiers are numeric. Assume them to represent read counts.";
	close SAM;
	if($all_identifiers_are_numeric==0)
		{
		foreach$seq(keys%extracted_sequences)
			{
			$read_counts=keys%{$extracted_sequences{$seq}};
			print EXTRACTED">$read_counts\n$seq\n";
			$input_map_seq{$seq}=$read_counts;
			$extracted_sequence_reads+=$read_counts;
			$extracted_sequences_nonid++;
			}
		}
	elsif($all_identifiers_are_numeric==1)
		{
		foreach$seq(keys%extracted_sequences)
			{
			foreach$val(keys%{$extracted_sequences{$seq}}) # there will be only one key!
				{
				$read_counts=$val;
				}
			print EXTRACTED">$read_counts\n$seq\n";
			$input_map_seq{$seq}=$read_counts;
			$extracted_sequence_reads+=$read_counts;
			$extracted_sequences_nonid++;
			}
		}
	print"\n  -> $extracted_sequence_reads reads. $extracted_sequences_nonid sequences.";
	
	# create a sorted ELAND3 map
	%lines_in_map=();
	open(INPUT_MAP,$map)||die print"\nERROR. Unable to open input map file $map.\n$!\n\n";
	$i=0;
	while(<INPUT_MAP>)
		{
		next if($_=~/^@/); # skip SAM header lines
		@sam=split("\t",$_);
		next if ($sam_mapped{$sam[1]}==0||$sam[3]==0); # skip unmapped reads
		$i++;
		# reverse complement if sequence was mapped to the minus strand (SAM gives sequence of the locus, not of the mapped sequence)
		if($sam_strand{$sam[1]}eq'-')
			{
			$sam[9]=reverse$sam[9];
			$sam[9]=~tr/ATGC/TACG/;
			}
		$lines_in_map{$sam[2]}{$sam[3]}{$input_map_seq{$sam[9]}}{$sam[9]}=$sam_strand{$sam[1]};
		}
	close INPUT_MAP;
	print"\n  -> Loaded $i genomic hits to memory.";
	
	open(SORTED_MAP,">$map.sorted.ELAND3")||die print"\nERROR. Unable to create sorted map file.\n$!\n\n";
	foreach$chr(sort{$a cmp $b}keys%lines_in_map)
		{
		foreach$pos(sort{$a<=>$b}keys%{$lines_in_map{$chr}})
			{
			foreach$reads(keys%{$lines_in_map{$chr}{$pos}})
				{
				foreach$seq(keys%{$lines_in_map{$chr}{$pos}{$reads}})
					{
					print SORTED_MAP"$chr\t$pos\t-\t$reads\t$seq\t-\t$lines_in_map{$chr}{$pos}{$reads}{$seq}\n";
					}
				}
			}
		}
	close SORTED_MAP;
	undef%lines_in_map;
	print"\n  -> Finished sorting. New map file: $map.sorted.ELAND3.";
	$map="$map.sorted.ELAND3";
	}

# read genome
print"\nread genome ($genome)...";
%list_titles_genome=();
%scaffolds_to_analyze=();
open(GENOME,$genome)||die print"Could not open genome file.\n$!\n";
$genome_size=0;
while(<GENOME>)
	{
	$_=~s/\s*$//;
	if($_=~s/^>//)
		{
		$title=$_;
		}
	else
		{
		$genome{$title}.=$_;
		if($_=~s/[NnXx-]//g)
			{
			$gap_size+=length$&;
			}
		$genome_size+=length$_;
		}
	}
close GENOME;
$number_of_scaffolds=keys%genome;
print" done.\ngenome size (without gaps): $genome_size bp\ngaps (N/X/-): $gap_size bp\nNumber of Chromosomes/Scaffolds: $number_of_scaffolds\n\n";
$proTRAC_runinfo_html.="Genome size (without gaps): ............ $genome_size bp<br>Gaps (N/X/-): .......................... $gap_size bp<br>";

# read map file or info file
$read_info=$map.".pTi";
if(-e$read_info)
	{
	print"read information from $read_info...";
	open(INFO,$read_info);
	while(<INFO>)
		{
		$_=~s/\s*$//;
		@d=split("\t",$_);
		$total_reads+=$d[1];
		$genomic_hits+=$d[2];
		
		if($normalize_by_hit_number==0)
			{
			$hits_per_seq{$d[0]}=1;
			}
		else
			{
			$hits_per_seq{$d[0]}=$d[2];
			}
		
		if($normalize_by_read_number==0)
			{
			$reads_per_seq{$d[0]}=1;
			}
		else
			{
			$reads_per_seq{$d[0]}=$d[1];
			}
		}
	close INFO;
	}
else
	{
	print"\nread map file ($map)...";
	open(MAP,$map)||die print$!;
	while(<MAP>)
		{
		next if($_=~/^trans_id/||$_=~/^[#\@]/); # skip SeqMap, sRNAmapper and SAM header lines
		$_=~s/\s*$//;
		@d=split("\t",$_);
		$hits_per_strand{$d[6]}++;
		if(!$total_reads{$d[4]})
			{
			$total_reads{$d[4]}=1;
			$total_reads+=$d[3];
			}
		
		if($normalize_by_hit_number==1)
			{
			$hits_per_seq{$d[4]}++;
			}
		else
			{
			$hits_per_seq{$d[4]}=1;
			}
		
		if($normalize_by_read_number==1)
			{
			$reads_per_seq{$d[4]}=$d[3];
			}
		else
			{
			$reads_per_seq{$d[4]}=1;
			}
		}
	close MAP;
	}
$non_ID=keys%total_reads;
$hits_plus_strand=$hits_per_strand{"+"};
$hits_minus_strand=$hits_per_strand{"-"};
$genomic_hits=$hits_plus_strand+$hits_minus_strand;
print" done.\nmapped reads: $total_reads\nnon-identical sequences: $non_ID\nGenomic hits: $genomic_hits (+:${hits_plus_strand} -:${hits_minus_strand})\n\n";
$proTRAC_runinfo_html.="Mapped reads: .......................... $total_reads<br>Non-identical sequences: ............... $non_ID<br>Genomic hits: .......................... $genomic_hits<br>";

# save read counts and hit counts to .pTi-file
if($save_info==1)
	{
	unless(-e$read_info)
		{
		$save_info_out=$map.".pTi";
		open(SAVEINFO,">$save_info_out")||die print$!;
		foreach(keys%hits_per_seq)
			{
			print SAVEINFO"$_\t$reads_per_seq{$_}\t$hits_per_seq{$_}\n";
			}
		close SAVEINFO;
		}
	else
		{
		print"INFORMATION: $read_info already exists. Will not overwrite it.\n";
		}
	}

# read RepeatMasker output file
$RM_index=0;
%list_titles_RM=();
if(-e$RMannotation)
	{
	if($html_files==1)
		{
		$RM_index=1;
		print"read RepeatMasker annotation $RMannotation...";
		open(RM,$RMannotation);
		$prev_RM_loc="";
		%masked_positions_loc=();
		$masked_positions_genome=0;
		$RM_file_lines=0;
		while(<RM>)
			{
			$RM_file_lines++;
			if($_=~s/^ *\d+ +//)
				{
				@d=split(/ +/,$_);
				$s=int($d[4]/1000000);
				$loc="$d[3]#$s";
				$list_titles_RM{$d[3]}=1;
				$repeats_on_loc{$loc}++;
				$div_to_cons=$d[0]+$d[1]+$d[2];
				$RM{$loc}{$repeats_on_loc{$loc}}=
					[
					$d[4],
					$d[5],
					$d[8],
					$d[7],
					$div_to_cons
					];
				
				# count overall repeat content of the genome
				foreach$p($d[4]..$d[5])
					{
					$masked_positions_loc{$p}=1;
					}
				if($d[3]ne$prev_RM_loc)
					{
					$masked_positions_genome+=keys%masked_positions_loc;
					undef%masked_positions_loc;
					%masked_positions_loc=();
					}
				}
			}
		close RM;
		$percent_repeats=(int((($masked_positions_genome/$genome_size)*10000)+0.5))/100;
		print" done. $masked_positions_genome bp masked ($percent_repeats% of the genome excl. Ns).\n";
		# check if RM titles are consistent with genome titles
		foreach$rmtitle(keys%list_titles_RM)
			{
			unless($genome{$rmtitle})
				{
				print"Location $rmtitle in RepeatMasker file is not part of $genome!\n";
				}
			}
		undef%list_titles_RM;
		}
	}

# read Gene Set (GTF) file
$GTF_index=0;
if(-e$GeneSet)
	{
	if($html_files==1)
		{
		$GTF_index=1;
		print"read Gene Set $GeneSet...";
		open(GTF,$GeneSet);
		while(<GTF>)
			{
			if($_!~/^#/)
				{
				@d=split("\t",$_);
				if($d[2]=~/exon/||$d[2]=~/UTR/)
					{
					$s=int($d[3]/1000000);
					$loc="$d[0]#$s";
					$list_titles_GTF{$d[0]}=1;
					$genes_on_loc{$loc}++;
					$gtf_info_text="";
					
					# extract information from gtf annotation
					if($d[8]=~s/gene_id "[^"]+"//)
						{
						$gene_id=$&;
						$gene_id=~s/^[^"]+"//;
						$gene_id=~s/"//;
						}
					else
						{
						$gene_id="unknown";
						}
					if($d[8]=~s/gene_name "[^"]+"//)
						{
						$gene_name=$&;
						$gene_name=~s/^[^"]+"//;
						$gene_name=~s/"//;	
						}
					else
						{
						$gene_name="unknown";
						}
					if($d[8]=~s/transcript_id "[^"]+"//)
						{
						$transcript_id=$&;
						$transcript_id=~s/\d+//;
						$transcript_id=$&;
						}
					else
						{
						$transcript_id="unknown";
						}
					if($d[8]=~s/biotype "[^"]+"//)
						{
						$biotype=$&;
						$biotype=~s/^[^"]+"//;
						$biotype=~s/"//;
						$biotype=$biotype.", ";	
						}
					else
						{
						$biotype="unknown";
						}
					$gtf_info_text="$gene_name ($biotype$gene_id) Tr:$transcript_id";
					if($d[8]=~s/exon_number "[^"]+"//)
						{
						$exon_number=$&;
						$exon_number=~s/\d+//;
						$gtf_info_text.=" Ex:$&";
						}
					elsif($d[2]=~/UTR/)
						{
						$gtf_info_text.=" UTR";
						}
					
					$GeneSet{$loc}{$genes_on_loc{$loc}}=
						[
						$d[3],
						$d[4],
						$gtf_info_text,
						$d[6]
						];
					}
				}
			}
		close GTF;
		print" done.\n";
		# check if GTF titles are consistent with genome titles
		foreach$gtftitle(keys%list_titles_GTF)
			{
			unless($genome{$gtftitle})
				{
				print"Location $gtftitle in GeneSet file is not part of $genome!\n";
				}
			}
		undef%list_titles_RM;
		}
	}
undef%list_titles_genome;


# load sequence motifs (selection of germline expressed transcription factor binding sites)
if($search_bindingsites==1)
	{
	#$binding_motifs{'TBP'}='TATAAA';
	#$binding_motifs{'TBP rc'}='TTTATA';
	#$binding_motifs{'RFX-palindrome'}='CCTAGG';				# RFX palindrome from Horvath et al. 2004
	$binding_motifs{'RFX4_1'}='C.T[AG]GCAAC';					# top k-mer in UniProbe database
	$binding_motifs{'RFX4_1 rc'}='GTTGC[TC]A.G';				# top k-mer in UniProbe database
	$binding_motifs{'RFX4_2'}='C.T[AG]G[TA]TAC';				# top k-mer in UniProbe database
	$binding_motifs{'RFX4_2 rc'}='GTA[AT]C[TC]A.G';			# top k-mer in UniProbe database
	$binding_motifs{'Mybl1_1'}='AACCGTTA';					# top k-mer in UniProbe database
	$binding_motifs{'Mybl1_1 rc'}='TAACGGTT';					# top k-mer in UniProbe database
	$binding_motifs{'A-MYB'}='[TA]G[GA]CAGTTGG';				# motif from Li et al. Mol Cell. 2013 50(1): 67–81
	$binding_motifs{'A-MYB rc'}='CCAACTG[CT]C[AT]';			# motif from Li et al. Mol Cell. 2013 50(1): 67–81
	$binding_motifs{'Gata4'}='[CG]TTATCT';					# JASPAR database
	$binding_motifs{'Gata4 rc'}='AGATAA[CG]';					# JASPAR database
	$binding_motifs{'SPZ1'}='[AG]GGGT[AT][AT][CG]AG';			# JASPAR database
	$binding_motifs{'SPZ1 rc'}='CT[CG][AT][AT]ACCC[CT]';			# JASPAR database
	$binding_motifs{'SOX9'}='[CT][CT]ATTGTT';					# JASPAR database
	$binding_motifs{'SOX9 rc'}='AACAAT[GA][GA]';				# JASPAR database
	$binding_motifs{'Nobox'}='TAATT[GA][GC][TC]';				# JASPAR database
	$binding_motifs{'Nobox rc'}='[AG][CG][CT]AATTA';			# JASPAR database
	$binding_motifs{'Lhx8'}='[CT]TAATTA[GA]';					# JASPAR database
	$binding_motifs{'Lhx8 rc'}='[CT]TAATTA[GA]';				# JASPAR database
	$binding_motifs{'FIGLA'}='[AT][CA]CA[CG][CG]TG[TG][TA]';		# JASPAR database
	$binding_motifs{'FIGLA rc'}='[AT][AC]CA[GC][GC]TG[GT][TA]';	# JASPAR database
	$binding_motifs{'POU4F2'}='TG[CA]ATA[AT]TTAAT[GT]A';			# JASPAR database
	$binding_motifs{'POU4F2 rc'}='T[CA]ATTAA[TA]TAT[GT]CA';		# JASPAR database
	$binding_motifs{'POU2F1'}='TAT[GT][CT][AT]AAT';			# JASPAR database
	$binding_motifs{'POU2F1 rc'}='ATT[TA][GA][CA]ATA';			# JASPAR database
	$binding_motifs{'POU5F1'}='ATGCAAA';					# JASPAR database
	$binding_motifs{'POU5F1 rc'}='TTTGCAT';					# JASPAR database
	$binding_motifs{'FOXO1'}='[CG][CT]TGTTT[AT][CT]';			# JASPAR database
	$binding_motifs{'FOXO1 rc'}='[GA][TA]AAACA[GA][GC]';			# JASPAR database
	$binding_motifs{'FOXO3_hsa'}='GTAAACA[AT]';				# JASPAR database
	$binding_motifs{'FOXO3_hsa rc'}='[TA]TGTTTAC';				# JASPAR database
	$binding_motifs{'FOXO3_mmu'}='[TG][GC][TA]AAACA';			# JASPAR database
	$binding_motifs{'FOXO3_mmu rc'}='TGTTT[AT][CG][AC]';			# JASPAR database
	$binding_motifs{'FOXP1'}='GTAAACA';						# JASPAR database
	$binding_motifs{'FOXP1 rc'}='TGTTTAC';					# JASPAR database
	$binding_motifs{'Sox5'}='ATTGTT';						# JASPAR database
	$binding_motifs{'Sox5 rc'}='AACAAT';					# JASPAR database
	$binding_motifs{'Rhox11'}='[CT]G[CG]TGT[AT][AT][AT]';		# JASPAR database
	$binding_motifs{'Rhox11 rc'}='[TA][TA][TA]ACA[GC]C[GA]';		# JASPAR database
	$binding_motifs{'RHOXF1'}='T[AG]A[TG]C[CT]';				# JASPAR database
	$binding_motifs{'RHOXF1 rc'}='[GA]G[AC]T[TC]A';			# JASPAR database
	}

if($normalize_by_read_number==0)
	{
	$total_reads=$non_ID;
	}
if($exp_or_obs==1)
	{
	# estimate significant hit density (normalized) per kb...
	print"estimating significant (p<=$threshold_density_pvalue) hit density...";
	$exp_hitdensity=$non_ID/$genome_size*1000*($total_reads/$non_ID);
	@values=values%hits_per_seq;
	foreach(1..$accuracy)
		{
		$i=0;
		foreach(1..1000)
			{
			if((rand)<=$genomic_hits/$genome_size)
				{
				$i+=$non_ID/$genomic_hits*($total_reads/$non_ID);
				}
			}
		push(@i,$i);
		}
	@i=sort{$b<=>$a}@i;
	foreach(1..$accuracy*$threshold_density_pvalue)
		{
		$last_value=shift@i;
		if($last_value>0)
			{
			$rescue_value=$last_value;
			}
		}
	$sig_hitdensity=shift@i;
	if($sig_hitdensity==0)
		{
		$sig_hitdensity=$rescue_value;
		}
	
	undef@i;
	print" done.\nexpectation: $exp_hitdensity reads/kb*\nsignificant density: $sig_hitdensity reads/kb*\n*=normalized\n\n";
	$proTRAC_runinfo_html.="Significant densitiy of mapped reads: .. $sig_hitdensity reads/kb<br><br>";
	}

elsif($exp_or_obs==2)
	{
	# or check density distribution
	print"\ncheck genomic distribution of mapped reads:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
	open(MAP,$map)||die print"Could not open map file: $map.\n$!\n";
	if($^O=~/MSWin/)
		{
		$seek=1; # -1 is necessary for windows systems because newline is represented by 2 characters.
		}
	else
		{
		$seek=0; # on linux or mac systems a newline is represented by 1 character.
		}
	$processed_genomic_hits=0;
	$printed_dots=0;
	@hit_densities=();
	@sw=('0')x$sliding_window_size;
	while(<MAP>)
		{
		next if($_=~/^trans_id/||$_=~/^[#\@]/); # skip SeqMap, sRNAmapper and SAM header lines
		@d=split("\t",$_);
		$processed_genomic_hits++;
		if($processed_genomic_hits>=$genomic_hits/71)
			{
			$processed_genomic_hits=0;
			print".";
			$printed_dots++;
			}
		
		# change of chromosome/scaffold/reference
		if($d[0]ne$prev_scaff)
			{
			# process previous sliding window 
			$sw_value=0;
			foreach$pos(0..$sliding_window_size-1)
				{
				$sw_value+=$sw[$pos];
				}
			push(@hit_densities,$sw_value);
			$sw_start=1;
			@sw=('0')x$sliding_window_size;
			seek(MAP,-(length$_)-$seek,1); # -1 is necessary for windows systems because newline is represented by 2 characters.
			$processed_genomic_hits=$processed_genomic_hits-1;
			}
		
		# hit is within the sliding window
		elsif($d[1]>=$sw_start&&$d[1]<$sw_start+$sliding_window_size)
			{
			# use standard normalized hit counts
			if(@d==7)
				{
				$sw[$d[1]-$sw_start]+=(1/$hits_per_seq{$d[4]})*$reads_per_seq{$d[4]};
				}
			# use weighted hit counts as calculated by map-arrange.pl
			elsif(@d==9)
				{
				$sw[$d[1]-$sw_start]+=$d[8];
				}
			}
		
		# process previous sliding window and move sliding window
		else
			{
			$sw_value=0;
			foreach$pos(0..$sliding_window_size-1)
				{
				$sw_value+=$sw[$pos];
				}
			push(@hit_densities,$sw_value);
			splice(@sw,0,$sliding_window_increament);
			$sw_start+=$sliding_window_increament;
			seek(MAP,-(length$_)-$seek,1); # -1 is necessary for windows systems because newline is represented by 2 characters.
			$processed_genomic_hits=$processed_genomic_hits-1;
			}
		$prev_scaff=$d[0];
		}
	close MAP;
	@hit_densities=sort{$b<=>$a}@hit_densities;
	$sig_hitdensity=$hit_densities[abs($threshold_density_pvalue*@hit_densities)];
	
	# check if significant hit density = 0
	while($sig_hitdensity==0)
		{
		$sig_hitdensity=pop@hit_densities;
		last if(@hit_densities==0);
		}
	if($sig_hitdensity==0)
		{
		die print"\nERROR: Hit density is 0 for all sliding windows!\nCheck map file.\n\n";
		}
	foreach($printed_dots..70)
		{
		print".";
		}
	print"\nsignificant (p<=$threshold_density_pvalue) density: $sig_hitdensity reads/kb*\n*=normalized\n\n";
	undef@hit_densities;
	$proTRAC_runinfo_html.="Significant densitiy of mapped reads: .. $sig_hitdensity reads/kb<br><br>";
	}
else
	{
	$sig_hitdensity=$threshold_density_absolute;
	}

# read map file with sliding window
$prev_trans_id="";
$sw_id=1-($sliding_window_increament/$sliding_window_size);
@sw=();
print"search and tag loci with significantly enriched number of mapped reads:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
open(MAP,$map)||die print$!;
$processed_genomic_hits=0;
$printed_dots=0;
while(<MAP>)
	{
	next if($_=~/^trans_id/||$_=~/^[#\@]/); # skip SeqMap, sRNAmapper and SAM header lines
	@d=split("\t",$_);
	$reads=$d[3];
	$seq=$d[4];
	$location=$d[0];
	$coordinate=$d[1];
	if(@d==7)
		{
		$read_counts=$reads_per_seq{$seq}/$hits_per_seq{$seq};
		}
	elsif(@d==9)
		{
		$read_counts=$d[8];
		}
	$processed_genomic_hits++;
	if($processed_genomic_hits>=$genomic_hits/71)
		{
		$processed_genomic_hits=0;
		print".";
		$printed_dots++;
		}
	
	# check if hit is on the same chromosome/scaffold
	if($location ne $prev_trans_id)
		{
		$sw_id=1-($sliding_window_increament/$sliding_window_size);
		if(@sw>0)
			{
			calc_and_check();
			}
		@sw=();
		}
	
	# add hit to array if coordinate is in the sliding window
	if($coordinate<=$sw_id*$sliding_window_size)
		{
		push(@sw,"$location\t$coordinate\t-\t$reads\t$seq\t-\t-\t-\t$read_counts");
		}
	
	elsif(@sw>0)
		{
		calc_and_check();
		sub calc_and_check
			{
			# calculate hit density in sliding window
			$sw_hits=0;
			%sw_coord=();
			$count=0;
			foreach$hit(@sw)
				{
				@hit=split("\t",$hit);
				$sw_coord{"$hit[0]-$hit[1]"}=1;
				$sw_hits+=$hit[8];
				}
			
			# check hit density
			if(($sw_hits/$sliding_window_size*1000)>=$sig_hitdensity&&($sw_hits/$sliding_window_size*1000)>=$threshold_density_absolute)
				{
				foreach$key(keys%sw_coord)
					{
					$all_coord{$key}=1;
					}
				}
			undef%sw_coord;
			}
		
		# move sliding window
		move_sw();
		sub move_sw
			{
			$sw_increment_steps=0;
			while($coordinate>$sw_id*$sliding_window_size)
				{
				$sw_id+=($sliding_window_increament/$sliding_window_size);
				$sw_increment_steps++;
				}
			if($sw_increment_steps==1)
				{
				@new_sw=();
				foreach$hit(@sw)
					{
					@hit=split("\t",$hit);
					if($hit[1]>($sw_id-1)*$sliding_window_size)
						{
						push(@new_sw,$hit);
						}
					}
				@sw=@new_sw;
				undef@new_sw;
				}
			else
				{
				@sw=();
				}
			push(@sw,"$location\t$coordinate\t-\t$reads\t$seq\t-\t-\t-\t$read_counts");
			}
		}
	else
		{
		move_sw();
		}
	$prev_trans_id=$location;
	}
close MAP;
foreach($printed_dots..70)
	{
	print".";
	}
print"\n\ncheck for compliance with defined piRNA cluster criteria:\n";

# print to output files
$id=0;
$stat=0;
open(MAP,$map)||die print$!;

if($fasta_cluster_file==1)
	{
	open(CLUSTER_FASTA,">$out_folder/clusters.fasta");
	}
if($results_table==1)
	{
	if($exp_or_obs==0)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb (user-defined): $threshold_density_absolute";
		}
	elsif($exp_or_obs==1)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb based on theoretically expected hit density. p=$threshold_density_pvalue -> $sig_hitdensity hits/kb\n1kb-window resampling steps to calculate p: $accuracy";
		}
	elsif($exp_or_obs==2)
		{
		$exp_or_obs_text="Threshold for normalized hits/kb based on observed hit density. p=$threshold_density_pvalue -> $sig_hitdensity hits/kb";
		}
	open(RESULTS_TEXT,">$out_folder/results.table");
	print RESULTS_TEXT$proTRAC_runinfo;
	}
if($gtf_for_clusters==1)
	{
	open(CLUSTER_GTF,">$out_folder/clusters.gtf");
	}


$total_size_of_all_clusters=0;
%total_clustered_sequences=();
@flank_buffer=();
$prev_cl_coordinate=0;
$prev_hit_location=0;
while(<MAP>)
	{
	if($_!~/^trans_id/&&$_!~/^[#\@]/) # skip SeqMap, sRNAmapper and SAM header lines
		{
		$_=~s/\s*$//;
		@hit_data=split("\t",$_);
		if(@hit_data==7)
			{
			$read_counts=$hit_data[3]/$hits_per_seq{$hit_data[4]};
			}
		elsif(@hit_data==9)
			{
			$read_counts=$hit_data[8];
			}
		$reads=$hit_data[3];
		$seq=$hit_data[4];
		$location=$hit_data[0];
		$coordinate=$hit_data[1];
		$strand=$hit_data[6];
		
		# check if chromosome/scaffold changed and last chromosome/scaffold hit is inside a cluster candidate
		if($location ne $prev_hit_location&&$stat==1)
			{
			cluster_candidate_check();
			$stat=0;
			}
		$prev_hit_location=$location;
		
		# store 5' hits if flank_size > 0
		if($flank_size>0)
			{
			push(@flank_buffer,"$location\t$coordinate\t-\t$reads\t$seq\t-\t$strand\t-\t$read_counts");
			if(@flank_buffer>$flank_size*10)
				{
				shift@flank_buffer;
				}
			}
		
		if($all_coord{"$location-$coordinate"})
			{
			# check if the gap between tagged hits is bigger than 10kb -> split
			if($hit_data[1]-$prev_cl_coordinate>=10000&&$stat==1)
				{
				cluster_candidate_check();
				$stat=0;
				}
			$prev_cl_coordinate=$coordinate;
			
			if($stat==0)
				{
				$hits_abs=0;
				$hits_norm=0;
				$cl_T1=0;
				$cl_A10=0;
				$cl_size=0;
				$id++;
				$start=$coordinate;
				@shape=();
				open(OUT,">$out_folder/$id.fasta");
				
				# check if there are 5' flanking hits if flank_size > 0
				if($flank_size>0)
					{
					pop@flank_buffer;
					foreach$flank_hit(@flank_buffer)
						{
						@d_flank=split("\t",$flank_hit);
						if($d_flank[0]eq$location&&$d_flank[1]>=$coordinate-$flank_size)
							{
							print OUT ">Location:$d_flank[0]\tCoordinate:$d_flank[1]\tReads:$reads_per_seq{$d_flank[4]}\tHits:$hits_per_seq{$d_flank[4]}\tAllocated_reads:$d_flank[8]\tStrand:$d_flank[6]\tFLANKING_SEQUENCE\n$d_flank[4]\n";
							}
						}
					}
				}
			$hits_abs++;
			$hits_norm+=$read_counts;
			
			if($seq=~/^T/)
				{
				$cl_T1+=$read_counts;
				}
			if($seq=~/^.{9}A/)
				{
				$cl_A10+=$read_counts;
				}
			if(length$seq>=$min_piRNAsize&&length$seq<=$max_piRNAsize)
				{
				$cl_size+=$read_counts;
				}
			$cl_location=$location;
			$end=$coordinate+(length$seq)-1;
			$seqs_in_candidate{$seq}=1;
			print OUT">Location:$location\tCoordinate:$coordinate\tReads:$reads\tHits:$hits_per_seq{$seq}\tAllocated_reads:$read_counts\tStrand:$strand\n$seq\n";
			push(@shape,$read_counts);
			$stat=1;
			}
		
		# check if there are 3' flanking hits if flank_size > 0
		elsif($location eq $cl_location&&$coordinate-$flank_size<=$end)
			{
			print OUT">Location:$location\tCoordinate:$coordinate\tReads:$reads\tHits:$hits_per_seq{$seq}\tAllocated_reads:$read_counts\tStrand:$strand\tFLANKING_SEQUENCE\n$seq\n";
			}
		
		else
			{
			if($stat==1)
				{
				cluster_candidate_check();
				sub cluster_candidate_check
					{
					close OUT;
					
					# check which value is lower, 1T or 10A | compare to $threshold_T1A10_both
					if($cl_A10<$cl_T1)
						{
						$low1T10A=$cl_A10/$hits_norm;
						}
					else
						{
						$low1T10A=$cl_T1/$hits_norm;
						}
					
					# calculate hit count for top fraction
					@shape=sort{$b<=>$a}@shape;
					$top_hits=0;
					foreach(0..int((@shape*($avoid_peaks[0]/100))-1))
						{
						$top_hits+=$shape[$_];
						}
					
					# check whether cluster results from one or few peaks
					if($top_hits>$hits_norm*($avoid_peaks[1]/100))
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# check cluster size
					elsif($end-$start+1<$threshold_clustersize)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# check number of hits (absolute)
					elsif($hits_abs<$threshold_hits)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# check number of hits (normalized)
					elsif($hits_norm<$threshold_hits_normalized)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# check number of hits (normalized) with typical piRNA size
					elsif($cl_size/$hits_norm<$threshold_size)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# check number of hits (normalized) with 1T or 10A
					elsif($cl_T1/$hits_norm<$threshold_T1A10&&$cl_A10/$hits_norm<$threshold_T1A10&&$low1T10A<$threshold_T1A10_both)
						{
						unlink"$out_folder/$id.fasta";
						$id=$id-1;
						}
					
					# calculate strand bias
					else
						{
						open(CHECKSTRAND,"$out_folder/$id.fasta");
						@values=();
						@coordinates_for_split=();
						$plus_strand=0;
						while(<CHECKSTRAND>)
							{
							$_=~s/\s*$//;
							if($_=~/^>/&&$_!~/FLANK/)
								{
								@head=split("\t",$_);
								foreach$element(0..5)
									{
									$head[$element]=~s/.+://;
									}
								push(@coordinates_for_split,$head[1]);
								if($head[5]eq'+')
									{
									$plus_strand+=$head[4];
									push(@values,$head[4]);
									}
								else
									{
									push(@values,($head[4])*-1);
									}
								}
							}
						close CHECKSTRAND;
						
						# calculate monodirectionality
						if($plus_strand/$hits_norm>=0.5)
							{
							$directionality="mono:plus";
							$strandbias=$plus_strand/$hits_norm;
							}
						else
							{
							$directionality="mono:minus";
							$strandbias=($hits_norm-$plus_strand)/$hits_norm;
							}
						
						# calculate bidirectionality
						@left_arm=();
						@left_fraction=();
						$fraction=0;
						$plus_strand=0;
						foreach$value(@values)
							{
							$replicate=$value;
							$fraction+=abs$replicate;
							push(@left_fraction,$fraction);
							if($value>0)
								{
								$plus_strand+=$value;
								}
							push(@left_arm,$plus_strand);
							}
							
						@values=reverse@values;
						
						@right_arm=();
						@right_fraction=();
						$fraction=0;
						$minus_strand=0;
						foreach$value(@values)
							{
							$replicate=$value;
							$fraction+=abs$replicate;
							push(@right_fraction,$fraction);
							if($value<0)
								{
								$minus_strand=$minus_strand-$value;
								}
							push(@right_arm,$minus_strand);
							}
						
						$best_split=-1;
						foreach$split_pos(0..@values-2)
							{
							# both arms with a minimum of normalized hits?
							if($left_fraction[$split_pos]/$hits_norm>=$bidirectionality_split&&$right_fraction[-2-$split_pos]/$hits_norm>=$bidirectionality_split)
								{
								# left arm with proper strans bias?
								if($left_arm[$split_pos]/$left_fraction[$split_pos]>=$threshold_strandbias||$left_arm[$split_pos]/$left_fraction[$split_pos]<=1-$threshold_strandbias)
									{
									# right arm with proper strans bias?
									if($right_arm[-2-$split_pos]/$right_fraction[-2-$split_pos]>=$threshold_strandbias||$right_arm[-2-$split_pos]/$right_fraction[-2-$split_pos]<=1-$threshold_strandbias)
										{
										$split_value=$left_arm[$split_pos]+$right_arm[-2-$split_pos];
										if($split_value/$hits_norm>0.5)
											{
											$plus_minus="plus-minus";
											$split_value=$split_value/$hits_norm;
											}
										else
											{
											$plus_minus="minus-plus";
											$split_value=($hits_norm-$split_value)/$hits_norm;
											}
										if($split_value>$strandbias)
											{
											$best_split=$split_pos;
											$strandbias=$split_value;
											$directionality="bi:$plus_minus";
											}
										}
									}
								}
							}
						
						# check strand bias
						if($strandbias<$threshold_strandbias)
							{
							unlink"$out_folder/$id.fasta";
							$id=$id-1;
							}
							
						# calculate spikiness
						else
							{
							$spike_reached=0;
							$spike_start=-1;
							$spike_end=0;
							$total_reads_in_this_cluster=0;
							open(CHECKSPIKE,"$out_folder/$id.fasta"); # get start end end coordinates
							while(<CHECKSPIKE>)
								{
								if($_=~/^>/&&$_!~/FLANK/)
									{
									$_=~s/Coordinate:\d+//;
									$spike_pos=$&;
									$spike_pos=~s/Coordinate://;
									$_=~s/Allocated_reads:[^\t]+//;
									$spike_reads=$&;
									$spike_reads=~s/Allocated_reads://;
									$total_reads_in_this_cluster+=$spike_reads;
									if($spike_start==-1)
										{
										$spike_start=$spike_pos;
										}
									if($spike_pos>$spike_end)
										{
										$spike_end=$spike_pos;
										}
									}
								}
							close CHECKSPIKE;
							
							$spike_start=$spike_start-int(($no_spikes[1]/2)+0.5);
							while(1) # calculate sliding windows
								{
								$spike_start=$spike_start+int(($no_spikes[1]/2)+0.5);
								$spike_sw_reads=0;
								open(CHECKSPIKE,"$out_folder/$id.fasta");
								while(<CHECKSPIKE>)
									{
									if($_=~/^>/)
										{
										$_=~s/Coordinate:\d+//;
										$spike_pos=$&;
										$spike_pos=~s/Coordinate://;
										$_=~s/Allocated_reads:[^\t]+//;
										$spike_reads=$&;
										$spike_reads=~s/Allocated_reads://;
										if($spike_pos>=$spike_start&&$spike_pos<=($spike_start+$no_spikes[1])-1)
											{
											$spike_sw_reads+=$spike_reads;
											}
										}
									}
								close CHECKSPIKE;
								if($spike_sw_reads>=$total_reads_in_this_cluster*($no_spikes[0]/100))
									{
									$spike_reached=1;
									}
								last if($spike_reached==1||$spike_start>$spike_end);
								}
						
							# check spikiness
							if($spike_reached==1)
								{
								unlink"$out_folder/$id.fasta";
								$id=$id-1;
								}
							else
								{
								print"Predicted Cluster: $cl_location $start-$end directionality: $directionality\n";
								if($fasta_cluster_file==1)
									{
									$cluster_sequence=substr($genome{$cl_location},$start-1-$flank_size,$end-$start+1+($flank_size*2));
									print CLUSTER_FASTA">$cl_location $start-$end (+-$flank_size bp) directionality: $directionality\n$cluster_sequence\n";
									}
								if($gtf_for_clusters==1)
									{
									if($directionality eq "mono:minus"){$gtf_directionality="-";}
									elsif($directionality eq "mono:plus"){$gtf_directionality="+";}
									else{$gtf_directionality=".";}
									$fraction_piRsized_reads=$cl_size/$hits_norm;
									$fraction_cluster_1U=$cl_T1/$hits_norm;
									$fraction_cluster_10A=$cl_A10/$hits_norm;
									print CLUSTER_GTF"$cl_location\tproTRAC\tpiRNA_cluster\t$start\t$end\t.\t$gtf_directionality\t.\tpiRNA cluster no: $id; directionality: $directionality; mapped sequence reads: $hits_norm; fraction 1U reads: $fraction_cluster_1U; fraction 10A reads: $fraction_cluster_10A; fraction piRNA-sized reads: $fraction_piRsized_reads\n";
									}
								foreach$cl_seq(keys%seqs_in_candidate)
									{
									$total_clustered_sequences{$cl_seq}=1;
									}
								if($results_table==1||$html_files==1)
									{
									# format some numbers for output
									$print_hits_normalized=(int(($hits_norm*10)+0.5))/10;
									$print_1T=(int((($cl_T1/$hits_norm)*1000)+0.5))/10;
									$print_10A=(int((($cl_A10/$hits_norm)*1000)+0.5))/10;
									$print_size=(int((($cl_size/$hits_norm)*1000)+0.5))/10;
									$print_strandbias=(int(($strandbias*1000)+0.5))/10;
									$print_clustersize=$end-$start+1;
									$print_density=(int((($hits_norm/$print_clustersize)*10000)+0.5))/10;
									$total_size_of_all_clusters+=$print_clustersize;
									$cl_location_for_image=$cl_location;
									
									# calculate some values for html output file
									$line_1Uor10A=110-($threshold_T1A10*100);
									$line_1Uand10A=110-($threshold_T1A10_both*100);
									$line_size=110-($threshold_size*100);
									$line_strand=110-($threshold_strandbias*100);
									
									$html_perc_1Uor10A=$threshold_T1A10*100;
									$html_perc_1Uand10A=$threshold_T1A10_both*100;
									$html_perc_size=$threshold_size*100;
									$html_perc_strand=$threshold_strandbias*100;
									
									$html='
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<title>piRNA cluster '."$id".'</title>
<style type="text/css"> #outer{ position:fixed; top:0; left:0; width:100%; height:100%; overflow: scroll; } #inner{width: 1000px; height: auto; margin: 20px auto; position: relative;} </style>
</head>
<body>
<div id="outer" style="background-color: white;">
<div id="inner" style="background-color: white; padding: 20px;">

<big>Predicted piRNA cluster no. '."$id".'</big>
<br>

<div id="block_runinfo" style="height:auto; position:relative; align:center">
<br>
<span style="display: block; cursor: pointer;" id="show_runinfo" onclick="document.getElementById(\'runinfo\').style.display=\'block\';document.getElementById(\'show_runinfo\').style.display=\'none\';document.getElementById(\'hide_runinfo\').style.display=\'block\';"><font color="#599C54" size="3"><u>Show proTRAC run info</u></font></span>
<span style="display: none;  cursor: pointer;" id="hide_runinfo" onclick="document.getElementById(\'runinfo\').style.display=\'none\' ;document.getElementById(\'hide_runinfo\').style.display=\'none\';document.getElementById(\'show_runinfo\').style.display=\'block\';"><font color="#C0696F" size="3"><u>Hide proTRAC run info</u></font></span>
<div id="runinfo" style="display: none; background-color:#EFF4FF; font-family:Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;">

<small>
'."$proTRAC_runinfo_html".'
</small>

</div>
</div>


<div id="block_clusterinfo" style="height:auto; position:relative; align:center">
<br>
<span style="display: block; cursor: pointer;" id="show_clusterinfo" onclick="document.getElementById(\'clusterinfo\').style.display=\'inline-block\';document.getElementById(\'show_clusterinfo\').style.display=\'none\';document.getElementById(\'hide_clusterinfo\').style.display=\'inline-block\';"><font color="#599C54" size="3"><u>Show proTRAC cluster info</u></font></span>
<span style="display: none;  cursor: pointer;" id="hide_clusterinfo" onclick="document.getElementById(\'clusterinfo\').style.display=\'none\' ;document.getElementById(\'hide_clusterinfo\').style.display=\'none\';document.getElementById(\'show_clusterinfo\').style.display=\'inline-block\';"><font color="#C0696F" size="3"><u>Hide proTRAC cluster info</u></font></span>
<div id="clusterinfo" style="display: none; background-color:#EFF4FF;">

<div style="float:left; position:relative; width:500px;">
 <table style="width:100%">
  <tr>
    <td>Location</td>
    <td>'."$cl_location_for_image".'</td>
  </tr>
  <tr>
    <td>Coordinates</td>
    <td>'."$start-$end".'</td>
  </tr>
    <tr>
    <td>Size [bp]</td>
    <td>'."$print_clustersize".'</td>
  </tr>
  <tr>
    <td>Sequence hit loci</td>
    <td>'."$hits_abs".'</td>
  </tr>
    <tr>
    <td>Mapped reads (normalized)</td>
    <td>'."$print_hits_normalized".'</td>
  </tr>
  <tr>
    <td>Mapped reads (normalized) per kb</td>
    <td>'."$print_density".'</td>
  </tr>
    <tr>
    <td>Normalized reads with 1T (1U)</td>
    <td>'."$print_1T".'%</td>
  </tr>
  <tr>
    <td>Normalized reads with 10A</td>
    <td>'."$print_10A".'%</td>
  </tr>
    <tr>
    <td>Normalized reads with length'." $min_piRNAsize-$max_piRNAsize".' nt</td>
    <td>'."$print_size".'%</td>
  </tr>
  <tr>
    <td>Normalized reads on the main strand(s)</td>
    <td>'."$print_strandbias".'%</td>
  </tr>
  <tr>
    <td>Predicted directionality</td>
    <td>'."#REPLACEDIRECTIONALITY#".'</td>
  </tr>
</table>
</div>


<div style="float:left; position:relative; height:300px; width:500px;">
<div style="position:absolute; height:100px; width:1px; background-color:#505050; left:60px; top: 10px;"></div>
<div style="position:absolute; height:1px; width:5px; background-color:#505050; left:58px; top: 10px;"></div>
<div style="position:absolute; height:1px; width:340px; background-color:#505050; left:60px; top: 110px;"></div>
<div style="position:absolute; height:auto; width:29px; text-align:right; left:20px; top:8px;">100%</div>
<div style="position:absolute; height:auto; width:29px; text-align:right; left:20px; top:96px;">0%</div>
<div style="position:absolute; text-align:center; height:auto; width:70px; left:100px; top:115px;"><small>1T (1U)<br>reads</small></div>
<div style="position:absolute; text-align:center; height:auto; width:70px; left:170px; top:115px;"><small>10A reads</small></div>
<div style="position:absolute; text-align:center; height:auto; width:70px; left:240px; top:115px;"><small>'."$min_piRNAsize-$max_piRNAsize nt<br>reads".'</small></div>
<div style="position:absolute; text-align:center; height:auto; width:70px; left:310px; top:115px;"><small>reads on mainstrand</small></div>
<div style="position:absolute; height:'."$print_1T".'px; width:35px; left:117px; bottom:190px; background-color:#EB8E8E;"></div>
<div style="position:absolute; height:'."$print_10A".'px; width:35px; left:187px; bottom:190px; background-color:#90D48D;"></div>
<div style="position:absolute; height:'."$print_size".'px; width:35px; left:257px; bottom:190px; background-color:#8DB3D4;"></div>
<div style="position:absolute; height:'."$print_strandbias".'px; width:35px; left:327px; bottom:190px; background-color:#8DB3D4;"></div>
<div style="position:absolute; height:2px; width:115px; background-color:#8648AF; left:112px; top:'."$line_1Uor10A".'px;"></div>
<div style="position:absolute; height:2px; width:115px; background-color:#FFA200; left:112px; top:'."$line_1Uand10A".'px;"></div>
<div style="position:absolute; height:2px; width:45px; background-color:#646464; left:252px; top:'."$line_size".'px;"></div>
<div style="position:absolute; height:2px; width:45px; background-color:#646464; left:322px; top:'."$line_strand".'px;"></div>
<div style="position:absolute; left:60px; top:150px;"><small><b>
<span style="color:#8648AF;">Either the amount of reads with 1T (1U) <u>OR</u> 10A has to exceed '."$html_perc_1Uor10A".'% (set with option: -1Tor10A)</span><br>
<span style="color:#FFA200;;">Alternatively the amount of reads with 1T (1U) <u>AND</u> 10A has to exceed '."$html_perc_1Uand10A".'% (set with option: -1Tand10A)</span><br>
<span style="color:#646464;">Minimum amount of reads with preferred size is '."$html_perc_size".'% (set with option: -pisize)</span><br>
<span style="color:#646464;">Minimum amount of reads on the main strand(s) is '."$html_perc_strand".'% (set with option: -clstrand)</span>
</b></small></div>

</div>
</div>


<div id="block_clustercov" style="height:auto; position:relative; align:center">
<br>
<span style="display: none; cursor: pointer;" id="show_clustercov" onclick="document.getElementById(\'clustercov\').style.display=\'inline-block\';document.getElementById(\'show_clustercov\').style.display=\'none\';document.getElementById(\'hide_clustercov\').style.display=\'inline-block\';"><font color="#599C54" size="3"><u>Show read coverage</u></font></span>
<span style="display: block;  cursor: pointer;" id="hide_clustercov" onclick="document.getElementById(\'clustercov\').style.display=\'none\' ;document.getElementById(\'hide_clustercov\').style.display=\'none\';document.getElementById(\'show_clustercov\').style.display=\'inline-block\';"><font color="#C0696F" size="3"><u>Hide read coverage</u></font></span>
<div id="clustercov" style="display:inline-block; width:1000px;">

#CONTAINER1#
#CONTAINER2#

</div>
</div>





<br>
</div>
</div>

';
								
									$html_annotation='
<div style="position:relative; float:right; width:400px;">
Gene Set Annotation<br>
<div style="width:400px; height:120px; padding:2px; overflow:scroll; border-style:solid; border-width:1px; border-color:#C9C9C9;">#REPLACEGENESET#</div>
<br>RepeatMasker Annotation<br>
<div style="width:400px; height:120px; padding:2px; overflow:scroll; border-style:solid; border-width:1px; border-color:#C9C9C9;">#REPLACEREPEATMASKER#</div>
<br>Transcription Factor Binding Sites<br>
<div style="width:400px; height:120px; padding:2px; overflow:scroll; border-style:solid; border-width:1px; border-color:#C9C9C9;">#REPLACETFBS#</div>
</div>
';
								
									# draw RepeatMasker annotation
									$html_RMline="";
									$html_RMannotation_box="";
									if($RM_index==1)
										{
										@html_RM_colors_plus=("#0035C6","#2656D8","#4F77E5","#7394EE","#90ABF6","#ADC2FD","#C0D1FF","#D0DDFF");
										@html_RM_colors_minus=("#D40000","#E22B2B","#E64A4A","#EC6A6A","#F58888","#FCA5A5","#FFBDBD","#FFD1D1");
										
										$s=int($start/1000000);
										$loc="$cl_location#$s";
										$repeat_id=0;
										$repeat_in_cl=0;
										draw_RM();
										sub draw_RM
											{
											while(1)
												{
												$repeat_id++;
												if($RM{$loc}{$repeat_id})
													{
													$ref=$RM{$loc}{$repeat_id};
													@rm=@$ref;
													$cluster_plus_flank=$end-$start+1+($flank_size*2);
													# repeat inside piRNA cluster
													if($rm[0]>=$start&&$rm[1]<=$end)
														{
														$rm_start=$rm[0];
														$rm_end=$rm[1];
														draw_RMelement();
														sub draw_RMelement
															{
															$repeat_in_cl++;
															if($rm[4]<2){$rm_color_code=0;}
															elsif($rm[4]<5){$rm_color_code=1;}
															elsif($rm[4]<10){$rm_color_code=2;}
															elsif($rm[4]<15){$rm_color_code=3;}
															elsif($rm[4]<20){$rm_color_code=4;}
															elsif($rm[4]<25){$rm_color_code=5;}
															elsif($rm[4]<30){$rm_color_code=6;}
															else{$rm_color_code=7;}
															
															# html
															$relative_position=$rm_start-$start+1+$flank_size;
															$pixel_pos=($relative_position/$cluster_plus_flank)*500+80;
															$pixel_width=((($rm_end-$start+1+$flank_size)/$cluster_plus_flank)*500+80)-$pixel_pos;
															if($pixel_width<=0.5)
																{
																$pixel_width=1;
																}
															if($pixel_pos+$pixel_width>581)
																{
																$pixel_width=581-$pixel_pos;
																}
															if($rm[3]=~/\+/)
																{
																$html_RMannotation_box.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\'RME'."$repeat_in_cl".'\').style.display=\'block\';" onmouseout="document.getElementById(\'RME'."$repeat_in_cl".'\').style.display=\'none\';"><small><b>'."$repeat_in_cl. $rm[2]".'</b>: '."$rm[0]-$rm[1] (+)".', Divergence to consensus: '."$rm[4]".'%</small></span><br>'."\n";
																$html_RMline.='<div style="width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:330px; background-color:'."$html_RM_colors_plus[$rm_color_code]".';"></div>'."\n";
																$html_RMline.='<div id="'."RME$repeat_in_cl".'" style="display:none; z-index:2; width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:330px; background-color:#FFE600;"></div>'."\n";
																}
															elsif($rm[3]=~/C/)
																{
																$html_RMannotation_box.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\'RME'."$repeat_in_cl".'\').style.display=\'block\';" onmouseout="document.getElementById(\'RME'."$repeat_in_cl".'\').style.display=\'none\';"><small><b>'."$repeat_in_cl. $rm[2]".'</b>: '."$rm[0]-$rm[1] (-)".', Divergence to consensus: '."$rm[4]".'%</small></span><br>'."\n";
																$html_RMline.='<div style="width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:330px; background-color:'."$html_RM_colors_minus[$rm_color_code]".';"></div>'."\n";
																$html_RMline.='<div id="'."RME$repeat_in_cl".'" style="display:none; z-index:2; width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:330px; background-color:#FFE600;"></div>'."\n";
																}
															}
														}
													# repeat overlaps the 5' end of the piRNA cluster
													elsif($rm[0]<=$start&&$rm[1]>=$start)
														{
														$rm_start=$start;
														$rm_end=$rm[1];
														draw_RMelement();
														}
													# repeat overlaps the 3' end of the piRNA cluster
													elsif($rm[0]<=$end&&$rm[1]>=$end)
														{
														$rm_start=$rm[0];
														$rm_end=$end;
														draw_RMelement();
														}
													}
												else
													{
													last;
													}
												}
											}
										if(int($start/1000000)<int($end/1000000))
											{
											$s=int($end/1000000);
											$loc="$cl_location#$s";
											draw_RM();
											}
										}
									
									# draw Gene Set (GTF) information
									$html_GTFline="";
									$html_GTFannotation_box="";
									if($GTF_index==1)
										{
										$s=int($start/1000000);
										$loc="$cl_location#$s";
										$gene_id=0;
										$gene_in_cl=0;
										draw_GTF();
										sub draw_GTF
											{
											while(1)
												{
												$gene_id++;
												if($GeneSet{$loc}{$gene_id})
													{
													$ref=$GeneSet{$loc}{$gene_id};
													@gtf=@$ref;
													$cluster_plus_flank=$end-$start+1+($flank_size*2);
													# gene inside piRNA cluster
													if($gtf[0]>=$start&&$gtf[1]<=$end)
														{
														$gtf_start=$gtf[0];
														$gtf_end=$gtf[1];
														draw_GTFelement();
														sub draw_GTFelement
															{
															$gene_in_cl++;
															
															# html
															$relative_position=$gtf_start-$start+1+$flank_size;
															$pixel_pos=($relative_position/$cluster_plus_flank)*500+80;
															$pixel_width=((($gtf_end-$start+1+$flank_size)/$cluster_plus_flank)*500+80)-$pixel_pos;
															if($pixel_width<=0.5)
																{
																$pixel_width=1;
																}
															if($pixel_pos+$pixel_width>581)
																{
																$pixel_width=581-$pixel_pos;
																}
															if($gtf[3]=~/\+/)
																{
																if($gtf[2]=~/pseudogene/||$gtf[2]=~/pseudo_gene/)
																	{
																	$html_color="#ADC2FD";
																	}
																elsif($gtf[2]=~/protein_coding/)
																	{
																	$html_color="#0035C6";
																	}
																else
																	{
																	$html_color="#162757";
																	}
																$html_GTFannotation_box.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\'GTFE'."$gene_in_cl".'\').style.display=\'block\';" onmouseout="document.getElementById(\'GTFE'."$gene_in_cl".'\').style.display=\'none\';"><small><b>'."$gene_in_cl. $gtf[2]".'</b>: '."$gtf[0]-$gtf[1] (+)".'</small></span><br>'."\n";
																$html_GTFline.='<div style="width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:310px; background-color:'."$html_color".';"></div>'."\n";
																$html_GTFline.='<div id="'."GTFE$gene_in_cl".'" style="display:none; z-index:2; width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:310px; background-color:#FFE600;"></div>'."\n";
																}
															elsif($gtf[3]=~/-/)
																{
																if($gtf[2]=~/pseudogene/||$gtf[2]=~/pseudo_gene/)
																	{
																	$html_color="#FCA5A5";
																	}
																elsif($gtf[2]=~/protein_coding/)
																	{
																	$html_color="#D40000";
																	}
																else
																	{
																	$html_color="#57163D";
																	}
																$html_GTFannotation_box.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\'GTFE'."$gene_in_cl".'\').style.display=\'block\';" onmouseout="document.getElementById(\'GTFE'."$gene_in_cl".'\').style.display=\'none\';"><small><b>'."$gene_in_cl. $gtf[2]".'</b>: '."$gtf[0]-$gtf[1] (-)".'</small></span><br>'."\n";
																$html_GTFline.='<div style="width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:310px; background-color:'."$html_color".';"></div>'."\n";
																$html_GTFline.='<div id="'."GTFE$gene_in_cl".'" style="display:none; z-index:2; width:'."$pixel_width".'px; height:12px; position:absolute; left:'."$pixel_pos".'px; top:310px; background-color:#FFE600;"></div>'."\n";
																}
															}
														}
													
													# gene overlaps the 5' end of the piRNA cluster
													elsif($gtf[0]<=$start&&$gtf[1]>=$start)
														{
														$gtf_start=$start;
														$gtf_end=$gtf[1];
														draw_GTFelement();
														}
													# gene overlaps the 3' end of the piRNA cluster
													elsif($gtf[0]<=$end&&$gtf[1]>=$end)
														{
														$gtf_start=$gtf[0];
														$gtf_end=$end;
														draw_GTFelement();
														}
													}
												else
													{
													last;
													}
												}
											}
										if(int($start/1000000)<int($end/1000000))
											{
											$s=int($end/1000000);
											$loc="$cl_location#$s";
											draw_GTF();
											}
										}
									
									# calculate values per million mapped reads (rpm)
									if($normalize_by_total_number_of_mapped_reads==1)
										{
										$hits_norm=($hits_norm/$total_reads)*1000000;
										$print_density=($print_density/$total_reads)*1000000;
										}
									
									if(length$cl_location>39)
										{
										$cl_location_for_image=substr($cl_location,0,34);
										$cl_location_for_image.='[...]';
										}
									
									if($best_split>-1)
										{
										$split1=$coordinates_for_split[$best_split];
										$split2=$coordinates_for_split[$best_split+1];
										$directionality.=" (split between $split1 and $split2)";
										}
									$html=~s/#REPLACEDIRECTIONALITY#/$directionality/;
									
									if($results_table==1)
										{
										print RESULTS_TEXT"Cluster $id\tLocation: $cl_location\tCoordinates: $start-$end\tSize [bp]: $print_clustersize\tHits (absolute): $hits_abs\tHits (normalized): $hits_norm\tHits (normalized) per kb: $print_density\tNormalized hits with 1T: $print_1T%\tNormalized hits with 10A: $print_10A%\tNormalized hits $min_piRNAsize-$max_piRNAsize nt: $print_size%\tNormalized hits on the main strand(s): $print_strandbias%\tPredicted directionality: $directionality";
										}
									
									$start_with_flank=$start-$flank_size;
									$end_with_flank=$end+$flank_size;
									$html_topo='
<div style="float:left; position:relative; width:600px;">

<div id=\'INFO1\' style="display:none; z-index:2; width:500px; height:71px; position:absolute; left:79px; top:20px; background-color:#000000; color:#FFFFFF; opacity:0.75; line-height:80%; font-size:11px; padding:2px; text-align:justify">WHAT DO I SEE HERE?<br>This chart shows the location of mapped sequence reads within a predicted piRNA cluster. The color refers to the number of genomic hits produced by the sequence read in question. A dark red bar indicates that this sequence read produces many other hits elsewhere in the genome. Many adjacent red or yellow bars can indicate the presence of a multi-copy element such as transposons or rRNA genes. A dark green bar indicates that this sequence read maps uniquely to this locus.</div>
<div style="opacity:0.0; z-index:2; width:600px; height:71px; position:absolute; left:0px; top:20px;" onmouseover="document.getElementById(\'INFO1\').style.display=\'inline-block\';" onmouseout="document.getElementById(\'INFO1\').style.display=\'none\';"></div>

<div style="width:10px; height:10px; position:absolute; left:2px; top:20px; background-color:#0D9011;"></div><div style="position:absolute; left:15px; top:20px; font-size:9px">1 hit</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:30px; background-color:#05B70B;"></div><div style="position:absolute; left:15px; top:30px; font-size:9px">2-5 hits</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:40px; background-color:#6AE102;"></div><div style="position:absolute; left:15px; top:40px; font-size:9px">6-10 hits</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:50px; background-color:#CAE102;"></div><div style="position:absolute; left:15px; top:50px; font-size:9px">11-20 hits</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:60px; background-color:#E1CA02;"></div><div style="position:absolute; left:15px; top:60px; font-size:9px">21-50 hits</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:70px; background-color:#E17902;"></div><div style="position:absolute; left:15px; top:70px; font-size:9px">51-100 hits</div>
<div style="width:10px; height:10px; position:absolute; left:2px; top:80px; background-color:#E10202;"></div><div style="position:absolute; left:15px; top:80px; font-size:9px">&gt; 100 hits</div>

<div style="position:absolute; left:0px; top:5px; width:80px; text-align:right; font-size:9px">'."$cl_location_for_image".'&nbsp;</div>
<div style="position:absolute; left:80px; top:5px; text-align:left; font-size:9px">'."$start_with_flank".'</div>
<div style="position:absolute; left:300px; top:5px; width:280px; text-align:right; font-size:9px">'."$end_with_flank".'</div>

<div style="width:500px; height:71px; position:absolute; left:80px; top:20px; background-color:#F5F5F5;"></div>
#REPLACEFLANK1#
<div style="width:500px; height:1px; position:absolute; left:80px; top:55px; background-color:#000000"></div>
<div style="width:1px; height:71px; position:absolute; left:79px; top:20px; background-color:#000000"></div>

<div style="width:500px; height:201px; position:absolute; left:80px; top:100px; background-color:#F5F5F5;"></div>
#REPLACEFLANK2#
<div style="width:500px; height:1px; position:absolute; left:80px; top:200px; background-color:#000000"></div>
<div style="width:1px; height:201px; position:absolute; left:79px; top:100px; background-color:#000000"></div>

<div style="width:500px; height:12px; position:absolute; left:80px; top:310px; background-color:#F5F5F5;"></div><div style="position:absolute; text-align:right; left:0px; width:75px; top:310px; font-size:10px">Gene Set</div>
<div style="width:500px; height:12px; position:absolute; left:80px; top:330px; background-color:#F5F5F5;"></div><div style="position:absolute; text-align:right; left:0px; width:75px; top:330px; font-size:10px">RepeatMasker</div>

#REPLACESPLITQUBES#
#REPLACESCALE#
#REPLACETOPO#
#REPLACECOVERAGE#

#REPLACERMLINE#
#REPLACEGTFLINE#
#REPLACETFBSPOINTS#

<div style="height:10px; position:absolute; left:80px; top:350px; font-size:10px;">RepeatMasker Color Code</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:360px; font-size:9px; text-align:center;"><b>+</b></div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:370px; background-color:#0035C6;"></div><div style="position:absolute; left:110px; top:370px; font-size:10px">100-98% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:380px; background-color:#2656D8;"></div><div style="position:absolute; left:110px; top:380px; font-size:10px">&lt;98-95% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:390px; background-color:#4F77E5;"></div><div style="position:absolute; left:110px; top:390px; font-size:10px">&lt;95-90% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:400px; background-color:#7394EE;"></div><div style="position:absolute; left:110px; top:400px; font-size:10px">&lt;90-85% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:410px; background-color:#90ABF6;"></div><div style="position:absolute; left:110px; top:410px; font-size:10px">&lt;85-80% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:420px; background-color:#ADC2FD;"></div><div style="position:absolute; left:110px; top:420px; font-size:10px">&lt;80-75% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:430px; background-color:#C0D1FF;"></div><div style="position:absolute; left:110px; top:430px; font-size:10px">&lt;75-70% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:80px; top:440px; background-color:#D0DDFF;"></div><div style="position:absolute; left:110px; top:440px; font-size:10px">&lt;70% Identity</div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:360px; font-size:9px; text-align:center;"><b>-</b></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:370px; background-color:#D40000;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:380px; background-color:#E22B2B;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:390px; background-color:#E64A4A;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:400px; background-color:#EC6A6A;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:410px; background-color:#F58888;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:420px; background-color:#FCA5A5;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:430px; background-color:#FFBDBD;"></div>
<div style="width:10px; height:10px; position:absolute; left:95px; top:440px; background-color:#FFD1D1;"></div>

<div style="height:10px; position:absolute; left:250px; top:350px; font-size:10px;">Gene Set Color Code</div>
<div style="width:10px; height:10px; position:absolute; left:250px; top:360px; font-size:9px; text-align:center;"><b>+</b></div>
<div style="width:10px; height:10px; position:absolute; left:250px; top:370px; background-color:#0035C6;"></div><div style="position:absolute; left:280px; top:370px; font-size:10px">Gene</div>
<div style="width:10px; height:10px; position:absolute; left:250px; top:380px; background-color:#ADC2FD;"></div><div style="position:absolute; left:280px; top:380px; font-size:10px">Pseudogene</div>
<div style="width:10px; height:10px; position:absolute; left:250px; top:390px; background-color:#162757;"></div><div style="position:absolute; left:280px; top:390px; font-size:10px">Other</div>
<div style="width:10px; height:10px; position:absolute; left:265px; top:360px; font-size:9px; text-align:center;"><b>-</b></div>
<div style="width:10px; height:10px; position:absolute; left:265px; top:370px; background-color:#D40000;"></div>
<div style="width:10px; height:10px; position:absolute; left:265px; top:380px; background-color:#FCA5A5;"></div>
<div style="width:10px; height:10px; position:absolute; left:265px; top:390px; background-color:#57163D;"></div>

<div style="height:10px; position:absolute; left:410px; top:350px; font-size:10px;">Topology/Coverage Color Code</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:370px; background-color:#1829BD;"></div><div style="position:absolute; left:425px; top:370px; font-size:10px">Coverage Plus Strand</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:380px; background-color:#BD1818;"></div><div style="position:absolute; left:425px; top:380px; font-size:10px">Coverage Minus Strand</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:395px; background-color:#E1F0FF;"></div><div style="position:absolute; left:425px; top:395px; font-size:10px">Mainstrand: Plus</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:405px; background-color:#FFD7D7;"></div><div style="position:absolute; left:425px; top:405px; font-size:10px">Mainstrand: Minus</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:420px; background-color:#F5F5F5;"></div><div style="position:absolute; left:425px; top:420px; font-size:10px">Complementary Strand</div>
<div style="width:10px; height:10px; position:absolute; left:410px; top:430px; background-color:#8E8E8E;"></div><div style="position:absolute; left:425px; top:430px; font-size:10px;">Flanking Region<br>(if option -flank &gt;0)</div>

</div>

';
								
								
								
									$cluster_plus_flank=$end-$start+1+($flank_size*2);
									if($best_split>-1)
										{
										$split_coordinate=(($split1+$split2)/2)-$start+$flank_size;
										$html_split_coordinate=($split_coordinate/$cluster_plus_flank)*500+80;
										if($directionality=~/plus-minus/)
											{
											# calculate some html values
											$html_start_blue=80+($flank_size/$cluster_plus_flank*500);
											$html_width_blue=$html_split_coordinate-80-($flank_size/$cluster_plus_flank*500);
											$html_width_red=580-$html_split_coordinate-($flank_size/$cluster_plus_flank*500);
											
											$html_splitqubes='
<div style="width:'."$html_width_blue".'px; height:35px; position:absolute; left:'."$html_start_blue".'px; top:20px; background-color:#E1F0FF;"></div>
<div style="width:'."$html_width_red".'px; height:35px; position:absolute; left:'."$html_split_coordinate".'px; top:56px; background-color:#FFD7D7;"></div>
<div style="width:1px; height:71px; position:absolute; left:'."$html_split_coordinate".'px; top:20px; background-color:#A4A4A4;"></div>

<div style="width:'."$html_width_blue".'px; height:100px; position:absolute; left:'."$html_start_blue".'px; top:100px; background-color:#E1F0FF;"></div>
<div style="width:'."$html_width_red".'px; height:100px; position:absolute; left:'."$html_split_coordinate".'px; top:201px; background-color:#FFD7D7;"></div>
<div style="width:1px; height:201px; position:absolute; left:'."$html_split_coordinate".'px; top:100px; background-color:#A4A4A4;"></div>
';
											$html_topo=~s/#REPLACESPLITQUBES#/$html_splitqubes/;
											}
										elsif($directionality=~/minus-plus/)
											{
											# calculate some html values
											$html_start_red=80+($flank_size/$cluster_plus_flank*500);
											$html_width_red=$html_split_coordinate-80-($flank_size/$cluster_plus_flank*500);
											$html_width_blue=580-$html_split_coordinate-($flank_size/$cluster_plus_flank*500);
											
											$html_splitqubes='
<div style="width:'."$html_width_blue".'px; height:35px; position:absolute; left:'."$html_split_coordinate".'px; top:20px; background-color:#E1F0FF;"></div>
<div style="width:'."$html_width_red".'px; height:35px; position:absolute; left:'."$html_start_red".'px; top:56px; background-color:#FFD7D7;"></div>
<div style="width:1px; height:71px; position:absolute; left:'."$html_split_coordinate".'px; top:20px; background-color:#A4A4A4;"></div>

<div style="width:'."$html_width_blue".'px; height:100px; position:absolute; left:'."$html_split_coordinate".'px; top:100px; background-color:#E1F0FF;"></div>
<div style="width:'."$html_width_red".'px; height:100px; position:absolute; left:'."$html_start_red".'px; top:201px; background-color:#FFD7D7;"></div>
<div style="width:1px; height:201px; position:absolute; left:'."$html_split_coordinate".'px; top:100px; background-color:#A4A4A4;"></div>
';
											$html_topo=~s/#REPLACESPLITQUBES#/$html_splitqubes/;
											}
										}
									elsif($directionality eq"mono:plus")
										{
										$html_start_blue=80+($flank_size/$cluster_plus_flank*500);
										$html_width_blue=500-(2*($flank_size/$cluster_plus_flank*500));
										
										$html_splitqubes='
<div style="width:'."$html_width_blue".'px; height:35px; position:absolute; left:'."$html_start_blue".'px; top:20px; background-color:#E1F0FF;"></div>
<div style="width:'."$html_width_blue".'px; height:100px; position:absolute; left:'."$html_start_blue".'px; top:100px; background-color:#E1F0FF;"></div>
';
										$html_topo=~s/#REPLACESPLITQUBES#/$html_splitqubes/;
										}
									elsif($directionality eq"mono:minus")
										{
										$html_start_red=80+($flank_size/$cluster_plus_flank*500);
										$html_width_red=500-(2*($flank_size/$cluster_plus_flank*500));
										
										$html_splitqubes='
<div style="width:'."$html_width_red".'px; height:35px; position:absolute; left:'."$html_start_red".'px; top:56px; background-color:#FFD7D7;"></div>
<div style="width:'."$html_width_red".'px; height:100px; position:absolute; left:'."$html_start_red".'px; top:201px; background-color:#FFD7D7;"></div>
';
										$html_topo=~s/#REPLACESPLITQUBES#/$html_splitqubes/;
										}
									
									$html_flank1="";
									$html_flank2="";
									if($flank_size>0)
										{
										$flank_block_up_width=($flank_size/$cluster_plus_flank*500);
										$flank_block_pos_down=580-($flank_size/$cluster_plus_flank*500);
										$html_flank1='
<div style="width:'."$flank_block_up_width".'px; height:71px; position:absolute; left:80px; top:20px; background-color:#8E8E8E;"></div>
<div style="width:'."$flank_block_up_width".'px; height:71px; position:absolute; left:'."$flank_block_pos_down".'px; top:20px; background-color:#8E8E8E;"></div>
';
									$html_flank2='
<div style="width:'."$flank_block_up_width".'px; height:201px; position:absolute; left:80px; top:100px; background-color:#8E8E8E;"></div>
<div style="width:'."$flank_block_up_width".'px; height:201px; position:absolute; left:'."$flank_block_pos_down".'px; top:100px; background-color:#8E8E8E;"></div>
';
										}
									$html_topo=~s/#REPLACEFLANK1#/$html_flank1/;
									$html_topo=~s/#REPLACEFLANK2#/$html_flank2/;
									
									open(IN,"$out_folder/$id.fasta");
									%transcription_plus=();
									%transcription_minus=();
									$extreme=0;
									while(<IN>)
										{
										if($_=~/^>/)
											{
											@d=split("\t",$_);
											foreach$element(0..5)
												{
												$d[$element]=~s/.+://;
												}
											if($d[5]=~/\+/)
												{
												if($normalize_by_total_number_of_mapped_reads==0)
													{
													foreach$pos($d[1]..($d[1]+length$d[1])-1)
														{
														$transcription_plus{$pos}+=$d[4];
														if($transcription_plus{$pos}>$extreme)
															{
															$extreme=$transcription_plus{$pos};
															}
														}
													}
												# calculate values per million mapped reads (rpm)
												else
													{
													foreach$pos($d[1]..($d[1]+length$d[1])-1)
														{
														$transcription_plus{$pos}+=($d[4]/$total_reads)*1000000;
														if($transcription_plus{$pos}>$extreme)
															{
															$extreme=$transcription_plus{$pos};
															}
														}
													}
													
												}
											else
												{
												if($normalize_by_total_number_of_mapped_reads==0)
													{
													foreach$pos($d[1]..($d[1]+length$d[1])-1)
														{
														$transcription_minus{$pos}+=$d[4];
														if($transcription_minus{$pos}>$extreme)
															{
															$extreme=$transcription_minus{$pos};
															}
														}
													}
												# calculate values per million mapped reads (rpm)
												else
													{
													foreach$pos($d[1]..($d[1]+length$d[1])-1)
														{
														$transcription_minus{$pos}+=(($d[4])/$total_reads)*1000000;
														if($transcription_minus{$pos}>$extreme)
															{
															$extreme=$transcription_minus{$pos};
															}
														}
													}
												}
											
											}
										}
									close IN;
									$print_extreme=(int(($extreme*100)+0.5))/100;
									$length_extremestring=length$print_extreme;
									$extreme_poscorr=$length_extremestring*5;
									
									$html_coverage_scale='
<div style="position:absolute; text-align:center; left:0px; width:75px; top:185px; font-size:12px">Mapped<br>Reads</div>
<div style="width:5px; height:1px; position:absolute; left:77px; top:100px; background-color:#000000;"></div>
<div style="width:5px; height:1px; position:absolute; left:77px; top:300px; background-color:#000000;"></div>
<div style="position:absolute; text-align:right; left:0px; width:75px; top:100px; font-size:9px">'."$print_extreme".'</div>
<div style="position:absolute; text-align:right; left:0px; width:75px; top:145px; font-size:9px">plus strand</div>
<div style="position:absolute; text-align:right; left:0px; width:75px; top:246px; font-size:9px">minus strand</div>
<div style="position:absolute; text-align:right; left:0px; width:75px; top:291px; font-size:9px">'."$print_extreme".'</div>
';
									$html_topo=~s/#REPLACESCALE#/$html_coverage_scale/;
									
									$prev_pos=-1;
									$html_coverage='';
									$html_max_plus=0;
									$html_max_minus=0;
									foreach$pos($start-$flank_size..$end+$flank_size)
										{
										$relative_position=$pos-$start+1+$flank_size;
										$html_pixel=int((($relative_position/$cluster_plus_flank)*500+80)+0.5);
										if($transcription_plus{$pos})
											{
											# save highest value per pixel for html
											if($transcription_plus{$pos}>$html_max_plus)
												{
												$html_max_plus=$transcription_plus{$pos};
												}
											}
										if($transcription_minus{$pos})
											{
											# save highest value per pixel for html
											if($transcription_minus{$pos}>$html_max_minus)
												{
												$html_max_minus=$transcription_minus{$pos};
												}
											}
										
										# check if pixel changes for html output
										if($html_pixel>$prev_pos&&$prev_pos!=-1)
											{
											$last_cluster_pos++;
											$html_max_plus_round=(int(($html_max_plus*100)+0.5))/100;
											$html_max_minus_round=(int(($html_max_minus*100)+0.5))/100;
											$html_coverage.='<div id=\''."$id-$html_pixel".'\' style="display:none; z-index:3; width:500px; height:15px; position:absolute; left:80px; top:100px; background-color:#000000; color:#FFFFFF; opacity:0.6; font-size:11px">'."&nbsp;Region: $cl_location_for_image $last_cluster_pos-$pos. Max. coverage (+): $html_max_plus_round. Max coverage (-): $html_max_minus_round".'</div>';
											$html_coverage.='<div style="opacity:0.0; z-index:2; width:1px; height:201px; position:absolute; left:'."$html_pixel".'px; top:100px;" onmouseover="document.getElementById(\''."$id-$html_pixel".'\').style.display=\'inline-block\';" onmouseout="document.getElementById(\''."$id-$html_pixel".'\').style.display=\'none\';"></div>';
												
											if($html_max_plus>0)
												{
												$html_height=int((($html_max_plus/$extreme)*100)+0.5);
												if($html_height>0)
													{
													$html_y_coord=201-$html_height;
													$html_coverage.='<div style="width:1px; height:'."$html_height".'px; position:absolute; left:'."$html_pixel".'px; top:'."$html_y_coord".'px; background-color:#1829BD;"></div>'."\n";
													}
												}
											if($html_max_minus>0)
												{
												$html_height=int((($html_max_minus/$extreme)*100)+0.5);
												if($html_height>0)
													{
													$html_coverage.='<div style="width:1px; height:'."$html_height".'px; position:absolute; left:'."$html_pixel".'px; top:201px; background-color:#BD1818;"></div>'."\n";
													}
												}
											$html_max_plus=0;
											$html_max_minus=0;
											$last_cluster_pos=$pos;
											}
										$prev_pos=$html_pixel;
										}
									# paint last html pixel
									$last_cluster_pos++;
									$html_max_plus_round=(int(($html_max_plus*100)+0.5))/100;
									$html_max_minus_round=(int(($html_max_minus*100)+0.5))/100;
									$html_coverage.='<div id=\''."$id-$html_pixel".'\' style="display:none; z-index:3; width:500px; height:15px; position:absolute; left:80px; top:100px; background-color:#000000; color:#FFFFFF; opacity:0.6; font-size:11px">'."&nbsp;Region: $cl_location_for_image $last_cluster_pos-$pos. Max. coverage (+): $html_max_plus_round. Max coverage (-): $html_max_minus_round".'</div>';
									$html_coverage.='<div style="opacity:0.0; z-index:2; width:1px; height:201px; position:absolute; left:'."$html_pixel".'px; top:100px;" onmouseover="document.getElementById(\''."$id-$html_pixel".'\').style.display=\'inline-block\';" onmouseout="document.getElementById(\''."$id-$html_pixel".'\').style.display=\'none\';"></div>';
									
									if($html_max_plus>0)
										{
										$html_height=int((($transcription_plus{$pos}/$extreme)*100)+0.5);
										$html_y_coord=201-$html_height;
										$html_coverage.='<div style="width:1px; height:'."$html_height".'px; position:absolute; left:'."$html_pixel".'px; top:'."$html_y_coord".'px; background-color:#1829BD;"></div>'."\n";
										}
									if($html_max_minus>0)
										{
										$html_height=int((($transcription_minus{$pos}/$extreme)*100)+0.5);
										$html_coverage.='<div style="width:1px; height:'."$html_height".'px; position:absolute; left:'."$html_pixel".'px; top:201px; background-color:#BD1818;"></div>'."\n";
										}
									$html_topo=~s/#REPLACECOVERAGE#/$html_coverage/;
									
									
									undef%transcription_plus;
									undef%transcription_minus;
									
									# search binding motifs
									$html_tfbs_list="";
									$html_tfbs_points="";
									if($search_bindingsites==1)
										{
										%sites_in_cluster=();
										$sites_in_cluster="";
										$found_motifs=0;
										$found_motifs_rc=0;
										foreach$motif_name(keys%binding_motifs)
											{
											$check_sequence=$cluster_sequence;
											while(1)
												{
												if($check_sequence=~/$binding_motifs{$motif_name}/)
													{
													$hit_motif=$&;
													$position=index($check_sequence,$hit_motif);
													$position++;
													$replace="";
													foreach(1..length$hit_motif)
														{
														$replace.="X";
														}
													$check_sequence=~s/$binding_motifs{$motif_name}//;
													$check_sequence=$replace.$check_sequence;
													$print_motiv_name=$motif_name;
													$print_motiv_name=~s/ rc$//;
													if($motif_name=~/ rc$/)
														{
														unless($sites_in_cluster{"$hit_motif-$position"})
															{
															$found_motifs_rc++;
															$sites_in_cluster{"$hit_motif-$position"}=1;
															$position=$position-$flank_size;
															$sites_in_cluster.="$print_motiv_name ($hit_motif: $position) ";
															$pixel_pos=int(((($position/$cluster_plus_flank)*500)+76)+0.5);
															$position+=$start-1;
															$html_tfbs_list.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\''."$print_motiv_name-$hit_motif-$position".'\').style.display=\'block\';" onmouseout="document.getElementById(\''."$print_motiv_name-$hit_motif-$position".'\').style.display=\'none\';"><small>'."<b>$print_motiv_name</b> (Sequence: $hit_motif (-): $position)".'</small></span><br>'."\n";
															$html_tfbs_points.='<div id=\''."$print_motiv_name-$hit_motif-$position".'\' style="white-space:nowrap; display:none; z-index=2; position:absolute; left:'."$pixel_pos".'px; top:196px; width:8px; height:8px; opacity:1; background-color:#FFBC00; border-style:solid; border-color:black; border-width:1px; webkit-border-radius:4px; moz-border-radius:4px; border-radius:4px;"></div>'."\n";
															}
														}
													else
														{
														unless($sites_in_cluster{"$hit_motif-$position"})
															{
															$found_motifs++;
															$sites_in_cluster{"$hit_motif-$position"}=1;
															$position=$position-$flank_size;
															$sites_in_cluster.="$print_motiv_name ($hit_motif: $position) ";
															$sites_in_cluster{"$hit_motif-$position"}=1;
															$pixel_pos=int(((($position/$cluster_plus_flank)*500)+76)+0.5);
															$position+=$start-1;
															$html_tfbs_list.='<span style="white-space:nowrap; cursor:pointer;" onmouseover="document.getElementById(\''."$print_motiv_name-$hit_motif-$position".'\').style.display=\'block\';" onmouseout="document.getElementById(\''."$print_motiv_name-$hit_motif-$position".'\').style.display=\'none\';"><small>'."<b>$print_motiv_name</b> (Sequence: $hit_motif (+): $position)".'</small></span><br>'."\n";
															$html_tfbs_points.='<div id=\''."$print_motiv_name-$hit_motif-$position".'\' style="white-space:nowrap; display:none; z-index=2; position:absolute; left:'."$pixel_pos".'px; top:196px; width:8px; height:8px; opacity:1; background-color:#FFBC00; border-style:solid; border-color:black; border-width:1px; webkit-border-radius:4px; moz-border-radius:4px; border-radius:4px;"></div>'."\n";
															}
														}
													}
												else
													{
													last;
													}
												}
											}
										if($results_table==1&&$found_motifs+$found_motifs_rc>0)
											{
											$sites_in_cluster=~s/ $//;
											print RESULTS_TEXT"\tBinding sites: $sites_in_cluster\n";
											}
										else
											{
											print RESULTS_TEXT"\n";
											}
										undef%sites_in_cluster;
										}
									elsif($results_table==1)
										{
										print RESULTS_TEXT"\n";
										}
									
									@redundancy_code=("1-1","2-5","6-10","11-20","21-50","51-100","101-999999999");
									$html_topo_bars="";
									@html_colors=("#0D9011","#05B70B","#6AE102","#CAE102","E1CA02","E17902","E10202");
									foreach$redundancy(0..6)
										{
										open(IN,"$out_folder/$id.fasta")||print"\nUnanble to open $out_folder/$id.fasta to generate html image file.\n$!\n\n";
										$color_id=6-$redundancy;
										$range=pop@redundancy_code;
										@range=split('-',$range);
										$html_prev_pixel=-1;
										while(<IN>)
											{
											if($_=~/^>/)
												{
												$_=~s/\s*$//;
												@d=split("\t",$_);
												$d[1]=~s/Coordinate://;
												$d[3]=~s/Hits://;
												$d[5]=~s/Strand://;
												if($d[3]>=$range[0]&&$d[3]<=$range[1])
													{
													$relative_position=$d[1]-$start+1+$flank_size;
													# html
													$html_pixel=int((($relative_position/$cluster_plus_flank)*500+80)+0.5);
													if($html_pixel>$html_prev_pixel)
														{
														if($d[5]=~/\+/)
															{
															$html_topo_bars.='<div style="width:1px; height:20px; position:absolute; left:'."$html_pixel".'px; top:35px; background-color:'."$html_colors[$color_id]".';"></div>'."\n";	
															}
														else
															{
															$html_topo_bars.='<div style="width:1px; height:20px; position:absolute; left:'."$html_pixel".'px; top:56px; background-color:'."$html_colors[$color_id]".';"></div>'."\n";	
															}
														}
													$html_prev_pixel=$html_pixel;
													}
												}
											}
										close IN;
										}
									$html_topo=~s/#REPLACETOPO#/$html_topo_bars/;
									if($html_files==1)
										{
										$html=~s/#CONTAINER1#/$html_topo/;
										$html=~s/#CONTAINER2#/$html_annotation/;
										$html=~s/#REPLACEREPEATMASKER#/$html_RMannotation_box/;
										$html=~s/#REPLACERMLINE#/$html_RMline/;
										$html=~s/#REPLACEGENESET#/$html_GTFannotation_box/;
										$html=~s/#REPLACEGTFLINE#/$html_GTFline/;
										$html=~s/#REPLACETFBSPOINTS#/$html_tfbs_points/;
										$html=~s/#REPLACETFBS#/$html_tfbs_list/;
										$html.='</div></div></body></html>';
										open(HTML,">$out_folder/$id.html");
										print HTML $html;
										close HTML;
										}
									}
								if($fasta_piRNA_files==0)
									{
									unlink"$out_folder/$id.fasta";
									}
								}
							}
						}
					undef%seqs_in_candidate;
					} 
				}
			$stat=0;
			}
		}
	}
close OUT;
close MAP;

#add buttons to html files
if($html_files==1&&$id>0)
	{
	print"\nFinalizing html output files...";
	$file_id=0;
	while(1)
		{
		$file_id++;
		last if(!-e"$out_folder/$file_id.html");
		open(HTML,"$out_folder/$file_id.html");
		@html=<HTML>;
		close HTML;
		
		$next=$file_id+1;
		$previous=$file_id-1;
		if($next>$id)
			{
			$next=1;
			}
		if($previous<1)
			{
			$previous=$id;
			}
		
		open(FINALHTML,">$out_folder/$file_id.html");
		$lines=@html;
		foreach(@html)
			{
			if($_=~/<big>Predicted piRNA cluster no. \d+<\/big>/)
				{
				$_=~s/\s*$//;
				$_.='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="'."$previous".'.html">previous</a>&nbsp;&nbsp;&nbsp;<a href="'."$next".'.html">next</a>'."\n";
				}
			print FINALHTML$_;
			}
		close FINALHTML;
		}
	print" done.";
	}

$percent_clusters_to_genome=(int((($total_size_of_all_clusters/$genome_size)*100000)+0.5))/1000;
$total_clustered_sequences=keys%total_clustered_sequences;
$total_clustered_reads=0;
foreach(keys%total_clustered_sequences)
	{
	$total_clustered_reads+=$reads_per_seq{$_};
	}
$percent_clustered_reads=(int((($total_clustered_reads/$total_reads)*100000)+0.5))/1000;
$percent_clustered_to_all_seq=(int((($total_clustered_sequences/$non_ID)*100000)+0.5))/1000;
if($results_table==1)
	{
	print RESULTS_TEXT"\nTotal size of $id predicted piRNA clusters: $total_size_of_all_clusters bp ($percent_clusters_to_genome%)\nNon identical sequences that can be assigned to clusters: $total_clustered_sequences ($percent_clustered_to_all_seq%)\nSequence reads that can be assigned to clusters: $total_clustered_reads ($percent_clustered_reads%)";
	close RESULTS_TEXT;
	$total_masked_in_clusters="n.a.";
	if(-e$RMannotation&&$id>0) # add information on transposon content in genome/piRNA clusters to summary table
		{
		print"\n\nCheck repeat content of predicted clusters:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
		%cluster_data=();
		open(RESULTS_TEXT,"$out_folder/results.table");
		while(<RESULTS_TEXT>)
			{
			next if($_!~/^Cluster \d+/);
			$_=~s/Cluster \d+//;
			$cluster_ID=$&;
			$_=~s/Location: [^\t]+//;
			$location=$&;
			$location=~s/Location: //;
			$_=~s/Coordinates: [^\t]+//;
			$coordinates=$&;
			$coordinates=~s/Coordinates: //;
			@coordinates=split('-',$coordinates);
			foreach$p($coordinates[0]..$coordinates[1])
				{
				$cluster_data{$location}{$p}=$cluster_ID;
				}
			}
		close RESULTS_TEXT;
		
		open(RM,$RMannotation);
		%masked_bases=();
		$printed_dots=0;
		$processed_lines=0;
		while(<RM>)
			{
			$processed_lines++;
			if($processed_lines>=$RM_file_lines/71)
				{
				print".";
				$printed_dots++;
				$processed_lines=0;
				}
			next if($_!~/^ *\d+/);
			$_=~s/^ *//;
			@d=split(/ +/,$_);
			foreach$p($d[5]..$d[6])
				{
				if($cluster_data{$d[4]}{$p})
					{
					$masked_bases{$cluster_data{$d[4]}{$p}}{$p}++;
					}
				}
			}
		close RM;
		foreach($printed_dots..70)
			{
			print".";
			}
		
		open(RESULTS_TEXT,"$out_folder/results.table");
		@results_table=();
		$total_masked_in_clusters=0;
		while(<RESULTS_TEXT>)
			{
			$_=~s/\s*$//;
			if($_=~/^Cluster \d+/)
				{
				$cluster_ID=$&;
				if($masked_bases{$cluster_ID})
					{
					$masked_bases_in_cluster=keys%{$masked_bases{$cluster_ID}};
					}
				else
					{
					$masked_bases_in_cluster=0;
					}
				$_.="\tmasked [bp]: $masked_bases_in_cluster";
				$total_masked_in_clusters+=$masked_bases_in_cluster;
				}
			push(@results_table,$_);
			}
		close RESULTS_TEXT;
		open(RESULTS_TEXT,">$out_folder/results.table");
		foreach$line(@results_table)
			{
			print RESULTS_TEXT"$line\n";
			}
		close RESULTS_TEXT;
		$percent_total_masked_in_clusters=(int((($total_masked_in_clusters/$total_size_of_all_clusters)*10000)+0.5))/100;
		print"\n";
		}
	}
print"\n
Total size of $id predicted piRNA clusters: $total_size_of_all_clusters bp ($percent_clusters_to_genome%)
Non identical sequences that can be assigned to clusters: $total_clustered_sequences ($percent_clustered_to_all_seq%)
Sequence reads that can be assigned to clusters: $total_clustered_reads ($percent_clustered_reads%)

Total repeat-masked bases in clusters: $total_masked_in_clusters ($percent_total_masked_in_clusters%)
Total repeat-masked bases in genome: $masked_positions_genome ($percent_repeats%)
";
exit;