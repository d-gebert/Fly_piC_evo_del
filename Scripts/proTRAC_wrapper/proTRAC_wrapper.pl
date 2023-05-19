#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $pp_script = "TBr2_pingpong.pl";
my $pT_script = "proTRAC_2.4.4.pl";
my $p_val = 0.05;
my $pT_opts = "-clsize 5000 -pimin 23 -pimax 29 -1Tor10A 0.3 -1Tand10A 0.3 -clstrand 0.0 -clsplit 1.0 -distr 1.0-99.0 -spike 90-1000 -nomotif -pdens $p_val";
# Global variables
my %fin_pics = ();
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <map_file> <genome_file> <repeatmasker_file>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $map_file = $ARGV[0];
my $gnm_file = $ARGV[1];
my $rep_file = $ARGV[2];
# Get existing proTRAC output directories
my $map_file_end = $map_file;
$map_file_end =~ s/.*\///;
my @pT_dirs = grep -e, glob "proTRAC_$map_file_end*";
# Run proTRAC
unless (@pT_dirs) {
	print("Run proTRAC on $map_file...");
	system("perl $pT_script -map $map_file -genome $gnm_file -repeatmasker $rep_file $pT_opts >>.pTlog 2>&1");
	@pT_dirs = grep -e, glob "proTRAC_$map_file_end*";
	system("mv .pTlog $pT_dirs[-1]");
	print("done.\n");
}
# Get proTRAC latest output directory
my $pT_dir = $pT_dirs[-1];
print("Processing $pT_dir...");
# Get proTRAC output files
my $pic_gtf_file = "$pT_dir/clusters.gtf";
my $pic_tbl_file = "$pT_dir/results.table";
# Get proTRAC output data
my $pic_gtf_props = get_pic_gtf_properties($pic_gtf_file);
my $pic_tbl_props = get_pic_tbl_properties($pic_tbl_file);
# Get total number of sample reads
my $pir_total_rds = get_total_pirna_reads($pic_tbl_file);

## Define groups of clusters that are in close proximity
# Initialize merge bins with first cluster id
my @pic_merge_bins = ([1]);
# Find close neighboring pics
foreach my $pic_c (sort {$a <=> $b} keys %{$pic_gtf_props}) {
	# Get properties of current pic
	my $chr_c = $pic_gtf_props->{$pic_c}->[0];
	my $beg_c = $pic_gtf_props->{$pic_c}->[1];
	my $end_c = $pic_gtf_props->{$pic_c}->[2];
	my $len_c = $end_c-$beg_c+1;
	# Get id of next pic
	my $pic_n = $pic_c+1;
	# Skip if no next pic exists
	unless ($pic_gtf_props->{$pic_n}) { last; }
	# Get properties of next pic
	my $chr_n = $pic_gtf_props->{$pic_n}->[0];
	my $beg_n = $pic_gtf_props->{$pic_n}->[1];
	my $end_n = $pic_gtf_props->{$pic_n}->[2];
	my $len_n = $end_n-$beg_n+1;
	# Check if current and next pic are on same contig
	if ($chr_c eq $chr_n) {
		# Check if gap between current and next pic is smaller than their combined length
		if (($beg_n-$end_c) < ($len_c+$len_n) && ($beg_n-$end_c) < 75_000) {
			# Save pair of pics to be combined
			my $included = 0;
			foreach my $merge_bin (@pic_merge_bins) {
				if (grep {$pic_c eq $_} @{$merge_bin} or grep {$pic_n eq $_} @{$merge_bin}) {
					push(@{$merge_bin},$pic_c) unless grep {$pic_c eq $_} @{$merge_bin};
					push(@{$merge_bin},$pic_n) unless grep {$pic_n eq $_} @{$merge_bin};
					$included = 1;
				}
			}
			if (not $included) {
				push(@pic_merge_bins,[$pic_c,$pic_n]);
			}
		} else {
			push(@pic_merge_bins,[$pic_n]);
		}
	} else {
		push(@pic_merge_bins,[$pic_n]);
	}
}

## Merge clusters and their properties in merge bins
# Go through each pic merge bin
foreach my $merge_bin (@pic_merge_bins) {
	# Get ids of first and last pic
	my $pic_fst = $merge_bin->[0];
	my $pic_lst = $merge_bin->[-1];
	# Get merged pic id
	my $pid = $pic_fst;
	foreach my $i (1..(scalar(@{$merge_bin})-1)) {
		$pid .= "_$merge_bin->[$i]";
	}
	# Get merged pic coordinates
	my $chr = $pic_gtf_props->{$pic_fst}->[0];
	my $beg = $pic_gtf_props->{$pic_fst}->[1];
	my $end = $pic_gtf_props->{$pic_lst}->[2];
	# Get merged pic length
	my $len = $end-$beg+1;
	# Get merged pic length (w/o gaps)
	my $pic_size_mrg = 0;
	foreach my $pic (@{$merge_bin}) {
		$pic_size_mrg += $pic_tbl_props->{$pic}->[3];
	}
	# Get merged pic reads
	my $rds = 0;
	foreach my $pic (@{$merge_bin}) {
		$rds += $pic_gtf_props->{$pic}->[4];
	}
	# Get rpm
	my $rpm = int(($rds/$pir_total_rds*1000000)+0.5);
	# Get rpkm
	my $rpkm = int($rpm/$pic_size_mrg*1000);
	# Get merged pic repeat share
	my $bps_mask_mrg = 0;
	foreach my $pic (@{$merge_bin}) {
		$bps_mask_mrg += $pic_tbl_props->{$pic}->[12];
	}
	my $rep_rate = int(($bps_mask_mrg/$pic_size_mrg*100)+0.5);
	# Get merged pic 1U and 10A rates
	my $T01_bps = 0;
	my $A10_bps = 0;
	foreach my $pic (@{$merge_bin}) {
		my $T01_rate = $pic_tbl_props->{$pic}->[7];
		my $A10_rate = $pic_tbl_props->{$pic}->[8];
		my $pic_size = $pic_tbl_props->{$pic}->[3];
		$T01_bps += $pic_size*($T01_rate/100);
		$A10_bps += $pic_size*($A10_rate/100);
	}
	my $T01_rate = int(($T01_bps/$pic_size_mrg*100)+0.5);
	my $A10_rate = int(($A10_bps/$pic_size_mrg*100)+0.5);
	# Get merged pic strandedness
	my %str_reads = ();
	foreach my $pic (@{$merge_bin}) {
		my $str_rate = $pic_tbl_props->{$pic}->[10];
		my $str_mdir = $pic_tbl_props->{$pic}->[11];
		my $pic_trds = $pic_gtf_props->{$pic}->[4];
		if ($str_mdir eq 'plus') {
			$str_reads{'plus'} += $pic_trds*($str_rate/100);
			$str_reads{'minus'} += $pic_trds*((100-$str_rate)/100);
		} elsif ($str_mdir eq 'minus') {
			$str_reads{'minus'} += $pic_trds*($str_rate/100);
			$str_reads{'plus'} += $pic_trds*((100-$str_rate)/100);
		}
	}
	my @str_reads = sort {$b <=> $a} values %str_reads;
	my @str_mdirs = sort {$str_reads{$b} <=> $str_reads{$a}} keys %str_reads;
	my $str_rate = int(($str_reads[0]/$rds*100)+0.5);
	my $str_mdir = $str_mdirs[0];
	# Get merged ping pong z-score
	my $pic_reads_file = "$pT_dir/$pid.fasta";
	if (scalar(@{$merge_bin}) > 1) {
		my $out = open_outfile($pic_reads_file);
		foreach my $pic (@{$merge_bin}) {
			my $infile = "$pT_dir/$pic.fasta";
			my $file_data = get_file_data_string($infile);
			print($out "$file_data\n");
		}
		close($out);
	}
	my $pp_z = get_ping_pong_z($pic_reads_file,$pp_script);
	$pp_z = (int($pp_z*100))/100;
	# Save results for each (merged) pic
	@{$fin_pics{$pid}} = ($chr,$beg,$end,$len,$rds,$rpm,$rpkm,$str_rate,$str_mdir,$rep_rate,$T01_rate,$A10_rate,$pp_z,$pid);
}

## Create final output
# Open output file
my $outfile = "$pT_dir/fin_pics.bed";
my $out = open_outfile($outfile);
# Output processed pics
print($out "chr\tbeg\tend\tlen\treads\trpm\trpkm\tstr_reads[%]\tmain_str\treps[%]\t1U[%]\t10A[%]\tpp_zscore\tpic_id\n");
foreach my $pid (sort {$fin_pics{$a}[1] cmp $fin_pics{$b}[1] || $fin_pics{$a}[2] <=> $fin_pics{$b}[2]} keys %fin_pics) {
	foreach my $f (@{$fin_pics{$pid}}) {
		print($out "$f\t");
	}
	print($out "\n");
}

print("done.\n");

exit;

################################# subroutines #################################

sub get_pic_gtf_properties {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
        # Get line data
        my @d = split(/\t/,$line);
        # Get relevant data fields
        my $chr = $d[0];
		my $beg = $d[3];
		my $end = $d[4];
		my $inf = $d[8];
		my($pid) = ($inf =~ /piRNA cluster no:\s([^;]+);/);
		my($rds) = ($inf =~ /mapped sequence reads:\s([^;]+);/);
		# Save pic properties
		@{$data_fields{$pid}} = ($chr,$beg,$end,$inf,$rds);
    }
	# Return data fields
	return \%data_fields;
}

sub get_pic_tbl_properties {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
		# Skip non cluster info line
		unless ($line =~ /Cluster/) { next; }
        # Get line data
        my @d = split(/\t/,$line);
        # Get relevant data fields
		foreach my $d (@d) {
			$d =~ s/.*://;
			$d =~ s/%//;
		}
		my($pid) = ($d[0] =~ /Cluster\s(\d+)/);
		# Save pic properties
		@{$data_fields{$pid}} = @d;
    }
	# Return data fields
	return \%data_fields;
}

sub get_total_pirna_reads {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my $total_pir_reads = 0;
	# Go through file data
	foreach my $line (@in_data) {
		# Skip non cluster info line
		if ($line =~ /Sequence reads that can be assigned to clusters: /) {
			$line =~ s/Sequence reads that can be assigned to clusters: //;
			my($pic_reads_bps) = ($line =~ /(\d+)\s/);
			my($pic_reads_per) = ($line =~ /\s\((.*)%\)/);
			$total_pir_reads = int(100*$pic_reads_bps/$pic_reads_per);
		}
    }
	# Return data fields
	return $total_pir_reads;
}

sub get_ping_pong_z {
	# Take name of pic reads file
	my($infile,$pp_script) = @_;
	# Get file data
	my @pic_reads_data = get_file_data_array($infile);
	# Name eland output file
	my $eland_file = "$infile.eland";
	my $eland = open_outfile($eland_file);
	# Convert pic reads file to eland format
	my @d = ();
	my $seq = '';
	foreach my $line (@pic_reads_data) {
		if ($line =~ /^>/) {
			@d = split(/\t/,$line);
			foreach my $d (@d) {
				$d =~ s/.*://;
			}
		} elsif ($line !~ /^$/) {
			$seq = $line;
			print($eland "$d[0]\t$d[1]\t$seq\t$d[4]\t$seq\t0\t$d[5]\n");
		}
	}
	# Name ping pong outfile
	my $pp_file = "$infile.pp.txt";
	# Call TBr2_pingpong.pl
	system("perl $pp_script -i $eland_file -o $pp_file >>.pplog 2>&1");
	#unlink($eland_file);
	unlink(".pplog");
	# Get ping pong outfile data
	my @pp_data = get_file_data_array($pp_file);
	# Extract z-score
	my $z_score = 0;
	foreach my $line (@pp_data) {
		if ($line =~ /Ping-Pong Z-Score:/) {
			($z_score) = ($line =~ /Ping-Pong Z-Score:\s(.*)/);
		}
	}
	if ($z_score eq 'infinite') { $z_score = 999 }
	return $z_score;
}

# Open input file
# Usage: my $in = open_infile($infile);
sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

# Open output file
# Usage: my $out = open_outfile($outfile);
sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

# Extract file data and save in array
# Usage: my @filedata = get_file_data_array($file);
sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}

# Extract file data and save in string
# Usage: my $filedata = get_file_data_string($file);
sub get_file_data_string {
    # Take input file name
	my($file,$ref_opt) = @_;
    $ref_opt = 0 unless $ref_opt;
	local $/=undef;
    # Open input file
	my $in = open_infile($file);
    # Save file data as string
	my $file_data = <$in>;
    # Return file data string/ref
    if ($ref_opt) {
    	return \$file_data;
    } else {
    	return $file_data;
    }
}