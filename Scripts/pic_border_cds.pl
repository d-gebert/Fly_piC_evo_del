#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $picfile_base = 'fin_pics.bed';
my $topfile_base = 'fin_pics.top.bed';
my $gnmfile_base = '_filtered_50k.fasta';
my $cdsfile_base = '_CDS.fasta';
my $bin_size = 500;
# Global variables
my @pic_locs = ();
my @pic_locs_top = ();
my @cds_locs = ();
my @chr_size = ();
# Options
$|=1; #Autoflush

# Species identifiers
my @species = ('Dmel','Dsec','Dsim','Dyak','Dana','Dpse','Dper','Dwil','Dmoj','Dvir');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <fin_pics_directory> <genomes_directory>\n";
unless ($ARGV[0]&&$ARGV[1]) {
   die("\nUsage: $USAGE\n");
}
# Input paths
my $fin_pics_dir = $ARGV[0];
my $genomes_dir = $ARGV[1];

# Get all pic loci
foreach my $sp (@species) {
	# Get pic file name
	my $pic_file = "$fin_pics_dir/$sp.$picfile_base";
	# Get pic locs
	$pic_locs[$spec_id{$sp}] = get_tab_fields($pic_file,1);
}

# Get all top pic loci
foreach my $sp (@species) {
	# Get pic file name
	my $pic_file = "$fin_pics_dir/$sp.$topfile_base";
	# Get pic locs
	$pic_locs_top[$spec_id{$sp}] = get_tab_fields($pic_file,1);
}

# Get pics that are not in top list
foreach my $i (0..$#species) {
	foreach my $pic (keys %{$pic_locs[$i]}) {
		# Get pic id
		my $pic_pid = $pic_locs[$i]->{$pic}->[14];
		# Check if pic is included in top pics
		my $in_top_list = 0;
		foreach my $top (keys %{$pic_locs_top[$i]}) {
			# Get pic id
			my $top_pid = $pic_locs_top[$i]->{$top}->[14];
			# Check if same pic id
			if ($pic_pid eq $top_pid) {
				$in_top_list = 1;
			}
		}
	}
}

# Get cds loci for each species
foreach my $sp (@species) {
	# Get cds, genome and blast file names
	my $cds_file = "$genomes_dir/${sp}$cdsfile_base";
	my $gnm_file = "$genomes_dir/${sp}$gnmfile_base";
	my $bln_file = "$genomes_dir/${sp}_cds.bln";
	# Get repeat loci
	$cds_locs[$spec_id{$sp}] = get_cds_blast_hits($cds_file,$gnm_file,$bln_file);
}

# Get genome size for each species
foreach my $sp (@species) {
	# Get genome file name
	my $gnm_file = "$genomes_dir/${sp}$gnmfile_base";
	# Get genome size
	$chr_size[$spec_id{$sp}] = get_chromosome_sizes($gnm_file);
}

my $pic_cds_distr = pic_flank_exon_distribution(\@pic_locs,$bin_size);
my $top_cds_distr = pic_flank_exon_distribution(\@pic_locs_top,$bin_size);

my $out = open_outfile('pic_flank_cds_distr.txt');
print($out "pos\tall_pics\ttop_pics\n");
foreach my $pos (sort {$a <=> $b} keys %{$pic_cds_distr}) {
	print($out "$pos\t$pic_cds_distr->{$pos}\t$top_cds_distr->{$pos}\n");
}

exit;

################################# subroutines #################################

sub get_tab_fields {
	# Take name of tab file
	my($infile,$skip_header) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	if ($skip_header) { shift(@in_data); }
	# Global tree hashes for each species
	my %data_fields = ();
	# Set 0 as start index
	my $id_i = 0;
	# Go through file data
	foreach my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        # Save data fields
        @{$data_fields{$id_i}} = @d;
		$id_i++;
    }
	# Return data fields
	return \%data_fields;
}

sub get_cds_blast_hits {
	# Take name of cds and genome file
	my($cdsfile,$gnmfile,$blnfile) = @_;
	# Blast output format
	my $outfmt = '7 sseqid sstart send qseqid score sstrand qstart qend length slen qlen pident evalue bitscore qcovs';
	# Blast genome
	unless (-e $blnfile) {
		system("blastn -query $cdsfile -subject $gnmfile -out $blnfile -outfmt \'$outfmt\'");
	}
	# Open blast file and bed output file
	my $bln = open_infile($blnfile);
	my $bed = open_outfile("$blnfile.bed");
	# Go through each blast file line
	while (my $line = <$bln>) {
		chomp($line);
		# Skip comments
		if ($line =~ /^#/) { next; }
		# Get line data
		my @d = split(/\t/,$line);
		if ($d[11] < 99) { next; }
		# Adjust fields
		$d[4] = '.';
		$d[5] = $d[5] eq 'plus' ? '+' : '-';
		print($bed "$d[0]\t$d[1]\t$d[2]\t$d[3]\t$d[4]\t$d[5]\n");
	}
	close($bed);
	system("sort -k1,1 -k2,2n <$blnfile.bed >$blnfile.sorted.bed");
	my $cds_locs = get_bed_data("$blnfile.sorted.bed");
	#unlink("$blnfile.bed");
	#unlink("$blnfile.sorted.bed");
	return $cds_locs;
}

sub get_bed_data {
	# Take repeatmasker file name
	my($bed_file) = @_;
	# Storage variable
	my %bed_data = ();
	# Get file data
	my @bed_data = get_file_data_array($bed_file);
	# Parse bed file
	foreach my $line (@bed_data) {
		# Get line data
		my @d = split(/\t/,$line);
		my $chr = $d[0];
		push(@{$bed_data{$chr}},\@d);
	}
	return \%bed_data;
}

sub get_chromosome_sizes {
	# Take fasta file name
	my($fasta_file) = @_;
	# Var init
	my %chr_sizes = ();
	my $chr = '';
	# Open input file
	my $in = open_infile($fasta_file);
	# Parse fasta file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Sequence
		if ($line =~ /^>/) {
			($chr) = ($line =~ /^>(\S+)/);
		}
		elsif ($line !~ /^>/) {
			$chr_sizes{$chr} += length($line);
		}
	}
	return \%chr_sizes;
}

sub pic_flank_exon_distribution {
	# Take name of tab file
	my($pic_locs,$bin_size) = @_;
	$bin_size = 1 unless $bin_size;
	# Parameters
	my $flank_len = 20_000;
	my $edge_len  = 5_000;
	# Get cds content for each relative position
	my %pos_rel = ();
	my $n_pics  = 0;
	my $n_no_flanks = 0;
	# Go through each species
	foreach my $i (0..$#species) {
		foreach my $pic (keys %{$pic_locs->[$i]}) {
			# Get pic coordinates
			my $pic_chr = $pic_locs->[$i]->{$pic}->[0];
			my $pic_beg = $pic_locs->[$i]->{$pic}->[1];
			my $pic_end = $pic_locs->[$i]->{$pic}->[2];
			# Get upstream flank coordinates
			my $ups_beg = $pic_beg-$flank_len;
			my $ups_end = $pic_beg-1;
			# Get downstream flank coordinates
			my $dos_beg = $pic_end+1;
			my $dos_end = $pic_end+$flank_len;
			# Get pic edge coordinates
			my $ups_edg = $pic_beg+$edge_len-1;
			my $dos_edg = $pic_end-$edge_len+1;
			# Skip if coordinates are out of chromosome bounds on both sides
			if ($ups_beg < 1 && $dos_end > $chr_size[$i]->{$pic_chr}) { $n_no_flanks++; next; }
			$n_pics++;
			# Check if only one flank is complete
			my $n_flanks = 2;
			my $up_flank = 1;
			my $do_flank = 1;
			if ($ups_beg < 1) {
				$up_flank = 0;
				$n_flanks = 1;
			}
			elsif ($dos_end > $chr_size[$i]->{$pic_chr}) {
				$do_flank = 0;
				$n_flanks = 1;
			}
			# Initialize position translator
			my %pos_trl = ();
			my %pos_count = ();
			# If upstream flank present
			if ($up_flank) {
				# Initialize relative start position
				my $rel_pos = 1;
				# Translate absolute to relative positions (upstream flank+edge)
				foreach my $pos ($ups_beg..$ups_edg) {
					$pos_trl{$pos} = $rel_pos;
					$pos_count{$pos}++;
					$rel_pos++;
				}
			}
			# If downstream flank present
			if ($do_flank) {
				# Initialize relative start position
				my $rel_pos = 1;
				# Translate absolute to relative positions (downstream flank+edge)
				foreach my $pos (reverse($dos_edg..$dos_end)) {
					$pos_trl{$pos} = $rel_pos;
					$pos_count{$pos}++;
					$rel_pos++;
				}
			}
			# Get cds content of pic loc and flanking region
			my %rep_poss = ();
			foreach my $rep (@{$cds_locs[$i]->{$pic_chr}}) {
				# Get cds coordinates
				my $rep_beg = $rep->[1];
				my $rep_end = $rep->[2];
				# Check for overlap with upstream flank
				if ($ups_beg <= $rep_end && $ups_edg >= $rep_beg) {
					# Save each position
					foreach my $pos ($rep_beg..$rep_end) {
						$rep_poss{$pos} = 1;
					}
				}
				# Check for overlap with downstream flank
				if ($dos_edg <= $rep_end && $dos_end >= $rep_beg) {
					# Save each position
					foreach my $pos ($rep_beg..$rep_end) {
						$rep_poss{$pos} = 1;
					}
				}
			}
			# Save each position
			foreach my $pos (keys %rep_poss) {
				$pos_rel{$pos_trl{$pos}} += 1/$n_flanks if $pos_trl{$pos};
			}
		}
	}
	# Adjust resolution
	my %bin_rel = ();
	if ($bin_size > 1) {
		# Go through each bin
		for (my $bin = 0; $bin < ($flank_len+$edge_len); $bin += $bin_size) {
			# Count all te positions for current bin
			my $bin_tes = 0;
			foreach my $pos ($bin..($bin+$bin_size)) {
				next unless $pos;
				next unless $pos_rel{$pos};
				$bin_tes += $pos_rel{$pos}/$bin_size;
			}
			$bin_rel{$bin+$bin_size} = $bin_tes;
		}
		%pos_rel = %bin_rel;
	}
	print("$n_pics \($n_no_flanks\)\n");
	# Get fractions instead of absolute cds position counts
	foreach my $pos (keys %pos_rel) {
		$pos_rel{$pos} = 0 unless $pos_rel{$pos};
		$pos_rel{$pos} = $pos_rel{$pos}/$n_pics;
	}
	return \%pos_rel;
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
