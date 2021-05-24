#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use lib '/home/dgebert/Dropbox/Perlmodules';
use FileIO;
use FastaIO;
use BioStat;

# Global constants
my $picfile_base = 'fin_pics.bed';
my $topfile_base = 'fin_pics.top.bed';
my $gnmfile_base = '_filtered_50k.fasta';
my $repfile_base = '_filtered_50k.fasta.out';
my $bin_size = 500;
# Global variables
my @pic_locs = ();
my @pic_locs_top = ();
my @pic_locs_not = ();
my @rep_locs = ();
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

# Get all pic loci for subsequent filtering top pics
foreach my $sp (@species) {
	# Get pic file name
	my $pic_file = "$fin_pics_dir/$sp.$picfile_base";
	# Get pic locs
	$pic_locs_not[$spec_id{$sp}] = get_tab_fields($pic_file,1);
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
		# Delete pic from not list if it is a top pic
		if ($in_top_list) {
			delete($pic_locs_not[$i]->{$pic});
		}
	}
}

# Get repeat loci for each species
foreach my $sp (@species) {
	# Get pic file name
	my $rep_file = "$genomes_dir/${sp}$repfile_base";
	# Get repeat loci
	$rep_locs[$spec_id{$sp}] = get_repeatmask_data($rep_file);
}

# Get genome size for each species
foreach my $sp (@species) {
	# Get pic file name
	my $gnm_file = "$genomes_dir/${sp}$gnmfile_base";
	# Get genome size
	$chr_size[$spec_id{$sp}] = get_chromosome_sizes($gnm_file);
}

my($pic_te_distr,$pic_te_distr_arr) = pic_flank_repeat_distribution(\@pic_locs,$bin_size);
my($top_te_distr,$top_te_distr_arr) = pic_flank_repeat_distribution(\@pic_locs_top,$bin_size);
my($not_te_distr,$not_te_distr_arr) = pic_flank_repeat_distribution(\@pic_locs_not,$bin_size);

my $out1 = FileIO::open_outfile('pic_flank_te_distr.txt');
print($out1 "pos\tall_pics\ttop_pics\tnot_pics\n");
foreach my $pos (sort {$a <=> $b} keys %{$pic_te_distr}) {
	print($out1 "$pos\t$pic_te_distr->{$pos}\t$top_te_distr->{$pos}\t$not_te_distr->{$pos}\n");
}

my($pic_te_div_distr,$pic_te_div_distr_arr) = pic_flank_rep_div_distribution(\@pic_locs,$bin_size);
my($top_te_div_distr,$top_te_div_distr_arr) = pic_flank_rep_div_distribution(\@pic_locs_top,$bin_size);
my($not_te_div_distr,$not_te_div_distr_arr) = pic_flank_rep_div_distribution(\@pic_locs_not,$bin_size);

my $out2 = FileIO::open_outfile('pic_flank_te_div_distr.txt');
print($out2 "pos\tall_pics\ttop_pics\tnot_pics\n");
foreach my $pos (sort {$a <=> $b} keys %{$pic_te_div_distr}) {
	print($out2 "$pos\t$pic_te_div_distr->{$pos}\t$top_te_div_distr->{$pos}\t$not_te_div_distr->{$pos}\n");
}

my $out3 = FileIO::open_outfile('pic_flank_te_div_distr_all.txt');
print($out3 "pos\tmean_div\tsder_div\n");
foreach my $pos (sort {$a <=> $b} keys %{$pic_te_div_distr_arr}) {
	my $mean_div = BioStat::get_mean($pic_te_div_distr_arr->{$pos});
	my $sder_div = BioStat::get_standard_error($pic_te_div_distr_arr->{$pos});
	print($out3 "$pos\t$mean_div\t$sder_div\n");
}

# Get pic te divergence for each species
my @pic_te_div_distr_sp = ();
my @pic_te_div_distr_arr_sp = ();
# Go through each species
foreach my $i (0..$#species) {
	($pic_te_div_distr_sp[$i],$pic_te_div_distr_arr_sp[$i]) = pic_flank_rep_div_distribution(\@pic_locs,$bin_size,1,$i);
}

# Prepare output
my $out4 = FileIO::open_outfile('pic_flank_te_div_distr_spec.txt');
print($out4 "pos");
foreach my $i (0..$#species) {
	print($out4 "\t$species[$i]");
}
print($out4 "\n");
# Go through each bin
foreach my $pos (sort {$a <=> $b} keys %{$pic_te_div_distr}) {
	print($out4 "$pos");
	# Go through each species
	foreach my $i (0..$#species) {
		print($out4 "\t$pic_te_div_distr_sp[$i]->{$pos}");
	}
	print($out4 "\n");
}

exit;

################################# subroutines #################################

sub get_tab_fields {
	# Take name of tab file
	my($infile,$skip_header) = @_;
	# Get file data
	my @in_data = FileIO::get_file_data_array($infile);
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

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file) = @_;
	# Storage variable
	my %rep_data = ();
	# Get file data
	my @repeatmask_data = FileIO::get_file_data_array($repeatmask_file);
	# Parse repeatmasker file
	foreach my $line (@repeatmask_data) {
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			$line =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line);
			my $chr = $d[4];
			my $cla = $d[10];
			if ($cla =~ /Satellite/ || $cla =~ /Simple/ || $cla =~ /Tandem/ || $cla =~ /Low/ || $cla =~ /Unknown/) { next; }
			push(@{$rep_data{$chr}},\@d);
		}
	}
	return \%rep_data;
}

sub get_chromosome_sizes {
	# Take fasta file name
	my($fasta_file) = @_;
	# Var init
	my %chr_sizes = ();
	my $chr = '';
	# Open input file
	my $in = FileIO::open_infile($fasta_file);
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

sub pic_flank_repeat_distribution {
	# Take name of tab file
	my($pic_locs,$bin_size) = @_;
	$bin_size = 1 unless $bin_size;
	# Parameters
	my $flank_len = 20_000;
	my $edge_len  = 5_000;
	# Get TE content for each relative position
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
			#if ($n_flanks == 1) { next; }
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
			# Get TE content of pic loc and flanking region
			my %rep_poss = ();
			foreach my $rep (@{$rep_locs[$i]->{$pic_chr}}) {
				# Get TE coordinates
				my $rep_beg = $rep->[5];
				my $rep_end = $rep->[6];
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
				$pos_rel{$pos_trl{$pos}} += $pos_count{$pos}/$n_flanks if $pos_trl{$pos};
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
				$bin_tes += $pos_rel{$pos}/$bin_size;
			}
			$bin_rel{$bin+$bin_size} = $bin_tes;
		}
		%pos_rel = %bin_rel;
	}
	print("$n_pics \($n_no_flanks\)\n");
	# Get fractions instead of absolute TE position counts
	foreach my $pos (keys %pos_rel) {
		$pos_rel{$pos} = 0 unless $pos_rel{$pos};
		$pos_rel{$pos} = $pos_rel{$pos}/$n_pics;
	}
	return \%pos_rel;
}

sub pic_flank_rep_div_distribution {
	# Take name of tab file
	my($pic_locs,$bin_size,$single,$spec_id) = @_;
	# Parameters
	my $flank_len = 20_000;
	my $edge_len  = 5_000;
	# Get TE content for each relative position
	my %pos_rel = ();
	my %pos_rel_div = ();
	my %pos_rel_cnt = ();
	my $n_pics  = 0;
	my $n_no_flanks = 0;
	my %pos_rel_div_arr = ();
	# Go through each species
	foreach my $i (0..$#species) {
		if ($single) {
			unless ($i == $spec_id) { next; }
		}
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
			#if ($n_flanks == 1) { next; }
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
			# Get TE content of pic loc and flanking region
			my %rep_poss = ();
			my %pos_divs = ();
			my %pos_reps = ();
			foreach my $rep (@{$rep_locs[$i]->{$pic_chr}}) {
				# Get TE coordinates
				my $rep_div = $rep->[1];
				my $rep_beg = $rep->[5];
				my $rep_end = $rep->[6];
				# Check for overlap with upstream flank
				if ($ups_beg <= $rep_end && $ups_edg >= $rep_beg) {
					# Save each position
					foreach my $pos ($rep_beg..$rep_end) {
						$rep_poss{$pos} = 1;
						$pos_divs{$pos} = $rep_div;
						$pos_reps{$pos}++;
					}
				}
				# Check for overlap with downstream flank
				if ($dos_edg <= $rep_end && $dos_end >= $rep_beg) {
					# Save each position
					foreach my $pos ($rep_beg..$rep_end) {
						$rep_poss{$pos} = 1;
						$pos_divs{$pos} = $rep_div;
						$pos_reps{$pos}++;
					}
				}
			}
			# Save each position
			foreach my $pos (keys %rep_poss) {
				$pos_rel{$pos_trl{$pos}} += $pos_count{$pos}/$n_flanks if $pos_trl{$pos};
				$pos_rel_div{$pos_trl{$pos}} += $pos_divs{$pos} if $pos_trl{$pos};
				$pos_rel_cnt{$pos_trl{$pos}} += $pos_reps{$pos} if $pos_trl{$pos};
				# Array
				push(@{$pos_rel_div_arr{$pos_trl{$pos}}},($pos_divs{$pos}/$pos_reps{$pos})) if $pos_trl{$pos};
			}
		}
	}
	foreach my $pos (keys %pos_rel) {
		$pos_rel_div{$pos} = $pos_rel_div{$pos}/$pos_rel_cnt{$pos};
	}
	# Adjust resolution
	my %bin_rel = ();
	my %bin_rel_div_arr = ();
	if ($bin_size > 1) {
		# Go through each bin
		for (my $bin = 0; $bin < ($flank_len+$edge_len); $bin += $bin_size) {
			# Count all te positions for current bin
			my $bin_tes = 0;
			foreach my $pos ($bin..($bin+$bin_size)) {
				next unless $pos;
				$bin_tes += $pos_rel_div{$pos}/$bin_size;
				# Array
				foreach my $te_i (0..scalar(@{$pos_rel_div_arr{$pos}})-1) {
					$bin_rel_div_arr{$bin+$bin_size}[$te_i] += $pos_rel_div_arr{$pos}[$te_i]/$bin_size;
				}
			}
			$bin_rel{$bin+$bin_size} = $bin_tes;
		}
		%pos_rel_div = %bin_rel;
		%pos_rel_div_arr = %bin_rel_div_arr;
	}
	print("$n_pics \($n_no_flanks\)\n");
	# Get fractions instead of absolute TE position counts
	foreach my $pos (keys %pos_rel_div) {
		$pos_rel_div{$pos} = 0 unless $pos_rel_div{$pos};
		$pos_rel_div{$pos} = $pos_rel_div{$pos};
	}
	return \(%pos_rel_div, %pos_rel_div_arr);
}
