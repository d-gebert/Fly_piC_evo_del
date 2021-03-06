#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $outfile = '';
my $iterations_max = 10;
# Global variables
my @files = ();
my @homlocs = ();
my @matrix = ();
my @new_locs = ();
my $br = 0;
my @branches = ();
my %deleted_br = ();
my $n_pics_all_prev = 0;
my $n_pics_uni_prev = 0;
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
# Usage message
my $USAGE = "perl $0 [-a (output all syntenic)]\n";
if ($ARGV[0] && $ARGV[0] ne '-a') {
	die("\nUsage: $USAGE\n");
}
my $all_opt = $ARGV[0] && $ARGV[0] eq '-a' ? '-a' : '';

if ($all_opt) {
	$outfile = 'homlocs_matrix_synt.txt';
} else {
	$outfile = 'homlocs_matrix_picp.txt';
}

# Go through each species
foreach my $i (0..$#species) {
	# Get species id
	my $spec_x = $species[$i];
	# Go through each species
	foreach my $j (0..$#species) {
		# Get species id
		my $spec_y = $species[$j];
		# Skip if same species
		next; if $spec_x eq $spec_y;
		#next; unless $spec_x eq 'Dmel' && $spec_y eq 'Dsec';
		# If no hompics file exists, run pipeline to get hompics
		unless (-e "${spec_x}_vs_${spec_y}.hompics.txt") {
			print("perl get_syntenic_loci.pl ${spec_x}.fin_pics.top.bed ${spec_x}_filtered_50k.fasta ${spec_x}_chromosome.fasta ${spec_x}_genes.gff dmel_orthologs_in_drosophila_species_fb_2019_04.tsv ${spec_y}_filtered_50k.fasta ${spec_y}_chromosome.fasta ${spec_y}_CDS.fasta ${spec_x}_filtered_50k.fasta.out ${spec_y}_filtered_50k.fasta.out ${spec_x}.filtered.23-29.f2.uni.rpm.bed ${spec_y}.filtered.23-29.f2.uni.rpm.bed");
			system("perl get_syntenic_loci.pl ${spec_x}.fin_pics.top.bed ${spec_x}_filtered_50k.fasta ${spec_x}_chromosome.fasta ${spec_x}_genes.gff dmel_orthologs_in_drosophila_species_fb_2019_04.tsv ${spec_y}_filtered_50k.fasta ${spec_y}_chromosome.fasta ${spec_y}_CDS.fasta ${spec_x}_filtered_50k.fasta.out ${spec_y}_filtered_50k.fasta.out ${spec_x}.filtered.23-29.f2.uni.rpm.bed ${spec_y}.filtered.23-29.f2.uni.rpm.bed");
		}
		push(@files,"${spec_x}_vs_${spec_y}.hompics.txt");
	}
}

exit;

## Extract homologous loci from homloc files
# Go through each homloc file
foreach my $file (@files) {
	# Get species and genome IDs
	my($spec_x) = ($file =~ /^(\w+)\_vs/);
	my($spec_y) = ($file =~ /vs\_(\w+)\./);
	my $x = $spec_id{$spec_x};
	my $y = $spec_id{$spec_y};
	# Get homologous loci
	if ($all_opt) {
		$homlocs[$x][$y] = get_homlocs_syntenic($file);
	} else {
		$homlocs[$x][$y] = get_homlocs_pc_present($file);
	}
}

## Build initial piC matrix
# Go through each x species
foreach my $x (0..$#species) {
	# Go through each y species
	foreach my $y (0..$#species) {
		# Go through each x species loc
		foreach my $loc_x (sort keys %{$homlocs[$x][$y]}) {
			# Increment branch id number
			$br++;
			push(@branches, $br);
			# Save x species loc in current matrix branch
			push(@{$matrix[$x]{$br}}, $loc_x);
			# Save y species homloc in current matrix branch
			my $loc_y = $homlocs[$x][$y]->{$loc_x};
			push(@{$matrix[$y]{$br}}, $loc_y) if $loc_y;
		}
	}
}

## Collapse matrix
# Get overall and non-redundant piC counts
my($n_pics_all, $n_pics_uni) = count_matrix_pics(\@matrix, \@branches);
# Repeat until overall and non-redundant counts are equal
while ($n_pics_all != $n_pics_uni) {
	# Go through each matrix branch (p)
	foreach my $br_p (@branches) {
		# Skip if branch is deleted
		if ($deleted_br{$br_p}) { next; }
		# Go through each following matrix branch (q)
		foreach my $br_q (($br_p+1)..scalar(@branches)) {
			# Skip if branch is deleted
			if ($deleted_br{$br_q}) { next; }
			# Check if loci in branches p and q overlap, set default: false
			my $loci_overlap = 0;
			my $p_len = 0;
			my $q_len = 0;
			# Go through each species
			foreach my $i (0..$#species) {
				# Get each loc of branch p
				foreach my $loc_p (@{$matrix[$i]{$br_p}}) {
					unless ($loc_p) { next; }
					if ($loc_p eq 'NA') { next; }
					# Get loc p coordinates
					my($chr_p,$beg_p,$end_p) = split(/-|:|\|/,$loc_p);
					# Get each loc of branch q
					foreach my $loc_q (@{$matrix[$i]{$br_q}}) {
						# Get loc q coordinates
						my($chr_q,$beg_q,$end_q) = split(/-|:|\|/,$loc_q);
						# Check if loci overlap
						if ($chr_p eq $chr_q && $beg_q <= $end_p && $end_q >= $beg_p) {
							# Save positions
							my %p_positions = ();
							foreach my $pos ($beg_p..$end_p) {
								$p_positions{$chr_p}{$pos} = 1;
							}
							# Compare positions
							foreach my $pos ($beg_q..$end_q) {
								$loci_overlap++ if $p_positions{$chr_q}{$pos};
							}
							$p_len = $end_p-$beg_p+1;
							$q_len = $end_q-$beg_q+1;
						}
					}
				}
			}
			# If loci overlap add loci in branch q to branch p
			if ($loci_overlap) {
				if ($loci_overlap/$p_len > 0.25 && $loci_overlap/$q_len > 0.25) {
					# Go through each species
					foreach my $i (0..$#species) {
						# Get each loc of branch q
						foreach my $loc_q (@{$matrix[$i]{$br_q}}) {
							# Get loc q coordinates
							my($chr_q,$beg_q,$end_q) = split(/-|:|\|/,$loc_q);
							# Check if q branch loc is included in p branch, set default: false
							my $loc_included = 0;
							# Get each loc of branch p
							foreach my $loc_p (@{$matrix[$i]{$br_p}}) {
								if ($loc_p eq 'NA') { next; }
								# Get loc p coordinates
								my($chr_p,$beg_p,$end_p) = split(/-|:|\|/,$loc_p);
								# Check if loci overlap
								if ($chr_p eq $chr_q && $beg_q <= $end_p && $end_q >= $beg_p) {
									# Locus included: true
									$loc_included = 1;
									# Locus might need update if loci not identical
									my $loc_needs_update = 0;
									# Loci overlap but are not identical
									if ($loc_p ne $loc_q) {
										# Get new start and end coordinates
										my $beg_new = $beg_p;
										my $end_new = $end_p;
										if ($beg_q < $beg_p) {
											$beg_new = $beg_q;
											# Locus needs update: true
											#$loc_needs_update = 1;
										}
										if ($end_q > $end_p) {
											$end_new = $end_q;
											# Locus needs update: true
											#$loc_needs_update = 1;
										}
										if ($loc_needs_update) {
											# Define new locus
											my $loc_new = $chr_p.':'.$beg_new.'-'.$end_new;
											# Get index of old locus to replace
											my $n_locs_p = scalar(@{$matrix[$i]{$br_p}})-1;
											my($i_old) = grep {$matrix[$i]{$br_p}[$_] eq $loc_p} 0..$n_locs_p;
											# Update locus
											$new_locs[$i]{$br_p}{$i_old} = $loc_new;
										}
									}
								}
							}
							# If not already included, add branch q locus to branch p
							unless ($loc_included) {
								push(@{$matrix[$i]{$br_p}}, $loc_q) if $matrix[$i]{$br_p};
							}
						}
						# Delete branch q
						@{$matrix[$i]{$br_q}} = ();
						$deleted_br{$br_q} = 1;
					}
				}
			}

		}
		foreach my $i (0..$#species) {
			foreach my $i_old (keys %{$new_locs[$i]{$br_p}}) {
				$matrix[$i]{$br_p}[$i_old] = $new_locs[$i]{$br_p}{$i_old};
			}
		}
	}
	# Get new overall and non-redundant piC counts
	$n_pics_all_prev = $n_pics_all;
	$n_pics_uni_prev = $n_pics_uni;
	($n_pics_all, $n_pics_uni) = count_matrix_pics(\@matrix, \@branches);
	# Abort if overall and non-redundant piC counts stay the same
	print("$n_pics_all\t$n_pics_uni\n");
	if ($n_pics_all == $n_pics_all_prev && $n_pics_uni == $n_pics_uni_prev) { last }
}

## Collapse overlapping loci within a branch and species
# Go through each matrix branch
foreach my $br_p (@branches) {
	# Skip if branch is deleted
	if ($deleted_br{$br_p}) { next; }
	# Go through each species
	foreach my $i (0..$#species) {
		# Skip if branch is not defined for species
		unless ($matrix[$i]{$br_p}) { next; }
		# Get locs
		my @update_locs = @{$matrix[$i]{$br_p}};
		# Compare all locs
		for (my $a = 0; $a < (scalar(@{$matrix[$i]{$br_p}})-1); $a++) {
			for (my $b = ($a+1); $b <= (scalar(@{$matrix[$i]{$br_p}})-1); $b++) {
				# Skip obsolete locs
				unless ($update_locs[$a] && $update_locs[$b]) { next; }
				if ($update_locs[$a] eq 'NA' && $update_locs[$b] eq 'NA') { next; }
				# Get each loc p coordinates (a and b)
				my($chr_a,$beg_a,$end_a) = split(/-|:|\|/,$update_locs[$a]);
				my($chr_b,$beg_b,$end_b) = split(/-|:|\|/,$update_locs[$b]);
				# Check overlap
				if ($chr_a eq $chr_b && $beg_b <= $end_a && $end_b >= $beg_a) {
					# Get new start and end coordinates
					my $beg_new = $beg_a;
					my $end_new = $end_a;
					if ($beg_b < $beg_a) { $beg_new = $beg_b }
					if ($end_b > $end_a) { $end_new = $end_b }
					# Update loc (a) and delete loc (b)
					$update_locs[$a] = $chr_a.' '.$beg_new.'-'.$end_new;
					splice(@update_locs, $b, 1);
				}
			}
		}
		@{$matrix[$i]{$br_p}} = @update_locs;
	}
}
($n_pics_all, $n_pics_uni) = count_matrix_pics(\@matrix, \@branches);
print("$n_pics_all\t$n_pics_uni\n");

## Output piC matrix
# Sort matrix branches by number of homologs
my @sorted_branches = sort_matrix_by_homologs(\@matrix, \@branches);
# Create plain text representation of piC matrix
draw_global_matrix_txt(\@matrix, \@sorted_branches, $outfile);
# Create html heatmap of piC matrix
draw_global_matrix_html_heatmap(\@matrix, \@sorted_branches, 'homlocs_matrix_hm.html');

exit;

################################# subroutines #################################

sub get_homlocs_syntenic {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = get_file_data_array($homlocs_file);
	# Initialize variables
	my %homlocs = ();
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^$/) {
			if ($line =~ /^pic/) { next; }
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[0];
			my $loc_x = $d[1];
			my $ups_x = $d[2];
			my $dos_x = $d[3];
			my $ups_y = $d[4];
			my $dos_y = $d[5];
			my $loc_y = $d[6];
			my $rpm_y = $d[7];
			my $rpk_y = $d[8];
			my $len_y = $d[9];
			my $lnp_y = $d[10];
			my $rep_y = $d[11];
			my $gen_y = $d[12];
			# Break mid-pic bool
			my $picbr = $d[15];
			# Link homologous loci
			if ($loc_y !~ /NA/) {
				$homlocs{$loc_x} = $loc_y;
				next;
			} else {
				$homlocs{$loc_x} = '';
				$homlocs{$loc_x} = 'NA' if $picbr;
				next;
			}
		}
	}
	return \%homlocs;
}

sub get_homlocs_pc_present {
	# File name
	my($homlocs_file) = @_;
	# Get homlocs data
	my @homlocs_data = get_file_data_array($homlocs_file);
	# Initialize variables
	my %homlocs = ();
	# Parse homlocs file
	foreach my $line (@homlocs_data) {
		# Result line
		if ($line !~ /^$/) {
			if ($line =~ /^pic/) { next; }
			my @d = split(/\t/,$line);
			# Get fields
			my $pic_x = $d[0];
			my $loc_x = $d[1];
			my $ups_x = $d[2];
			my $dos_x = $d[3];
			my $ups_y = $d[4];
			my $dos_y = $d[5];
			my $loc_y = $d[6];
			my $rpm_y = $d[7];
			my $rpk_y = $d[8];
			my $len_y = $d[9];
			my $lnp_y = $d[10];
			my $rep_y = $d[11];
			my $gen_y = $d[12];
			# Undef if locus is not present
			if ($loc_y =~ /NA/ || $rpk_y =~ /NA/ || $rep_y =~ /NA/) {
				$homlocs{$loc_x} = '';
				next;
			}
			$gen_y = 0 if $gen_y =~ /NA/;
			# Link homologous loci
			if ($rpm_y > 100 && $rpk_y > 9 && $len_y < 999_999 && $len_y > 5000 && $lnp_y < 500 && $lnp_y > 20 && $rep_y > 25 && $gen_y < 25) {
				$homlocs{$loc_x} = $loc_y;
			} else {
				$homlocs{$loc_x} = '';
			}
		}
	}
	return \%homlocs;
}

sub count_matrix_pics {
	# Take piC matrix
	my($picmatrix, $matrixbranches) = @_;
	# Count variables
	my $n_pics_all = 0;
	my $n_pics_uni = 0;
	my %uni_matrix_pics = ();
	# Go through each matrix branch
	foreach my $br (@{$matrixbranches}) {
		# Go through each species
		foreach my $i (0..(scalar(@{$picmatrix})-1)) {
			# Skip if no piCs
			if ($picmatrix->[$i]->{$br}) {
				# Add piC counts
				$n_pics_all += scalar(@{$picmatrix->[$i]->{$br}});
				foreach my $pic (@{$picmatrix->[$i]->{$br}}) {
					$uni_matrix_pics{$pic} = 1;
				}
			}
		}
		# Get non-redundant piC counts
		my @uni_matrix_pics = keys %uni_matrix_pics;
		$n_pics_uni = @uni_matrix_pics;
	}
	return ($n_pics_all, $n_pics_uni);
}

sub sort_matrix_by_homologs {
	# Take piC matrix and its branches
	my($matrix, $branches) = @_;
	# Variables initialization
	my %n_homs = ();
	my @sorted_brs = ();
	# Go through each matrix branch
	foreach my $br (@{$branches}) {
		# Get number of branch homologs
		my $homs = 0;
		foreach my $i (0..(scalar(@{$matrix})-1)) {
			if ($matrix->[$i]->{$br}) {
				if (scalar(@{$matrix->[$i]->{$br}})>0) { $homs++ }
			}
		}
		$n_homs{$br} = $homs;
	}
	# Get branch list sorted by number of homologs
	foreach my $br (sort { $n_homs{$b} <=> $n_homs{$a} } keys %n_homs) {
		push(@sorted_brs,$br) if $n_homs{$br};
	}
	return @sorted_brs;
}

sub draw_global_matrix_txt {
	# Take piC matrix and its branches
	my($matrix, $branches, $outfile) = @_;
	# Open output file
	$outfile = find_unique_filename($outfile);
	my $out = open_outfile($outfile);
	# Go through each global matrix branch
	foreach my $br (@{$branches}) {
		# Get max array size for the groups in this branch
		my @sizes = ();
		foreach my $i (0..(scalar(@{$matrix})-1)) {
			push(@sizes,scalar(@{$matrix->[$i]->{$br}})) if $matrix->[$i]->{$br};
		}
		my $size_max = get_maximum(\@sizes);
		# Skip empty branch
		if ($size_max == 0) { next; }
		# Go through each putative array position
		foreach my $j (0..($size_max-1)) {
			# If this position in the group array exists print the element
			foreach my $i (0..(scalar(@{$matrix})-1)) {
				if (${$matrix->[$i]->{$br}}[$j]) {
					print($out "${$matrix->[$i]->{$br}}[$j]\t");
				} else {
					print($out "\t");
				}
			}
			print($out "\n");
		}
		print($out "\n");
	}
}

sub draw_global_matrix_html_heatmap {
	# Take piC matrix and its branches
	my($matrix, $branches, $outfile) = @_;
	# Option for skipping single piCs
	my $skip_single_clusters = 0;
	# Open output file
	$outfile = find_unique_filename($outfile);
	my $out = open_outfile($outfile);
	# Parameters
	my $width = 200;
	my $height = 500;
	my $cell_spacing = 0;
	my $cell_border = 0;
	# RGB colors
	my @no_pics = (0,0,0);
	my @one_pic = (3,101,192);
	my @many_pics = (81,167,249);
	# Print html header
	print($out "\n<table width=\"${width}px\" height=\"${height}px\" border=\"${cell_border}px\" cellspacing=\"${cell_spacing}px\">");

	# Go through each global matrix branch
	foreach my $br (@{$branches}) {
		# Skip single clusters if option true
		if ($skip_single_clusters == 1) {
			# Get number of branch homologs
			my $homs = 0;
			foreach my $i (0..(scalar(@{$matrix})-1)) {
				if ($matrix->[$i]->{$br}) {
					if (scalar(@{$matrix->[$i]->{$br}})>0) { $homs++ }
				}
			}
			if ($homs <= 1) { next; }
		}
		# Initialize RGB
		my @rgb = (0,0,0);
		# Draw heatmap fields for current branch
		print($out "\n<tr>");
		foreach my $i (0..(scalar(@{$matrix})-1)) {
			# Choose field color according to number of grouped intra-species pics
			if (scalar(@{$matrix->[$i]->{$br}}) == 1) { @rgb = @one_pic }
			if (scalar(@{$matrix->[$i]->{$br}}) >= 2) { @rgb = @many_pics }
			if (scalar(@{$matrix->[$i]->{$br}}) == 0) { @rgb = @no_pics }
			print($out "<td style=\"background-color:rgb($rgb[0],$rgb[1],$rgb[2]);\"></td>");
		}
		print($out "</tr>");
	}
	print($out "</table>");
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

# Find file name that does not already exist
# Usage: my $filename = find_unique_filename($filename);
sub find_unique_filename {
    # Take file name
    my($filename) = @_;
    # Check if file already exists
    if (-e $filename) {
    	# Filename contains a suffix
    	if ($filename =~ /.+\..+/) {
			my($prefix) = ($filename =~ /(.*)\.[^.]*$/);
			my($suffix) = ($filename =~ /.*(\.[^.]*$)/);
			my $orig_prefix = $prefix;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$prefix = $orig_prefix.'_'.$count;
				$filename = $prefix.$suffix;
			}
		}
		# Filename does not contain a suffix
		else {
			my $orig_name = $filename;
			# Find a non-existent name
			my $count = 0;
			while (-e $filename) {
				$count++;
				$filename = $orig_name;
				$filename = $orig_name.'_'.$count;
			}
		}
    }
    # Return unique filename
    return $filename;
}

# Get maximum of array values
sub get_maximum {
	# Take array list
	my($array) = @_;
	# Sort array list
	my @ordered_list = sort { $a <=> $b } @{$array};
	# Get maximum value
	my $max = $ordered_list[-1];
	# Return maximum value
	return $max;
}
