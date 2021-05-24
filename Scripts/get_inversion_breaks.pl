#!/usr/bin/perl
use strict;
use warnings;

# Global constants
# Global variables
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

if ($0 =~ /\//) {
	my($dir) = ($0 =~ /(.*)\//);
	chdir($dir);
}

my $num_locs = 0;
my @rel_breakpoints = ();
my %possible_flank_poss = ();
my @zero_pt_lens = ();
my @zero_pt_lens_hom = ();
# Go through each species
foreach my $i (0..$#species) {
	# Get species id
	my $spec_x = $species[$i];
	# Go through each species
	foreach my $j (0..$#species) {
		# Get species id
		my $spec_y = $species[$j];
		# Skip if same species
		next if $spec_x eq $spec_y;
		#next unless $spec_x eq 'Dmel' && $spec_y eq 'Dyak';
        # Get name of flankgenes file
        my $file = "${spec_x}_vs_${spec_y}.flankgenes.txt";
        # Get synteny file data
        my @flankgenes_data = get_file_data_array($file);
        # Initialize synteny data variables
        my %pic_loc = ();
        my $pic_idn = 0;
        my($ups,$dos) = (0,0);
        my %ups_gcs = ();
        my %dos_gcs = ();
        my %ups_hcs = ();
        my %dos_hcs = ();
        # Extract synteny data
        foreach my $line (@flankgenes_data) {
            if ($line =~ /_pic_/) {
        		($pic_idn) = ($line =~ /\S+_(\d+)\s/);
                my($pic_chr) = ($line =~ /(\S+):/);
                my($pic_beg) = ($line =~ /:(\d+)-/);
                my($pic_end) = ($line =~ /-(\d+)$/);
        		@{$pic_loc{$pic_idn}} = ($pic_chr,$pic_beg,$pic_end);
            } elsif ($line =~ /^Upstream/) {
                my($ups_beg) = ($line =~ /:(\d+)-/);
                my($ups_end) = ($line =~ /-(\d+)$/);
        		push(@{$pic_loc{$pic_idn}},($ups_beg,$ups_end));
        		($ups,$dos) = (1,0);
            } elsif ($line =~ /^Downstream/) {
                my($dos_beg) = ($line =~ /:(\d+)-/);
                my($dos_end) = ($line =~ /-(\d+)$/);
        		push(@{$pic_loc{$pic_idn}},($dos_beg,$dos_end));
        		($ups,$dos) = (0,1);
            } elsif ($line =~ /^$/) {
                ($ups,$dos) = (0,0);
            } elsif ($ups) {
                my @gen_locs = split(/\t/,$line);
        		unless ($gen_locs[3]) { next; }
        		my($os1_gen) = ($gen_locs[0] =~ /^(\S+)\s/);
                my($ns1_chr) = ($gen_locs[0] =~ /(\S+):/);
                my($ns1_beg) = ($gen_locs[0] =~ /:(\S+)-/);
                my($ns1_end) = ($gen_locs[0] =~ /-(\S+)/);
                my($ns1_str) = ($gen_locs[0] =~ /(\S+)$/);
                my($ns2_chr) = ($gen_locs[3] =~ /(\S+):/);
                my($ns2_beg) = ($gen_locs[3] =~ /:(\S+)-/);
                my($ns2_end) = ($gen_locs[3] =~ /-(\S+)/);
                my($ns2_str) = ($gen_locs[3] =~ /(\S+)$/);
                push(@{$ups_gcs{$pic_idn}},[$ns1_chr,$ns1_beg,$ns1_end,$ns1_str,$os1_gen]);
                push(@{$ups_hcs{$pic_idn}},[$ns2_chr,$ns2_beg,$ns2_end,$ns2_str,$os1_gen]);
            } elsif ($dos) {
                my @gen_locs = split(/\t/,$line);
        		unless ($gen_locs[3]) { next; }
        		my($os1_gen) = ($gen_locs[0] =~ /^(\S+)\s/);
                my($ns1_chr) = ($gen_locs[0] =~ /(\S+):/);
                my($ns1_beg) = ($gen_locs[0] =~ /:(\S+)-/);
                my($ns1_end) = ($gen_locs[0] =~ /-(\S+)/);
                my($ns1_str) = ($gen_locs[0] =~ /(\S+)$/);
                my($ns2_chr) = ($gen_locs[3] =~ /(\S+):/);
                my($ns2_beg) = ($gen_locs[3] =~ /:(\S+)-/);
                my($ns2_end) = ($gen_locs[3] =~ /-(\S+)/);
                my($ns2_str) = ($gen_locs[3] =~ /(\S+)$/);
                push(@{$dos_gcs{$pic_idn}},[$ns1_chr,$ns1_beg,$ns1_end,$ns1_str,$os1_gen]);
                push(@{$dos_hcs{$pic_idn}},[$ns2_chr,$ns2_beg,$ns2_end,$ns2_str,$os1_gen]);
            }
        }
		# Go through each piRNA cluster
        foreach my $pic_idn (sort keys %pic_loc) {
			# Initialize variables for signed permutation
			my @sig_per_1 = ();
			my %sig_per_locs_1 = ();
	        my %sig_per_2 = ();
	        my %hom_glocs = ();
			my %hom_loc_borders = ();
            # Gene permuation number
        	my $i_gen = 0;
			my $i_gen_ups_last_end;
			my $i_gen_dos_first_beg;
			my $i_gen_ups_first_beg;
			my $i_gen_dos_last_end;
			# Go through each upstream gene and save permuation number
            foreach my $i (0..scalar(@{$ups_gcs{$pic_idn}})-1) {
				# Increment permuation number
    			$i_gen++;
				push(@sig_per_1,$i_gen);
    			# Get gene coordinates
                my $gen_chr = $ups_gcs{$pic_idn}[$i]->[0];
    			my $gen_beg = $ups_gcs{$pic_idn}[$i]->[1];
    			my $gen_end = $ups_gcs{$pic_idn}[$i]->[2];
    			my $gen_str = $ups_gcs{$pic_idn}[$i]->[3];
				# Save original coordinates for permuation number
                $sig_per_locs_1{$i_gen} = [$gen_beg,$gen_end];
				# Get last gene end of upstream flank
				$i_gen_ups_last_end = $gen_end;
				# Get first gene start of upstream flank
				$i_gen_ups_first_beg = $gen_beg unless $i_gen_ups_first_beg;
    			# Get ortholog coordinates
                my $hom_chr = $ups_hcs{$pic_idn}[$i]->[0];
                my $hom_beg = $ups_hcs{$pic_idn}[$i]->[1];
    			my $hom_end = $ups_hcs{$pic_idn}[$i]->[2];
    			my $hom_str = $ups_hcs{$pic_idn}[$i]->[3];
                my $hom_nam = $ups_hcs{$pic_idn}[$i]->[4];
				# Get signed permutation number for species 2
                my $i_per = $gen_str eq $hom_str ? $i_gen : -$i_gen;
                # Save signed permutation number for species 2 gene id
                $sig_per_2{$hom_nam} = $i_per;
				# Save coordinates for gene id
                $hom_glocs{$hom_chr}{$hom_nam} = [$hom_beg,$hom_end];
				# Get homloc ups border
				my $hom_loc_ctr = $hom_end-(($hom_end-$hom_beg+1)/2);
				undef(%hom_loc_borders);
				$hom_loc_borders{0} = [$hom_chr,$hom_beg,$hom_end,$hom_loc_ctr];
            }
			# Go through each downstream gene and save permuation number
            foreach my $i (0..scalar(@{$dos_gcs{$pic_idn}})-1) {
				# Increment permuation number
    			$i_gen++;
				push(@sig_per_1,$i_gen);
    			# Get gene coordinates
                my $gen_chr = $dos_gcs{$pic_idn}[$i]->[0];
    			my $gen_beg = $dos_gcs{$pic_idn}[$i]->[1];
    			my $gen_end = $dos_gcs{$pic_idn}[$i]->[2];
    			my $gen_str = $dos_gcs{$pic_idn}[$i]->[3];
				# Save original coordinates for permuation number
                $sig_per_locs_1{$i_gen} = [$gen_beg,$gen_end];
				# Get first gene start of downstream flank
				$i_gen_dos_first_beg = $gen_beg unless $i_gen_dos_first_beg;
				# Get last gene end of downstream flank
				$i_gen_dos_last_end = $gen_end;
    			# Get ortholog coordinates
                my $hom_chr = $dos_hcs{$pic_idn}[$i]->[0];
                my $hom_beg = $dos_hcs{$pic_idn}[$i]->[1];
    			my $hom_end = $dos_hcs{$pic_idn}[$i]->[2];
    			my $hom_str = $dos_hcs{$pic_idn}[$i]->[3];
                my $hom_nam = $dos_hcs{$pic_idn}[$i]->[4];
				# Get signed permutation number for species 2
                my $i_per = $gen_str eq $hom_str ? $i_gen : -$i_gen;
                # Save signed permutation number for species 2 gene id
                $sig_per_2{$hom_nam} = $i_per;
				# Save coordinates for gene id
                $hom_glocs{$hom_chr}{$hom_nam} = [$hom_beg,$hom_end];
				# Get homloc dos border
				my $hom_loc_ctr = $hom_end-(($hom_end-$hom_beg+1)/2);
				$hom_loc_borders{1} = [$hom_chr,$hom_beg,$hom_end,$hom_loc_ctr] unless $hom_loc_borders{1};
            }
			# Get homloc length
			my @hom_chrs = ();
			my @hom_crds = ();
			foreach my $f (sort {$hom_loc_borders{$a}[3] <=> $hom_loc_borders{$b}[3]} keys %hom_loc_borders) {
				push(@hom_chrs,$hom_loc_borders{$f}[0]);
				push(@hom_crds,[$hom_loc_borders{$f}[1],$hom_loc_borders{$f}[2]]);
			}
			if ($hom_chrs[0] eq $hom_chrs[1]) {
				my $hl_beg = $hom_crds[0]->[1];
				my $hl_end = $hom_crds[1]->[0];
				my $hl_len = $hl_end-$hl_beg+1;
				#print("$hl_beg\t$hl_end\t$hl_len\n");
				push(@zero_pt_lens_hom,$hl_len) unless $hl_len > 200_000;
			}
			# Get signed permutation array (in order of species 2)
			my @sig_per_2 = ();
			my @sig_per_2_each_chr = ();
			# Go through each chromosome in syntenic region of species 2
	        foreach my $hom_chr (sort keys %hom_glocs) {
				my @sig_per_2_curr_chr = ();
				# Go through each homologous flank gene, ordered by their coordinates on the chromosome
	            foreach my $hom_nam (sort {$hom_glocs{$hom_chr}{$a}->[0] <=> $hom_glocs{$hom_chr}{$b}->[0]} keys %{$hom_glocs{$hom_chr}}) {
					# Save permutation number to signed permutation array
	                push(@sig_per_2,$sig_per_2{$hom_nam});
					push(@sig_per_2_curr_chr,$sig_per_2{$hom_nam});
					# Get plus sign for positive permutation numbers (for output)
	                my $sign = $sig_per_2{$hom_nam} > 0 ? '+' : '';
	                print("${sign}$sig_per_2{$hom_nam} ");
	            }
				push(@sig_per_2_each_chr,\@sig_per_2_curr_chr);
	        }
			if (scalar(@sig_per_2) <= 2) { next; }
			print("\n");
			# Get break point permutation number pairs for permuation splits
			my %br_pairs_splits = ();
			# Check if signed permuation spans multiple chromosomes (treat as breaks)
			if (scalar(@sig_per_2_each_chr) > 1) {
				# Go through each split permutation
				foreach my $x (0..$#sig_per_2_each_chr) {
					# Go through each split permutation
					foreach my $y (0..$#sig_per_2_each_chr) {
						# Skip identical split permutation
						next if $x == $y;
						# Compare each first and last number
						foreach my $r (0,-1) {
							foreach my $s (0,-1) {
								# Get potential break pair
								my @br_per_pair = sort { $a <=> $b } (abs($sig_per_2_each_chr[$x]->[$r]), abs($sig_per_2_each_chr[$y]->[$s]));
								# Save break pair if within bounds and adjacent in original species 1 order
								if ($br_per_pair[0] >= 1 && $br_per_pair[1] <= $i_gen && $br_per_pair[1]-$br_per_pair[0] == 1) {
									$br_pairs_splits{"$br_per_pair[0]|$br_per_pair[1]"} = 1;
								}
							}
						}
					}
				}
			}
			# Get break point permutation number pairs (initialize with permuation split points)
			my %br_pairs = %br_pairs_splits;
	        # Go through signed permutation
	    	foreach my $i (0..($#sig_per_2-1)) {
	    	    # Get signed permutation distance
	    	    my $dist = $sig_per_2[$i+1]-$sig_per_2[$i];
				# Get previous and next signed permutation distance if possible
	            my $dist_next = $sig_per_2[$i+2] && $sig_per_2[$i+1] ? $sig_per_2[$i+2]-$sig_per_2[$i+1] : 0;
	    		my $dist_prev = $sig_per_2[$i] && $sig_per_2[$i-1] && $i ? $sig_per_2[$i]-$sig_per_2[$i-1] : 0;
				my $dist_next2 = $sig_per_2[$i+3] && $sig_per_2[$i+2] ? $sig_per_2[$i+3]-$sig_per_2[$i+2] : 0;
	    		my $dist_prev2 = $sig_per_2[$i-1] && $sig_per_2[$i-2] && $i ? $sig_per_2[$i-1]-$sig_per_2[$i-2] : 0;
				# Get plus sign for positive distance (for output)
				my $sign = $dist > 0 ? '+' : '';
	            print(" ${sign}$dist ");
	            # Ignore two adjacent breakpoints
				if ($dist_next) { if ($dist_next > 3 || $dist_next < -2) { next; } }
				if ($dist_prev) { if ($dist_prev > 3 || $dist_prev < -2) { next; } }
				if ($dist_next2) { if ($dist_next2 > 3 || $dist_next2 < -2) { next; } }
				if ($dist_prev2) { if ($dist_prev2 > 3 || $dist_prev2 < -2) { next; } }
	            # Possible break point if the distance is greater than 3 or smaller than -2 (room for order error)
	    		if ($dist > 3 || $dist < -2) {
					# Check for inversion
					my $inversion_flag = 0;
					if ($sig_per_2[$i] > 0 && $sig_per_2[$i+1] < 0) {
						$inversion_flag = 1;
					}
					if ($sig_per_2[$i] < 0 && $sig_per_2[$i+1] > 0) {
						$inversion_flag = 1;
					}
					# Get break point permutation number pairs as ordered absolute numbers
					my @br_per_pair_a = sort { $a <=> $b } (abs($sig_per_2[$i]), abs($sig_per_2[$i]+1));
					my @br_per_pair_b = sort { $a <=> $b } (abs($sig_per_2[$i+1]), abs($sig_per_2[$i+1]-1));
					# Save break point permutation number pairs if not out of bounds
					if ($br_per_pair_a[0] >= 1 && $br_per_pair_a[1] <= $i_gen) {
						$br_pairs{"$br_per_pair_a[0]|$br_per_pair_a[1]"} = 1 if $inversion_flag;
					}
					if ($br_per_pair_b[0] >= 1 && $br_per_pair_b[1] <= $i_gen) {
						$br_pairs{"$br_per_pair_b[0]|$br_per_pair_b[1]"} = 1 if $inversion_flag;
					}
	            }
	        }
			print("\n");
			# Go through each permutation break gene pair for this piC
			foreach my $br_pair (sort keys %br_pairs) {
				# Get gene ids up- and downstream of break
				my @br_i = split(/\|/,$br_pair);
				# Get absolute coordinates of the space between genes in species 1
				my $br_loc_ups = $sig_per_locs_1{$br_i[0]}->[1];
				my $br_loc_dos = $sig_per_locs_1{$br_i[1]}->[0];
				# Get absolute coordinate of center between genes in species 1
				my $br_abs = int($br_loc_dos-(($br_loc_dos-$br_loc_ups+1)/2));
				# Get the coordinate of center between genes in species 1 relative to the piC
				my $br_rel;
				# Break in upstream flank
				if ($br_abs < $i_gen_ups_last_end) {
					$br_rel = $i_gen_ups_last_end - $br_abs;
				}
				# Break in downstream flank
				elsif ($br_abs > $i_gen_dos_first_beg) {
					$br_rel = $br_abs - $i_gen_dos_first_beg;
				}
				# Break in between flanks
				elsif ($br_abs >= $i_gen_ups_last_end && $br_abs <= $i_gen_dos_first_beg) {
					$br_rel = 0;
				}
				print("$br_pair\n");
				print("$br_i[0] $br_i[1]\n");
				print("$br_loc_ups $br_loc_dos\n");
				print("$br_abs $br_rel\n");
				push(@rel_breakpoints,$br_rel);
			}
	        print("\n");
			print("\n");
			# Count possible positions for normalization of breakpoint counts
			my $max_ups_flank_len = $i_gen_ups_last_end - $i_gen_ups_first_beg;
			for (my $p = 1; $p < $max_ups_flank_len; $p++) {
				$possible_flank_poss{$p} += 1;
			}
			my $max_dos_flank_len = $i_gen_dos_last_end - $i_gen_dos_first_beg;
			for (my $p = 1; $p < $max_dos_flank_len; $p++) {
				$possible_flank_poss{$p} += 1;
			}
			$possible_flank_poss{0}++;

			my $zero_pt_len = $i_gen_dos_first_beg-$i_gen_ups_last_end+1;
			push(@zero_pt_lens,$zero_pt_len);
			# Count processed loci
			$num_locs++;
			#last if $num_locs == 535;
        }
		#last if $num_locs == 535;
	}
	#last if $num_locs == 535;
}

#print(STDERR "$num_locs\n");
my $median_loc_len = get_median(\@zero_pt_lens);
my $averag_loc_len = get_mean(\@zero_pt_lens);
my $median_loc_len_hom = get_median(\@zero_pt_lens_hom);
my $averag_loc_len_hom = get_mean(\@zero_pt_lens_hom);

printf(STDERR "%d\t%d\t%d\t%d\n", $median_loc_len,$averag_loc_len,$median_loc_len_hom,$averag_loc_len_hom);

# Open output file
my $outfile = "BREAKPOINTS.txt";
my $out = open_outfile($outfile);

my $bin_size = 15_000;
my %brp_bins = ();
foreach my $brp (sort {$a <=> $b} @rel_breakpoints) {
	for (my $i = 0; $i <= 210_000; $i += $bin_size) {
		my $min = $i+1;
		my $max = $i+$bin_size;
        if ($brp >= $min && $brp <= $max) {
            $brp_bins{$max} += 1/$possible_flank_poss{$brp};
        }
    }
	if ($brp == 0) {
		my $j = 0;
		$brp_bins{$j} += 1/$possible_flank_poss{$brp};
	}
}

foreach my $i (sort {$a <=> $b} keys %brp_bins) {
	print($out "$i\t$brp_bins{$i}\n");
}

exit;

################################# subroutines #################################

# Save fasta data as hash
# Usage: my $sequences = get_fasta_seqs($infile);
sub get_fasta_seqs {
   # Take fasta file name
   my($file,$short) = @_;
   # Variables
   my $name = '';
   my %sequences = ();
   # Open file
   my $in = open_infile($file);
   # Extract sequence and save with name
   while (my $line = <$in>) {
       $line =~ s/\s+$//; #better chomp
       if ($line =~ /^>/) {
           ($name) = ($line =~ />(.*)$/);
           if ($short) { $name =~ s/\s.*// }
       } else {
           $sequences{$name} .= $line;
       }
   }
   return \%sequences;
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

# Calculate mean of array values
sub get_mean {
	# Take array list
	my($array) = @_;
	# Number of values
	my $N = scalar(@{$array});
	if ($N == 0) { return 0 }
	# Sum values
	my $sum = get_sum($array);
	# Calculate mean value
	my $mean = $sum/$N;
	# Return mean
	return $mean;
}

# Calculate meadian of array values
sub get_median {
	# Take array list
	my($array) = @_;
	# Get median, 50th percentile
	my $median = get_percentile_value($array, 0.5);
	$median = get_mean($array) if scalar(@{$array})<3;
	# Return median
	return $median;
}
