#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $flank_rel = 0.5;
my $res = 100;
# Global variables
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 genome.fa.out pic_locs.bed\n";
unless ($ARGV[0]&&$ARGV[1]) {
   die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $rep_file = $ARGV[0];
my $pic_file = $ARGV[1];

# Get repeat loci
my $rep_locs = get_repeatmask_data($rep_file);

# Get TAD loci
my $pic_locs = get_tab_fields($pic_file);

# Open output file 1
my $outfile1 = "$pic_file.tes.tbl";
my $out1 = open_outfile($outfile1);

my %div_pre_pos = ();
my %ins_pre_pos = ();
# Go through each TAD
foreach my $pic (sort {$a <=> $b} keys %{$pic_locs}) {
    # Get TAD coordinates
    my $pic_chr = $pic_locs->{$pic}->[0];
    my $pic_beg = $pic_locs->{$pic}->[1];
    my $pic_end = $pic_locs->{$pic}->[2];
    my $pic_len = $pic_end-$pic_beg+1;
    my $ups_beg = int($pic_beg-($pic_len*$flank_rel)+1);
    my $dos_end = int($pic_end+($pic_len*$flank_rel));
    #unless ($pic_len > 20_000) { next; }
    # Get conversion for contig positions into graphic positions for sp_a
    my %conv_pos = ();
	# Calculate relative coordinate for each position
	foreach my $pos_abs ($ups_beg..$dos_end) {
	    my $pos_rel = ($pos_abs-$ups_beg+1)/(($dos_end-$ups_beg+1)/$res);
	    $pos_rel = int($pos_rel+0.5);
	    $conv_pos{$pos_abs} = $pos_rel;
	}
    print($out1 "#$pic_locs->{$pic}->[3]\n$pic_chr\t$pic_beg\t$pic_end\t$conv_pos{$pic_beg}\t$conv_pos{$pic_end}\t$ups_beg\t$dos_end\n");
    # Go through each repeat
    foreach my $rep (@{$rep_locs->{$pic_chr}}) {
        # Get repeat coordinates
        my $rep_chr = $rep->[4];
        my $rep_beg = $rep->[5];
        my $rep_end = $rep->[6];
        my $rep_len = $rep_end-$rep_beg+1;
        # Get repeat divergence
        my $rep_div = $rep->[1];
        # Check if repeat is located in/around the current TAD
        if ($ups_beg < $rep_end && $dos_end > $rep_beg) {
            # Save divergence values for relative position bin
            my $rel_rep_pos1 = $conv_pos{$rep_beg} ? $conv_pos{$rep_beg} : $conv_pos{$ups_beg};
            my $rel_rep_pos2 = $conv_pos{$rep_end} ? $conv_pos{$rep_end} : $conv_pos{$dos_end};
            print($out1 "$rep_chr\t$rep_beg\t$rep_end\t$rel_rep_pos1\t$rel_rep_pos2\t$rep->[9]\t$rep->[10]\t$rep->[1]\n");
            foreach my $rel_rep_pos ($rel_rep_pos1..$rel_rep_pos2) {
                push(@{$div_pre_pos{$rel_rep_pos}},$rep_div);
                push(@{$ins_pre_pos{$rel_rep_pos}},1);
            }
        } elsif ($dos_end < $rep_beg) {
            last;
        }
    }
}
close($out1);

# Open output file 2
my $outfile2 = "$pic_file.tedivs.tbl";
my $out2 = open_outfile($outfile2);

my $perc_pos = -50;
foreach my $pos (sort {$a <=> $b} keys %div_pre_pos) {
    my $mean_div = get_mean($div_pre_pos{$pos});
    my $totl_ins = get_sum($ins_pre_pos{$pos});
    printf($out2 "%d\t%.2f\t%d\n", $perc_pos,$mean_div,$totl_ins);
    $perc_pos += 200/$res;
}
close($out2);

my %div_flipped = ();
my %ins_flipped = ();
$perc_pos = -50;
foreach my $pos (sort {$a <=> $b} keys %div_pre_pos) {
    my $mean_div = get_mean($div_pre_pos{$pos});
    my $totl_ins = get_sum($ins_pre_pos{$pos});

    $div_flipped{$perc_pos} += $mean_div/2;
    $ins_flipped{$perc_pos} += $totl_ins/2;

    if ($pos < $res/2) {
        $perc_pos += 200/$res;
    } else {
        $perc_pos -= 200/$res;
    }
}

foreach my $perc_pos (sort {$b <=> $a} keys %div_flipped) {
    $div_flipped{$perc_pos} = $div_flipped{$perc_pos}*2;
    $ins_flipped{$perc_pos} = $ins_flipped{$perc_pos}*2;
    last;
}

my $outfile3 = "$pic_file.tedivs.flip.tbl";
my $out3 = open_outfile($outfile3);

foreach my $perc_pos (sort {$a <=> $b} keys %div_flipped) {
    printf($out3 "%d\t%.2f\t%d\n", $perc_pos,$div_flipped{$perc_pos},$ins_flipped{$perc_pos});
}

close($out3);

exit;

################################# subroutines #################################

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file) = @_;
	# Storage variable
	my %rep_data = ();
	# Get file data
	my @repeatmask_data = get_file_data_array($repeatmask_file);
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
        next unless $d[1] =~ /\d+/;
        # Save data fields
        @{$data_fields{$id_i}} = @d;
		$id_i++;
    }
	# Return data fields
	return \%data_fields;
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

# Calculate sum of array values
sub get_sum {
	# Take array list
	my($array) = @_;
	# Sum values
	my $sum = 0;
	grep { $sum += $_ } @{$array};
	# Return sum
	return $sum;
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
