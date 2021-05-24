#!/usr/bin/perl
use strict;
use warnings;
use lib '/Users/dgebert/Dropbox/Perlmodules';
use lib '/home/dgebert/Dropbox/Perlmodules';
use FileIO;
use FastaIO;
use SeqMan;
use BioStat;

# Constants
my $res = 100_000;
my $genomes_dir = 'Genomes';
my $fin_pic_dir = 'Final_pics';
my $fin_pic_te_dir = $fin_pic_dir.'/piC_TEs';
my $picfile_base = '.fin_pics.bed';
my $repfile_base = '_filtered_50k.fasta.out';
$|=1; #Autoflush
# Variables

# Species identifiers
my @species = ('Dmel','Dsec','Dsim','Dyak','Dana','Dpse','Dper','Dwil','Dmoj','Dvir');
my %spec_id = ();
foreach my $i (0..$#species) {
	$spec_id{$species[$i]} = $i;
}

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0\n";

foreach my $i (0..$#species) {
    #next unless $i == 0;
    # Get repeatmasker file name
	my $rep_file = "$genomes_dir/$species[$i]$repfile_base";
	# Get final pics file name
	my $pic_file = "$fin_pic_dir/$species[$i]$picfile_base";
	# Get pic loci
	my $pic_locs = get_locs_per_chr($pic_file);
	# Get repeat loci
	my $rep_data = get_repeatmask_data($rep_file);
	# Get repeat density
	my %rep_content_num = ();
	my %rep_content_per = ();
	# Go through each chromosome
	foreach my $chr (sort keys %{$rep_data}) {
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			# Get repeat length
			my $rep_len = $rep_end-$rep_beg+1;
			# Get interval
			my $interval = int($rep_beg/$res);
			# Save data for each interval
			$rep_content_num{$chr}{$interval}++;
			$rep_content_per{$chr}{$interval} += $rep_len/$res*100;
		}
	}


	### Testing ###
	my $out = FileIO::open_outfile("$species[$i].rep_per.tbl");
	my @rep_content_num = ();
	foreach my $chr (sort keys %rep_content_per) {
		foreach my $int (sort {$a <=> $b} keys %{$rep_content_per{$chr}}) {
			push(@rep_content_num,$rep_content_per{$chr}{$int});
			print($out "$rep_content_per{$chr}{$int}\n");
		}
	}
	my $mean = BioStat::get_mean(\@rep_content_num);
	my $gmin = BioStat::get_minimum(\@rep_content_num);
	my $gmax = BioStat::get_maximum(\@rep_content_num);
	my $q1st = BioStat::get_first_quartile(\@rep_content_num);
	my $medi = BioStat::get_median(\@rep_content_num);
	my $q3rd = BioStat::get_third_quartile(\@rep_content_num);
	my $hcth = $mean*0.75;
	printf("%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", $species[$i],$gmin,$q1st,$medi,$q3rd,$gmax,$mean,$hcth);
	### Testing ###


	my $out1 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/ec.rm.out.tbl");
	my $out2 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/hc.rm.out.tbl");
	my $out3 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/allpic.rm.out.tbl");
	my $out4 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/nonpic.rm.out.tbl");
	my $out5 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/hc_nonpic.rm.out.tbl");
	my %pic_reps = ();
	my $pic_rep_i = 0;
	# Go through each chromosome
	foreach my $chr (sort keys %{$rep_data}) {
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			# Get repeat length
			my $rep_len = $rep_end-$rep_beg+1;
			# Get interval
			my $interval = int($rep_beg/$res);
			# Allocate to eu- or heterochromatin according to TE density
			if ($rep_content_per{$chr}{$interval} >= $hcth) {
				foreach my $idx (0..14) {
					print($out2 "$rep->[$idx]\t");
				}
				print($out2 "\n");
			} else {
				foreach my $idx (0..14) {
					print($out1 "$rep->[$idx]\t");
				}
				print($out1 "\n");
			}
			# Allocate to pic or non-pic regions
			foreach my $pic (@{$pic_locs->{$chr}}) {
				# Get pic coordinates
				my $pic_beg = $pic->[1];
				my $pic_end = $pic->[2];
				# Check if repeat lies in pic
				if ($rep_end >= $pic_beg && $rep_beg <= $pic_end) {
					$pic_rep_i++;
					$pic_reps{$chr}{$rep} = $pic_rep_i;
				}
			}
		}
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			# Get interval
			my $interval = int($rep_beg/$res);
			# Check if repeat is within any pic
			if ($pic_reps{$chr}{$rep}) {
				foreach my $idx (0..14) {
					print($out3 "$rep->[$idx]\t");
				}
				print($out3 "\n");
			} else {
				foreach my $idx (0..14) {
					print($out4 "$rep->[$idx]\t");
				}
				print($out4 "\n");
				# Repeat not in pic but in heterochromatin
				if ($rep_content_per{$chr}{$interval} >= $hcth) {
					foreach my $idx (0..14) {
						print($out5 "$rep->[$idx]\t");
					}
					print($out5 "\n");
				}
			}
		}
	}
	close($out1);
	close($out2);
	close($out3);
	close($out4);
	close($out5);

	#my $out3b = FileIO::open_outfile("$fin_pic_te_dir/$species[$i]/toppic.rm.out.tbl");
	# Extract final pics file data
	#my $pics_data = get_tab_fields($pic_file);

	# Get repeat loci
	my $rep_data_ec = get_repeatmask_data("$fin_pic_te_dir/$species[$i]/ec.rm.out.tbl");
	my $rep_data_hc = get_repeatmask_data("$fin_pic_te_dir/$species[$i]/hc_nonpic.rm.out.tbl");
	my $rep_data_pc = get_repeatmask_data("$fin_pic_te_dir/$species[$i]/allpic.rm.out.tbl");

	my $ec_rep_landscape = get_gnom_total_te_landscape($rep_data_ec);
	my $hc_rep_landscape = get_gnom_total_te_landscape($rep_data_hc);
	my $pc_rep_landscape = get_gnom_total_te_landscape($rep_data_pc);
	my $out6 = FileIO::open_outfile("$fin_pic_te_dir/$species[$i].te_total_ls.phe.txt");
	print($out6 "div\tpic\thet\teuc\n");
	foreach my $div (0..40) {
		$pc_rep_landscape->{$div} = 0 unless $pc_rep_landscape->{$div};
		$hc_rep_landscape->{$div} = 0 unless $hc_rep_landscape->{$div};
		$ec_rep_landscape->{$div} = 0 unless $ec_rep_landscape->{$div};
		printf($out6 "%d\t%.6f\t%.6f\t%.6f\n",$div,$pc_rep_landscape->{$div},$hc_rep_landscape->{$div},$ec_rep_landscape->{$div});
	}
	close($out6);
}

exit;

################################# subroutines #################################

sub get_locs_per_chr {
	# Take name of tab file
	my($infile,$skip_header) = @_;
	# Get file data
	my @in_data = FileIO::get_file_data_array($infile);
	if ($skip_header) { shift(@in_data); }
	# Global tree hashes for each species
	my %data_fields = ();
	# Go through file data
	foreach my $line (@in_data) {
		# Get line data
        my @d = split(/\t/,$line);
        # Save data fields
        push(@{$data_fields{$d[0]}},[@d]);
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
			my $line_s = $line;
			$line_s =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line_s);
			my $chr = $d[4];
			my $cla = $d[10];
			if ($cla =~ /Satellite/ || $cla =~ /Simple/ || $cla =~ /Tandem/ || $cla =~ /Low/ || $cla =~ /Unknown/) { next; }
			push(@{$rep_data{$chr}},[@d,$line]);
		}
	}
	return \%rep_data;
}

sub get_gnom_total_te_landscape {
	# Take pic and repeat data
	my($rep_data) = @_;
	# Storage variable
	my %rep_div = ();
	my $n_rep_total = 0;
	foreach my $chr (sort keys %{$rep_data}) {
		# Go through each repeat
		foreach my $rep (@{$rep_data->{$chr}}) {
			# Get repeat coordinates
			my $rep_div = $rep->[1];
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			my $rep_len = $rep_end-$rep_beg+1;
			# Save divergence
			my $div = int($rep_div+0.5);
			$rep_div{$div}++;
			$n_rep_total++;
		}
	}
	# Calculate relative repeats shares
	foreach my $div (keys %rep_div) {
		$rep_div{$div} = $rep_div{$div}/$n_rep_total*100;
	}
	return \%rep_div;
}
