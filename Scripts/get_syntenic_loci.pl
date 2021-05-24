#!/usr/bin/perl
use strict;
use warnings;
# Modules for graphical output
use LWP::Simple;
use GD::Simple;

# Global constants
my %specs_inc = ('Dmel'=>1,'Dsec'=>1,'Dsim'=>1,'Dyak'=>1,'Dana'=>1,'Dpse'=>1,'Dper'=>1,'Dwil'=>1,'Dmoj'=>1,'Dvir'=>1);
my $flank_len = 200_000;
my $min_flank_len = 50_000;
my $outfmt = '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore';
# Global variables
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <fin_pics1.bed> <genome1_new.fas> <genome1_old.fas> <genes1_old.gff> <FB_orthologs.tsv> <genome2_new.fas> <genome2_old.fas> <cds2.fas>\n";
unless ($ARGV[0]&&$ARGV[1]&&$ARGV[2]&&$ARGV[3]&&$ARGV[4]&&$ARGV[5]&&$ARGV[6]&&$ARGV[7]&&$ARGV[8]&&$ARGV[9]&&$ARGV[10]&&$ARGV[11]) {
	die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $pc_loc_s1_new_file = $ARGV[0];
my $genome_s1_new_file = $ARGV[1];
my $genome_s1_old_file = $ARGV[2];
my $genset_s1_old_file = $ARGV[3];
my $orthologs_old_file = $ARGV[4];
my $genome_s2_new_file = $ARGV[5];
my $genome_s2_old_file = $ARGV[6];
my $cdseqs_s2_old_file = $ARGV[7];
my $spec1_rep_out_file = $ARGV[8];
my $spec2_rep_out_file = $ARGV[9];
my $spec1_bin_rds_file = $ARGV[10];
my $spec2_bin_rds_file = $ARGV[11];

# Get species ids
my($spec1_id) = ($genome_s1_new_file =~ /^([^_]+)_/);
my($spec2_id) = ($genome_s2_new_file =~ /^([^_]+)_/);

## Create directories for output files
# Make directory for flank blast
my $flank_dir = "${spec1_id}_pic_flanks";
mkdir($flank_dir) unless (-e $flank_dir);
# Make directory for orthologs blast
my $ortho_dir = "$flank_dir/${spec2_id}_orthologs";
mkdir($ortho_dir);
# Make directory for orthologs blast
my $synrg_dir = "$flank_dir/${spec2_id}_syn_regions";
mkdir($synrg_dir);

# Get pic coordinates for new genome assembly (spec1)
my $pic_locs_s1_new = get_tab_fields($pc_loc_s1_new_file,1);
# Extract genomic sequences from new genome assembly (spec1)
my $chr_seqs_s1_new = get_fasta_seqs($genome_s1_new_file,1);
# Extract pic flank sequences from new genome assembly (spec1)
my $flank_seqs_s1_new = get_pic_flank_seqs($pic_locs_s1_new,$chr_seqs_s1_new,$flank_len,$min_flank_len);
# Search corresponding pic flank sequences in old assembly (spec1)
my $flank_locs_s1_old = get_flank_locs_old($flank_seqs_s1_new,$genome_s1_old_file,$spec1_id,$flank_dir);
# Get gene locations in old assembly (spec1)
my $gen_locs_s1_old = get_gff_gene_locs($genset_s1_old_file);
# Get gene locations in new assembly (spec1)
my $gen_locs_s1_new = get_gene_locs_new($pic_locs_s1_new,$flank_seqs_s1_new,$flank_locs_s1_old,$gen_locs_s1_old,$genome_s1_old_file);
# Get orthologs for each species
my $orthologs = get_orthologs($orthologs_old_file);
# Get ortholog locations
my $ortho_locs = get_ortholog_locs($orthologs_old_file);
# Get pic flank ortholog sequences from old assembly (spec2)
my $flank_ortho_seqs = get_flank_ortho_seqs($flank_locs_s1_old,$gen_locs_s1_old,$orthologs,$ortho_locs,$genome_s2_old_file,$spec1_id,$ortho_dir);
# Get locations of flank ortholog sequences in new assembly (spec2)
my $flank_ortho_locs_s2_new = get_flank_orthos_new($flank_ortho_seqs,$flank_locs_s1_old,$spec1_id,$spec2_id,$ortho_dir);

# Extract genomic sequences from genome assemblies
my $chr_seqs_s1 = get_fasta_seqs($genome_s1_new_file,1);
my $chr_seqs_s2 = get_fasta_seqs($genome_s2_new_file,1);

# Extract repeat data
my $repeats_s2 = get_repeatmask_data($spec2_rep_out_file);

# Get genome-wide read counts
my($global_rds1,$bin_rds_s1,$pls_rds_s1,$mns_rds_s1) = get_bed_reads($spec1_bin_rds_file);
my($global_rds2,$bin_rds_s2,$pls_rds_s2,$mns_rds_s2) = get_bed_reads($spec2_bin_rds_file);

## Extract definite flanking genes
my %pic_loc = ();
my %flk_gens_ns1 = ();
my %flk_gens_ns2 = ();
my %flk_gens_os1 = ();
my %flk_gens_os2 = ();
# Go through each pic
foreach my $pic (sort {$a <=> $b} keys %{$flank_ortho_locs_s2_new}) {
	# Get pic coordinates
	my $pic_chr = $pic_locs_s1_new->{$pic}->[0];
	my $pic_beg = $pic_locs_s1_new->{$pic}->[1];
	my $pic_end = $pic_locs_s1_new->{$pic}->[2];
	@{$pic_loc{$pic}} = ($pic_chr,$pic_beg,$pic_end);
	# Get upstream coordinates of pic
	my $ups_beg = $flank_seqs_s1_new->{$pic}->[0]->[2];
	my $ups_end = $flank_seqs_s1_new->{$pic}->[0]->[3];
	push(@{$pic_loc{$pic}},($ups_beg,$ups_end));
	# Get downstream coordinates of pic
	my $dos_beg = $flank_seqs_s1_new->{$pic}->[1]->[2];
	my $dos_end = $flank_seqs_s1_new->{$pic}->[1]->[3];
	push(@{$pic_loc{$pic}},($dos_beg,$dos_end));
	# Go through each flank
	foreach my $i (0..1) {
		# Sort flank contig regions of old assembly by position on new assembly
		my @chrs_os1 = sort {$flank_locs_s1_old->{$pic}->[$i]->{$a}->[0] <=> $flank_locs_s1_old->{$pic}->[$i]->{$b}->[0]} keys %{$flank_locs_s1_old->{$pic}->[$i]};
		# Go through each flank contig of old assembly
		foreach my $chr (@chrs_os1) {
			# Get direction of old assembly contig relative to new assembly contig
			my $chr_str = $flank_locs_s1_old->{$pic}->[$i]->{$chr}->[6];
			# Sort genes by position on contig
			my @genes_os1 = sort {$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$a}->[0] <=> $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$b}->[0]} keys %{$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}};
			# Go through each flanking gene on contig and search ortholog sequences in new assembly
			foreach my $gene (@genes_os1) {
				# Get gene coordinates
				my $gen_ns1_beg = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[0];
				my $gen_ns1_end = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[1];
				my $gen_ns1_str = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[2] eq 'plus' ? '+' : '-';
				my $gen_os1_beg = $gen_locs_s1_old->{$chr}->{$gene}->[3];
				my $gen_os1_end = $gen_locs_s1_old->{$chr}->{$gene}->[4];
				my $gen_os1_str = $gen_locs_s1_old->{$chr}->{$gene}->[6];
				# Save all flanking genes
				$flk_gens_os1{$pic}{$i}{$gene} = [$chr,$gen_os1_beg,$gen_os1_end,$gen_os1_str,$gene];
				$flk_gens_ns1{$pic}{$i}{$gene} = [$pic_chr,$gen_ns1_beg,$gen_ns1_end,$gen_ns1_str,$gene];
				# Skip if no orthologs for this gene exists
				unless ($flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene} && $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}) { next; }
				# Get ortholog coordinates
				my $ort_os2_gid = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[0];
				my $ort_ns2_chr = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[1];
				my $ort_ns2_beg = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[2];
				my $ort_ns2_end = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[3];
				my $ort_ns2_str = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[4] eq 'plus' ? '+' : '-';
				my $ort_os2_chr = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[5];
				my $ort_os2_beg = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[6];
				my $ort_os2_end = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[7];
				my $ort_os2_str = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[8];
				# Save all flanking genes
			    $flk_gens_os2{$pic}{$i}{$gene} = [$ort_os2_chr,$ort_os2_beg,$ort_os2_end,$ort_os2_str,$ort_os2_gid];
				$flk_gens_ns2{$pic}{$i}{$gene} = [$ort_ns2_chr,$ort_ns2_beg,$ort_ns2_end,$ort_ns2_str,$ort_os2_gid];
			}
		}
	}
}

## Get coordinates for each syntenic contig
my %hom_reg = ();
foreach my $pic (sort {$a <=> $b} keys %pic_loc) {
	# Skip if no flank gene homologs on both sides
	unless ($flk_gens_ns2{$pic}{0} || $flk_gens_ns2{$pic}{1}) { next; }
	# Combine upstream and downstream orthologous genes
	my @ups_hcs = ();
	my @dos_hcs = ();
	foreach my $gen (sort {$flk_gens_ns2{$pic}{0}{$a}[1] <=> $flk_gens_ns2{$pic}{0}{$b}[1]} keys %{$flk_gens_ns2{$pic}{0}}) {
		push(@ups_hcs,$flk_gens_ns2{$pic}{0}{$gen});
	}
	foreach my $gen (sort {$flk_gens_ns2{$pic}{1}{$a}[1] <=> $flk_gens_ns2{$pic}{1}{$b}[1]} keys %{$flk_gens_ns2{$pic}{1}}) {
		push(@dos_hcs,$flk_gens_ns2{$pic}{1}{$gen});
	}
	my @all_hcs = (@ups_hcs,@dos_hcs);
	# Find regional start and stop of orthologs for each syntenic contig
	my $i_g = 0;
	foreach my $gen (@all_hcs) {
		my $chr = $gen->[0];
		my $beg = $gen->[1];
		my $end = $gen->[2];
		$hom_reg{$pic}{$chr}[0] = $beg unless $hom_reg{$pic}{$chr}[0];
		$hom_reg{$pic}{$chr}[0] = $beg if $beg < $hom_reg{$pic}{$chr}[0];
		$hom_reg{$pic}{$chr}[1] = $end unless $hom_reg{$pic}{$chr}[1];
		$hom_reg{$pic}{$chr}[1] = $end if $end > $hom_reg{$pic}{$chr}[1];
		my $score = 1;
		if ($i_g == $#ups_hcs || $i_g == ($#ups_hcs-1) || $i_g == ($#ups_hcs-2)) {
			$score = 2;
		}
		if ($i_g == ($#ups_hcs+1) || $i_g == ($#ups_hcs+2) || $i_g == ($#ups_hcs+3)) {
			$score = 2;
		}
		$hom_reg{$pic}{$chr}[2] += $score;
		$i_g++;
	}
	foreach my $chr (sort keys %{$hom_reg{$pic}}) {
		my $len = length($chr_seqs_s2->{$chr});
		$hom_reg{$pic}{$chr}[0] = ($hom_reg{$pic}{$chr}[0]-200_000) > 1 ? ($hom_reg{$pic}{$chr}[0]-200_000) : 1;
		$hom_reg{$pic}{$chr}[1] = ($hom_reg{$pic}{$chr}[1]+200_000) <= $len ? ($hom_reg{$pic}{$chr}[1]+200_000) : $len;
	}
	delete($hom_reg{$pic}{'NA'});
	# Get region length
	foreach my $chr (sort keys %{$hom_reg{$pic}}) {
		$hom_reg{$pic}{$chr}[3] = $hom_reg{$pic}{$chr}[1]-$hom_reg{$pic}{$chr}[0]+1;
	}
}

## Get all gene locations in new species 2 assembly
my %syn_gns = ();
# Go through each pic
foreach my $pic (sort {$a <=> $b} keys %pic_loc) {
	# Go through each contig
	foreach my $chr (sort keys %{$hom_reg{$pic}}) {
		# Get coordinates of syntenic region
		my $hom_reg_beg = $hom_reg{$pic}{$chr}[0];
		my $hom_reg_end = $hom_reg{$pic}{$chr}[1];
		my $hom_reg_len = $hom_reg{$pic}{$chr}[3];
		# Extract sequence of syntenic region and print to file
		my $syn_seq = substr($chr_seqs_s2->{$chr},$hom_reg_beg-1,$hom_reg_len);
		# Open output file for sequence of syntenic region
		my $synseqfile = $synrg_dir.'/'.$spec1_id.'.'.$pic.'_'.$spec2_id.'.syn_region.'.$chr.'.fas';
		my $syn_out = open_outfile($synseqfile);
		# Print sequence of syntenic region and to file
		print($syn_out ">${spec1_id}_${pic}_${spec2_id}_hom_reg\n$syn_seq\n");
		close($syn_out);
		# Blast gene sequence to genome of species 2
		my $blast_syn_out = $synrg_dir.'/'.$spec1_id.'.'.$pic.'_'.$spec2_id.'.syn_region.'.$chr.'.blast';
		unless (-e $blast_syn_out) {
			system("blastn -query $cdseqs_s2_old_file -subject $synseqfile -out $blast_syn_out -outfmt \'$outfmt\' >>.log 2>&1");
		}
		# Extract filtered blast hits
		my $gsyn_hits = get_blastn_gene_hits($blast_syn_out,90);
		# Group blast hits
		my $gsyn_hit_locs = get_gene_hit_group_locs($gsyn_hits);
		# Save unique gene locations
		foreach my $reg (keys %{$gsyn_hit_locs}) {
			# Go through each gene
			foreach my $gen (sort {$gsyn_hit_locs->{$reg}->{$a}->[0] <=> $gsyn_hit_locs->{$reg}->{$b}->[0]} keys %{$gsyn_hit_locs->{$reg}}) {
				# Get gene coordinates and properties
				my ($gid) = ($gen =~ /(.*)-/);
				my $g_beg = $gsyn_hit_locs->{$reg}->{$gen}->[0]+$hom_reg{$pic}{$chr}[0]-1;
				my $g_end = $gsyn_hit_locs->{$reg}->{$gen}->[1]+$hom_reg{$pic}{$chr}[0]-1;
				my $g_str = $gsyn_hit_locs->{$reg}->{$gen}->[3] eq 'plus' ? '+' : '-';
				# Save gene if it does not overlap with other genes
				if ($syn_gns{$pic}) {
					# Check overlap with saved genes
					my $overlap = 0;
					foreach my $g (keys %{$syn_gns{$pic}{$chr}}) {
						if ($g_str eq $syn_gns{$pic}{$chr}{$g}[2] && $syn_gns{$pic}{$chr}{$g}[0] < $g_end && $syn_gns{$pic}{$chr}{$g}[1] > $g_beg) {
							$overlap = 1;
						}
					}
					# Save if no overlaps found
					if (not $overlap) {
						@{$syn_gns{$pic}{$chr}{$gid}} = ($g_beg,$g_end,$g_str);
					}
				} else {
					@{$syn_gns{$pic}{$chr}{$gid}} = ($g_beg,$g_end,$g_str);
				}
			}
		}
	}
}

# Find and delete pics for which orthologs do not exist for both flanks
my %flank_genes_present = ();
# Go through each pic
foreach my $pic (sort {$a <=> $b} keys %pic_loc) {
	# Go through each flank
	foreach my $i (0..1) {
		# Go through each gene
		foreach my $gen (sort {$flk_gens_ns1{$pic}{$i}{$a}[1] <=> $flk_gens_ns1{$pic}{$i}{$b}[1]} keys %{$flk_gens_ns1{$pic}{$i}}) {
			# Skip genes without orthologs
			unless ($flk_gens_ns2{$pic}{$i}{$gen}) { next; }
			# Orthologs are present
			$flank_genes_present{$pic}{$i} = 1;
		}
	}
	# Delete pic if orthologs do not exist for both flanks
	unless ($flank_genes_present{$pic}{0} && $flank_genes_present{$pic}{1}) {
		delete $pic_loc{$pic};
	}
}

## Output flank genes and their orthologs for each pic
my $outfile1 = "${spec1_id}_vs_${spec2_id}.flankgenes.txt";
my $out1 = open_outfile($outfile1);
# Go through each pic
foreach my $pic (sort {$a <=> $b} keys %pic_loc) {
	# Get pic coordinates
	my $pic_chr = $pic_loc{$pic}[0];
	my $pic_beg = $pic_loc{$pic}[1];
	my $pic_end = $pic_loc{$pic}[2];
	# Print pic coordinates
	print($out1 "${spec1_id}_pic_$pic $pic_chr:$pic_beg-$pic_end\n\n");
	# Get upstream coordinates of pic
	my $ups_beg = $pic_loc{$pic}[3];
	my $ups_end = $pic_loc{$pic}[4];
	# Get downstream coordinates of pic
	my $dos_beg = $pic_loc{$pic}[5];
	my $dos_end = $pic_loc{$pic}[6];
	# Go through each flank
	foreach my $i (0..1) {
		# Print flank coordinates
		print($out1 "Upstream    $pic_chr:$ups_beg-$ups_end\n") if $i == 0;
		print($out1 "Downstream  $pic_chr:$dos_beg-$dos_end\n") if $i == 1;
		# Go through each gene
		foreach my $gen (sort {$flk_gens_ns1{$pic}{$i}{$a}[1] <=> $flk_gens_ns1{$pic}{$i}{$b}[1]} keys %{$flk_gens_ns1{$pic}{$i}}) {
			# Skip genes without orthologs
			unless ($flk_gens_ns2{$pic}{$i}{$gen}) { next; }
			# Get ns1 gene coordinates
			my $ns1_chr = $flk_gens_ns1{$pic}{$i}{$gen}[0];
			my $ns1_beg = $flk_gens_ns1{$pic}{$i}{$gen}[1];
			my $ns1_end = $flk_gens_ns1{$pic}{$i}{$gen}[2];
			my $ns1_str = $flk_gens_ns1{$pic}{$i}{$gen}[3];
			my $sp1_gen = $flk_gens_ns1{$pic}{$i}{$gen}[4];
			# Get os1 gene coordinates
			my $os1_chr = $flk_gens_os1{$pic}{$i}{$gen}[0];
			my $os1_beg = $flk_gens_os1{$pic}{$i}{$gen}[1];
			my $os1_end = $flk_gens_os1{$pic}{$i}{$gen}[2];
			my $os1_str = $flk_gens_os1{$pic}{$i}{$gen}[3];
			# Get os2 gene coordinates
			my $os2_chr = $flk_gens_os2{$pic}{$i}{$gen}[0];
			my $os2_beg = $flk_gens_os2{$pic}{$i}{$gen}[1];
			my $os2_end = $flk_gens_os2{$pic}{$i}{$gen}[2];
			my $os2_str = $flk_gens_os2{$pic}{$i}{$gen}[3];
			my $sp2_gen = $flk_gens_os2{$pic}{$i}{$gen}[4];
			# Get ns2 gene coordinates
			my $ns2_chr = $flk_gens_ns2{$pic}{$i}{$gen}[0];
			my $ns2_beg = $flk_gens_ns2{$pic}{$i}{$gen}[1];
			my $ns2_end = $flk_gens_ns2{$pic}{$i}{$gen}[2];
			my $ns2_str = $flk_gens_ns2{$pic}{$i}{$gen}[3];
			# Print gene coordinates
			unless ($os1_beg) {
				$os1_chr = 'NA';
				$os1_beg = 'NA';
				$os1_end = 'NA';
				$os1_str = 'NA';
			}
			print($out1 "$sp1_gen $ns1_chr:$ns1_beg-$ns1_end $ns1_str\t$os1_chr:$os1_beg-$os1_end $os1_str\t");
			print($out1 "$sp2_gen $os2_chr:$os2_beg-$os2_end $os2_str\t$ns2_chr:$ns2_beg-$ns2_end $ns2_str\n");
		}
		print($out1 "\n");
	}
	print($out1 "\n");
}
close($out1);

my $outfile2 = "${spec1_id}_vs_${spec2_id}.hompics.txt";
my $out2 = open_outfile($outfile2);
print($out2 "pic_idn\tpic_loc\tpic_ups\tpic_dos\thom_ups\thom_dos\thom_pic\thom_rpm\thom_rpkm\thom_len\thom_lenp\trep_perc\tgen_perc\tass_dir\n");

# Output syntenic region and check for inversions
my %ort_locs_s2_new = ();
my %syn_flanks = ();
my %break_pos_rel = ();
my %break_pos_sp2 = ();
my %hom_pic_locs = ();
my %hpc_rpms = ();
my %hpc_rpkms = ();
my %rep_percs = ();
my %gen_percs = ();
my %hpc_lens = ();
my %hpc_lenp = ();
# Go through each pic
foreach my $pic (sort {$a <=> $b} keys %pic_loc) {
	# Get pic coordinates
	my $pic_chr = $pic_locs_s1_new->{$pic}->[0];
	my $pic_beg = $pic_locs_s1_new->{$pic}->[1];
	my $pic_end = $pic_locs_s1_new->{$pic}->[2];
	my $pic_len = $pic_end-$pic_beg+1;
	# Get upstream coordinates of pic
	my $ups_beg = $flank_seqs_s1_new->{$pic}->[0]->[2];
	my $ups_end = $flank_seqs_s1_new->{$pic}->[0]->[3];
	# Get downstream coordinates of pic
	my $dos_beg = $flank_seqs_s1_new->{$pic}->[1]->[2];
	my $dos_end = $flank_seqs_s1_new->{$pic}->[1]->[3];
	# Save possible syntenic contigs
	my @ns2_chrs = ();
	my %gen_poss = ();
	# Go through each flank
	foreach my $i (0..1) {
		# Collect flank genes
		my @genes = ();
		foreach my $gene (sort {$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$a}->[0] <=> $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$b}->[0]} keys %{$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}}) {
			my $gen_inc = 0;
			foreach my $chr (sort keys %{$flank_ortho_locs_s2_new->{$pic}->[$i]}) {
	 			if ($flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}) { $gen_inc = 1; }
			}
			if ($gen_inc) { push(@genes,$gene); }
		}
		# Go through each flanking gene on contig and search ortholog sequences in new assembly
		foreach my $gene (@genes) {
			# Get gene coordinates
			my $gen_ns1_beg = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[0];
			my $gen_ns1_end = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[1];
			my $gen_ns1_str = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[2] eq 'plus' ? '+' : '-';
			# Get ortholog coordinates
			my $ort_ns2_chr = '';
			my $ort_ns2_beg = 0;
			my $ort_ns2_end = 0;
			my $ort_ns2_str = '';
			foreach my $chr (sort keys %{$flank_ortho_locs_s2_new->{$pic}->[$i]}) {
				# Skip if no orthologs for this gene exist
	 			unless ($flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}) { next; }
	 			# Get ortholog coordinates
	 			$ort_ns2_chr = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[1];
	 			$ort_ns2_beg = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[2];
	 			$ort_ns2_end = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[3];
	 			$ort_ns2_str = $flank_ortho_locs_s2_new->{$pic}->[$i]->{$chr}->{$gene}->[4] eq 'plus' ? '+' : '-';
				my $ort_ns2_len = $ort_ns2_end-$ort_ns2_beg+1;
				# Get possible overlap with already saved genes
				my $ort_ns2_ovr = 0;
				foreach my $pos ($ort_ns2_beg..$ort_ns2_end) {
					if ($gen_poss{$ort_ns2_chr}{$pos}) { $ort_ns2_ovr++; }
				}
				# Do not save genes that overlap with already saved genes
				unless ($ort_ns2_ovr/$ort_ns2_len > 0.1) {
					# Save for each new assembly contig
					@{$ort_locs_s2_new{$pic}{$i}{$ort_ns2_chr}{$gene}} = ($ort_ns2_chr,$ort_ns2_beg,$ort_ns2_end,$ort_ns2_str);
				}
				# Save gene positions on ns2 contig
				foreach my $pos ($ort_ns2_beg..$ort_ns2_end) {
					$gen_poss{$ort_ns2_chr}{$pos} = 1;
				}
				last;
			}
			unless ($ort_ns2_chr) { next; }
			# Save possible syntenic contigs
			foreach my $ns2_chr (sort keys %{$ort_locs_s2_new{$pic}{$i}}) {
				# Award points for each flank gene ortholog with genes close to pic yielding more points
				unless ($ort_locs_s2_new{$pic}{$i}{$ns2_chr}{$gene}) { next; }
				my $points = 1;
				if ($i == 0) {
					if ($gene eq $genes[-1] || $gene eq $genes[-2] || $gene eq $genes[-3] || $gene eq $genes[-4]) {
						$points = 2;
					}
				} elsif ($i == 1) {
					if ($gene eq $genes[0] || $gene eq $genes[1] || $gene eq $genes[2] || $gene eq $genes[3]) {
						$points = 2;
					}
				}
				$ns2_chrs[$i]{$ns2_chr} += $points;
			}
	   	}
	}
	# Check for existence of contig for each flank
	my @up_do_chr = ();
	my @up_chr = ();
	my @do_chr = ();
	# Sort contigs by number of flank gene orthologs upstream
	foreach my $ns2_chr (sort {$ns2_chrs[0]{$b} <=> $ns2_chrs[0]{$a}} keys %{$ns2_chrs[0]}) {
		# Save contigs present for both flanks
		if ($ns2_chrs[1]{$ns2_chr}) {
			push(@up_do_chr,$ns2_chr);
		} else {
			push(@up_chr,$ns2_chr);
		}
	}
	# Sort contigs by number of flank gene orthologs downstream
	foreach my $ns2_chr (sort {$ns2_chrs[1]{$b} <=> $ns2_chrs[1]{$a}} keys %{$ns2_chrs[1]}) {
		unless (grep {$_ eq $ns2_chr} @up_do_chr) {
			push(@do_chr,$ns2_chr);
		}
	}

	## Choose syntenic contig or contigs
	my @syn_chr = ();
	# Case 1: At least one contig covering both flanks
	if (scalar(@up_do_chr) > 0) {
		# One contig contains likely the whole syntenic region
		my $best_chr = '';
		my $max_ngen = 0;
		# Find best fitting contig
		foreach my $chr (@up_do_chr) {
			# Get contig with largest number of flank gene orthologs
			my $n_gen = $ns2_chrs[0]{$chr}+$ns2_chrs[1]{$chr};
			if ($n_gen > $max_ngen) {
				$max_ngen = $n_gen;
				$best_chr = $chr;
			}
		}
		# Save syntenic contigs up- and downstream
		@syn_chr = ($best_chr,$best_chr);
	}
	# Case 2: No contigs covering both flanks
	elsif (scalar(@up_do_chr) == 0) {
		# Case 2.1: Only contigs covering upstream flank found
		if (scalar(@up_chr) > 0 && scalar(@do_chr) == 0) {
			# Contig contains likely the whole syntenic region
			@syn_chr = ($up_chr[0],'');
		}
		# Case 2.2: Only contigs covering downstream flank found
		elsif (scalar(@do_chr) > 0 && scalar(@up_chr) == 0) {
			# Contig contains likely the whole syntenic region
			@syn_chr = ('',$do_chr[0]);
		}
		# Case 2.3: Both, contigs covering either up- or downstream flank found
		elsif (scalar(@do_chr) > 0 && scalar(@up_chr) > 0) {
			# Likely broken synteny or assembly break
			@syn_chr = ($up_chr[0],$do_chr[0]);
		}
	}

	## Transform gene ortholog list to permutation
	my %syn_lst_ns1 = ();
	my %syn_lst_ns2 = ();
	my @syn_lst_ns1 = ();
	my @syn_lst_ns2 = ();
	my %all_ort_loc = ();
	my %all_ort_loc_chr = ();
	my %permut_gens = ();
	# Gene order id
	my $i_gen = 0;
	# Go through each flank
	foreach my $i (0..1) {
		# Collect flank genes
		my @genes = ();
		foreach my $gene (sort {$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$a}->[0] <=> $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$b}->[0]} keys %{$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}}) {
			my $gen_inc = 0;
			if ($syn_chr[$i] && $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}) { $gen_inc = 1; }
			if ($gen_inc) { push(@genes,$gene); }
		}
		# Go through each flanking gene on contig and search ortholog sequences in new assembly
		foreach my $gene (@genes) {
			$i_gen++;
			# Get gene coordinates
			my $gen_ns1_beg = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[0];
			my $gen_ns1_end = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[1];
			my $gen_ns1_str = $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$gene}->[2] eq 'plus' ? '+' : '-';
			# Get ortholog coordinates
			unless ($ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}) { next; }
			my $ort_ns2_beg = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}[1];
			my $ort_ns2_end = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}[2];
			my $ort_ns2_str = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}[3];
			@{$all_ort_loc{$gene}} = ($syn_chr[$i],$ort_ns2_beg,$ort_ns2_end,$ort_ns2_str);
			@{$all_ort_loc_chr{$syn_chr[$i]}{$gene}} = ($syn_chr[$i],$ort_ns2_beg,$ort_ns2_end,$ort_ns2_str);
			# Assign number to genes according to order in species 1
			$syn_lst_ns1{$gene} = $i_gen;
			push(@syn_lst_ns1,$i_gen);
			# Assign negative number if ortholog is in reverse orientation
			if ($gen_ns1_str eq $ort_ns2_str) {
				$syn_lst_ns2{$gene} = $i_gen;
			} else {
				$syn_lst_ns2{$gene} = -($i_gen);
			}
			$permut_gens{$i_gen} = $gene;
		}
	}
	# Is the assembly broken
	my $ass_broken = 0;
	if ($syn_chr[0] && $syn_chr[1] && $syn_chr[0] ne $syn_chr[1]) {
		$ass_broken = 1;
	}
	# Is the assembly reverse
	my $neg_per = 0;
	my $pos_per = 0;
	foreach my $gene (keys %syn_lst_ns2) {
		if ($syn_lst_ns2{$gene} < 0) { $neg_per++; }
		if ($syn_lst_ns2{$gene} > 0) { $pos_per++; }
	}
	my $ass_reverse = $neg_per >= $pos_per ? 1 : 0;
	my @flank_order = (0);
	# If assembly broken and reverse, change flank chr order
	if ($ass_broken && $ass_reverse) {
		@flank_order = (1,0);
	} elsif ($ass_broken) {
		@flank_order = (0,1);
	}
	# Order list of species 2 orthologs by order in species 1
	foreach my $i (@flank_order) {
		foreach my $gene (sort {$all_ort_loc_chr{$syn_chr[$i]}{$a}[1] <=> $all_ort_loc_chr{$syn_chr[$i]}{$b}[1]} keys %{$all_ort_loc_chr{$syn_chr[$i]}}) {
			push(@syn_lst_ns2,$syn_lst_ns2{$gene});
		}
	}

	## Get homologous flank coordinates
	my @flk_beg = ();
	my @flk_end = ();
	my @pic_neigbrs = ('NA','NA');
	my @pic_distals = ('NA','NA');
	foreach my $i (0..$#syn_chr) {
		if ($syn_chr[$i]) {
			# Collect flank genes
			my @genes = ();
			foreach my $gene (sort {$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$a}->[0] <=> $gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}->{$b}->[0]} keys %{$gen_locs_s1_new->{$pic}->[$i]->{$pic_chr}}) {
				my $gen_inc = 0;
				if ($ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}) { $gen_inc = 1; }
				if ($gen_inc) { push(@genes,$gene); }
			}
			# Get pic neighboring and distal genes
			if ($i == 0) {
				$pic_neigbrs[$i] = $genes[-1]; # Last gene of upstream flank
				$pic_distals[$i] = $genes[0]; # First gene of upstream flank
			} elsif ($i == 1) {
				$pic_neigbrs[$i] = $genes[0]; # First gene of downstream flank
				$pic_distals[$i] = $genes[-1]; # Last gene of downstream flank
			}
			# Get coordinates of homologous flank genes (first and last genes of original flank)
			my $gen_beg_nbr = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_neigbrs[$i]}[1];
			my $gen_end_nbr = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_neigbrs[$i]}[2];
			my $gen_beg_dst = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_distals[$i]}[1];
			my $gen_end_dst = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_distals[$i]}[2];
			# Get border coordinates of homologous flanks
			foreach my $gene (sort keys %{$ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}}) {
				my $gen_beg = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}[1];
				my $gen_end = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$gene}[2];
				$flk_beg[$i] = $flk_beg[$i] && $gen_beg > $flk_beg[$i] ? $flk_beg[$i] : $gen_beg;
				$flk_end[$i] = $flk_end[$i] && $gen_end < $flk_end[$i] ? $flk_end[$i] : $gen_end;
			}
		}
	}
	if ($flk_beg[0] && $flk_beg[1] && $flk_end[0] > $flk_beg[1] && $flk_beg[0] < $flk_end[1]) {
		foreach my $i (0..$#syn_chr) {
			if ($syn_chr[$i]) {
				# Get coordinates of homologous flank genes (first and last genes of original flank)
				my $gen_beg_nbr = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_neigbrs[$i]}[1];
				my $gen_end_nbr = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_neigbrs[$i]}[2];
				my $gen_beg_dst = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_distals[$i]}[1];
				my $gen_end_dst = $ort_locs_s2_new{$pic}{$i}{$syn_chr[$i]}{$pic_distals[$i]}[2];
				# Get border coordinates of homologous flanks (first and last genes of original flank)
				$flk_beg[$i] = $gen_beg_dst < $gen_beg_nbr ? $gen_beg_dst : $gen_beg_nbr;
				$flk_end[$i] = $gen_end_dst > $gen_end_nbr ? $gen_end_dst : $gen_end_nbr;
			}
		}
	}
	if ($syn_chr[0]) {
		@{$syn_flanks{$pic}{0}} = ($syn_chr[0],$flk_beg[0],$flk_end[0]);
	} else {
		@{$syn_flanks{$pic}{0}} = ("NA","NA","NA");
	}
	if ($syn_chr[1]) {
		@{$syn_flanks{$pic}{1}} = ($syn_chr[1],$flk_beg[1],$flk_end[1]);
	} else {
		@{$syn_flanks{$pic}{1}} = ("NA","NA","NA");
	}

	## Find most likely putative homologous pic region
	my $hom_pic_loc = "NA:NA-NA";
	# Check if neighboring flank gene orthologs have unbroken synteny
	if ($syn_chr[0] && $syn_chr[1] && $syn_chr[0] eq $syn_chr[1]) {
		# Neighboring flank gene identity permutation numbers
		my $ups_nbr = $syn_lst_ns1{$pic_neigbrs[0]};
		my $dos_nbr = $syn_lst_ns1{$pic_neigbrs[1]};
		# Search for neighboring flank gene orthologs
		foreach my $i (0..($#syn_lst_ns2-1)) {
			foreach my $j (1..4) {
				if ($syn_lst_ns2[$i+$j]) {
					if ($syn_lst_ns2[$i] eq $ups_nbr && $syn_lst_ns2[$i+$j] eq $dos_nbr) {
						# If homologous flanks overlap
						if ($flk_beg[0] < $flk_end[1] && $flk_end[0] > $flk_beg[1]) {
							$hom_pic_loc = "$syn_chr[0]:$all_ort_loc{$permut_gens{$ups_nbr}}[2]-$all_ort_loc{$permut_gens{$dos_nbr}}[1]";
						} else {
							$hom_pic_loc = "$syn_chr[0]:$flk_end[0]-$flk_beg[1]";
						}
					}
					if ($syn_lst_ns2[$i] eq ($dos_nbr*-1) && $syn_lst_ns2[$i+$j] eq ($ups_nbr*-1)) {
						# If homologous flanks overlap
						if ($flk_beg[0] < $flk_end[1] && $flk_end[0] > $flk_beg[1]) {
							$hom_pic_loc = "$syn_chr[0]:$all_ort_loc{$permut_gens{$dos_nbr}}[2]-$all_ort_loc{$permut_gens{$ups_nbr}}[1]";
						} else {
							$hom_pic_loc = "$syn_chr[0]:$flk_end[1]-$flk_beg[0]";
						}
					}
				}
			}
		}
		# Re-evaluate homologous pic loc if it overlaps with flanks
		unless ($hom_pic_loc eq "NA:NA-NA") {
			# Get homologous pic loc coordinates
			my @loc = split(/:|-/,$hom_pic_loc);
			# Look for overlaps with flanks
			my $flank_overlaps = 0;
			if ($loc[1] < $flk_end[0] && $loc[2] > $flk_beg[0]) {
				$flank_overlaps++;
			}
			if ($loc[1] < $flk_end[1] && $loc[2] > $flk_beg[1]) {
				$flank_overlaps++;
			}
			# If homologous pic loc overlaps with flanks, check to what extent
			if ($flank_overlaps) {
				my $overlap_bps = 0;
				# Save positions
				my %homloc_poss = ();
				foreach my $pos ($loc[1]..$loc[2]) { $homloc_poss{$pos} = 1; }
				my %flanks_poss = ();
				foreach my $pos ($flk_beg[0]..$flk_end[0]) { $flanks_poss{$pos} = 1; }
				foreach my $pos ($flk_beg[1]..$flk_end[1]) { $flanks_poss{$pos} = 1; }
				# Count overlapping positions
				foreach my $pos (keys %flanks_poss) {
					if ($homloc_poss{$pos}) { $overlap_bps++; }
				}
				my $flanks_bps = keys %flanks_poss;
				my $homloc_bps = keys %homloc_poss;
				# If overlap exceeds threshold (33%) undefine homologous pic loc
				if ($overlap_bps/$flanks_bps > 0.33) {
					$hom_pic_loc = "NA:NA-NA";
				}
				if ($homloc_bps > 999_999) {
					$hom_pic_loc = "NA:NA-NA";
				}
				if ($homloc_bps > ($pic_len*10)) {
					$hom_pic_loc = "NA:NA-NA";
				}
			}
			if (($loc[2]-$loc[1]+1) > ($pic_len*15)) {
				$hom_pic_loc = "NA:NA-NA";
			}
		}
	} elsif ($syn_chr[0] && $syn_chr[1]) {
		# Homologous pic locus cannot be determined
		$hom_pic_loc = "NA:NA-NA";
		# Check if assembly brake is within expected homologous pic region
		my $len_upc = length($chr_seqs_s2->{$syn_chr[0]});
		my $len_doc = length($chr_seqs_s2->{$syn_chr[1]});
		if ($len_upc-$flk_end[0] < 250_000 && $flk_beg[1] < 250_000) {
			$hom_pic_loc = "$syn_chr[0]:$flk_end[0]-$len_upc|$syn_chr[1]:1-$flk_end[1]";
		} elsif ($flk_beg[0] < 250_000 && $len_doc-$flk_end[1] < 250_000) {
			$hom_pic_loc = "$syn_chr[0]:1-$flk_beg[0]|$syn_chr[1]:$flk_end[1]-$len_doc";
		}
	} elsif ($syn_chr[0]) {
		# Neighboring upstream flank gene identity permutation number
		my $ups_nbr = $syn_lst_ns1{$pic_neigbrs[0]};
		# Search for neighboring flank gene ortholog at beginning or end
		if ($syn_lst_ns2[-1] eq $ups_nbr) {
			$hom_pic_loc = "$syn_chr[0]:$flk_end[0]-NA";
		}
		if ($syn_lst_ns2[0] eq ($ups_nbr*-1)) {
			$hom_pic_loc = "$syn_chr[0]:NA-$flk_beg[0]"
		}
	} elsif ($syn_chr[1]) {
		# Neighboring downstream flank gene identity permutation number
		my $dos_nbr = $syn_lst_ns1{$pic_neigbrs[1]};
		# Search for neighboring flank gene ortholog at beginning or end
		if ($syn_lst_ns2[0] eq $dos_nbr) {
			$hom_pic_loc = "$syn_chr[1]:NA-$flk_beg[1]";
		}
		if ($syn_lst_ns2[-1] eq ($dos_nbr*-1)) {
			$hom_pic_loc = "$syn_chr[1]:$flk_end[1]-NA";
		}
	}
	# Save putative homologous pic locus
	my @hom_pic_loc = split(/:|-|\|/,$hom_pic_loc);
	@{$hom_pic_locs{$pic}} = @hom_pic_loc;

	## Determine assembly orientation
	my $assembly_status = "";
	# Check if signed permutation includes minus genes
	my @for_minus_gens = grep {$_ < 0} @syn_lst_ns2;
	# Check if reversed permutation includes minus genes
	my @syn_lst_ns2_rev = map {$_ * -1} reverse(@syn_lst_ns2);
	my @rev_minus_gens = grep {$_ < 0} @syn_lst_ns2_rev;
	# Compare shares of minus genes in forward and reverse signed permutations
	if (scalar(@rev_minus_gens) >= scalar(@for_minus_gens)) {
		$assembly_status = "Parallel assembly";
	} else {
		$assembly_status = "Reversed assembly";
	}

	# Get rpm for putative homologous pic region
	if ($hom_pic_locs{$pic}[0] ne 'NA' && $global_rds2) {
		# Get pic coordinates
		my $hpc_chr = $hom_pic_locs{$pic}[0];
		my $hpc_beg = $hom_pic_locs{$pic}[1] eq 'NA' ? $hom_reg{$pic}{$hpc_chr}[0] : $hom_pic_locs{$pic}[1];
		my $hpc_end = $hom_pic_locs{$pic}[2] eq 'NA' ? $hom_reg{$pic}{$hpc_chr}[1] : $hom_pic_locs{$pic}[2];
		# Sum up read number falling in homologous pic region
		my $hpc_rds = 0;
		# Sum up size of putative homologous pic region
		my $hpc_len = $hpc_end-$hpc_beg+1;
		# Sum up repeat base pairs
		my $rep_bps = 0;
		# Sum up gene base pairs
		my $gen_bps = 0;
		# Go through each bin on plus strand
		foreach my $bin (sort {$a <=> $b} keys %{$pls_rds_s2}) {
			# Get pic coordinates
			my $bin_chr = $pls_rds_s2->{$bin}->[0];
			my $bin_beg = $pls_rds_s2->{$bin}->[1];
			my $bin_end = $pls_rds_s2->{$bin}->[2];
			my $bin_rds = $pls_rds_s2->{$bin}->[3];
			# Add bin reads to pic if bin is located in pic
			if ($hpc_chr eq $bin_chr) {
				if ($bin_beg >= $hpc_beg && $bin_end <= $hpc_end) {
					$hpc_rds += $bin_rds;
				}
			}
		}
		# Go through each bin on minus strand
		foreach my $bin (sort {$a <=> $b} keys %{$mns_rds_s2}) {
			# Get pic coordinates
			my $bin_chr = $mns_rds_s2->{$bin}->[0];
			my $bin_beg = $mns_rds_s2->{$bin}->[1];
			my $bin_end = $mns_rds_s2->{$bin}->[2];
			my $bin_rds = $mns_rds_s2->{$bin}->[3];
			# Add bin reads to pic if bin is located in pic
			if ($hpc_chr eq $bin_chr) {
				if ($bin_beg >= $hpc_beg && $bin_end <= $hpc_end) {
					$hpc_rds += $bin_rds;
				}
			}
		}
		# Get reapeat share for putative homologous pic region
		foreach my $rep (@{$repeats_s2->{$hpc_chr}}) {
			# Repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			# Check if repeat is in region
			if ($hpc_beg < $rep_beg && $hpc_end > $rep_end) {
				$rep_bps += $rep_end-$rep_beg+1;
			}
		}
		# Get gene share for putative homologous pic region
		foreach my $gen (sort {$syn_gns{$pic}{$hpc_chr}{$a}[0] <=> $syn_gns{$pic}{$hpc_chr}{$b}[0]} keys %{$syn_gns{$pic}{$hpc_chr}}) {
			my $gen_beg = $syn_gns{$pic}{$hpc_chr}{$gen}[0];
			my $gen_end = $syn_gns{$pic}{$hpc_chr}{$gen}[1];
			# Check if gene is in region
			if ($hpc_beg < $gen_beg && $hpc_end > $gen_end) {
				$gen_bps += $gen_end-$gen_beg+1;
			}
		}
		if ($hom_pic_locs{$pic}[3]) {
			# Get pic coordinates
			$hpc_chr = $hom_pic_locs{$pic}[3];
			$hpc_beg = $hom_pic_locs{$pic}[4] eq 'NA' ? $hom_reg{$pic}{$hpc_chr}[0] : $hom_pic_locs{$pic}[4];
			$hpc_end = $hom_pic_locs{$pic}[5] eq 'NA' ? $hom_reg{$pic}{$hpc_chr}[1] : $hom_pic_locs{$pic}[5];
			# Sum up size of putative homologous pic region
			$hpc_len += $hpc_end-$hpc_beg+1;
			# Go through each bin on plus strand
			foreach my $bin (sort {$a <=> $b} keys %{$pls_rds_s2}) {
				# Get pic coordinates
				my $bin_chr = $pls_rds_s2->{$bin}->[0];
				my $bin_beg = $pls_rds_s2->{$bin}->[1];
				my $bin_end = $pls_rds_s2->{$bin}->[2];
				my $bin_rds = $pls_rds_s2->{$bin}->[3];
				# Add bin reads to pic if bin is located in pic
				if ($hpc_chr eq $bin_chr) {
					if ($bin_beg >= $hpc_beg && $bin_end <= $hpc_end) {
						$hpc_rds += $bin_rds;
					}
				}
			}
			# Go through each bin on minus strand
			foreach my $bin (sort {$a <=> $b} keys %{$mns_rds_s2}) {
				# Get pic coordinates
				my $bin_chr = $mns_rds_s2->{$bin}->[0];
				my $bin_beg = $mns_rds_s2->{$bin}->[1];
				my $bin_end = $mns_rds_s2->{$bin}->[2];
				my $bin_rds = $mns_rds_s2->{$bin}->[3];
				# Add bin reads to pic if bin is located in pic
				if ($hpc_chr eq $bin_chr) {
					if ($bin_beg >= $hpc_beg && $bin_end <= $hpc_end) {
						$hpc_rds += $bin_rds;
					}
				}
			}
			# Get reapeat share for putative homologous pic region
			foreach my $rep (@{$repeats_s2->{$hpc_chr}}) {
				# Repeat coordinates
				my $rep_beg = $rep->[5];
				my $rep_end = $rep->[6];
				# Check if repeat is in region
				if ($hpc_beg < $rep_beg && $hpc_end > $rep_end) {
					$rep_bps += $rep_end-$rep_beg+1;
				}
			}
			# Get gene share for putative homologous pic region
			foreach my $gen (sort {$syn_gns{$pic}{$hpc_chr}{$a}[0] <=> $syn_gns{$pic}{$hpc_chr}{$b}[0]} keys %{$syn_gns{$pic}{$hpc_chr}}) {
				my $gen_beg = $syn_gns{$pic}{$hpc_chr}{$gen}[0];
				my $gen_end = $syn_gns{$pic}{$hpc_chr}{$gen}[1];
				# Check if gene is in region
				if ($hpc_beg < $gen_beg && $hpc_end > $gen_end) {
					$gen_bps += $gen_end-$gen_beg+1;
				}
			}
		}
		# Calculate rpm for putative homologous pic region
		my $hpc_rpm = (int((($hpc_rds*1_000_000)/$global_rds2)*100)/100);
		$hpc_rpms{$pic} = $hpc_rpm;
		# Calculate rpkm for putative homologous pic region
		my $hpc_rpkm = (int(($hpc_rpm/($hpc_len/1000))*100)/100);
		$hpc_rpkms{$pic} = $hpc_rpkm;
		# Calculate repeat share for putative homologous pic region
		my $rep_perc = (int(($rep_bps/$hpc_len*100)*100)/100);
		$rep_percs{$pic} = $rep_perc;
		# Calculate repeat share for putative homologous pic region
		my $gen_perc = (int(($gen_bps/$hpc_len*100)*100)/100);
		$gen_percs{$pic} = $gen_perc;
		# Save length of putative homologous pic region
		$hpc_lens{$pic} = $hpc_len;
		# Save length proportion of putative homologous pic region to original pic
		$hpc_lenp{$pic} = (int(($hpc_len/$pic_len*100)*100)/100);
	}
	$hpc_rpms{$pic} = 'NA' unless $hpc_rpms{$pic};
	$hpc_rpkms{$pic} = 'NA' unless $hpc_rpkms{$pic};
	$rep_percs{$pic} = 'NA' unless $rep_percs{$pic};
	$gen_percs{$pic} = 'NA' unless $gen_percs{$pic};
	$hpc_lens{$pic} = 'NA' unless $hpc_lens{$pic};
	$hpc_lenp{$pic} = 'NA' unless $hpc_lenp{$pic};

	### Output 2: pic_idn pic_loc pic_ups pic_dos hom_ups hom_dos hom_pic ass_dir breakps
	# Output original pic information
	print($out2 "${spec1_id}_pic_$pic\t$pic_chr:$pic_beg-$pic_end\t");
	print($out2 "$pic_chr:$ups_beg-$ups_end\t$pic_chr:$dos_beg-$dos_end\t");
	# Output identified homologous flank regions
	print($out2 "$syn_flanks{$pic}{0}[0]:$syn_flanks{$pic}{0}[1]-$syn_flanks{$pic}{0}[2]\t");
	print($out2 "$syn_flanks{$pic}{1}[0]:$syn_flanks{$pic}{1}[1]-$syn_flanks{$pic}{1}[2]\t");
	# Output identified homologous pic region
	print($out2 "$hom_pic_loc\t");
	# Output rpm for identified homologous pic region
	print($out2 "$hpc_rpms{$pic}\t");
	# Output rpkm for identified homologous pic region
	print($out2 "$hpc_rpkms{$pic}\t");
	# Output length for identified homologous pic region
	print($out2 "$hpc_lens{$pic}\t");
	# Output length proportion of putative homologous pic region to original pic
	print($out2 "$hpc_lenp{$pic}\t");
	# Output repeat share for identified homologous pic region
	print($out2 "$rep_percs{$pic}\t");
	# Output gene share for identified homologous pic region
	print($out2 "$gen_percs{$pic}\t");
	# Output assembly orientation
	print($out2 "$assembly_status\n");
	# Output number of breakpoints
	#print($out2 "$n_breakpoints\t");
	# Output number of breakpoints
	#print($out2 "$pic_break\n");
}
close($out2);

############### Graphics output ###############

# Make directory for graphical output
my $gd_dir = "${spec1_id}_vs_${spec2_id}_graphics";
mkdir($gd_dir);

# Collect command line arguments
my $pic_synteny_file = $outfile1;
my $spec1_rmout_file = $spec1_rep_out_file ? $spec1_rep_out_file : '';
my $spec2_rmout_file = $spec2_rep_out_file ? $spec2_rep_out_file : '';

# Prefix for image file name
my($img_prefix) = ($pic_synteny_file =~ /^(\w+)\./);

# Get synteny file data
my @pic_synteny_data = get_file_data_array($pic_synteny_file);

# Get repeats
my $rep_data_s1 = get_repeatmask_data($spec1_rmout_file) if $spec1_rmout_file;
my $rep_data_s2 = get_repeatmask_data($spec2_rmout_file) if $spec2_rmout_file;

# Initialize synteny data variables
%pic_loc = ();
my $pic_idn = 0;
my($ups,$dos) = (0,0);
my %ups_gcs = ();
my %dos_gcs = ();
my %ups_hcs = ();
my %dos_hcs = ();
my %gen_cols = ();
# Extract synteny data
foreach my $line (@pic_synteny_data) {
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
		$gen_cols{$pic_idn}{$os1_gen} = [0,134,43] if $ns1_str eq '+';
		$gen_cols{$pic_idn}{$os1_gen} = [78,202,101] if $ns1_str eq '-';
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
		$gen_cols{$pic_idn}{$os1_gen} = [43,39,182] if $ns1_str eq '+';
		$gen_cols{$pic_idn}{$os1_gen} = [120,125,216] if $ns1_str eq '-';
    }
}

# Get all homologous genes (sorted)
my %all_hcs = ();
foreach my $pid (sort {$a <=> $b} keys %pic_loc) {
	foreach my $chr (sort {$hom_reg{$pid}{$b}[2] <=> $hom_reg{$pid}{$a}[2]} keys %{$hom_reg{$pid}}) {
		# Upstream orthologous flank genes
	    foreach my $gen (@{$ups_hcs{$pid}}) {
	        unless ($gen->[0] eq $chr) { next; }
			push(@{$all_hcs{$pid}{$chr}},$gen);
	    }
		# Downstream orthologous flank genes
	    foreach my $gen (@{$dos_hcs{$pid}}) {
	        unless ($gen->[0] eq $chr) { next; }
			push(@{$all_hcs{$pid}{$chr}},$gen);
	    }
		# Sort homologous genes by position
		@{$all_hcs{$pid}{$chr}} = sort {$a->[1] <=> $b->[1]} @{$all_hcs{$pid}{$chr}};
	}
}

# Get conversion for contig positions into graphic positions for sp_a
my %convx_a = ();
foreach my $pid (sort keys %pic_loc) {
	# Get regional start and stop positions
	my $reg_beg = $pic_loc{$pid}[3];
	my $reg_end = $pic_loc{$pid}[6];
	# Calculate graphic coordinate for each position
	foreach my $xo ($reg_beg..$reg_end) {
	    my $xn = ($xo-$reg_beg+1)/(($reg_end-$reg_beg+1)/1000);
	    $xn = int($xn+0.5);
	    $convx_a{$pid}{$xo} = $xn;
	}
}
# Get conversion for contig positions into graphic positions for sp_b
my %slices = ();
my %region_lens = ();
my %convx_b = ();
foreach my $pid (sort keys %pic_loc) {
	foreach my $chr (sort keys %{$hom_reg{$pid}}) {
		# Initialize slice number
		my $s_i = 1;
		# Beginning of first slice
		$slices{$pid}{$chr}{1}[0] = $hom_reg{$pid}{$chr}[0];
		# Identify genes in homologous regions that are far apart (>200_000)
		my $prev_gen = '';
		foreach my $gen (@{$all_hcs{$pid}{$chr}}) {
			# Get current gene coordinates
			my $gc_beg = $gen->[1];
			my $gc_end = $gen->[2];
			# Skip if previous gene is not yet defined
			if ($prev_gen) {
				# Get previous gene coordinates
				my $gp_beg = $prev_gen->[1];
				my $gp_end = $prev_gen->[2];
				# Check if distance between current and previous gene is above threshold
				if ($gc_beg-$gp_end > 400_000) {
					$slices{$pid}{$chr}{$s_i}[1] = $gp_beg+(400_000/2);
					$slices{$pid}{$chr}{$s_i+1}[0] = $gc_end-(400_000/2);
					$s_i++;
				}
			}
			# Save current genes as next previous gene
			$prev_gen = $gen;
		}
		# End of last slice
		$slices{$pid}{$chr}{$s_i}[1] = $hom_reg{$pid}{$chr}[1];
		# Get combined length of all slices
		foreach my $slice (sort {$a <=> $b} keys %{$slices{$pid}{$chr}}) {
			# Get slice start and stop positions
			my $reg_beg = $slices{$pid}{$chr}{$slice}[0];
			my $reg_end = $slices{$pid}{$chr}{$slice}[1];
			# Get slice lengths
			$region_lens{$pid}{$chr} += $reg_end-$reg_beg+1;
		}
		my $pre_prc = 0;
		# Go through each slice
		foreach my $slice (sort {$a <=> $b} keys %{$slices{$pid}{$chr}}) {
			# Get slice start and stop positions
			my $reg_beg = $slices{$pid}{$chr}{$slice}[0];
			my $reg_end = $slices{$pid}{$chr}{$slice}[1];
			# Get slice length
			my $reg_len = $reg_end-$reg_beg+1;
			# Get slice length percentage
			my $reg_prc = $reg_len/$region_lens{$pid}{$chr};
			# Calculate graphic coordinate for each position
			foreach my $xo ($reg_beg..$reg_end) {
				my $xn = ($xo-$reg_beg+1)/(($reg_end-$reg_beg+1)/(1000*$reg_prc))+((1000*$pre_prc));
				$xn = int($xn+0.5);
				$convx_b{$pid}{$chr}{$xo} = $xn;
			}
			$pre_prc += $reg_prc;
		}
	}
}

# Create graphical output on synteny for each pic
# Global graphic coordinates
my $l_space = 20;
foreach my $pid (sort keys %pic_loc) {
	# Skip if no flank gene homologs on either side
	unless ($ups_hcs{$pid} || $dos_hcs{$pid}) { next; }
	my $pic_chr = $pic_loc{$pid}[0];
	# Create a new image (width, height)
	my($width,$height) = (1050,1000);
	my $img = GD::Simple->new($width,$height);
	## Species 1 contig region
	# Caption
	$img->bgcolor('black');
	$img->fgcolor('black');
	$img->moveTo($l_space, 30);
	$img->fontsize(10);
	my $sp1_chr = $pic_loc{$pid}[0];
	my $sp1_beg = add_thousands_separators($pic_loc{$pid}[3]);
	my $sp1_end = add_thousands_separators($pic_loc{$pid}[6]);
	$img->string("$sp1_chr:$sp1_beg-$sp1_end");
	# Scale
	my $scale_len = 100;
	$img->moveTo(900+$l_space, 30);
	$img->lineTo(900+$l_space+$scale_len, 30);
	$img->moveTo(900+$l_space+(($scale_len/2)-10), 30);
	my $scale_bps = (($pic_loc{$pid}[6]-$pic_loc{$pid}[3]+1)/1000)*$scale_len;
	my $scale_kbs = int($scale_bps/1000);
	$img->string("${scale_kbs}kb");
	# Contig of pic and flanks
	$img->bgcolor(255,255,255);
	$img->fgcolor('black');
	$img->rectangle($l_space, 10+30, 1000+$l_space, 40+30);
	# Y-axis labels
	$img->fontsize(10);
	$img->moveTo($l_space-10, 40+10);
	$img->string("+");
	$img->moveTo($l_space-10, 70+2);
	$img->string("-");
	# Contig of pic
	$img->bgcolor('gray');
	$img->fgcolor('gray');
	# Horizontal line
	$img->moveTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[1]}, 50-15);
	$img->lineTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[2]}, 50-15);
	# Vertical line
	$img->moveTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[1]}, 50-15);
	$img->lineTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[1]}, 50-13);
	# Vertical line
	$img->moveTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[2]}, 50-15);
	$img->lineTo($l_space+$convx_a{$pid}{$pic_loc{$pid}[2]}, 50-13);

	# Load pdf module
	use PDF::Create;
	use PDF::Create::Page;
	# Create pdf file
	my $pdf = PDF::Create->new(
	    'filename'     => "$gd_dir/$img_prefix.pic_$pid.pdf",
	    'Author'       => 'Daniel Gebert',
	    'Title'        => "$img_prefix.pic_$pid.pdf",
	    'CreationDate' => [ localtime ]
	);
	# Add a A4 sized page
	my $root = $pdf->new_page('MediaBox' => $pdf->get_page_size('A2'));
	# Add a page which inherits its attributes from $root
	my $page1 = $root->new_page;
	# Prepare a font
	my $font = $pdf->font('BaseFont' => 'Helvetica');
	## Species 1 contig region
	# Caption
	$page1->string($font, 15, 40, 1650-20, "$sp1_chr:$sp1_beg-$sp1_end");
	# Scale pdf
	$page1->setrgbcolorstroke(0.1,0.1,0.1);
	$page1->line(900+40, 1650-20, 900+40+$scale_len, 1650-20);
	$page1->stringc($font, 15, 900+$scale_len-10, 1650-16, "${scale_kbs}kb");
	# piC loc
	$page1->setrgbcolorstroke(0.6,0.6,0.6);
	$page1->line(40+$convx_a{$pid}{$pic_loc{$pid}[1]}, 1650-(20)-8, 40+$convx_a{$pid}{$pic_loc{$pid}[2]}, 1650-(20)-8);
	$page1->line(40+$convx_a{$pid}{$pic_loc{$pid}[1]}, 1650-(20)-8, 40+$convx_a{$pid}{$pic_loc{$pid}[1]}, 1650-(3+20)-8);
	$page1->line(40+$convx_a{$pid}{$pic_loc{$pid}[2]}, 1650-(20)-8, 40+$convx_a{$pid}{$pic_loc{$pid}[2]}, 1650-(3+20)-8);

	# Positions of other genes in region
	foreach my $i (0..1) {
	    # Sort genes by position on contig
	    my @genes_os1 = sort {$gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}->{$a}->[0] <=> $gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}->{$b}->[0]} keys %{$gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}};
	    # Go through each flanking gene on contig and search ortholog sequences in new assembly
	    foreach my $gene (@genes_os1) {
	        # Get gene coordinates
	        my $g_beg = $gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}->{$gene}->[0];
	        my $g_end = $gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}->{$gene}->[1];
	        my $g_str = $gen_locs_s1_new->{$pid}->[$i]->{$pic_chr}->{$gene}->[2] eq 'plus' ? '+' : '-';
			# Check overlap with flank genes
			my $overlap = 0;
			foreach my $g (@{$ups_gcs{$pid}}) {
				if ($g_str eq $g->[3] && $g->[1] < $g_end && $g->[2] > $g_beg) {
					$overlap = 1;
				}
			}
			foreach my $g (@{$dos_gcs{$pid}}) {
				if ($g_str eq $g->[3] && $g->[1] < $g_end && $g->[2] > $g_beg) {
					$overlap = 1;
				}
			}
			# Skip if overlaps found
			if ($overlap) { next; }
			$img->bgcolor(210,210,210);
	        $img->fgcolor(210,210,210);
			if ($convx_a{$pid}{$g_end}-$convx_a{$pid}{$g_beg} <= 1) { $img->fgcolor(210,210,210); }
	        if ($g_str eq '+') {
	            $img->rectangle($l_space+$convx_a{$pid}{$g_beg}, 10+30, $l_space+$convx_a{$pid}{$g_end}, 24+30);
				$page1->setrgbcolor(240/255,240/255,240/255);
		        $page1->rectangle(40+$convx_a{$pid}{$g_beg}, 1650-(20+30), ($convx_a{$pid}{$g_end}-$convx_a{$pid}{$g_beg}+1), 15);
		        $page1->fill();
	        } elsif ($g_str eq '-') {
	            $img->rectangle($l_space+$convx_a{$pid}{$g_beg}, 26+30, $l_space+$convx_a{$pid}{$g_end}, 40+30);
				$page1->setrgbcolor(240/255,240/255,240/255);
		        $page1->rectangle(40+$convx_a{$pid}{$g_beg}, 1650-(35+30), ($convx_a{$pid}{$g_end}-$convx_a{$pid}{$g_beg}+1), 15);
		        $page1->fill();
	        }
	    }
	}
	# Positions of upstream flank genes
	foreach my $gen (@{$ups_gcs{$pid}}) {
		$img->bgcolor(@{$gen_cols{$pid}{$gen->[4]}});
	    $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}});
		$page1->setrgbcolor($gen_cols{$pid}{$gen->[4]}->[0]/255,$gen_cols{$pid}{$gen->[4]}->[1]/255,$gen_cols{$pid}{$gen->[4]}->[2]/255);
		#if ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]} <= 1) { $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}}); }
	    if ($gen->[3] eq '+') {
	        $img->rectangle($l_space+$convx_a{$pid}{$gen->[1]}, 10+30, $l_space+$convx_a{$pid}{$gen->[2]}, 24+30);
	        $page1->rectangle(40+$convx_a{$pid}{$gen->[1]}, 1650-(20+30), ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]}+1), 15);
	        $page1->fill();
	    } elsif ($gen->[3] eq '-') {
	        $img->rectangle($l_space+$convx_a{$pid}{$gen->[1]}, 26+30, $l_space+$convx_a{$pid}{$gen->[2]}, 40+30);
	        $page1->rectangle(40+$convx_a{$pid}{$gen->[1]}, 1650-(35+30), ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]}+1), 15);
	        $page1->fill();
	    }
	}
	# Positions of downstream flank genes
	foreach my $gen (@{$dos_gcs{$pid}}) {
		$img->bgcolor(@{$gen_cols{$pid}{$gen->[4]}});
	    $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}});
		$page1->setrgbcolor($gen_cols{$pid}{$gen->[4]}->[0]/255,$gen_cols{$pid}{$gen->[4]}->[1]/255,$gen_cols{$pid}{$gen->[4]}->[2]/255);
		#if ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]} <= 1) { $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}}); }
	    if ($gen->[3] eq '+') {
	        $img->rectangle($l_space+$convx_a{$pid}{$gen->[1]}, 10+30, $l_space+$convx_a{$pid}{$gen->[2]}, 24+30);
	        $page1->rectangle(40+$convx_a{$pid}{$gen->[1]}, 1650-(20+30), ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]}+1), 15);
	        $page1->fill();
	    } elsif ($gen->[3] eq '-') {
	        $img->rectangle($l_space+$convx_a{$pid}{$gen->[1]}, 26+30, $l_space+$convx_a{$pid}{$gen->[2]}, 40+30);
	        $page1->rectangle(40+$convx_a{$pid}{$gen->[1]}, 1650-(35+30), ($convx_a{$pid}{$gen->[2]}-$convx_a{$pid}{$gen->[1]}+1), 15);
	        $page1->fill();
	    }
	}
	# Contig rectangle
	$page1->setrgbcolorstroke(0.1,0.1,0.1);
	$page1->rectangle(40,1650-65,1000,30);
	$page1->stroke();
	# Contig of pic and flanks
	$img->bgcolor(undef);
	$img->fgcolor('black');
	$img->rectangle($l_space, 10+30, 1000+$l_space, 40+30);
	# Repeat annotation
	# Y-axis line
	$img->bgcolor('black');
	$img->fgcolor('black');
	$img->moveTo($l_space, 80);
	$img->lineTo($l_space, 110);
	# Y-axis labels
	$img->fontsize(10);
	$img->moveTo($l_space-10, 80+10);
	$img->string("+");
	$img->moveTo($l_space-10, 110+2);
	$img->string("-");
	if ($spec1_rmout_file) {
		# Repeat positions
		foreach my $rep (@{$rep_data_s1->{$pic_loc{$pid}[0]}}) {
			# Repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_str = $rep->[8];
			my $rep_div = $rep->[1];
			# Check if repeat is in region
			if ($pic_loc{$pid}[3] < $rep_beg && $pic_loc{$pid}[6] > $rep_end) {
				# Calculate gray shade according to sequence divergence
				my @rgb = (255*($rep_div/40),255*($rep_div/40),255*($rep_div/40));
				# Get repeat position to consensus sequence
				my @con_pos = $rep_str eq '+' ? ($rep->[11],$rep->[12],$rep->[13]) : ($rep->[13],$rep->[12],$rep->[11]);
				# Color full length repeats in red
				#if ($con_pos[0] == 1 && $con_pos[2] eq '(0)' && $rep_div < 3) { @rgb = (255,0,0); }
				# Color repeat
				$img->bgcolor(@rgb);
				$img->fgcolor(@rgb);
				$page1->setrgbcolor($rgb[0]/255,$rgb[1]/255,$rgb[2]/255);
				# Draw repeat
			    if ($rep_str eq '+') {
			        $img->rectangle($l_space+$convx_a{$pid}{$rep_beg}, 10+70, $l_space+$convx_a{$pid}{$rep_end}, 24+70);
			        $page1->rectangle(40+$convx_a{$pid}{$rep_beg}, 1650-(20+30+40), ($convx_a{$pid}{$rep_end}-$convx_a{$pid}{$rep_beg}+1), 15);
			        $page1->fill();
			    } elsif ($rep_str eq 'C') {
			        $img->rectangle($l_space+$convx_a{$pid}{$rep_beg}, 26+70, $l_space+$convx_a{$pid}{$rep_end}, 40+70);
					$page1->rectangle(40+$convx_a{$pid}{$rep_beg}, 1650-(35+30+40), ($convx_a{$pid}{$rep_end}-$convx_a{$pid}{$rep_beg}+1), 15);
			        $page1->fill();
			    }
			}
		}
	}

	# piRNA expression
	my $max_rpm = 40;
	# Y-axis line
	$img->bgcolor('black');
	$img->fgcolor('black');
	$img->moveTo($l_space, 120);
	$img->lineTo($l_space, 160);
	# Y-axis labels
	$img->fontsize(10);
	$img->moveTo($l_space-15, 120+10);
	$img->string("$max_rpm");
	$img->moveTo($l_space-15, 161+2);
	$img->string("$max_rpm");
	# Y-axis line pdf
	$page1->setrgbcolorstroke(0.1,0.1,0.1);
	$page1->line(40, 1650-(20+30+40+50-20), 40, 1650-(20+30+40+50+20)-8);
	$page1->line(37, 1650-(20+30+40+50-20), 40, 1650-(20+30+40+50-20));
	$page1->line(37, 1650-(20+30+40+50+20)-8, 40, 1650-(20+30+40+50+20)-8);
	# Y-axis labels pdf
	$page1->setrgbcolor(0,0,0);
	$page1->stringc($font, 15, 40-15, 1650-(20+30+40+50-10), "$max_rpm");
	$page1->stringc($font, 15, 40-15, 1650-(20+30+40+50+25), "$max_rpm");
	# Go through each bin on plus strand
	foreach my $bin (sort {$a <=> $b} keys %{$pls_rds_s1}) {
		# Get pic coordinates
		my $bin_chr = $pls_rds_s1->{$bin}->[0];
		my $bin_beg = $pls_rds_s1->{$bin}->[1];
		my $bin_end = $pls_rds_s1->{$bin}->[2];
		my $bin_rds = $pls_rds_s1->{$bin}->[3];
		# Skip if not same contig
		unless ($pic_chr eq $bin_chr) { next; }
		# Check if reads bin is in region
		if ($pic_loc{$pid}[3] < $bin_beg && $pic_loc{$pid}[6] > $bin_end) {
			# Skip if in slice brake
			unless ($convx_a{$pid}{$bin_beg} && $convx_a{$pid}{$bin_end}) { next; }
			# Calculate rpm for bin
			my $bin_rpm = (int($bin_rds/$global_rds1*1_000_000))/($max_rpm/20);
			my $rpm_nrm = $bin_rpm <= 20 ? $bin_rpm : 20;
			# Draw rpm
			if ($rpm_nrm) {
				# Color rpm
				$img->bgcolor(255,148,0);
				$img->fgcolor(255,148,0);
				$img->rectangle($l_space+$convx_a{$pid}{$bin_beg}, 140-$rpm_nrm, $l_space+$convx_a{$pid}{$bin_end}, 140);
				$page1->setrgbcolor(255/255,148/255,0/255);
				$page1->rectangle(40+$convx_a{$pid}{$bin_beg}, 1650-(20+30+40+50), ($convx_a{$pid}{$bin_end}-$convx_a{$pid}{$bin_beg}+1), $rpm_nrm);
				$page1->fill();
			}
			if ($bin_rpm > 20) {
				# Color rpm
				$img->bgcolor(255,0,0);
				$img->fgcolor(255,0,0);
				$img->rectangle($l_space+$convx_a{$pid}{$bin_beg}, 120, $l_space+$convx_a{$pid}{$bin_end}, 120);
				$page1->setrgbcolor(255/255,0/255,0/255);
				$page1->rectangle(40+$convx_a{$pid}{$bin_beg}, 1650-(20+30+40+50)+20, ($convx_a{$pid}{$bin_end}-$convx_a{$pid}{$bin_beg}+1), 1);
				$page1->fill();
			}
		}
	}
	# Go through each bin on minus strand
	foreach my $bin (sort {$a <=> $b} keys %{$mns_rds_s1}) {
		# Get pic coordinates
		my $bin_chr = $mns_rds_s1->{$bin}->[0];
		my $bin_beg = $mns_rds_s1->{$bin}->[1];
		my $bin_end = $mns_rds_s1->{$bin}->[2];
		my $bin_rds = $mns_rds_s1->{$bin}->[3];
		# Skip if not same contig
		unless ($pic_chr eq $bin_chr) { next; }
		# Check if reads bin is in region
		if ($pic_loc{$pid}[3] < $bin_beg && $pic_loc{$pid}[6] > $bin_end) {
			# Skip if in slice brake
			unless ($convx_a{$pid}{$bin_beg} && $convx_a{$pid}{$bin_end}) { next; }
			# Calculate rpm for bin
			my $bin_rpm = (int($bin_rds/$global_rds1*1_000_000))/($max_rpm/20);
			my $rpm_nrm = $bin_rpm <= 20 ? $bin_rpm : 20;
			# Draw rpm
			if ($rpm_nrm) {
				# Color rpm
				$img->bgcolor(255,192,0);
				$img->fgcolor(255,192,0);
				$img->rectangle($l_space+$convx_a{$pid}{$bin_beg}, 141, $l_space+$convx_a{$pid}{$bin_end}, 141+$rpm_nrm);
				$page1->setrgbcolor(255/255,192/255,0/255);
				$page1->rectangle(40+$convx_a{$pid}{$bin_beg}, 1650-(20+30+40+50)-$rpm_nrm, ($convx_a{$pid}{$bin_end}-$convx_a{$pid}{$bin_beg}+1), $rpm_nrm);
				$page1->fill();
			}
			if ($bin_rpm > 20) {
				# Color rpm
				$img->bgcolor(255,0,0);
				$img->fgcolor(255,0,0);
				$img->rectangle($l_space+$convx_a{$pid}{$bin_beg}, 161, $l_space+$convx_a{$pid}{$bin_end}, 161);
				$page1->setrgbcolor(255/255,0/255,0/255);
				$page1->rectangle(40+$convx_a{$pid}{$bin_beg}, 1650-(20+30+40+50)-20, ($convx_a{$pid}{$bin_end}-$convx_a{$pid}{$bin_beg}+1), 1);
				$page1->fill();
			}
		}
	}

	## Species 2 syntenic contig regions
	# Number of contigs and their space in the image
	my $chr_i = 0;
	my $chr_n = keys %{$hom_reg{$pid}} <= 5 ? keys %{$hom_reg{$pid}} : 5;
	my $chr_shift = $chr_n == 5 ? $height/$chr_n : $height/($chr_n+1);
	# Output graphical representations for each syntenic contig
	foreach my $chr (sort {$hom_reg{$pid}{$b}[2] <=> $hom_reg{$pid}{$a}[2]} keys %{$hom_reg{$pid}}) {
		# Current contig number
		$chr_i++;
		if ($chr_i > 4) { last; }
		# Caption
	    $img->bgcolor('black');
	    $img->fgcolor('black');
	    $img->moveTo($l_space, 20+($chr_shift*$chr_i)-20);
	    $img->fontsize(10);
		my $sp2_beg = add_thousands_separators($hom_reg{$pid}{$chr}[0]);
		my $sp2_end = add_thousands_separators($hom_reg{$pid}{$chr}[1]);
	    $img->string("$chr:$sp2_beg-$sp2_end");
		# Caption pdf
		$page1->setrgbcolor(0/255,0/255,0/255);
	    $page1->string($font, 15, 40, 1650-20-($chr_shift*$chr_i), "$chr:$sp2_beg-$sp2_end");
		# Scale
		$scale_len = 100;
		$img->moveTo(900+$l_space, 20+($chr_shift*$chr_i)-20);
		$img->lineTo(900+$l_space+$scale_len, 20+($chr_shift*$chr_i)-20);
		$img->moveTo(900+$l_space+(($scale_len/2)-10), 20+($chr_shift*$chr_i)-20);
		$scale_bps = (($region_lens{$pid}{$chr}+1)/1000)*$scale_len;
		$scale_kbs = int($scale_bps/1000);
		$img->string("${scale_kbs}kb");
		# Scale pdf
		$page1->setrgbcolorstroke(0.1,0.1,0.1);
		$page1->line(900+40, 1650-20-($chr_shift*$chr_i), 900+40+$scale_len, 1650-20-($chr_shift*$chr_i));
		$page1->stringc($font, 15, 900+$scale_len-10, 1650-16-($chr_shift*$chr_i), "${scale_kbs}kb");
		# Contig region
	    $img->bgcolor(255,255,255);
	    $img->fgcolor('black');
	    $img->rectangle($l_space, 10+($chr_shift*$chr_i), 1000+$l_space, 40+($chr_shift*$chr_i));
		# Y-axis labels
		$img->fontsize(10);
		$img->moveTo($l_space-10, 10+($chr_shift*$chr_i)+10);
		$img->string("+");
		$img->moveTo($l_space-10, 40+($chr_shift*$chr_i)+2);
		$img->string("-");
		# Positions of other genes in region
		foreach my $gen (sort {$syn_gns{$pid}{$chr}{$a}[0] <=> $syn_gns{$pid}{$chr}{$b}[0]} keys %{$syn_gns{$pid}{$chr}}) {
			my $g_beg = $syn_gns{$pid}{$chr}{$gen}[0];
			my $g_end = $syn_gns{$pid}{$chr}{$gen}[1];
			my $g_str = $syn_gns{$pid}{$chr}{$gen}[2];
			unless ($convx_b{$pid}{$chr}{$g_beg} && $convx_b{$pid}{$chr}{$g_end}) { next; }
			# Check overlap with flank genes
			my $overlap = 0;
			foreach my $g (@{$ups_hcs{$pid}}) {
				if ($g_str eq $g->[3] && $g->[1] < $g_end && $g->[2] > $g_beg) {
					$overlap = 1;
				}
			}
			foreach my $g (@{$dos_hcs{$pid}}) {
				if ($g_str eq $g->[3] && $g->[1] < $g_end && $g->[2] > $g_beg) {
					$overlap = 1;
				}
			}
			# Skip if overlaps found
			if ($overlap) { next; }
			unless ($convx_b{$pid}{$chr}{$g_beg} && $convx_b{$pid}{$chr}{$g_end}) { next; }
			$img->bgcolor(210,210,210);
	        $img->fgcolor(210,210,210);
			if ($convx_b{$pid}{$chr}{$g_end}-$convx_b{$pid}{$chr}{$g_beg} <= 1) { $img->fgcolor(210,210,210); }
	        if ($g_str eq '+') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$g_beg}, 10+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$g_end}, 24+($chr_shift*$chr_i));
				$page1->setrgbcolor(240/255,240/255,240/255);
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$g_beg}, 1650-(20+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$g_end}-$convx_b{$pid}{$chr}{$g_beg}+1), 15);
	            $page1->fill();
	        } elsif ($g_str eq '-') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$g_beg}, 26+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$g_end}, 40+($chr_shift*$chr_i));
				$page1->setrgbcolor(240/255,240/255,240/255);
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$g_beg}, 1650-(35+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$g_end}-$convx_b{$pid}{$chr}{$g_beg}+1), 15);
	            $page1->fill();
	        }
	    }
		# Upstream orthologous flank genes in green
	    foreach my $gen (@{$ups_hcs{$pid}}) {
	        unless ($gen->[0] eq $chr) { next; }
			$img->bgcolor(@{$gen_cols{$pid}{$gen->[4]}});
	        $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}});
			$page1->setrgbcolor($gen_cols{$pid}{$gen->[4]}->[0]/255,$gen_cols{$pid}{$gen->[4]}->[1]/255,$gen_cols{$pid}{$gen->[4]}->[2]/255);
			if ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]} <= 1) { $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}}); }
	        if ($gen->[3] eq '+') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$gen->[1]}, 10+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$gen->[2]}, 24+($chr_shift*$chr_i));
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$gen->[1]}, 1650-(20+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]}+1), 15);
	            $page1->fill();
	        } elsif ($gen->[3] eq '-') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$gen->[1]}, 26+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$gen->[2]}, 40+($chr_shift*$chr_i));
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$gen->[1]}, 1650-(35+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]}+1), 15);
	            $page1->fill();
	        }
	    }
		# Downstream orthologous flank genes in blue
	    foreach my $gen (@{$dos_hcs{$pid}}) {
	        unless ($gen->[0] eq $chr) { next; }
			$img->bgcolor(@{$gen_cols{$pid}{$gen->[4]}});
	        $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}});
			$page1->setrgbcolor($gen_cols{$pid}{$gen->[4]}->[0]/255,$gen_cols{$pid}{$gen->[4]}->[1]/255,$gen_cols{$pid}{$gen->[4]}->[2]/255);
			if ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]} <= 1) { $img->fgcolor(@{$gen_cols{$pid}{$gen->[4]}}); }
	        if ($gen->[3] eq '+') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$gen->[1]}, 10+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$gen->[2]}, 24+($chr_shift*$chr_i));
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$gen->[1]}, 1650-(20+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]}+1), 15);
	            $page1->fill();
	        } elsif ($gen->[3] eq '-') {
	            $img->rectangle($l_space+$convx_b{$pid}{$chr}{$gen->[1]}, 26+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$gen->[2]}, 40+($chr_shift*$chr_i));
	            $page1->rectangle(40+$convx_b{$pid}{$chr}{$gen->[1]}, 1650-(35+30)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$gen->[2]}-$convx_b{$pid}{$chr}{$gen->[1]}+1), 15);
	            $page1->fill();
	        }
	    }
		# pdf contig rectangle
		$page1->setrgbcolorstroke(0.1,0.1,0.1);
	    $page1->rectangle(40,1650-65-($chr_shift*$chr_i),1000,30);
	    $page1->stroke();

		# Putative homologous pic locus
		my $hom_chr = $hom_pic_locs{$pid}[0];
		my $hom_beg = $hom_pic_locs{$pid}[1] eq 'NA' ? $hom_reg{$pid}{$chr}[0] : $hom_pic_locs{$pid}[1];
		my $hom_end = $hom_pic_locs{$pid}[2] eq 'NA' ? $hom_reg{$pid}{$chr}[1] : $hom_pic_locs{$pid}[2];
		if ($hom_chr eq $chr) {
			$img->bgcolor('gray');
	        $img->fgcolor('gray');
			# Horizontal line
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$hom_beg}, 20+($chr_shift*$chr_i)-15);
			$img->lineTo($l_space+$convx_b{$pid}{$chr}{$hom_end}, 20+($chr_shift*$chr_i)-15);
			# Vertical line
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$hom_beg}, 20+($chr_shift*$chr_i)-15);
			$img->lineTo($l_space+$convx_b{$pid}{$chr}{$hom_beg}, 20+($chr_shift*$chr_i)-13);
			# Vertical line
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$hom_end}, 20+($chr_shift*$chr_i)-15);
			$img->lineTo($l_space+$convx_b{$pid}{$chr}{$hom_end}, 20+($chr_shift*$chr_i)-13);
			# pdf
			$page1->setrgbcolorstroke(0.6,0.6,0.6);
			$page1->line(40+$convx_b{$pid}{$chr}{$hom_beg}, 1650-(20)-($chr_shift*$chr_i)-8, 40+$convx_b{$pid}{$chr}{$hom_end}, 1650-(20)-($chr_shift*$chr_i)-8);
			$page1->line(40+$convx_b{$pid}{$chr}{$hom_beg}, 1650-(20)-($chr_shift*$chr_i)-8, 40+$convx_b{$pid}{$chr}{$hom_beg}, 1650-(3+20)-($chr_shift*$chr_i)-8);
			$page1->line(40+$convx_b{$pid}{$chr}{$hom_end}, 1650-(20)-($chr_shift*$chr_i)-8, 40+$convx_b{$pid}{$chr}{$hom_end}, 1650-(3+20)-($chr_shift*$chr_i)-8);
		}

		# Breakpoints
		#foreach my $br_pair (@{$break_pos_sp2{$pid}{$chr}}) {
		#	foreach my $br (@{$br_pair}) {
		#		$img->bgcolor('red');
		#        $img->fgcolor('red');
		#        $img->rectangle(10+$convx_b{$pid}{$chr}{$br}, 10+($chr_shift*$chr_i)-5, 10+$convx_b{$pid}{$chr}{$br}, 24+($chr_shift*$chr_i)-15);
		#	}
		#}

		# Contig region
	    $img->bgcolor(undef);
	    $img->fgcolor('black');
	    $img->rectangle($l_space, 10+($chr_shift*$chr_i), 1000+$l_space, 40+($chr_shift*$chr_i));
		# Repeat annotation
		# Y-axis line
		$img->moveTo($l_space, 10+($chr_shift*$chr_i)+40);
		$img->lineTo($l_space, 40+($chr_shift*$chr_i)+40);
		# Y-axis labels
		$img->fontsize(10);
		$img->moveTo($l_space-10, 10+($chr_shift*$chr_i)+40+10);
		$img->string("+");
		$img->moveTo($l_space-10, 40+($chr_shift*$chr_i)+40+2);
		$img->string("-");
		if ($spec2_rmout_file) {
			# Repeat positions
			foreach my $rep (@{$rep_data_s2->{$chr}}) {
				# Repeat coordinates
				my $rep_beg = $rep->[5];
				my $rep_end = $rep->[6];
				my $rep_str = $rep->[8];
				my $rep_div = $rep->[1];
				# Check if repeat is in region
				if ($hom_reg{$pid}{$chr}[0] < $rep_beg && $hom_reg{$pid}{$chr}[1] > $rep_end) {
					unless ($convx_b{$pid}{$chr}{$rep_beg} && $convx_b{$pid}{$chr}{$rep_end}) { next; }
					# Calculate gray shade according to sequence divergence
					my @rgb = (255*($rep_div/40),255*($rep_div/40),255*($rep_div/40));
					# Get repeat position to consensus sequence
					my @con_pos = $rep_str eq '+' ? ($rep->[11],$rep->[12],$rep->[13]) : ($rep->[13],$rep->[12],$rep->[11]);
					# Color full length repeats in red
					#if ($con_pos[0] == 1 && $con_pos[2] eq '(0)' && $rep_div < 3) { @rgb = (255,0,0); }
					# Color repeat
					$img->bgcolor(@rgb);
				    $img->fgcolor(@rgb);
					$page1->setrgbcolor($rgb[0]/255,$rgb[1]/255,$rgb[2]/255);
					# Draw repeat
				    if ($rep_str eq '+') {
				        $img->rectangle($l_space+$convx_b{$pid}{$chr}{$rep_beg}, 10+($chr_shift*$chr_i)+40, $l_space+$convx_b{$pid}{$chr}{$rep_end}, 24+($chr_shift*$chr_i)+40);
			            $page1->rectangle(40+$convx_b{$pid}{$chr}{$rep_beg}, 1650-(20+30+40)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$rep_end}-$convx_b{$pid}{$chr}{$rep_beg}+1), 15);
			            $page1->fill();
				    } elsif ($rep_str eq 'C') {
				        $img->rectangle($l_space+$convx_b{$pid}{$chr}{$rep_beg}, 26+($chr_shift*$chr_i)+40, $l_space+$convx_b{$pid}{$chr}{$rep_end}, 40+($chr_shift*$chr_i)+40);
						$page1->rectangle(40+$convx_b{$pid}{$chr}{$rep_beg}, 1650-(35+30+40)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$rep_end}-$convx_b{$pid}{$chr}{$rep_beg}+1), 15);
			            $page1->fill();
				    }
				}
			}
		}

		# piRNA expression
		# Y-axis line
		$img->bgcolor('black');
		$img->fgcolor('black');
		$img->moveTo($l_space, 50+($chr_shift*$chr_i)+40);
		$img->lineTo($l_space, 50+($chr_shift*$chr_i)+80);
		# Y-axis labels
		$img->fontsize(10);
		$img->moveTo($l_space-15, 50+($chr_shift*$chr_i)+40+10);
		$img->string("$max_rpm");
		$img->moveTo($l_space-15, 50+($chr_shift*$chr_i)+81+2);
		$img->string("$max_rpm");
		# Y-axis line pdf
		$page1->setrgbcolorstroke(0.1,0.1,0.1);
		$page1->line(40, 1650-(20+30+40+50-20)-($chr_shift*$chr_i), 40, 1650-(20+30+40+50+20)-($chr_shift*$chr_i)-8);
		$page1->line(37, 1650-(20+30+40+50-20)-($chr_shift*$chr_i), 40, 1650-(20+30+40+50-20)-($chr_shift*$chr_i));
		$page1->line(37, 1650-(20+30+40+50+20)-($chr_shift*$chr_i)-8, 40, 1650-(20+30+40+50+20)-($chr_shift*$chr_i)-8);
		# Y-axis labels pdf
		$page1->setrgbcolor(0,0,0);
		$page1->stringc($font, 15, 40-15, 1650-(20+30+40+50-10)-($chr_shift*$chr_i), "$max_rpm");
		$page1->stringc($font, 15, 40-15, 1650-(20+30+40+50+25)-($chr_shift*$chr_i), "$max_rpm");
		# Go through each bin on plus strand
		foreach my $bin (sort {$a <=> $b} keys %{$pls_rds_s2}) {
			# Get pic coordinates
			my $bin_chr = $pls_rds_s2->{$bin}->[0];
			my $bin_beg = $pls_rds_s2->{$bin}->[1];
			my $bin_end = $pls_rds_s2->{$bin}->[2];
			my $bin_rds = $pls_rds_s2->{$bin}->[3];
			# Skip if not same contig
			unless ($chr eq $bin_chr) { next; }
			# Check if reads bin is in region
			if ($hom_reg{$pid}{$chr}[0] < $bin_beg && $hom_reg{$pid}{$chr}[1] > $bin_end) {
				# Skip if in slice brake
				unless ($convx_b{$pid}{$chr}{$bin_beg} && $convx_b{$pid}{$chr}{$bin_end}) { next; }
				# Calculate rpm for bin
				my $bin_rpm = (int($bin_rds/$global_rds2*1_000_000))/($max_rpm/20);
				my $rpm_nrm = $bin_rpm <= 20 ? $bin_rpm : 20;
				# Draw rpm
				if ($rpm_nrm) {
					# Color rpm
					$img->bgcolor(255,148,0);
					$img->fgcolor(255,148,0);
					$img->rectangle($l_space+$convx_b{$pid}{$chr}{$bin_beg}, 50+($chr_shift*$chr_i)+60-$rpm_nrm, $l_space+$convx_b{$pid}{$chr}{$bin_end}, 50+($chr_shift*$chr_i)+60);
					$page1->setrgbcolor(255/255,148/255,0/255);
					$page1->rectangle(40+$convx_b{$pid}{$chr}{$bin_beg}, 1650-(20+30+40+50)-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$bin_end}-$convx_b{$pid}{$chr}{$bin_beg}+1), $rpm_nrm);
					$page1->fill();
				}
				if ($bin_rpm > 20) {
					# Color rpm
					$img->bgcolor(255,0,0);
					$img->fgcolor(255,0,0);
					$img->rectangle($l_space+$convx_b{$pid}{$chr}{$bin_beg}, 50+($chr_shift*$chr_i)+40, $l_space+$convx_b{$pid}{$chr}{$bin_end}, 50+($chr_shift*$chr_i)+40);
					$page1->setrgbcolor(255/255,0/255,0/255);
					$page1->rectangle(40+$convx_b{$pid}{$chr}{$bin_beg}, 1650-(20+30+40+50)+20-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$bin_end}-$convx_b{$pid}{$chr}{$bin_beg}+1), 1);
					$page1->fill();
				}
			}
		}
		# Go through each bin on minus strand
		foreach my $bin (sort {$a <=> $b} keys %{$mns_rds_s2}) {
			# Get pic coordinates
			my $bin_chr = $mns_rds_s2->{$bin}->[0];
			my $bin_beg = $mns_rds_s2->{$bin}->[1];
			my $bin_end = $mns_rds_s2->{$bin}->[2];
			my $bin_rds = $mns_rds_s2->{$bin}->[3];
			# Skip if not same contig
			unless ($chr eq $bin_chr) { next; }
			# Check if reads bin is in region
			if ($hom_reg{$pid}{$chr}[0] < $bin_beg && $hom_reg{$pid}{$chr}[1] > $bin_end) {
				# Skip if in slice brake
				unless ($convx_b{$pid}{$chr}{$bin_beg} && $convx_b{$pid}{$chr}{$bin_end}) { next; }
				# Calculate rpm for bin
				my $bin_rpm = (int($bin_rds/$global_rds2*1_000_000))/($max_rpm/20);
				my $rpm_nrm = $bin_rpm <= 20 ? $bin_rpm : 20;
				# Draw rpm
				if ($rpm_nrm) {
					# Color rpm
					$img->bgcolor(255,192,0);
					$img->fgcolor(255,192,0);
					$img->rectangle($l_space+$convx_b{$pid}{$chr}{$bin_beg}, 50+($chr_shift*$chr_i)+61, $l_space+$convx_b{$pid}{$chr}{$bin_end}, 50+($chr_shift*$chr_i)+61+$rpm_nrm);
					$page1->setrgbcolor(255/255,192/255,0/255);
					$page1->rectangle(40+$convx_b{$pid}{$chr}{$bin_beg}, 1650-(20+30+40+50)-($chr_shift*$chr_i)-$rpm_nrm, ($convx_b{$pid}{$chr}{$bin_end}-$convx_b{$pid}{$chr}{$bin_beg}+1), $rpm_nrm);
					$page1->fill();
				}
				if ($bin_rpm > 20) {
					# Color rpm
					$img->bgcolor(255,0,0);
					$img->fgcolor(255,0,0);
					$img->rectangle($l_space+$convx_b{$pid}{$chr}{$bin_beg}, 50+($chr_shift*$chr_i)+81, $l_space+$convx_b{$pid}{$chr}{$bin_end}, 50+($chr_shift*$chr_i)+81);
					$page1->setrgbcolor(255/255,0/255,0/255);
					$page1->rectangle(40+$convx_b{$pid}{$chr}{$bin_beg}, 1650-(20+30+40+50)-20-($chr_shift*$chr_i), ($convx_b{$pid}{$chr}{$bin_end}-$convx_b{$pid}{$chr}{$bin_beg}+1), 1);
					$page1->fill();
				}
			}
		}

		# Insert positions of slice breaks
		foreach my $slice (sort {$a <=> $b} keys %{$slices{$pid}{$chr}}) {
			if ($slice == 1) { next; }
			# White color for slice point rectangle
			$img->bgcolor(255,255,255);
			$img->fgcolor(255,255,255);
			# Get slice point coordinate
			my $slice_pt = $slices{$pid}{$chr}{$slice}[0];
			my $sgap_len = int(($slices{$pid}{$chr}{$slice}[0]-$slices{$pid}{$chr}{$slice-1}[1]+500)/1000);
			# Draw white rectangle around slice point coordinate
			$img->rectangle($l_space+$convx_b{$pid}{$chr}{$slice_pt}-5, 5+($chr_shift*$chr_i), $l_space+$convx_b{$pid}{$chr}{$slice_pt}+5, 50+($chr_shift*$chr_i)+80);
			# Black slice point lines
			$img->bgcolor(0,0,0);
			$img->fgcolor(0,0,0);
			# Left slice point line
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$slice_pt}-5, 5+($chr_shift*$chr_i));
			$img->lineTo($l_space+$convx_b{$pid}{$chr}{$slice_pt}-5, 45+($chr_shift*$chr_i));
			# Right slice point line
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$slice_pt}+5, 5+($chr_shift*$chr_i));
			$img->lineTo($l_space+$convx_b{$pid}{$chr}{$slice_pt}+5, 45+($chr_shift*$chr_i));
			# Gap length label
			$img->fontsize(10);
			$img->moveTo($l_space+$convx_b{$pid}{$chr}{$slice_pt}-15, 5+($chr_shift*$chr_i));
			$img->string("${sgap_len}kb");

			# Draw white rectangle around slice point coordinate pdf
			$page1->setrgbcolor(1,1,1);
			$page1->setrgbcolorstroke(0,1,0);
			$page1->rectangle(40+$convx_b{$pid}{$chr}{$slice_pt}-5, 1650-69-($chr_shift*$chr_i), 10, 40);
			$page1->fill();
			# Black slice point lines pdf
			$page1->setrgbcolor(0,0,0);
			$page1->setrgbcolorstroke(0,0,0);
			$page1->line(40+$convx_b{$pid}{$chr}{$slice_pt}-5, 1650-29-($chr_shift*$chr_i), 40+$convx_b{$pid}{$chr}{$slice_pt}-5, 1650-29-($chr_shift*$chr_i)-41);
			$page1->line(40+$convx_b{$pid}{$chr}{$slice_pt}+5, 1650-29-($chr_shift*$chr_i), 40+$convx_b{$pid}{$chr}{$slice_pt}+5, 1650-29-($chr_shift*$chr_i)-41);
		}
	}
	# Convert into png data
	open my $out, '>', "$gd_dir/$img_prefix.pic_$pid.png" or die;
	binmode $out;
	print $out $img->png;
	# Close the file and write the PDF
	$pdf->close;
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

sub get_pic_flank_seqs {
	# Take pic locs
	my($pc_locs,$chr_seqs1_new,$flank_len,$min_flank_len) = @_;
	# Storage variable
	my %pic_flank_seqs = ();
	# Go through each pic
	foreach my $pic (sort {$a <=> $b} keys %{$pc_locs}) {
		# Get pic coordinates
		my $pic_chr = $pc_locs->{$pic}->[0];
		my $pic_beg = $pc_locs->{$pic}->[1];
		my $pic_end = $pc_locs->{$pic}->[2];
		# Get pic length
		my $pic_len = $pic_end-$pic_beg+1;
		# Get contig sequence and length
		my $chr_seq = $chr_seqs1_new->{$pic_chr};
		my $chr_len = length($chr_seq);
		# Get upstream sequence of pic
		my $ups_beg = $pic_beg > $flank_len ? $pic_beg-$flank_len : 0;
		my $ups_len = $pic_beg > $flank_len ? $flank_len : $pic_beg;
		my $ups_end = $ups_beg+$ups_len;
		my $ups_seq = substr($chr_seq,$ups_beg-1,$ups_len);
		# Get downstream sequence of pic
		my $dos_beg = $pic_end;
		my $dos_len = $pic_end+$flank_len < $chr_len ? $flank_len : $chr_len-$pic_end;
		my $dos_end = $dos_beg+$dos_len;
		my $dos_seq = substr($chr_seq,$dos_beg-1,$dos_len);
		# Save pic flank sequences if larger than minimum length
		#if ($ups_len >= $min_flank_len) {
			@{$pic_flank_seqs{$pic}[0]} = ($ups_seq,$pic_chr,$ups_beg,$ups_end);
		#} else {
			#@{$pic_flank_seqs{$pic}[0]} = ('','',0,0,);
		#}
		#if ($dos_len >= $min_flank_len) {
			@{$pic_flank_seqs{$pic}[1]} = ($dos_seq,$pic_chr,$dos_beg,$dos_end);
		#} else {
			#@{$pic_flank_seqs{$pic}[1]} = ('','',0,0,);
		#}
	}
	# Return pic flank sequences
	return \%pic_flank_seqs;
}

sub get_flank_locs_old {
	# Take name of tab file
	my($pic_flank_seqs,$genome1_old_file,$spec1_id,$flank_dir) = @_;
	# Storage variable
	my %flank_locs_old = ();
	# Go through each pic
	foreach my $pic (sort {$a <=> $b} keys %{$pic_flank_seqs}) {
		# Go through each flank
		foreach my $i (0..1) {
			# If flank sequence present
			if ($pic_flank_seqs->{$pic}->[$i]->[0]) {
				# Flank name
				my $flank = $i == 0 ? 'ups' : 'dos';
				# Print upstream and downstream sequences to files
				my $outfile = $flank_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'.fas';
				unless (-e $outfile) {
					my $out = open_outfile($outfile);
					print($out ">${pic}_$flank\n$pic_flank_seqs->{$pic}->[$i]->[0]\n");
					close($out);
				}
				# Call blast
				my $blast_out = $flank_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'.blast';
				unless (-e $blast_out) {
					system("blastn -query $outfile -subject $genome1_old_file -out $blast_out -outfmt \'$outfmt\' >>.log 2>&1");
				}
				# Extract filtered blast hits
				my $blast_hits = get_blastn_hits($blast_out);
				# Group blast hits
				my $chr_hit_locs = get_hit_group_locs($blast_hits);
				# Group homologous contig regions
				my @hom_chrs = ();
				foreach my $chr (sort {$chr_hit_locs->{$b}->[5] <=> $chr_hit_locs->{$a}->[5]} keys %{$chr_hit_locs}) {
					# Get query coverage
					my $qcov = $chr_hit_locs->{$chr}->[5];
					# Subject contig has query coverage >80
					if ($qcov > 75) {
						push(@hom_chrs,$chr);
						last;
					} elsif ($qcov > 45 && $qcov <= 75) {
						# If split over two contigs
						push(@hom_chrs,$chr);
					}
				}
				# Choose non-overlapping contigs
				my @chrs = delete_overlapping_chrs(\@hom_chrs,$chr_hit_locs);
				# Filter corresponding chromosome regions
				foreach my $chr (@chrs) {
					@{$flank_locs_old{$pic}[$i]{$chr}} = @{$chr_hit_locs->{$chr}};
				}
			}
		}
	}
	# Return pic flank locs of old assembly
	return \%flank_locs_old;
}

sub get_blastn_hits {
	# Take infile name
	my($file) = @_;
	# Get file data
	my @file_data = get_file_data_array($file);
	# Initialize variables
	my %blast_hits = ();
	# Parse tblastn output
	foreach my $line (@file_data) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr_id = $d[1];
			# Filter hits
			if ($d[2] < 95) { next; }
			if ($d[4] < 20) { next; }
			if ($d[12] < 2000) { next; }
			# Save hits
			push(@{$blast_hits{$chr_id}},\@d);
		}
	}
	return \%blast_hits;
}

sub get_hit_group_locs {
	# Take gene blast hits
	my($blast_hits) = @_;
	# Storage variable
	my %chr_hit_locs = ();
	# Go through each contig
	foreach my $chr (sort keys %{$blast_hits}) {
		# Initialize variables
		my $group = 0;
		my $prev_hit_end = 0;
		my $max_distance = 10_000;
		my %grouped_hits = ();
		my %group_score = ();
		my %group_locs = ();
		my %group_locq = ();
		my %group_ident = ();
		my %group_qcov = ();
		my %group_match = ();
		my %group_mis_m = ();
		my %hit_strands = ();
		my %group_strand = ();
		# Go through each hit
		foreach my $hit (sort {$a->[9] <=> $b->[9]} @{$blast_hits->{$chr}}) {
			# Hit positions
			my $identity = $hit->[2];
			my $querycov = $hit->[4];
			my $qry_beg  = $hit->[6];
			my $qry_end  = $hit->[7];
			my $hit_beg  = $hit->[9];
			my $hit_end  = $hit->[10];
			my $sstrand  = $hit->[11];
			my $hit_len  = $hit->[12];
			my $bitscore = $hit->[13];
			# Adjust symbols and coordinates according to subject strand
			if ($hit_beg > $hit_end) {
				$hit_beg = $hit->[10];
				$hit_end = $hit->[9];
			}
			# Save strand
			$hit_strands{$sstrand}++;
			# Get distance
			my $distance = $hit_beg-$prev_hit_end;
			# If distance above threshold open new group
			if ($distance > $max_distance) { $group++ }
			# Allocate hit to group
			push(@{$grouped_hits{$group}},$hit);
			# Save group properties
			$group_score{$group} += $bitscore;
			$group_qcov{$group} = $querycov;
			$group_locs{$group}[0] = $hit_beg unless $group_locs{$group}[0];
			$group_locs{$group}[0] = $hit_beg if $hit_beg < $group_locs{$group}[0];
			$group_locs{$group}[1] = $hit_end unless $group_locs{$group}[1];
			$group_locs{$group}[1] = $hit_end if $hit_end > $group_locs{$group}[1];
			$group_locq{$group}[0] = $qry_beg unless $group_locq{$group}[0];
			$group_locq{$group}[0] = $qry_beg if $qry_beg < $group_locq{$group}[0];
			$group_locq{$group}[1] = $qry_end unless $group_locq{$group}[1];
			$group_locq{$group}[1] = $qry_end if $qry_end > $group_locq{$group}[1];
			# Get matches and mismatches
			my $n_match = int(($hit_len*$identity/100)+0.5);
			my $n_mis_m = $hit_len-$n_match;
			$group_match{$group} += $n_match;
			$group_mis_m{$group} += $n_mis_m;
			# Save end position for next hit
			$prev_hit_end = $hit_end;
		}
		# Get group strand
		foreach my $gr (keys %group_match) {
			foreach my $str (sort {$hit_strands{$b} <=> $hit_strands{$a}} keys %hit_strands) {
				$group_strand{$gr} = $str and last;
			}
		}
		# Calculate identities
		foreach my $gr (keys %group_match) {
			$group_ident{$gr} = $group_match{$gr}/($group_match{$gr}+$group_mis_m{$gr})*100;
		}
		# Choose best group
		my $best_gr = 0;
		foreach my $gr (sort {$group_score{$b} <=> $group_score{$a}} keys %group_score) {
			$best_gr = $gr and last;
		}
		@{$chr_hit_locs{$chr}} = (
			@{$group_locq{$best_gr}},
			@{$group_locs{$best_gr}},
			$group_ident{$best_gr},
			$group_qcov{$best_gr},
			$group_strand{$best_gr}
		);
	}
	return \%chr_hit_locs;
}

sub delete_overlapping_chrs {
	# Take array ref of contigs
	my($all_chrs,$hit_locs) = @_;
	# Array for deletion
	my @del_chrs = ();
	# Final array of contigs
	my @fin_chrs = ();
	# If more than one contig
	if (scalar(@$all_chrs) > 1) {
		my $prev_chr = '';
		foreach my $chr (sort {$hit_locs->{$a}->[0] <=> $hit_locs->{$b}->[0]} @$all_chrs) {
			unless ($prev_chr) {
				$prev_chr = $chr;
				next;
			}
			# Get properties
			my $qbeg = $hit_locs->{$chr}->[0];
			my $qend = $hit_locs->{$chr}->[1];
			my $qlen = $qend-$qbeg+1;
			my $prev_qbeg = $hit_locs->{$prev_chr}->[0];
			my $prev_qend = $hit_locs->{$prev_chr}->[1];
			my $prev_qlen = $prev_qend-$prev_qbeg+1;
			# Check if regions overlap
			if ($prev_qend > $qbeg) {
				# Get overlap length
				my $overlap_len = $prev_qend-$qbeg+1;
				# If substential overlap (>25%), discard shorter region
				if ($overlap_len > ($qlen/4) || $overlap_len > ($prev_qlen/4)) {
					if ($prev_qlen > $qlen) {
						push(@del_chrs,$chr);
					} else {
						push(@del_chrs,$prev_chr);
					}
				}
			}
			# Save contig name
			$prev_chr = $chr;
		}
	}
	# Delete bad contigs from list
	foreach my $chr (sort {$hit_locs->{$a}->[0] <=> $hit_locs->{$b}->[0]} @$all_chrs) {
		unless (grep {$_ eq $chr} @del_chrs) {
			push(@fin_chrs,$chr);
		}
	}
	return @fin_chrs;
}

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = open_infile($gff_file);
	# Initialize variables
	my %gff_gen_locs = ();
	my $gene = '';
	my $name = '';
	my $tran = '';
	# Parse gff file
	while (my $line = <$gff>) {
		# Exit loop if fasta part begins
		if ($line =~ /^>/) { last; }
		# Results line
		if ($line =~ /^[^#]+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $beg  = $d[3];
			my $end  = $d[4];
			my $info = $d[8];
			# Gene line
			if ($type =~ /gene/i) {
				# Get gene id
				if ($info =~ /ID=[^;:]+;/) {
					($gene) = ($info =~ /ID=([^;:]+);/);
				}
				# Save entries
				$gff_gen_locs{$chr}{$gene} = \@d;
			}
		}
	}
	return \%gff_gen_locs;
}

sub get_gene_locs_new {
	# Take name of tab file
	my($pic_locs_s1_new,$flank_seqs_s1_new,$flank_locs_s1_old,$gen_locs_s1_old,$genome_s1_old_file) = @_;
	# Storage variable
	my %gen_locs_s1_new = ();
	# Extract genomic sequences from old assembly (spec2)
	my $chr_seqs_s1_old = get_fasta_seqs($genome_s1_old_file,1);
	# Go through each pic
	foreach my $pic (sort {$a <=> $b} keys %{$pic_locs_s1_new}) {
		# Get pic coordinates
		my $pic_chr = $pic_locs_s1_new->{$pic}->[0];
		my $pic_beg = $pic_locs_s1_new->{$pic}->[1];
		my $pic_end = $pic_locs_s1_new->{$pic}->[2];
		# Go through each flank
		foreach my $i (0..1) {
			# Flank name
			my $flank = $i == 0 ? 'ups' : 'dos';
			# Print flank sequence to file
			my $foutfile = $flank_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'.fas';
			unless (-e $foutfile) {
				my $fout = open_outfile($foutfile);
				print($fout ">${pic}_$flank\n$flank_seqs_s1_new->{$pic}->[0]\n");
				close($fout);
			}
			# Open output file for ortholog gene sequences
			my $gseqfile = $flank_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'.genes.fas';
			my $gout = open_outfile($gseqfile);
			# Get flank sequence
			my $flank_seq = $flank_seqs_s1_new->{$pic}->[$i]->[0];
			my $flank_lng = length($flank_seq);
			# Go through each flank contig in old assembly corresponding to pic flank of new assembly
			foreach my $chr_os1 (sort keys %{$flank_locs_s1_old->{$pic}->[$i]}) {
				# Go through each gene
				foreach my $gene (sort keys %{$gen_locs_s1_old->{$chr_os1}}) {
					# Get gene coordinates
					my $gen_beg = $gen_locs_s1_old->{$chr_os1}->{$gene}->[3];
					my $gen_end = $gen_locs_s1_old->{$chr_os1}->{$gene}->[4];
					my $gen_str = $gen_locs_s1_old->{$chr_os1}->{$gene}->[6];
					my $gen_len = $gen_end-$gen_beg+1;
					# Get gene sequence from old assembly sequence
					my $gen_seq = substr($chr_seqs_s1_old->{$chr_os1},$gen_beg-1,$gen_len);
					# RevCom gene sequence if on minus strand
					$gen_seq = rev_com($gen_seq) if $gen_str eq '-';
					# Print sequence to output file
					print($gout ">$gene\n$gen_seq\n");
				}
			}
			# Blast gene sequence to genome of species 2
			my $blast_gout = $flank_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'.genes.blast';
			unless (-e $blast_gout) {
				system("blastn -query $gseqfile -subject $foutfile -out $blast_gout -outfmt \'$outfmt\' >>.log 2>&1");
			}
			unless (-e $blast_gout) { system("touch $blast_gout"); }
			# Extract filtered blast hits
			my $gene_hits = get_blastn_gene_hits($blast_gout,90);
			# Group blast hits
			my $gene_hit_locs = get_gene_hit_group_locs($gene_hits);
			# Go through each flank contig of old assembly
			foreach my $chr_os1 (sort keys %{$flank_locs_s1_old->{$pic}->[$i]}) {
				# Print ortholog locs in new assembly
				foreach my $flank_ns1 (sort keys %{$gene_hit_locs}) {
					foreach my $gen_os1 (sort keys %{$gene_hit_locs->{$flank_ns1}}) {
						# Get ortholog coordinates in new assembly
						if ($gen_locs_s1_new{$pic}[$i]{$pic_chr}{$gen_os1}) { next; }
						my $gen_ns1_beg = $gene_hit_locs->{$flank_ns1}->{$gen_os1}->[0]+$pic_beg-$flank_lng;
						my $gen_ns1_end = $gene_hit_locs->{$flank_ns1}->{$gen_os1}->[1]+$pic_beg-$flank_lng;
						if ($i == 1) {
							$gen_ns1_beg = $gene_hit_locs->{$flank_ns1}->{$gen_os1}->[0]+$pic_end;
							$gen_ns1_end = $gene_hit_locs->{$flank_ns1}->{$gen_os1}->[1]+$pic_end;
						}
						my $gen_ns1_str = $gene_hit_locs->{$flank_ns1}->{$gen_os1}->[3];
						# Save coordinates of flank gene orthologs
						@{$gen_locs_s1_new{$pic}[$i]{$pic_chr}{$gen_os1}} = ($gen_ns1_beg,$gen_ns1_end,$gen_ns1_str);
					}
				}
			}
		}
	}
	return \%gen_locs_s1_new;
}

sub get_blastn_gene_hits {
	# Take infile name
	my($file,$cov_th) = @_;
	# Get file data
	my @file_data = get_file_data_array($file);
	# Initialize variables
	my %blast_hits = ();
	# Parse tblastn output
	foreach my $line (@file_data) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $gen_id = $d[0];
			my $chr_id = $d[1];
			my $coverg = $d[4];
			# Save hits
			if ($coverg < $cov_th) { next }
			push(@{$blast_hits{$chr_id}{$gen_id}},\@d);
		}
	}
	# Filter repetitive hits
	foreach my $chr_id (keys %blast_hits) {
		foreach my $gen_id (keys %{$blast_hits{$chr_id}}) {
			$blast_hits{$chr_id}{$gen_id} = filter_repetitive_hits($blast_hits{$chr_id}{$gen_id});
			unless (@{$blast_hits{$chr_id}{$gen_id}}) { delete $blast_hits{$chr_id}{$gen_id} }
		}
	}
	# Return blast hits
	return \%blast_hits;
}

sub get_gene_hit_group_locs {
	# Take gene blast hits
	my($gen_hits) = @_;
	# Storage variable
	my %gen_hit_locs = ();
	# Go through each subject hit chromosome
	foreach my $chr (sort keys %{$gen_hits}) {
		# Go through each gene
		foreach my $gen (sort keys %{$gen_hits->{$chr}}) {
			# Group gene hits
			my $group = 0;
			my $prev_hit_end = 0;
			my $max_distance = 10_000;
			my %grouped_hits = ();
			my %group_score = ();
			my %group_locs = ();
			my %group_ident = ();
			my %group_match = ();
			my %group_mis_m = ();
			my %hit_strands = ();
			my %group_strand = ();
			# Go through each hit
			foreach my $hit (sort {$a->[9] <=> $b->[9]} @{$gen_hits->{$chr}->{$gen}}) {
				# Hit positions
				my $identity = $hit->[2];
				my $hit_beg  = $hit->[9];
				my $hit_end  = $hit->[10];
				my $sstrand  = $hit->[11];
				my $hit_len  = $hit->[12];
				my $bitscore = $hit->[13];
				if ($bitscore =~ /e/) {
					my($b) = ($bitscore =~ /^(.*)e/);
					my($p) = ($bitscore =~ /e\+(.*)/);
					$bitscore = $b*(10**$p);
				}
				# Adjust symbols and coordinates according to subject strand
				if ($hit_beg > $hit_end) {
					$hit_beg = $hit->[10];
					$hit_end = $hit->[9];
				}
				# Save strand
				$hit_strands{$sstrand} += $bitscore;
				# Get distance
				my $distance = $hit_beg-$prev_hit_end;
				# If distance above threshold open new group
				if ($prev_hit_end > $hit_end) { next; }
				if ($distance > $max_distance) { $group++ }
				# Allocate hit to group
				push(@{$grouped_hits{$group}},$hit);
				# Save group properties
				$group_score{$group} += $bitscore;
				$group_locs{$group}[0] = $hit_beg unless $group_locs{$group}[0];
				$group_locs{$group}[0] = $hit_beg if $hit_beg < $group_locs{$group}[0];
				$group_locs{$group}[1] = $hit_end unless $group_locs{$group}[1];
				$group_locs{$group}[1] = $hit_end if $hit_end > $group_locs{$group}[1];
				# Get matches and mismatches
				my $n_match = int(($hit_len*$identity/100)+0.5);
				my $n_mis_m = $hit_len-$n_match;
				$group_match{$group} += $n_match;
				$group_mis_m{$group} += $n_mis_m;
				# Save end position for next hit
				$prev_hit_end = $hit_end;
			}
			# Get group strand
			foreach my $gr (keys %group_match) {
				foreach my $str (sort {$hit_strands{$b} <=> $hit_strands{$a}} keys %hit_strands) {
					$group_strand{$gr} = $str and last;
				}
			}
			# Calculate identities
			foreach my $gr (keys %group_match) {
				$group_ident{$gr} = $group_match{$gr}/($group_match{$gr}+$group_mis_m{$gr})*100;
			}
			# Choose best group
			my $best_gr = 0;
			foreach my $gr (sort {$group_score{$b} <=> $group_score{$a}} keys %group_score) {
				$best_gr = $gr and last;
			}
			@{$gen_hit_locs{$chr}{$gen}} = (@{$group_locs{$best_gr}},$group_ident{$best_gr},$group_strand{$best_gr});
		}
	}
	return \%gen_hit_locs;
}

sub get_orthologs {
	# Take name of tab file
	my($dmel_orthologs_file) = @_;
	# Get Dmel orthologs
	my $dmel_orthos = get_dmel_orthologs($dmel_orthologs_file);
	# Save orthologs for each species
	my %orthologs = ();
	# Go through each Dmel gene id
	foreach my $dmel_id (sort keys %{$dmel_orthos}) {
		# Go through each ortholog
		foreach my $ortho (@{$dmel_orthos->{$dmel_id}}) {
			# Get ortholog id and species id
			my $ortho_id = $ortho->[5];
			my($ortho_sp) = ($ortho->[6] =~ /([^\\]+)\\/);
			# Go through each ortholog
			foreach my $ortho_b (@{$dmel_orthos->{$dmel_id}}) {
				# Get ortholog id and species id
				my $ortho_id_b = $ortho_b->[5];
				my($ortho_sp_b) = ($ortho_b->[6] =~ /([^\\]+)\\/);
				# Skip same ortholog or species
				if ($ortho eq $ortho_b) { next; }
				if ($ortho_sp eq $ortho_sp_b) { next; }
				# Save orthologs
				unless (grep {$_ eq $ortho_id_b} @{$orthologs{$ortho_sp}{$ortho_id}{$ortho_sp_b}}) {
					push(@{$orthologs{$ortho_sp}{$ortho_id}{$ortho_sp_b}},$ortho_id_b);
				}
			}
			# Save Dmel ortholog
			unless ($ortho_sp eq 'Dmel') {
				unless (grep {$_ eq $dmel_id} @{$orthologs{$ortho_sp}{$ortho_id}{'Dmel'}}) {
					push(@{$orthologs{$ortho_sp}{$ortho_id}{'Dmel'}},$dmel_id);
				}
			}
		}
		# Go through each ortholog and save Dmel orthologs
		foreach my $ortho_b (@{$dmel_orthos->{$dmel_id}}) {
			# Get ortholog id and species id
			my $ortho_id_b = $ortho_b->[5];
			my($ortho_sp_b) = ($ortho_b->[6] =~ /([^\\]+)\\/);
			# Save orthologs
			unless (grep {$_ eq $ortho_id_b} @{$orthologs{'Dmel'}{$dmel_id}{$ortho_sp_b}}) {
				push(@{$orthologs{'Dmel'}{$dmel_id}{$ortho_sp_b}},$ortho_id_b);
			}
		}
	}
	return \%orthologs;
}

sub get_dmel_orthologs {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
	my @in_data = get_file_data_array($infile);
	# Global tree hashes for each species
	my %dmel_orthos = ();
	# Go through file data
	foreach my $line (@in_data) {
		if ($line =~ /^#/) { next }
		# Get line data
		my @d = split(/\t/,$line);
		my $dmel_id = $d[0];
		# Save all entries linked to this dmel gn_id
		push(@{$dmel_orthos{$dmel_id}},\@d);
	}
	return \%dmel_orthos;
}

sub get_ortholog_locs {
	# Take name of tab file
	my($dmel_orthologs_file) = @_;
	# Get file data
	my @file_data = get_file_data_array($dmel_orthologs_file);
	# Storage variable
	my %ortho_locs = ();
	# Go through file data
	foreach my $line (@file_data) {
		if ($line =~ /^#/) { next; }
        # Get line data
        my @d = split(/\t/,$line);
		# Dmel genes
		my $m_gid = $d[0];
		my $m_chr = $d[2];
		my $m_loc = $d[3];
		my($m_beg) = ($m_loc =~ /(\d+)\.\./);
		my($m_end) = ($m_loc =~ /\.\.(\d+)/);
		my $m_str = $d[4] eq '1' ? '+' : '-';
		# Save coordinates
		@{$ortho_locs{$m_gid}} = ($m_chr,$m_beg,$m_end,$m_str);
		# Orthologs
		unless ($d[5]) { next; }
		my $o_gid = $d[5];
		my $o_chr = $d[7];
		my $o_loc = $d[8];
		my($o_beg) = ($o_loc =~ /(\d+)\.\./);
		my($o_end) = ($o_loc =~ /\.\.(\d+)/);
		my $o_str = $d[9] eq '1' ? '+' : '-';
		# Save coordinates
		@{$ortho_locs{$o_gid}} = ($o_chr,$o_beg,$o_end,$o_str);
	}
	return \%ortho_locs;
}

sub get_flank_ortho_seqs {
	# Take name of tab file
	my($flank_locs_s1_old,$gen_locs_s1_old,$orthologs,$ortho_locs,$genome_s2_old_file,$spec1_id,$ortho_dir) = @_;
	# Storage variable
	my %flank_orthologs = ();
	# Extract genomic sequences from old assembly (spec2)
	my $chr_seqs_s2_old = get_fasta_seqs($genome_s2_old_file,1);
	# Go through each pic
	foreach my $pic (sort {$a <=> $b} keys %{$flank_locs_s1_old}) {
		# Go through each flank
		foreach my $i (0..1) {
			# Flank name
			my $flank = $i == 0 ? 'ups' : 'dos';
			# Sort flank contig regions of old assembly by position on new assembly
			my @chrs_os1 = sort {$flank_locs_s1_old->{$pic}->[$i]->{$a}->[0] <=> $flank_locs_s1_old->{$pic}->[$i]->{$b}->[0]} keys %{$flank_locs_s1_old->{$pic}->[$i]};
			# Go through each flank contig of old assembly
			foreach my $chr (@chrs_os1) {
				# Get direction of old assembly contig relative to new assembly contig
				my $chr_str = $flank_locs_s1_old->{$pic}->[$i]->{$chr}->[6];
				# Sort genes by position on contig depending on direction of old assembly contig relative to new assembly contig
				my @genes_os1 = ();
				if ($chr_str eq 'plus') {
					@genes_os1 = sort {$gen_locs_s1_old->{$chr}->{$a}->[3] <=> $gen_locs_s1_old->{$chr}->{$b}->[3]} keys %{$gen_locs_s1_old->{$chr}};
				} elsif ($chr_str eq 'minus') {
					@genes_os1 = sort {$gen_locs_s1_old->{$chr}->{$b}->[4] <=> $gen_locs_s1_old->{$chr}->{$a}->[4]} keys %{$gen_locs_s1_old->{$chr}};
				}
				my %ort_chr_prev = ();
				# Go through each gene on same contig and get flanking genes
				foreach my $gene (@genes_os1) {
					# Get gene coordinates
					my $gen_beg = $gen_locs_s1_old->{$chr}->{$gene}->[3];
					my $gen_end = $gen_locs_s1_old->{$chr}->{$gene}->[4];
					my $gen_str = $gen_locs_s1_old->{$chr}->{$gene}->[6];
					# Save upstream flanking genes
					if ($gen_beg >= $flank_locs_s1_old->{$pic}->[$i]->{$chr}->[2] && $gen_end <= $flank_locs_s1_old->{$pic}->[$i]->{$chr}->[3]) {
						# Go through each species for which orthologs exist
						foreach my $ortho_sp_b (sort keys %{$orthologs->{$spec1_id}->{$gene}}) {
							# Skip species that are not included in analysis
							unless ($ortho_sp_b eq $spec2_id) { next; }
							# Print orthologs
							my $ort_seq = '';
							my $ort_gid = '';
							my $ort_chr = '';
							my $ort_beg = 0;
							my $ort_end = 0;
							my $ort_str = '';
							# Initialize best fit ortho_id as first available ortholog
							my $ortho_id_b_bfit = '';
							foreach my $ortho_id_b (@{$orthologs->{$spec1_id}->{$gene}->{$ortho_sp_b}}) {
								$ortho_id_b_bfit = $ortho_id_b;
								last;
							}
							# See if other orthologs are on the same contig as the previous ortholog
							foreach my $ortho_id_b (@{$orthologs->{$spec1_id}->{$gene}->{$ortho_sp_b}}) {
								my $chr = $ortho_locs->{$ortho_id_b}->[0];
								if ($ort_chr_prev{$ortho_sp_b} && $chr eq $ort_chr_prev{$ortho_sp_b}) {
									$ortho_id_b_bfit = $ortho_id_b;
									last;
								}
							}
							# Get coordinates and sequence of best fitting ortholog
							$ort_chr = $ortho_locs->{$ortho_id_b_bfit}->[0];
							$ort_beg = $ortho_locs->{$ortho_id_b_bfit}->[1];
							$ort_end = $ortho_locs->{$ortho_id_b_bfit}->[2];
							$ort_str = $ortho_locs->{$ortho_id_b_bfit}->[3];
							my $ort_len = $ort_end-$ort_beg+1;
							# Get ortholog sequence
							$ort_seq = substr($chr_seqs_s2_old->{$ort_chr},$ort_beg-1,$ort_len) unless $ort_seq;
							# RevCom ortholog sequence if on minus strand
							$ort_seq = rev_com($ort_seq) if $ort_str eq '-';
							# Get gene id
							$ort_gid = $ortho_id_b_bfit unless $ort_gid;
							# Save chromosome location of ortholog
							$ort_chr_prev{$ortho_sp_b} = $ort_chr;
							# Save ortholog sequence
							@{$flank_orthologs{$pic}[$i]{$chr}{$gene}} = ($ort_gid,$ort_seq,$gen_beg,$gen_end,$gen_str,$ort_chr,$ort_beg,$ort_end,$ort_str);
						}
					}
				}
			}
		}
	}
	# Return pic flank locs of old assembly
	return \%flank_orthologs;
}

sub get_flank_orthos_new {
	# Take pic locs
	my($flank_ortho_seqs,$flank_locs_s1_old,$spec1_id,$spec2_id,$ortho_dir) = @_;
	# Storage variable
	my %flank_ortho_locs_s2_new = ();
	# Go through each pic
	foreach my $pic (sort {$a <=> $b} keys %{$flank_ortho_seqs}) {
		# Go through each flank
		foreach my $i (0..1) {
			# Flank name
			my $flank = $i == 0 ? 'ups' : 'dos';
			# Open output file for ortholog gene sequences
			my $gseqfile = $ortho_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'_'.$spec2_id.'.orthologs.fas';
			my $gout = open_outfile($gseqfile);
			# Sort flank contig regions of old assembly by position on new assembly
			my @chrs_os1 = sort {$flank_locs_s1_old->{$pic}->[$i]->{$a}->[0] <=> $flank_locs_s1_old->{$pic}->[$i]->{$b}->[0]} keys %{$flank_locs_s1_old->{$pic}->[$i]};
			# Go through each flank contig of old assembly
			foreach my $chr (@chrs_os1) {
				# Get direction of old assembly contig relative to new assembly contig
				my $chr_str = $flank_locs_s1_old->{$pic}->[$i]->{$chr}->[6];
				# Sort genes by position on contig depending on direction of old assembly contig relative to new assembly contig
				my @genes_os1 = ();
				if ($chr_str eq 'plus') {
					@genes_os1 = sort {$flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$a}->[2] <=> $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$b}->[2]} keys %{$flank_ortho_seqs->{$pic}->[$i]->{$chr}};
				} elsif ($chr_str eq 'minus') {
					@genes_os1 = sort {$flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$b}->[3] <=> $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$a}->[3]} keys %{$flank_ortho_seqs->{$pic}->[$i]->{$chr}};
				}
				# Go through each flanking gene on contig and search ortholog sequences in new assembly
				foreach my $gene (@genes_os1) {
					# Get ortholog id and sequence
					my $ort_gid = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[0];
					my $ort_seq = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gene}->[1];
					# Print sequence to output file
					print($gout ">$gene\t$ort_gid\n$ort_seq\n");
				}
			}
			# Blast gene sequence to genome of species 2
			my $blast_gout = $ortho_dir.'/'.$spec1_id.'.'.$pic.'_'.$flank.'_'.$spec2_id.'.orthologs.blast';
			unless (-e $blast_gout) {
				system("blastn -query $gseqfile -subject $genome_s2_new_file -out $blast_gout -outfmt \'$outfmt\' >>.log 2>&1");
			}
			# Extract filtered blast hits
			my $gene_hits = get_blastn_gene_hits($blast_gout,90);
			# Group blast hits
			my $gene_hit_locs = get_gene_hit_group_locs($gene_hits);
			# Go through each flank contig of old assembly
			foreach my $chr (@chrs_os1) {
				# Print ortholog locs in new assembly
				foreach my $chr_ns2 (sort keys %{$gene_hit_locs}) {
					foreach my $gen_os1 (sort keys %{$gene_hit_locs->{$chr_ns2}}) {
						# Get ortholog coordinates in new assembly
						unless ($flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gen_os1}) { next; }
						my $ort_os2_gid = $flank_ortho_seqs->{$pic}->[$i]->{$chr}->{$gen_os1}->[0];
						my $ort_ns2_beg = $gene_hit_locs->{$chr_ns2}->{$gen_os1}->[0];
						my $ort_ns2_end = $gene_hit_locs->{$chr_ns2}->{$gen_os1}->[1];
						my $ort_ns2_str = $gene_hit_locs->{$chr_ns2}->{$gen_os1}->[3];
						# Save coordinates of flank gene orthologs
						@{$flank_ortho_locs_s2_new{$pic}[$i]{$chr}{$gen_os1}} = ($ort_os2_gid,$chr_ns2,$ort_ns2_beg,$ort_ns2_end,$ort_ns2_str);
					}
				}
			}
		}
	}
	return \%flank_ortho_locs_s2_new;
}

sub filter_repetitive_hits {
	# Take blast hits
	my($hits) = @_;
	# Constants
	my $max_reps_hits = 5;
	my $max_rep_share = 0.5;
	# Variables
	my @hits_norep = ();
	# Get pos repetitions
	my $pos_repetitions = get_hits_per_pos($hits);
	# Go through each hit
	foreach my $hit (@{$hits}) {
		# Get pic hit coordinates
		my $pic_beg = $hit->[6];
		my $pic_end = $hit->[7];
		my $hit_len = $hit->[12];
		## Discard gene hits located in a repeat
		my $pos_in_rep = 0; #count
		# Check if region is repetitive (not annotated repeats)
		foreach my $pos ($pic_beg..$pic_end) {
			if ($pos_repetitions->{$pos} >= $max_reps_hits) {
				$pos_in_rep++ unless $hit_len > 1000;
			}
		}
		my $rep_share = $pos_in_rep/$hit_len;
		# Save hit if share in repeat location below threshold
		if ($rep_share < $max_rep_share) {
			push(@hits_norep,$hit);
		}
	}
	return \@hits_norep;
}
## Get hits per position to find repetitive regions
sub get_hits_per_pos {
	# Take blast hits
	my($hits) = @_;
	# Variables
	my %pos_repetitions = ();
	# Go through each hit
	foreach my $hit (@{$hits}) {
		# Get pic hit coordinates
		my $pic_beg = $hit->[6];
		my $pic_end = $hit->[7];
		# Count hits on each position
		foreach my $pos ($pic_beg..$pic_end) {
			$pos_repetitions{$pos}++;
		}
	}
	return \%pos_repetitions;
}

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file,$locs) = @_;
	# Storage variable
	my %rep_data = ();
	# If slice locs are given, get positions
	my %chr_poss = ();
	if ($locs) {
		foreach my $loc (sort keys %{$locs}) {
			# Get loc coordinates
			my $chr = $locs->{$loc}->[1];
			my $beg = $locs->{$loc}->[2];
			my $end = $locs->{$loc}->[3];
			# Save positions
			foreach my $pos ($beg..$end) {
				$chr_poss{$chr}{$pos} = 1;
			}
		}
	}
	# Open input file
	my $in = open_infile($repeatmask_file);
	# Parse repeatmasker file
	while (my $line = <$in>) {
		$line =~ s/\s+$//; #better chomp
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			$line =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line);
			my $chr = $d[4];
			my $beg = $d[5];
			my $end = $d[6];
			my $str = $d[8];
			my $cla = $d[10];
			# Ignore non-TE repeats
			if ($cla =~ /Simple_repeat/ || $cla =~ /Low_complexity/ || $cla =~ /Satellite/) { next }
			# Ignore very short TE family repeats
			my @con_pos = $str eq '+' ? ($d[11],$d[12],$d[13]) : ($d[13],$d[12],$d[11]);
			($con_pos[2]) = ($con_pos[2] =~ /(\d+)/);
			my $len = $con_pos[1]+$con_pos[2];
			if ($len < 100) { next }
			# If slice locs are given, skip repeats outside locs
			my $rep_in = 0;
			if ($chr_poss{$chr}{$beg}) { $rep_in++ }
			if ($chr_poss{$chr}{$end}) { $rep_in++ }
			if ($locs && not $rep_in) { next }
			# If slice locs are given, update repeat coordinates
			if ($locs) {
				# If repeat not completely in pic, update
				if ($rep_in < 2) {
					foreach my $loc (sort keys %{$locs}) {
						# Get loc coordinates
						my $loc_chr = $locs->{$loc}->[1];
						my $loc_beg = $locs->{$loc}->[2];
						my $loc_end = $locs->{$loc}->[3];
						# Same chr/scaff
						if ($loc_chr eq $chr) {
							# Same location
							if ($loc_beg <= $end && $loc_end >= $beg) {
								# Update repeat coordinates if needed
								if ($loc_beg > $beg) { $beg = $loc_beg }
								if ($loc_end < $end) { $end = $loc_end }
							}
						}
					}
				}
				$d[5] = $beg;
				$d[6] = $end;
			}
			# Save repeat data
			push(@{$rep_data{$chr}},\@d);
		}
	}
	return \%rep_data;
}

sub get_bed_reads {
	# Take name of tab file
	my($bin_rds_file) = @_;
	# Initialize variables
	my $global_rds = 0;
	my $pls_rds;
	my $mns_rds;
	# Get bin coordinates and read counts
	my $gnm_bin_rds = get_tab_fields($bin_rds_file);
	# Get total genome-wide number of reads
	foreach my $bin (sort {$a <=> $b} keys %{$gnm_bin_rds}) {
		# Get reads and strand
		my $bin_rds = $gnm_bin_rds->{$bin}->[3];
		my $bin_str = $gnm_bin_rds->{$bin}->[5];
		# Save plus/minus reads
		@{$pls_rds->{$bin}} = @{$gnm_bin_rds->{$bin}} if $bin_str eq '+';
		@{$mns_rds->{$bin}} = @{$gnm_bin_rds->{$bin}} if $bin_str eq '-';
		# Count global reads
		$global_rds += $bin_rds;
	}
	# Return global read count, all read bins and plus and minus read bins
	return ($global_rds,$bin_rds_file,$pls_rds,$mns_rds);
}

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

sub rev_com {
	my($seq) = @_;

	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
	$revcom =~ s/\s//g;

	return $revcom;
}

sub rev {
    my($seq) = @_;
    my $revseq = reverse $seq;
    $revseq =~ s/\s//g;
    return $revseq;
}

sub com {
    my($seq) = @_;
    $seq =~ tr/ACGTUNacgtun/TGCAANtgcaan/;
    return $seq;
    $seq =~ s/\s//g;
}

# Add commas as thounds separators
sub add_thousands_separators {
   # Take number
   my($number) = @_;
   my($predec,$decimal) = split(/\./,$number);
   # Add commas as thousands separators to number
   $predec = reverse join ",", (reverse $predec) =~ /(\d{1,3})/g;
   # Return number
   $number = $decimal ? $predec.'.'.$decimal : $predec;
   return $number;
}
