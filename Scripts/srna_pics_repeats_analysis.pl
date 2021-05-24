#!/usr/bin/perl
use strict;
use warnings;

# Global constants
my $min_len = 23;
my $max_len = 29;
my $mir_seq_file = "Dmel_miRNA";
my $pic_loc_file = "dm3.piRNAcluster.bed";
my $rep_loc_file = "dm3.transposon.rm.bed";
my @pic_names = ("42AB","20A","38C","flam");
# Global variables
# Options
$|=1; #Autoflush

# Program name
print("\n--- $0 ---\n");
# Usage messaging
my $USAGE = "perl $0 <Index (Integer)>\n";
unless ($ARGV[0] || $ARGV[0] == 0) {
   die("\nUsage: $USAGE\n");
}
# Collect command line arguments
my $index = $ARGV[0];

my @srna_file_prefixes = (
    "NEBNext03.trimmed.cat",
    "NEBNext13.trimmed.cat",
    "w1118_trimmed",
    "Cl2pD_trimmed",
    "w_38CDf1_trimmed",
    "w_42ABDf1_trimmed",
    "Cl2pD_38CDf1_trimmed",
    "Cl2pD_42ABDf1_trimmed",
    "w_38CDf1_42ABDf1_trimmed",
);

foreach my $i (0..$#srna_file_prefixes) {

    unless ($i == $index) { next; }

    my $file_prefix = $srna_file_prefixes[$i];
    print("$i $file_prefix\n");

    # Count trimmed reads
    my $trimmed_count = `wc -l $file_prefix.fq`;
    ($trimmed_count) = ($trimmed_count =~ /\s+(\d+)\s+/);
    $trimmed_count = $trimmed_count/4;

    ## Filtering of processed reads
    # Map trimmed reads to microRNAs, extract matching reads
    unless (-e "$file_prefix.miR.bam") {
        system("bowtie -a --best -n 0 -k 1 $mir_seq_file -q $file_prefix.fq -S | samtools view -S -b -F 4 | samtools sort > $file_prefix.miR.bam");
    }
    # Count microRNA-matching reads
    my $mir_count = `samtools view -c $file_prefix.miR.bam`;
    chomp($mir_count);
    # Map trimmed reads to ncRNAs, extract NON-matching (i.e. filtered) reads and save as bam file
    unless (-e "$file_prefix.filt.bam") {
        system("bowtie -a --best -n 2 -k 1 refRNA -q $file_prefix.fq -S | samtools view -S -b -f 4 | samtools sort > $file_prefix.filt.bam");
    }
    # Extract filtered reads from bam file and save as fasta file
    unless (-e "$file_prefix.filt.fa") {
        system("samtools fasta $file_prefix.filt.bam > $file_prefix.filt.fa");
    }
    # Remove filtered reads that are outside of piRNA size range
    unless (-e "$file_prefix.filt.$min_len-$max_len.fa") {
        my $outfile1 = filter_read_length_fasta("$file_prefix.filt.fa",$min_len,$max_len);
    }

    ## Genome mapping of filtered reads
    # Map piRNA-sized filtered reads to genome (unique mappers) and save as bam file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.bam") {
        system("bowtie -m 1 -n 0 dm3 -f $file_prefix.filt.$min_len-$max_len.fa -S | samtools view -S -b -F 4 | samtools sort > $file_prefix.filt.$min_len-$max_len.dm3.uni.bam");
    }
    # Convert bam file of genome-mapped piRNA-sized filtered reads (unique mappers) to sam file and annotate with unitas
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.sam") {
        system("samtools view $file_prefix.filt.$min_len-$max_len.dm3.uni.bam > $file_prefix.filt.$min_len-$max_len.dm3.uni.sam");
        system("perl unitas_1.7.7.pl -i $file_prefix.filt.$min_len-$max_len.dm3.uni.sam -s drosophila_melanogaster -skip_mapping -skip_dust -pp");
    }
    # Count genome-mapped piRNA-sized filtered reads (unique mappers)
    my $pir_uni_count = 0;
    if (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.sam") {
        $pir_uni_count = get_total_sam_reads("$file_prefix.filt.$min_len-$max_len.dm3.uni.sam");
        # Normlize read counts
        $pir_uni_count = $pir_uni_count/$mir_count*1_000_000;
    }
    # Convert bam file of genome-mapped piRNA-sized filtered reads (unique mappers) to bed file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.bed") {
        system("bedtools bamtobed -i $file_prefix.filt.$min_len-$max_len.dm3.uni.bam > $file_prefix.filt.$min_len-$max_len.dm3.uni.bed");
    }
    # Convert bam file of genome-mapped piRNA-sized filtered reads (unique mappers) to bedgraph file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.bedgraph") {
        system("bedtools genomecov -bg -ibam $file_prefix.filt.$min_len-$max_len.dm3.uni.bam > $file_prefix.filt.$min_len-$max_len.dm3.uni.bedgraph");
    }

    # Map piRNA-sized filtered reads to genome (all mappers) and save as bam file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.bam") {
        system("bowtie -a --best -n 0 -k 1 dm3 -f $file_prefix.filt.$min_len-$max_len.fa -S | samtools view -S -b -F 4 | samtools sort > $file_prefix.filt.$min_len-$max_len.dm3.all.bam");
    }
    # Convert bam file of genome-mapped piRNA-sized filtered reads (all mappers) to sam file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.sam") {
        system("samtools view $file_prefix.filt.$min_len-$max_len.dm3.all.bam > $file_prefix.filt.$min_len-$max_len.dm3.all.sam");
    }
    # Count genome-mapped piRNA-sized filtered reads (all mappers)
    my $pir_all_count = 0;
    if (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.sam") {
        $pir_all_count = get_total_sam_reads("$file_prefix.filt.$min_len-$max_len.dm3.all.sam");
        # Normlize read counts
        $pir_all_count = $pir_all_count/$mir_count*1_000_000;
    }
    # Convert bam file of genome-mapped piRNA-sized filtered reads (unique mappers) to bed file
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.bed") {
        system("bedtools bamtobed -i $file_prefix.filt.$min_len-$max_len.dm3.all.bam > $file_prefix.filt.$min_len-$max_len.dm3.all.bed");
    }

    ## Print read counts
    #$raw_count = add_thousands_separators($raw_count);
    $trimmed_count = add_thousands_separators($trimmed_count);
    $mir_count = add_thousands_separators($mir_count);
    $pir_uni_count = add_thousands_separators(int($pir_uni_count));
    $pir_all_count = add_thousands_separators(int($pir_all_count));
    print("$trimmed_count\t$mir_count\t$pir_uni_count\t$pir_all_count\n");
    next;

    ## Get piRNA expression per piRNA cluster
    # Clustered piRNAs count (unique mappers)
    my $pic_uni_count = 0;
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.pics.txt") {
        # Determine expression of known piRNA clusters from unique mappers
        my $pic_expr_uni = get_loc_expression($pic_loc_file,"$file_prefix.filt.$min_len-$max_len.dm3.uni.bed");
        # Open output file
        my $outfile1 = "$file_prefix.filt.$min_len-$max_len.dm3.uni.pics.txt";
        my $out1 = open_outfile($outfile1);
        # Go through each line of locus data
        foreach my $pic (sort { $a cmp $b } keys %{$pic_expr_uni}) {
        	printf($out1 "%s\t%.2f\t%.2f\n", $pic,($pic_expr_uni->{$pic}->[0]/$mir_count*1_000_000),($pic_expr_uni->{$pic}->[1]/$mir_count*1_000_000));
            $pic_uni_count += $pic_expr_uni->{$pic}->[0]+$pic_expr_uni->{$pic}->[1];
        }
        # Normlize read counts
        $pic_uni_count = $pic_uni_count/$mir_count*1_000_000;
        # Close file and delete storage variable
        close($out1);
        undef(%{$pic_expr_uni});
    } else {
        $pic_uni_count = get_table_sum("$file_prefix.filt.$min_len-$max_len.dm3.uni.pics.txt");
    }

    # Clustered piRNAs count (all mappers)
    my $pic_all_count = 0;
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.pics.txt") {
        # Determine expression of known piRNA clusters from all mappers
        my $pic_expr_all = get_loc_expression($pic_loc_file,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
        # Open output file
        my $outfile2 = "$file_prefix.filt.$min_len-$max_len.dm3.all.pics.txt";
        my $out2 = open_outfile($outfile2);

        # Go through each line of locus data
        foreach my $pic (sort { $a cmp $b } keys %{$pic_expr_all}) {
        	printf($out2 "%s\t%.2f\t%.2f\n", $pic,($pic_expr_all->{$pic}->[0]/$mir_count*1_000_000),($pic_expr_all->{$pic}->[1]/$mir_count*1_000_000));
            $pic_all_count += $pic_expr_all->{$pic}->[0]+$pic_expr_all->{$pic}->[1];
        }
        # Normlize read counts
        $pic_all_count = $pic_all_count/$mir_count*1_000_000;
        # Close file and delete storage variable
        close($out2);
        undef(%{$pic_expr_all});
    } else {
        $pic_all_count = get_table_sum("$file_prefix.filt.$min_len-$max_len.dm3.all.pics.txt");
    }

    ## Get piRNA reads per transposon family
    # Transposon piRNAs count (unique mappers)
    my $rep_uni_count = 0;
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.txt") {
        # Determine expression of known piRNA clusters from unique mappers
        my $rep_expr_uni = get_loc_expression($rep_loc_file,"$file_prefix.filt.$min_len-$max_len.dm3.uni.bed");
        # Open output file
        my $outfile3 = "$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.txt";
        my $out3 = open_outfile($outfile3);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_uni}) {
        	printf($out3 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_uni->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_uni->{$rep}->[1]/$mir_count*1_000_000));
            $rep_uni_count += $rep_expr_uni->{$rep}->[0]+$rep_expr_uni->{$rep}->[1];
        }
        # Normlize read counts
        $rep_uni_count = $rep_uni_count/$mir_count*1_000_000;
        # Close file and delete storage variable
        close($out3);
        undef(%{$rep_expr_uni});
    } else {
        $rep_uni_count = get_table_sum("$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.txt");
    }

    # Transposon piRNAs count (all mappers)
    my $rep_all_count = 0;
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.txt") {
        # Determine expression of known piRNA clusters from all mappers
        my $rep_expr_all = get_loc_expression($rep_loc_file,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
        # Open output file
        my $outfile4 = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.txt";
        my $out4 = open_outfile($outfile4);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all}) {
        	printf($out4 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all->{$rep}->[1]/$mir_count*1_000_000));
            $rep_all_count += $rep_expr_all->{$rep}->[0]+$rep_expr_all->{$rep}->[1];
        }
        # Normlize read counts
        $rep_all_count = $rep_all_count/$mir_count*1_000_000;
        # Close file and delete storage variable
        close($out4);
        undef(%{$rep_expr_all});
    } else {
        $rep_all_count = get_table_sum("$file_prefix.filt.$min_len-$max_len.dm3.all.reps.txt");
    }

    # Get piRNA expression for each transposon family for selected clusters
    foreach my $pic_name (@pic_names) {
        my $pic_rep_file_c1 = get_reps_in_pic($pic_loc_file,$rep_loc_file,$pic_name);
        my $rep_expr_uni_c1 = get_loc_expression($pic_rep_file_c1,"$file_prefix.filt.$min_len-$max_len.dm3.uni.bed");
        # Open output file
        my $outfile_c1 = "$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.$pic_name.txt";
        my $out_c1 = open_outfile($outfile_c1);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_uni_c1}) {
            printf($out_c1 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_uni_c1->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_uni_c1->{$rep}->[1]/$mir_count*1_000_000));
        }
        close($out_c1);
    }
    # Get piRNA expression for each transposon family for non-selected clusters
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.otherpics.txt") {
        # Get piRNA expression for each transposon family for non-selected clusters
        my $pic_rep_file_other = get_reps_outside_given_pics($pic_loc_file,$rep_loc_file,\@pic_names);
        my $rep_expr_uni_other = get_loc_expression($pic_rep_file_other,"$file_prefix.filt.$min_len-$max_len.dm3.uni.bed");
        # Open output file
        my $outfile_other = "$file_prefix.filt.$min_len-$max_len.dm3.uni.reps.otherpics.txt";
        my $out_other = open_outfile($outfile_other);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_uni_other}) {
            printf($out_other "%s\t%.2f\t%.2f\n", $rep,($rep_expr_uni_other->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_uni_other->{$rep}->[1]/$mir_count*1_000_000));
        }
        close($out_other);
    }

    # Get piRNA expression for each transposon family for selected clusters
    foreach my $pic_name (@pic_names) {
        unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.$pic_name.txt") {
            my $pic_rep_file_c1 = get_reps_in_pic($pic_loc_file,$rep_loc_file,$pic_name);
            my $rep_expr_all_c1 = get_loc_expression($pic_rep_file_c1,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
            # Open output file
            my $outfile_c1 = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.$pic_name.txt";
            my $out_c1 = open_outfile($outfile_c1);
            # Go through each line of locus data
            foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all_c1}) {
                printf($out_c1 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all_c1->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all_c1->{$rep}->[1]/$mir_count*1_000_000));
            }
            close($out_c1);
        }
    }
    # Get piRNA expression for each transposon family for non-selected clusters
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.otherpics.txt") {
        # Get piRNA expression for each transposon family for non-selected clusters
        my $pic_rep_file_other = get_reps_outside_given_pics($pic_loc_file,$rep_loc_file,\@pic_names);
        my $rep_expr_all_other = get_loc_expression($pic_rep_file_other,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
        # Open output file
        my $outfile_other = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.otherpics.txt";
        my $out_other = open_outfile($outfile_other);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all_other}) {
            printf($out_other "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all_other->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all_other->{$rep}->[1]/$mir_count*1_000_000));
        }
        close($out_other);
    }

    # Get piRNA expression for each transposon family for selected clusters - potential
    foreach my $pic_name (@pic_names) {
        unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.$pic_name.pot.txt") {
            my $pic_rep_file_c1 = get_reps_in_pic($pic_loc_file,$rep_loc_file,$pic_name);
            my $rep_expr_all_c1 = get_loc_expression_potential($pic_rep_file_c1,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
            # Open output file
            my $outfile_c1 = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.$pic_name.pot.txt";
            my $out_c1 = open_outfile($outfile_c1);
            # Go through each line of locus data
            foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all_c1}) {
                printf($out_c1 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all_c1->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all_c1->{$rep}->[1]/$mir_count*1_000_000));
            }
            close($out_c1);
        }
    }

    # Get piRNA expression for each transposon family for non-selected clusters - potential
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.otherpics.pot.txt") {
        # Get piRNA expression for each transposon family for non-selected clusters
        my $pic_rep_file_other = get_reps_outside_given_pics($pic_loc_file,$rep_loc_file,\@pic_names);
        my $rep_expr_all_other = get_loc_expression_potential($pic_rep_file_other,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
        # Open output file
        my $outfile_other = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.otherpics.pot.txt";
        my $out_other = open_outfile($outfile_other);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all_other}) {
            printf($out_other "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all_other->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all_other->{$rep}->[1]/$mir_count*1_000_000));
        }
        close($out_other);
    }

    # Get piRNA expression for each transposon family for for genome - potential
    unless (-e "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.pot.txt") {
        my $rep_expr_all_c1 = get_loc_expression_potential($rep_loc_file,"$file_prefix.filt.$min_len-$max_len.dm3.all.bed");
        # Open output file
        my $outfile_c1 = "$file_prefix.filt.$min_len-$max_len.dm3.all.reps.pot.txt";
        my $out_c1 = open_outfile($outfile_c1);
        # Go through each line of locus data
        foreach my $rep (sort { $a cmp $b } keys %{$rep_expr_all_c1}) {
            printf($out_c1 "%s\t%.2f\t%.2f\n", $rep,($rep_expr_all_c1->{$rep}->[0]/$mir_count*1_000_000),($rep_expr_all_c1->{$rep}->[1]/$mir_count*1_000_000));
        }
        close($out_c1);
    }

    ## Print read counts
    #printf("%d\t", $mir_count);
    printf("%.2f\t", $pir_uni_count);
    printf("%.2f\t", $pir_all_count);
    printf("%.2f\t", $pic_uni_count);
    printf("%.2f\t", $pic_all_count);
    printf("%.2f\t", $rep_uni_count);
    printf("%.2f\n", $rep_all_count);

}

exit;

################################# subroutines #################################

sub fastq_read_number {
    # Take name of tab file
	my($infile) = @_;
    # Initialize read count
    my $reads_n = 0;
    # Open input file
	my $in = open_infile($infile);
    # Go through each line
    while (my $line = <$in>) {
        $line =~ s/\s+$//; #better chomp
        # Non-empty line
        if ($line !~ /^$/) {
            # Four lines per read
            $reads_n += 0.25;
        }
    }
    return $reads_n;
}

sub filter_read_length_fasta {
    # Take name of tab file
	my($infile,$min,$max) = @_;
    # Get file data
	my @in_data = get_file_data_array($infile);
    # Open output file
    my $outfile = $infile;
    $outfile =~ s/\.fa/\.$min-$max\.fa/;
    my $out = open_outfile($outfile);
    # Parse file content
    my $head = '';
    foreach my $line (@in_data) {
        if ($line =~ /^>/) {
            $head = $line;
        } else {
            if (length($line) >= $min && length($line) <= $max) {
                print($out "$head\n$line\n")
            }
        }
    }
    return $outfile;
}

sub get_total_sam_reads {
    # Take name of tab file
    my($infile) = @_;
    # Read count variable
    my $reads_count = 0;
    # Get file data
    my $in = open_infile($infile);
    # Storage variable
    my %rna_count_per_seq = ();
    # Parse sam file
    while (my $line = <$in>) {
        # Split line
        my @d = split(/\t/,$line);
        # Count unique read ids
        $rna_count_per_seq{$d[0]}++;
    }
    # Get total read count
    $reads_count = keys %rna_count_per_seq;
    # Return total read count
    return $reads_count
}

sub get_reps_in_pic {
    # Collect command line arguments
    my($gnmloci_bed_file,$repeats_bed_file,$pic_name) = @_;
    # Get bed file data
    my($genome_loci,$loc_cnts) = get_bed_data($gnmloci_bed_file);
    # Get bedgraph file data
    my($repeat_loci,$rep_cnts) = get_bed_data($repeats_bed_file);
    # Open output file
    my $outfile = "$repeats_bed_file.$pic_name.bed";
    my $out = open_outfile($outfile);
    # Go through each chromosome
    foreach my $chr (sort { $a cmp $b } keys %{$genome_loci}) {
        # Go through each locus
        foreach my $gloc (sort { $a <=> $b } keys %{$genome_loci->{$chr}}) {
        	# Get loc data
        	my $loc_chr = $genome_loci->{$chr}->{$gloc}->[0];
        	my $loc_beg = $genome_loci->{$chr}->{$gloc}->[1];
        	my $loc_end = $genome_loci->{$chr}->{$gloc}->[2];
        	my $loc_nam = $genome_loci->{$chr}->{$gloc}->[3];
            # Skip if not pic of interest
            if ($loc_nam ne $pic_name) { next;}
            # Go through each line of repeat data
        	foreach my $rloc (sort { $a <=> $b } keys %{$repeat_loci->{$chr}}) {
                # Get seq data
        		my $rep_chr = $repeat_loci->{$chr}->{$rloc}->[0];
        		my $rep_beg = $repeat_loci->{$chr}->{$rloc}->[1];
        		my $rep_end = $repeat_loci->{$chr}->{$rloc}->[2];
        		my $rep_sid = $repeat_loci->{$chr}->{$rloc}->[3];
                my $rep_scr = $repeat_loci->{$chr}->{$rloc}->[4];
                my $rep_str = $repeat_loci->{$chr}->{$rloc}->[5];
        		# Check if rep is in locus
        		if ($loc_beg < $rep_beg && $loc_end > $rep_end) {
        			print($out "$rep_chr\t$rep_beg\t$rep_end\t$rep_sid\t$rep_scr\t$rep_str\n");
        		} elsif ($loc_beg > $rep_end) {
                    delete($repeat_loci->{$chr}->{$rloc});
                } elsif ($loc_end < $rep_beg) {
                    last;
                }
            }
        }
    }
    return $outfile;
}
sub get_reps_outside_given_pics {
    # Collect command line arguments
    my($gnmloci_bed_file,$repeats_bed_file,$pic_names) = @_;
    # Get bed file data
    my($genome_loci,$loc_cnts) = get_bed_data($gnmloci_bed_file);
    # Get bedgraph file data
    my($repeat_loci,$rep_cnts) = get_bed_data($repeats_bed_file);
    my %pic_names = map { $_ => 1 } @{$pic_names};
    # Open output file
    my $outfile = "$repeats_bed_file.otherpics.bed";
    my $out = open_outfile($outfile);
    # Go through each chromosome
    foreach my $chr (sort { $a cmp $b } keys %{$genome_loci}) {
        # Go through each locus
        foreach my $gloc (sort { $a <=> $b } keys %{$genome_loci->{$chr}}) {
        	# Get loc data
        	my $loc_chr = $genome_loci->{$chr}->{$gloc}->[0];
        	my $loc_beg = $genome_loci->{$chr}->{$gloc}->[1];
        	my $loc_end = $genome_loci->{$chr}->{$gloc}->[2];
        	my $loc_nam = $genome_loci->{$chr}->{$gloc}->[3];
            # Skip if pic included in given pics
            if ($pic_names{$loc_nam}) { next;}
            # Go through each line of repeat data
        	foreach my $rloc (sort { $a <=> $b } keys %{$repeat_loci->{$chr}}) {
                # Get seq data
        		my $rep_chr = $repeat_loci->{$chr}->{$rloc}->[0];
        		my $rep_beg = $repeat_loci->{$chr}->{$rloc}->[1];
        		my $rep_end = $repeat_loci->{$chr}->{$rloc}->[2];
        		my $rep_sid = $repeat_loci->{$chr}->{$rloc}->[3];
                my $rep_scr = $repeat_loci->{$chr}->{$rloc}->[4];
                my $rep_str = $repeat_loci->{$chr}->{$rloc}->[5];
        		# Check if rep is in locus
        		if ($loc_beg < $rep_beg && $loc_end > $rep_end) {
        			print($out "$rep_chr\t$rep_beg\t$rep_end\t$rep_sid\t$rep_scr\t$rep_str\n");
        		} elsif ($loc_beg > $rep_end) {
                    delete($repeat_loci->{$chr}->{$rloc});
                } elsif ($loc_end < $rep_beg) {
                    last;
                }
            }
        }
    }
    return $outfile;
}

sub get_loc_expression {
    # Collect command line arguments
    my($gnmloci_bed_file,$srnaseq_bed_file) = @_;
    # Get bed file data
    my($genome_loci,$loc_cnts) = get_bed_data($gnmloci_bed_file);
    # Get bedgraph file data
    my($srnaseq_bed,$seq_hits) = get_bed_data($srnaseq_bed_file);
    # Storage variable
    my %gloc_expr = ();
    my $outfile = "test.$index.txt";
    my $out = open_outfile($outfile);
    # Go through each chromosome
    foreach my $chr (sort { $a cmp $b } keys %{$genome_loci}) {
        # Go through each locus
        foreach my $gloc (sort { $a <=> $b } keys %{$genome_loci->{$chr}}) {
        	# Get loc data
        	my $loc_chr = $genome_loci->{$chr}->{$gloc}->[0];
        	my $loc_beg = $genome_loci->{$chr}->{$gloc}->[1];
        	my $loc_end = $genome_loci->{$chr}->{$gloc}->[2];
        	my $loc_nam = $genome_loci->{$chr}->{$gloc}->[3];
            print($out "$loc_chr\t$loc_beg\t$loc_end\t$loc_nam\n");
        	@{$gloc_expr{$loc_nam}} = (0,0) unless $gloc_expr{$loc_nam};
        	# Go through each line of rnaseq data
        	foreach my $srna (sort { $a <=> $b } keys %{$srnaseq_bed->{$chr}}) {
        		# Get seq data
        		my $seq_chr = $srnaseq_bed->{$chr}->{$srna}->[0];
        		my $seq_beg = $srnaseq_bed->{$chr}->{$srna}->[1];
        		my $seq_end = $srnaseq_bed->{$chr}->{$srna}->[2];
        		my $seq_sid = $srnaseq_bed->{$chr}->{$srna}->[3];
                my $seq_str = $srnaseq_bed->{$chr}->{$srna}->[5];
        		# Check if seq is in locus
        		if ($loc_beg < $seq_beg && $loc_end > $seq_end) {
        			$gloc_expr{$loc_nam}[0] += 1/$seq_hits->{$seq_sid} if $seq_str eq '+';
                    $gloc_expr{$loc_nam}[1] += 1/$seq_hits->{$seq_sid} if $seq_str eq '-';
                    #print("$seq_hits->{$seq_sid}\n");
        		} elsif ($loc_beg > $seq_end) {
                    delete($srnaseq_bed->{$chr}->{$srna});
                } elsif ($loc_end < $seq_beg) {
                    last;
                }
        	}
        }
    }

    return \%gloc_expr;
}

sub get_loc_expression_potential {
    # Collect command line arguments
    my($gnmloci_bed_file,$srnaseq_bed_file) = @_;
    # Count hits per piC for potential read counts
    my %pic_hits = ();
    # Seperate counting
    my $potential = 1;
    if ($potential) {
        # Get bed file data
        my($genome_loci,$loc_cnts) = get_bed_data($gnmloci_bed_file);
        # Get bed file data
        my($srnaseq_bed,$seq_hits) = get_bed_data($srnaseq_bed_file);
        # Go through each chromosome
        foreach my $chr (sort { $a cmp $b } keys %{$genome_loci}) {
            # Go through each locus
            foreach my $gloc (sort { $a <=> $b } keys %{$genome_loci->{$chr}}) {
                # Get loc data
                my $loc_chr = $genome_loci->{$chr}->{$gloc}->[0];
                my $loc_beg = $genome_loci->{$chr}->{$gloc}->[1];
                my $loc_end = $genome_loci->{$chr}->{$gloc}->[2];
                my $loc_nam = $genome_loci->{$chr}->{$gloc}->[3];
                # Go through each line of rnaseq data
                foreach my $srna (sort { $a <=> $b } keys %{$srnaseq_bed->{$chr}}) {
                    # Get seq data
                    my $seq_chr = $srnaseq_bed->{$chr}->{$srna}->[0];
                    my $seq_beg = $srnaseq_bed->{$chr}->{$srna}->[1];
                    my $seq_end = $srnaseq_bed->{$chr}->{$srna}->[2];
                    my $seq_sid = $srnaseq_bed->{$chr}->{$srna}->[3];
                    my $seq_str = $srnaseq_bed->{$chr}->{$srna}->[5];
                    # Check if seq is in locus
                    if ($loc_beg < $seq_beg && $loc_end > $seq_end) {
                        $pic_hits{$seq_sid}++;
                    } elsif ($loc_beg > $seq_end) {
                        delete($srnaseq_bed->{$chr}->{$srna});
                    } elsif ($loc_end < $seq_beg) {
                        last;
                    }
                }
            }
        }
    }

    # Get bed file data
    my($genome_loci,$loc_cnts) = get_bed_data($gnmloci_bed_file);
    # Get bed file data
    my($srnaseq_bed,$seq_hits) = get_bed_data($srnaseq_bed_file);
    # Get sequences and read ids from sam file
    #my $sam_file = $srnaseq_bed_file;
    #$sam_file =~ s/\.bed/\.sam/;
    #my($sam_seqs,$sam_rids) = get_sam_ids_seqs($sam_file);
    my %incl_seqs = ();
    # Storage variable
    my %gloc_expr = ();
    my $outfile = "test.$index.txt";
    my $out = open_outfile($outfile);
    # Go through each chromosome
    foreach my $chr (sort { $a cmp $b } keys %{$genome_loci}) {
        # Go through each locus
        foreach my $gloc (sort { $a <=> $b } keys %{$genome_loci->{$chr}}) {
        	# Get loc data
        	my $loc_chr = $genome_loci->{$chr}->{$gloc}->[0];
        	my $loc_beg = $genome_loci->{$chr}->{$gloc}->[1];
        	my $loc_end = $genome_loci->{$chr}->{$gloc}->[2];
        	my $loc_nam = $genome_loci->{$chr}->{$gloc}->[3];
            print($out "$loc_chr\t$loc_beg\t$loc_end\t$loc_nam\n");
        	@{$gloc_expr{$loc_nam}} = (0,0) unless $gloc_expr{$loc_nam};
        	# Go through each line of rnaseq data
        	foreach my $srna (sort { $a <=> $b } keys %{$srnaseq_bed->{$chr}}) {
        		# Get seq data
        		my $seq_chr = $srnaseq_bed->{$chr}->{$srna}->[0];
        		my $seq_beg = $srnaseq_bed->{$chr}->{$srna}->[1];
        		my $seq_end = $srnaseq_bed->{$chr}->{$srna}->[2];
        		my $seq_sid = $srnaseq_bed->{$chr}->{$srna}->[3];
                my $seq_str = $srnaseq_bed->{$chr}->{$srna}->[5];
                #my $seq_seq = $sam_rids->{$seq_sid};
                #my $seq_ids = keys %{$sam_seqs->{$seq_seq}};
        		# Check if seq is in locus
        		if ($loc_beg < $seq_beg && $loc_end > $seq_end) {
                    $gloc_expr{$loc_nam}[0] += 1/$pic_hits{$seq_sid} if $seq_str eq '+';
                    $gloc_expr{$loc_nam}[1] += 1/$pic_hits{$seq_sid} if $seq_str eq '-';
        			#$gloc_expr{$loc_nam}[0] += 1 if $seq_str eq '+';
                    #$gloc_expr{$loc_nam}[1] += 1 if $seq_str eq '-';
                    #print("$seq_hits->{$seq_sid}\n");
        		} elsif ($loc_beg > $seq_end) {
                    delete($srnaseq_bed->{$chr}->{$srna});
                } elsif ($loc_end < $seq_beg) {
                    last;
                }
        	}
        }
    }

    return \%gloc_expr;
}

sub get_sam_ids_seqs {
	# Take name of sam file
	my($samfile) = @_;
	# Initialize variables
	my %sam_seqs = ();
    my %sam_rids = ();
	# Open all mappers sam file
	my $sam_a = open_infile($samfile);
	# Go through each line
	while (my $line = <$sam_a>) {
		$line =~ s/\s+$//; #better chomp
		# Get map data
		my @d = split(/\t/,$line);
		my $rid = $d[0];
		my $flg = $d[1];
		my $str = $flg == 16 ? '-' : '+';
		my $seq = $d[9];
		# Reverse sequence if on minus strand
		if ($str eq '-') {
			$seq = reverse($seq);
			$seq =~ tr/ATGC/TACG/;
		}
		# Save sequence and read id
		$sam_seqs{$seq}{$rid} = 1;
        $sam_rids{$rid} = $seq;
	}
	return(\%sam_seqs,\%sam_rids);
}

sub get_bed_data {
	# Take name of tab file
	my($infile) = @_;
	# Get file data
    my $in = open_infile($infile);
	# Global tree hashes for each species
	my %data_fields = ();
    my %hits_per_id = ();
	# Set 0 as start index
	my $id_i = 0;
	# Go through file data
	while (my $line = <$in>) {
        $line =~ s/\s+$//; #better chomp
        # Get line data
        my @d = split(/\t/,$line);
        # Save data fields
        @{$data_fields{$d[0]}{$id_i}} = @d;
        $hits_per_id{$d[3]}++;
		$id_i++;
    }
	# Return data fields
	return(\%data_fields,\%hits_per_id);
}

sub get_table_sum {
    # Take name of tab file
    my($infile) = @_;
    # Count
    my $count = 0;
    # Get file data
    my @file_data = get_file_data_array($infile);
    # Parse file data
    foreach my $line (@file_data) {
        # Get data fields
        my @d = split(/\t/,$line);
        # Get sum of column 1 and 2
        $count += ($d[1]+$d[2]);
    }
    return $count;
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
