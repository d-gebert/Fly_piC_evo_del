#!/usr/bin/perl
use strict;
use warnings;

# Load modules
use Getopt::Long;
# Check if PDF::Create is installed
my $pdf_installed = 0;
eval { require PDF::Create };
unless ($@) { $pdf_installed = 1; }

# Options
my $query_file = '';
my $subject_file = '';
my $int_loc_file = '';
my $sensitive = 0;
my $ident_th = 0;
my $help = 0;
# Program names/paths (dependencies)
my $blastn = 'blastn';
# Global constants
my $outfmt = '7 qseqid sseqid pident evalue qcovs qlen qstart qend slen sstart send sstrand length bitscore';
# Global variables
my @argv = @ARGV;
# Autoflush
$|=1;

# Program name
print("\n--- $0 ---\n");

# Test dependencies
my $depends = '';
foreach my $prog ($blastn) {
	system("$prog -h 2> .progtest_err > .progtest_log");
	open(my $fh, '<', ".progtest_err");
	my $frst_ln = <$fh>;
	if ($frst_ln && $frst_ln =~ /not found/) { $depends .= "$prog\n" }
	unlink('.progtest_err','.progtest_log');
}
if ($depends) { die("\nERROR: Cannot start program.\nFollowing dependencies need to be installed:\n$depends\n"); }

# Collect command line arguments
GetOptions (
	"repeat|r=s"	=> \$query_file,	# Repeat fasta file
	"genome|g=s"	=> \$subject_file,	# Genome fasta file
	"target|t=s"	=> \$int_loc_file,	# Loci of interest bed file
	"sensitive|s"	=> \$sensitive,		# Flag: Use sensitive search
    "min_ident|m=i"	=> \$ident_th,		# Minimum %-identity cut-off
	"help|h"		=> \$help)			# Help
or die(usage());

## Handle options / error messaging
# Help message
if ($help) { help() };
# Usage message if run without options
if (!@argv) { die(usage()) }
# Input files
if (!$query_file || !$subject_file) { die(usage("\nError! Input files not specified.")) }
unless (-e $query_file) { die(usage("\nError! File \'$query_file\' does not exist.")) }
unless (-e $subject_file) { die(usage("\nError! File \'$subject_file\' does not exist.")) }

# Set sensitive search algorithm (dc-megablast) if -s option is given
my $dcmopt = $sensitive ? '-task dc-megablast' : '';

# Call blast
my $blast_out = "$query_file.$subject_file.bln";
if ($sensitive) { $blast_out =~ s/bln/dcm/; }
unless (-e $blast_out) {
	system("$blastn $dcmopt -query $query_file -subject $subject_file -out $blast_out -outfmt \'$outfmt\' >>.log 2>&1");
    system("$blastn $dcmopt -query $query_file -subject $subject_file -out $blast_out.aln >>.log 2>&1");
}
# Extract filtered blast hits
my $blast_hits = get_blastn_hits($blast_out);
# Group blast hits
my $chr_hit_locs = get_hit_group_locs($blast_hits);

# Get alignment mismatches
my $aln_mms = get_blastn_alignment_mismatches("$blast_out.aln");

# Get loci of interest coordinates
my $int_locs = get_int_loc_poss($int_loc_file) if -e "$int_loc_file";

# Hashes for hit group sorting
my %group_scores = ();
my %group_scores_keys = ();
# Sort hits groups with new ids
my $id = 1;
# Go through each chromosome in hits group loc list
foreach my $chr (sort {$chr_hit_locs->{$a} <=> $chr_hit_locs->{$b}} keys %{$chr_hit_locs}) {
    if ($chr eq 'chrY') { next; }
    # Go through each group in hits group loc list on same chromosome
    foreach my $grp (sort {$chr_hit_locs->{$chr}->{$a} <=> $chr_hit_locs->{$chr}->{$b}} keys %{$chr_hit_locs->{$chr}}) {
        # Go through each strand in hits group loc list
        foreach my $str (sort keys %{$chr_hit_locs->{$chr}->{$grp}}) {
            # Go through each hit in group
            foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
                # Get hit identity, coordinates and alignment length
                my $iden = $hit->[2];
                my $hit_beg = $hit->[9];
        		my $hit_end = $hit->[10];
                my $alen = $hit->[12];
                # Set multiplying factor (default 1) high if the hit is located in a target locus of interest
                my $factor = 1;
                if ($int_locs->{$chr}->{$hit_beg} || $int_locs->{$chr}->{$hit_end}) { $factor = 1_000_000; }
                # Calculate a total ordering score for hit group (alignment lengths weighted by identity)
                $group_scores{$id} += ($alen*($iden/100)*$factor);
                # Save chromosome, group and strand
                @{$group_scores_keys{$id}} = ($chr,$grp,$str);
            }
            # Increment id number
            $id++;
        }
    }
}

# Open output bed file
my $outfile = "$query_file.$subject_file.bln.bed";
if ($sensitive) { $outfile =~ s/bln/dcm/; }
my $out = open_outfile($outfile);
# Title line
print($out "query\tqbeg\tqend\tins_id\talen\tqlen\tstrand\tchr\thit_beg\thit_end\tident\tint_ins\n");
# Go through each hit group, sorted by their ordering score
foreach my $order_id (sort {$group_scores{$b} <=> $group_scores{$a}} keys %group_scores) {
    # Get chromosome, group and strand
    my $chr = $group_scores_keys{$order_id}[0];
    my $grp = $group_scores_keys{$order_id}[1];
    my $str = $group_scores_keys{$order_id}[2];
    # Get mean identity
    my $mean_ident = 100;
    my $totl_alen = 0;
    foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
        my $iden = $hit->[2];
        my $alen = $hit->[12];
        $totl_alen += $alen;
        $mean_ident += ($iden*$alen);
    }
    $mean_ident = $mean_ident/$totl_alen;
    # Minimum identity cutoff
    if ($mean_ident < $ident_th) { next; }
    # Go through each hit of the current group
    foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
        # Adjust symbols and coordinates according to subject strand
        my $strandl = $hit->[11];
        my $strands = $hit->[11] eq 'plus' ? '+' : '-';
        my $hit_beg = $hit->[9];
		my $hit_end = $hit->[10];
		if ($hit_beg > $hit_end) {
			$hit_beg = $hit->[10];
			$hit_end = $hit->[9];
		}
        my $name = $hit->[0];
        my $iden = $hit->[2];
        my $qlen = $hit->[5];
        my $qbeg = $hit->[6];
        my $qend = $hit->[7];
        my $alen = $hit->[12];
        # Check if hit is in target locus of interest
        my $int_ins = 'no';
        if ($int_locs->{$chr}->{$hit_beg} || $int_locs->{$chr}->{$hit_end}) { $int_ins = 'int_loc'; }
        # Print hit data to output bed file
        print($out "$name\t$qbeg\t$qend\t$chr.$grp.$strandl\t$alen\t$qlen\t$strands\t$chr\t$hit_beg\t$hit_end\t$iden\t$int_ins\n");
    }
}
close($out);

### Graphics Output ###
# Ask if PDF::Create is installed
if ($pdf_installed) {
    # Get query sequence length and name
    my $query_seq = get_fasta_seqs($query_file);
    my $query_len = 0;
    my $query_name = '';
    foreach my $seq_name (keys %{$query_seq}) {
        $query_len = length($query_seq->{$seq_name});
        $query_name = $seq_name;
    }
    # Get conversion for contig positions into graphic positions
    my %convx_a = ();
    # Get regional start and stop positions
    my $reg_beg = 1;
    my $reg_end = $query_len;
    # Calculate graphic coordinate for each position
    foreach my $xo ($reg_beg..$reg_end) {
        my $xn = ($xo-$reg_beg+1)/(($reg_end-$reg_beg+1)/1000);
        $xn = int($xn+0.5);
        $convx_a{$xo} = $xn;
    }
    # Prepare graphics output file
    my $pdf_name = "$blast_out.pdf";
    # Create pdf file
    my $pdf = PDF::Create->new(
        'filename'     => "$pdf_name",
        'Author'       => 'Perl',
        'Title'        => "$pdf_name",
        'CreationDate' => [ localtime ]
    );
    # Add a A2 sized page
    my $root = $pdf->new_page('MediaBox' => $pdf->get_page_size('A2'));
    # Add a page which inherits its attributes from $root
    my $page1 = $root->new_page;
    # Prepare a font
    my $font = $pdf->font('BaseFont' => 'Helvetica');
    # Query name
    $page1->setrgbcolorstroke(0.1,0.1,0.1);
    $page1->string($font, 12, 20+10, 1600+30, "$query_name");
    # Query rectangle
    $page1->setrgbcolorstroke(0.1,0.1,0.1);
    $page1->rectangle(20,1600,1000,10);
    $page1->stroke();
    # Line number
    my $i = 1;
    # Go through each hit group, sorted by their ordering score
    foreach my $order_id (sort {$group_scores{$b} <=> $group_scores{$a}} keys %group_scores) {
        # Get chromosome, group and strand
        my $chr = $group_scores_keys{$order_id}[0];
        my $grp = $group_scores_keys{$order_id}[1];
        my $str = $group_scores_keys{$order_id}[2];
        # Get mean identity
        my $mean_ident = 100;
        my $totl_alen = 0;
        foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
            my $iden = $hit->[2];
            my $alen = $hit->[12];
            $totl_alen += $alen;
            $mean_ident += ($iden*$alen);
        }
        $mean_ident = $mean_ident/$totl_alen;
        # Minimum identity cutoff
        if ($mean_ident < $ident_th) { next; }
        # Coordinates for subject and query sequence
        my $min_beg = 1_000_000_000;
        my $max_end = 0;
        my $min_beg_rep = 1_000_000_000;
        my $max_end_rep = 0;
        foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
            # Adjust symbols and coordinates according to subject strand
            my $strand = $hit->[11] eq 'plus' ? '+' : '-';
            my $hit_beg = $hit->[9];
    		my $hit_end = $hit->[10];
    		if ($hit_beg > $hit_end) {
    			$hit_beg = $hit->[10];
    			$hit_end = $hit->[9];
    		}
            $min_beg = $hit_beg < $min_beg ? $hit_beg : $min_beg;
            $max_end = $hit_end > $max_end ? $hit_end : $max_end;
            $min_beg_rep = $hit->[6] < $min_beg_rep ? $hit->[6] : $min_beg_rep;
            $max_end_rep = $hit->[7] > $max_end_rep ? $hit->[7] : $max_end_rep;
        }
        # Insert line connecting grouped hits
        $page1->setrgbcolorstroke(0.5,0.5,0.5);
        $page1->line(20+$convx_a{$min_beg_rep}, 1600+5-(12*$i), (20+$convx_a{$max_end_rep}), 1600+5-(12*$i));
        # Insert group coordinates in subject sequence
        $page1->setrgbcolor(0.1,0.1,0.1);
        $page1->string($font, 12, 20+1000+10, 1600+2-(12*$i), "$chr:$min_beg-$max_end");
        # Go through each hit of the current group
        foreach my $hit (sort {$a->[6] <=> $b->[6]} @{$chr_hit_locs->{$chr}->{$grp}->{$str}}) {
            # Get strand and hit coordinates
            my $strand = $hit->[11] eq 'plus' ? '+' : '-';
            my $hit_beg = $hit->[9];
    		my $hit_end = $hit->[10];
    		if ($hit_beg > $hit_end) {
    			$hit_beg = $hit->[10];
    			$hit_end = $hit->[9];
    		}
            # Insert the hit as an alignment block
            $page1->setrgbcolor(0,0,0);
            $page1->setrgbcolorstroke(0,0,0);
            $page1->rectangle(20+$convx_a{$hit->[6]}, 1600-(12*$i), ($convx_a{$hit->[7]}-$convx_a{$hit->[6]}+1), 10);
            $page1->fill();
            # Insert alignment information on mismatches, small deletions and insertions
            foreach my $pos (sort keys %{$aln_mms->{$hit->[-1]}}) {
                # Skip if position is unknown
                next unless $convx_a{$pos};
                # Set RGB according to type: mismatch, deletion or insertion
                my @rgb = ();
                if ($aln_mms->{$hit->[-1]}->{$pos} eq 'mis') {
                    @rgb = (1,0.2,0.2);
                } elsif ($aln_mms->{$hit->[-1]}->{$pos} eq 'ins') {
                    @rgb = (0.2,0.2,1);
                } elsif ($aln_mms->{$hit->[-1]}->{$pos} eq 'del') {
                    @rgb = (1,1,1);
                }
                $page1->setrgbcolor(@rgb);
                $page1->setrgbcolorstroke(@rgb);
                # Insert misaligned position to alignment graphic
                $page1->line(20+$convx_a{$pos}, 1600+9-(12*$i), (20+$convx_a{$pos}, 1600+1-(12*$i)));
            }
        }
        # Increment line number
        $i++;
    }
    # Close the file and write the PDF
    $pdf->close;
} else {
    # Say that graphical output will be skipped
    print("Perl module \'PDF::Create\' is not installed: Graphical output skipped.\n");
    print("Install Perl modules with CPAN. See http://www.cpan.org/modules/INSTALL.html\n");
}
# Delete log file
unlink('.log');

exit;

################################# subroutines #################################

sub get_blastn_alignment_mismatches {
	# Take infile name
	my($file) = @_;
	# Get file data
	my @file_data = get_file_data_array($file);
    # Initialize variables
    my $hid = 0;
    my %aln_seq_q = ();
    my %aln_seq_s = ();
    my %begs_q = ();
    my %ends_q = ();
    my %mms_q = ();
    # Parse tblastn output
    foreach my $line (@file_data) {
        # Results line
        if ($line =~ /^ Score = / || $line =~ /^Matrix/) {
            if ($hid) {
                my $pos_q = $begs_q{$hid}[0];
                for (my $i = 0; $i < length($aln_seq_q{$hid}); $i++) {
                    my $base_q = substr($aln_seq_q{$hid}, $i, 1);
                    my $base_s = substr($aln_seq_s{$hid}, $i, 1);
                    if ($base_q ne '-') {
                        if ($base_q ne $base_s) {
                            $mms_q{$hid}{$pos_q} = 'mis';
                        }
                    } elsif ($base_q eq '-') {
                        $mms_q{$hid}{$pos_q} = 'ins';
                    }
                    if ($base_s eq '-') {
                        $mms_q{$hid}{$pos_q} = 'del';
                    }
                    $pos_q++ unless $base_q eq '-';
                }
            }
            $hid++;
        }
        if ($line =~ /^Query  /) {
            my @d = split(/\s+/, $line);
            my $beg = $d[1];
            my $seq = $d[2];
            my $end = $d[3];
            $aln_seq_q{$hid} .= uc($seq);
            push(@{$begs_q{$hid}},$beg);
            push(@{$ends_q{$hid}},$end);
        }
        if ($line =~ /^Sbjct  /) {
            my @d = split(/\s+/, $line);
            my $seq = $d[2];
            $aln_seq_s{$hid} .= uc($seq);
        }
    }
	return \%mms_q;
}

sub get_blastn_hits {
	# Take infile name
	my($file) = @_;
	# Get file data
	my @file_data = get_file_data_array($file);
	# Initialize variables
	my %blast_hits = ();
    my $i = 0;
	# Parse tblastn output
	foreach my $line (@file_data) {
		# Results line
		if ($line !~ /^#/) {
            $i++;
			# Get info
			my @d = split(/\t/, $line);
            push(@d,$i);
			my $chr_id = $d[1];
			# Filter hits
			#if ($d[2] < 95) { next; }
			#if ($d[4] < 20) { next; }
			#if ($d[12] < 2000) { next; }
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
		my $max_distance = 5_000;
		my %grouped_hits = ();
        my %grouped_hits_clean = ();
        my %grouped_hits_clean2 = ();
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
			# Get distance
			my $distance = $hit_beg-$prev_hit_end;
			# If distance above threshold open new group
			if ($distance > $max_distance) { $group++ }
			# Allocate hit to group
			push(@{$grouped_hits{$group}{$sstrand}},$hit);
			# Save end position for next hit
			$prev_hit_end = $hit_end;
		}
        # Delete overlapping hits
        foreach my $gr (sort keys %grouped_hits) {
            foreach my $str (sort keys %{$grouped_hits{$gr}}) {
                my %del_hits = ();
                foreach my $i_x ((0 .. @{$grouped_hits{$gr}{$str}}-1)) {
                    # Skip if hit to be deleted
                    if ($del_hits{$i_x}) { next; }
                    # Get properties
                    my $qbeg_x = $grouped_hits{$gr}{$str}[$i_x]->[6];
                    my $qend_x = $grouped_hits{$gr}{$str}[$i_x]->[7];
                    my $qlen_x = $qend_x-$qbeg_x+1;

                    foreach my $i_y ((0 .. @{$grouped_hits{$gr}{$str}}-1)) {
                        # Get actual index of element
                        if ($i_x eq $i_y) { next; }
                        # Skip if hit to be deleted
                        if ($del_hits{$i_y}) { next; }
                        # Get properties
                        my $qbeg_y = $grouped_hits{$gr}{$str}[$i_y]->[6];
                        my $qend_y = $grouped_hits{$gr}{$str}[$i_y]->[7];
                        my $qlen_y = $qend_y-$qbeg_y+1;
                        # Check if regions overlap
                        if ($qend_y >= $qbeg_x && $qbeg_y <= $qend_x) {
                            # Get overlap length
                            my $overlap_len = $qend_y-$qbeg_x+1;
                            # If substential overlap discard shorter region
                            if ($overlap_len > 0) {
                                if ($qlen_y > $qlen_x) {
                                    $del_hits{$i_x} = 1;
                                } else {
                                    $del_hits{$i_y} = 1;
                                }
                            }
                        }
                    }
                }
                # Save non-overlapping hits
                foreach my $i ((0 .. @{$grouped_hits{$gr}{$str}}-1)) {
                    unless ($del_hits{$i}) {
                        push(@{$grouped_hits_clean{$gr}{$str}},$grouped_hits{$gr}{$str}[$i]);
                    }
                }

            }
        }
		# Save hit groups
		foreach my $gr (sort keys %grouped_hits_clean) {
            @{$chr_hit_locs{$chr}{$gr}{'plus'}} = @{$grouped_hits_clean{$gr}{'plus'}} if $grouped_hits_clean{$gr}{'plus'};
            @{$chr_hit_locs{$chr}{$gr}{'minus'}} = @{$grouped_hits_clean{$gr}{'minus'}} if $grouped_hits_clean{$gr}{'minus'};
		}

	}
	return \%chr_hit_locs;
}

sub get_int_loc_poss {
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
        my $chr = $d[0];
        my $beg = $d[1];
        my $end = $d[2];
        # Save data fields
        foreach my $pos ($beg..$end) {
            $data_fields{$chr}{$pos} = 1;
        }
    }
	# Return data fields
	return \%data_fields;
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

sub usage {
	# Take optional message
	my($message) = @_;
	# Usage message
	my $USAGE = "perl $0 -r repeat.fas -g genome.fas [options]\n";
	# Handle message
	if ($message) { $message .= "\n" unless $message =~ /\n$/ }
	if (!$message) { $message = '' }
	# Print usage to stderr
	return ($message, "\nUsage: $USAGE\n");
}

sub help {
	# Print usage and help message
	die(usage(),
		"Options:\n",
		"  -r, --repeat [file]     repeat (query) fasta file\n",
        "  -g, --genome [file]     genome (subject) fasta file\n",
		"  -t, --target [file]     target loci of interest bed file\n",
		"  -s, --sensitive         search with more sensitive algorithm\n",
        "  -m, --min_ident [int]   minimum sequence %-identity cut-off\n",
		"  -h, --help              show help message\n",
		"\n",
        "  int = integer\n",
        "\n",
	);
}

################################ documentation ################################

=pod

=head1 DESCRIPTION

This script requires at least two parameters. The name of a query
sequences fasta file (-r; e.g. a repeat) and the name of a corresponding
subject sequence fasta file (-g; e.g. a genome). Additional options:
A bed file with loci of interest in the subject sequence (-t), for which
hits will appear at the top of the output. The -s flag determines that a
more sensitive search algorithm, namely discontiguous megablast, should
be used. With the -m option, you can choose a minimum threshold value
for the percent sequence identity of hits to be reported. Output files
include tabular and alignmnt blast files and a bed file of the final
query hit groups. If the Perl module PDF::Create is installed, an
additional graphical output is created in a pdf format.

You can install Perl modules with CPAN. See
http://www.cpan.org/modules/INSTALL.html

=cut
