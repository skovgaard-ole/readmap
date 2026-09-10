#!/usr/bin/perl
# Printing ONP reads:
# #\x1F#	The complementary fasta files are mapped to a reference with blastn fx:
# #\x1F#	"blastn -num_threads 8 -db P19A_chr_plsm -query D:\P19A_Nanopore_reads\P19A_all.fa -out P19A_2.btop -outfmt "6 delim=	 qseqid sframe qstart qend sstart send qlen sseqid btop"
#
# C:\scripts\read_print_23.pl AP027148.gb  DRR325755.btop -x 10
#
use strict;
use warnings;
use warnings FATAL => 'uninitialized';
use POSIX qw(strftime);
use IO::Handle;
use Getopt::Long qw(GetOptions);
use File::Basename;

# "Global variables:"
my ($read_tot,$single_read_count, $long_div_read, $short_div_read, $arrow_count) = (0, 0, 0, 0, 0);
my ($btop, $sample,  $rs_name, $ref_seq_type, $cov);
my (%seqs, %seq_length, $seq_length_total, %mapstring_plus, %mapstring_minus,%cont_cov);
my (%bt_r_single, %bt_r_multi_short, %bt_r_multi_long, %map_pages);
my (@seqs,@extract_anno);
our (%read_extract, %btop_acc);

# Define parameters (only numeric parameters)
my @PARAM_META = (    
    { key => 'min_read_length',    cli => '-r  --min-read-length',    help => 'Minimum read length',    fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'min_match_length',   cli => '-m  --min-match-length',   help => 'Minimum match length',   fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'micro_match_length', cli => '-y  --micro-match-length', help => 'Micro match length',     fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'cov_max',            cli => '-x  --cov-max',            help => 'Max coverage',           fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'map_scale',          cli => '-u  --map-scale',          help => 'Map scale',              fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'page_scale',         cli => '-p  --page-scale',         help => 'Page scale',             fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'space',              cli => '-s  --space',              help => 'Spacing',                fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'line_space',         cli => '-c  --line-space',         help => 'Line spacing',           fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'space2reads',        cli => '-v  --space2reads',        help => 'Axis spacing',           fmt => sub { sprintf "%8d", $_[0] } },
    { key => 'l_w',                cli => '-l  --l-w',                help => 'Line width',             fmt => sub { sprintf "%8d", $_[0] } },
);

# GenBank Variables:
my ($gb_seq, %gb_annot, $seq_name, $organism, $bio_project, $bio_sample, $sra, $asm_meth, $seq_tech);

# extracted annotations
our @ext;

my $platform;
if ($^O eq 'MSWin32') {$platform = "Windows";} else {$platform = "Unix-like";}  
my @settings = @ARGV;

# Hash of PostScript FileHandles:
our %PSFH;

#	Define settings
our %color;
our %arrow;

define_PS_settings();

logmsg("Begin mapping");
my $star = "********************************************************************************";
print "\n$star\n***   Read_Print.pl - Visualizing Long Sequence Reads Mapped by BLAST BTOP   ***\n$star\n\n";
print "System: $^O - $platform\n";

our (%opt, %opt_default);

READ_PARAMETERS();

# Fixed PS variables, depends on READ_PARAMETERS:
my $axis_width = 10000 * $opt{page_scale};  # length of x-axis in "PS_pixels"; used in many SUBs
my $k_max = 3700/$opt{line_space} - 8;	# setting max lines up or down with reads; used in many SUBs
my $front = "front";

# global debug handles:
our %debug;
OPEN_DEBUG() if $opt{debug};

my $txt = logmsg("Begin mapping");
print {$debug{log}} "System: $^O - $platform;\n$txt\nCommand line: $0 @settings\n" if $debug{log};

# Use IO::Uncompress::Gunzip if installed, else fallback to gzip -dc
my $have_gunzip = 0;
eval {
	require IO::Uncompress::Gunzip;
	IO::Uncompress::Gunzip->import();
	$have_gunzip = 1;
	1;
} or do { $have_gunzip = 0 };

OPEN_FILES();

PS_init();

Init_mapstrings();

Print_annot() if $ref_seq_type eq "gb";

Sort_reads($btop);

## --> if extract reads:
if ( @{ $opt{extract_reads} } ) {
	for my $args ( @{ $opt{extract_reads} } ) {
		extract_reads(@$args);
	}
}
## <--if extract reads:

Print_reads();

PS_finish();

Report_Analysis() if $opt{report_analysis};

logmsg("Finish");


# ============================================================================
# SUBROUTINES
# ============================================================================


# READ_PARAMETERS begin
sub READ_PARAMETERS {
	
	use constant NO_LIMIT          => 0;

	%opt = (
	min_read_length    => 10000,   # no minimum
	min_match_length   => 2000,
	micro_match_length => 100,
	cov_max            => NO_LIMIT,   # unlimited
	map_scale          => 100,   # no scaling
	page_scale         => 1,
	space              => 100,
	line_space         => 15,
	space2reads        => 70,
	l_w                => 2,   # default width
	extract_reads      => [],   # off
	report_analysis    => [],   # off
	debug              => 0,   # off
	help               => 0,   # off
	out                => '',  # auto
	);

	%opt_default = %opt;
	$opt_default{extract_reads} = [];   # break shared reference

	GetOptions(
	"out|o=s"               => \$opt{out},
	"min-read-length|r=i"   => \$opt{min_read_length},
	"min-match-length|m=i"  => \$opt{min_match_length},
	"micro-match-length|y=i"=> \$opt{micro_match_length},
	"cov-max|x=i"           => \$opt{cov_max},
	"map-scale|u=f"         => \$opt{map_scale},
	"page-scale|p=f"        => \$opt{page_scale},
	"space|s=i"             => \$opt{space},
	"line-space|c=i"        => \$opt{line_space},   
	"space-to-reads|v=i"    => \$opt{space2reads},
	"l-w|l=i"               => \$opt{l_w},
	
	"extract-reads|e=s" => sub {
	    my ($opt_name, $val) = @_;
	    my @vals = split /,/, $val;
	    push @{ $opt{extract_reads} }, \@vals;
	},
	"report_analysis|a"     => \$opt{report_analysis},
	"debug|d"               => \$opt{debug},
	"help|h"                => \$opt{help},
	) or usage("Invalid options");

	for my $k (keys %opt_default) {
		$opt{$k} = $opt_default{$k} unless defined $opt{$k};
	}
	usage() if $opt{help};

	# ---------- positional args ----------

	($opt{seq_file}, $opt{btop_file}) = @ARGV;

	unless ($opt{seq_file} && $opt{btop_file}) {
		usage("Missing required arguments for sequence file and for btop file");
	}

	unless (-e $opt{seq_file}) {
		die "ERROR: seq_file not found: $opt{seq_file}\n";
	}

	unless (-e $opt{btop_file}) {
		die "ERROR: btop_file not found: $opt{btop_file}\n";
	}

	# ---------- defaults ----------

	if (!$opt{out}) {
		my $base = basename($opt{btop_file});
		$base =~ s/\.\w+$//;
		$opt{out} = $base;
	}

	# ---------- debug ----------

	if ($opt{debug}) {
		print "\n--- PARAMETERS ---\n";
		printf "%-20s  %15s  :  %-15s\n", "Parameter name","Default value","Set value";

	for my $p (@PARAM_META) {
	    my $k = $p->{key};

	    printf "%-20s  %15s  :  %-15s\n",
	        $k,
	        fmt_opt($opt_default{$k}),
	        fmt_opt($opt{$k});
	}
		print "File with reference sequence           :  $opt{seq_file}\n";
		print "File with btop mapped reads            :  $opt{btop_file}\n";
		print "Output prefix                          :  $opt{out}\n";
		print "------------------\n";
		# $debug = "debug";
	}
	return;
	#return ($seq_file, $opt{btop_file}, \%opt);
}

sub fmt_opt { #helper for de-referencing arrays, used in READ_PARAMETERS 
    my ($v) = @_;
    return "" unless defined $v;
    if (ref $v eq 'ARRAY') {
         # return  map { join(",", @$_) } @$v; 
         return join(" | ", map { join(",", @$_) } @$v);
         # Known bug here: value in $opt_default{$extract_reads} is changed to the value value in $opt{$extract_reads}
    }
    return $v;
}
# READ_PARAMETERS end

#	OPEN_FILES begin
sub OPEN_FILES {
	if ($ARGV[0]) {
		open(my $rs, "<", $opt{seq_file}) or die "Unable to open the Reference Sequence: $!";
		$rs_name = <$rs>;
		print {$debug{log}} "\nFirst line reads:\n$rs_name" if $debug{log};				
		if (substr($rs_name, 0, 1) eq ">") {
			close($rs);
			readfa($opt{seq_file});
			$ref_seq_type = "fa";
		}
		if (substr($rs_name, 0, 6) eq "LOCUS ") {
			close($rs);
			readgb($opt{seq_file});
			annotation_extract();
			$ref_seq_type = "gb";
		}
	}

	while ($opt{btop_file} !~ /btop/) { # fall-back open of .btop file
		print "Input file.btop with blast mapped reads to print \n\\tFile extension must be \".btop\"\n";
		$opt{btop_file} = <STDIN>;
		chomp $opt{btop_file};
	}
	$opt{btop_file} =~ /(\w+)\.btop/;
	$sample = $1;
	unless($opt{out}) { $opt{out} = $sample }

	@seqs = sort { length($seqs{$b}) <=> length($seqs{$a}) } (keys %seqs);

	foreach (@seqs) {	
		$map_pages{$_} = int(	($seq_length{$_} / $opt{map_scale} + $opt{map_scale}) / $axis_width	);
		print "\tContig: $sample; $_ is ".sprintf("%10d",$seq_length{$_})." bases long \n";
	}

	open($btop, "<", $opt{btop_file}) or die "\n*****\n**   The btop file is not accessible\n*****\n";
}
#	OPEN_FILES end

# Read reference sequence begin
sub readgb
{
	my ($file) = @_;
	open(my $gb, "<", $file) or die "Cannot open GenBank file: $!";
	my @gb_all = <$gb>;
	close($gb);

	# Extract header metadata (first 100 lines)
	for my $i (0 .. 100) {
		if ($i >= @gb_all) { last; }
		if ($gb_all[$i] =~ /ORGANISM/) { $organism = substr($gb_all[$i], 12); }
		if ($gb_all[$i] =~ /BioProject:/) { $bio_project = substr($gb_all[$i], 12); $bio_project =~ s/:/:             /; }
		if ($gb_all[$i] =~ /BioSample:/) { $bio_sample = substr($gb_all[$i], 12); $bio_sample =~ s/:/:              /; }
		if ($gb_all[$i] =~ /Sequence Read Archive:/) { $sra = substr($gb_all[$i], 12); $sra =~ s/:/:  /; }
		if ($gb_all[$i] =~ /Assembly Method/) { $asm_meth = substr($gb_all[$i], 12); }
		if ($gb_all[$i] =~ /Sequencing Technology/) { $seq_tech = substr($gb_all[$i], 12); $seq_tech =~ s/[()]/ /;     }
	}

	print "\n--- Reference sequence information: ---\nORG: $organism";
	print "B_P: $bio_project";
	print "B_S: $bio_sample";
	print "SRA: $sra";
	print "ASM: $asm_meth";
	print "SEQ: $seq_tech\n";

	# Parse sequences and annotations (multiple sequences per file)
	for my $i (0 .. $#gb_all) {
		if (substr($gb_all[$i], 0, 6) eq "LOCUS ") {
			# New sequence record
			$gb_all[$i] =~ /\s+(\S+)/;
			$seq_name = $1;
			$gb_seq = "";
			next;
		}

		if ($gb_all[$i] =~ /^origin/i) {
			$gb_seq = 1;
			next;
		}

		if (!$gb_seq) {
			$gb_annot{$seq_name} .= $gb_all[$i];
		}
		
		if ($gb_seq) {
			$seqs{$seq_name} .= $gb_all[$i];
		}
	}

	# Process all sequences: clean and calculate lengths
	for (keys %seqs) {
		$seqs{$_} =~ s/[^acgtACGTnN]//g;
		$seq_length{$_} = length($seqs{$_});
		$seq_length_total += length($seqs{$_});
	}

	# Sort sequences by length (descending)
	@seqs = sort { length($seqs{$b}) <=> length($seqs{$a}) } (keys %seqs);
	my $seqs = join "\n\t", @seqs;
	print {$debug{log}} "\nThese Accession numbers were found in the GenBank file:\n\t$seqs\n" if $debug{log};		

	return;
}

sub annotation_extract {
	@extract_anno = ();
	foreach my $key (@seqs) {
		my $text = $gb_annot{$key} // next;
		my @lines = split /\n/, $text;
		my @anno;
		for my $line (@lines) {
			# start of a new gene feature
			if ($line =~ /^\s{5}gene\b/) {
				# process previous feature
				if (@anno) {
					my $ext = extanno(@anno);
					push @extract_anno, "$key\t$ext" if defined $ext;
				}
				@anno = ($line);
				next;
			}
			# accumulate feature block
			push @anno, $line if @anno;
		}
		# last feature
		if (@anno) {
			my $ext = extanno(@anno);
			push @extract_anno, "$key\t$ext" if defined $ext;
		}
	}
	warn "Extracted ".scalar(@extract_anno)." annotations from $opt{seq_file}:\n" if $opt{debug};
	return;
}

sub extanno {

	my @anno = @_;
	my @result;
	my $l_t;

	# Merge feature block into one string (handles multiline qualifiers)
	my $qual_text = join "\n", @anno;

	# continuation lines (GenBank column ≥21)
	$qual_text =~ s/\n\s{21,}/ /g;

	# Strand
	$result[0] = ($qual_text =~ /complement/) ? 'c' : '';

	# Coordinates (handles join/complement/< > etc.)
	my @coords = ($qual_text =~ /[<>]?\s*(\d+)\s*\.\.\s*[<>]?\s*(\d+)/g);

	if (@coords) {
		$result[1] = min(@coords);
		$result[2] = max(@coords);
	}

	# Feature type
	if ($qual_text =~ /^\s*(tRNA|rRNA|RNA)/m) {
		$result[3] = $1;
		$result[5] = 'RNA';
	}
	elsif ($qual_text =~ /^\s*(CDS)/m) {
		$result[3] = $1;
	}

	# Gene / locus_tag
	$result[4] = $1
	if $qual_text =~ m{/gene\s*=\s*"([^"]+)"};

	$l_t = $1
	if $qual_text =~ m{/locus_tag\s*=\s*"([^"]+)"};

	# fallback
	$result[4] //= $l_t // '';

	# Extract PRODUCT only (avoid translation noise)
	my $product = '';

	$product = $1
	if $qual_text =~ m{/product\s*=\s*"([^"]+)"}s;

	# Classification (case-sensitive, safer)
	$result[5] = 'phage'
	if $product =~ /\b(phage|tail)\b/;

	$result[5] = 'hyp'
	if $product =~ /\bhypothetical\b/;

	$result[5] = 'IS'
	if $product =~ /\b(transposase|recombinase)\b/;

	$result[5] = 'pseudo'
	if $qual_text =~ /\/pseudo\b/;

	$_ //= '' for @result;
	# --------------------------------------------------
	return join "\t", @result;
}

sub readfa {
	my ($file) = @_;

	open(my $fasta, "<", $file) or die "Cannot open FASTA file: $!";
	my @fasta = <$fasta>;
	close($fasta);

	for my $i (0 .. $#fasta) {
		if (substr($fasta[$i], 0, 1) eq ">") {
			if ($fasta[$i] =~ /\|/) {
				$fasta[$i] =~ /\|(\S+)/;
				$seq_name = $1;
			} else {
				$fasta[$i] =~ />(\S+)/;
				$seq_name = $1;
			}
			next;
		}
		$seqs{$seq_name} .= $fasta[$i];
	}

	for (keys %seqs) {
		$seqs{$_} =~ s/\s//g;
		$seq_length{$_} = length($seqs{$_});
		$seq_length_total += length($seqs{$_});
	}

	@seqs = sort { length($seqs{$b}) <=> length($seqs{$a}) } (keys %seqs);
	my $seqs = join "\n\t", @seqs;
	print {$debug{log}} "\nThese Accession numbers were found in the Fasta file:\n\t$seqs\n" if $debug{log};	
	
	return;
}
# Read reference sequence end

# Processing btop reads begin
sub Sort_reads {
	my ($btop_fh) = @_;
	my $cm = 0;
	my $bt_prev = "";
	my @bt_read = ();
	my $short = 0;
	$cov = "";

	while (my $line = <$btop_fh>) {
		chomp $line;
		my @btop_read = split("\t", $line);

		# eliminate short reads
		if ($btop_read[6] < $opt{min_read_length}) { $short++; next; }

		#index for extraction of btop reads:
		my $win1 = int(min($btop_read[4], $btop_read[5])/1000)+1;
		my $win2 = int(max($btop_read[4], $btop_read[5])/1000);

		if ($win2 - $win1 > $opt{min_match_length}/1000) {
			for my $i ($win1 .. $win2) {
				push @{ $read_extract{$btop_read[7]}{$i} }, {
					read  => $btop_read[0],
					start => $btop_read[4],
					end   => $btop_read[5],
				}
			}
		}

		$btop_read[0] =~ /\w+[\._](\d+)/;
		$btop_read[0] = "R_" . $1;

		if ($btop_read[3] - $btop_read[2] > 100) {    # $opt{micro_match_length}) { 	# "micro" bc. also small inversions must pass here
			# New read encountered -> process accumulated @bt_read
			if ($bt_prev ne $btop_read[0] && @bt_read) {
				if ($opt{cov_max}) {
					$cm += ($btop_read[6] / $seq_length_total);
					if ($cm > $opt{cov_max}) { $cov = " (limited by cov_max)"; return; }
				}
				if ($debug{bt_read}) {
					for my $i (0 .. $#bt_read) {
						print {$debug{bt_read}} "bt_$i  \t$bt_read[$i]\n";
					}
				}
				
				_process_bt_read(@bt_read);
				@bt_read = ();
			}

			# push current match into @bt_read
			my $bt_read = join "\t", @btop_read;       
		  push @bt_read, $bt_read;
		}
		
		$btop_acc{$btop_read[7]} = 1;
		$bt_prev = $btop_read[0];
	}
	
	my @accs = sort keys %btop_acc;
	my $seqs = join "\n\t", @accs;
	print {$debug{log}} "\nThese Accession numbers were found in the BTOP file:\nThey must match the Accession numbers from the GB/Fasta file\n\t$seqs\n\n" if $debug{log};	
	
	# process any remaining accumulated @bt_read after EOF
	if (@bt_read) {
		_process_bt_read(@bt_read);
	}

	return;
}

sub _process_bt_read {
	my (@bt_read) = @_;
	return unless @bt_read;

	if ($#bt_read == 0) {
		my @single_read = split "\t", $bt_read[0];
		push @{$bt_r_single{$single_read[7]}}, $bt_read[0];
		print {$debug{single_read}} "SR: \t$bt_read[0];\n"if $debug{single_read};
		return;
	}

	# Multi-read processing
	my @taken = ();
	my $m_rd_pr = 0;
	my $m_rd_start = 0;
	my @multi_rd_pr = ();

	for my $j (0 .. $#bt_read) {
		my @multi_rd = split "\t", $bt_read[$j];
		if ($j == 0) { $m_rd_start = $multi_rd[4] }
		$multi_rd[7] = abs($multi_rd[2] - $multi_rd[3]);
		if ($j >= 2 && abs($m_rd_pr - $multi_rd[7]) < 50) {
			if (abs($multi_rd_pr[4] - $m_rd_start) > abs($multi_rd[4] - $m_rd_start + $multi_rd[6])) {
				$bt_read[$j-1] = $bt_read[$j];
			}
		}
		$m_rd_pr = $multi_rd[7];
		@multi_rd_pr = @multi_rd;
	}

	my @multi_quest = ();
	my $tol = $opt{micro_match_length} / 4; # tolerance for repeats allowed before taken; there is no exact way to solve this.
	for my $j (0 .. $#bt_read) {
		my @multi_rd = split "\t", $bt_read[$j];
		my $overlap = 0;
		for my $k ($multi_rd[2] .. $multi_rd[3]) {
			if ($taken[$k]) { $overlap = 1; last; }
		}
		unless ($overlap) {
			for my $k ($multi_rd[2]+$tol .. $multi_rd[3]-$tol) { $taken[$k] = 1; }
			push @multi_quest, $bt_read[$j];
			print {$debug{to_multiquest}} "$bt_read[$j]\n" if $debug{to_multiquest};
		}
	}
	_process_multi_quest(@multi_quest);
}

sub _process_multi_quest {
	my (@multi_quest) = @_;

	print {$debug{into_multiquest}} "In: $#multi_quest\n@multi_quest\n\n" if $debug{into_multiquest};

	if ($#multi_quest == 0) {
		my @multi_quest_read = split "\t", $multi_quest[0];
		push @{$bt_r_single{$multi_quest_read[7]}}, $multi_quest[0];
		print {$debug{single_read}} "M: $#multi_quest\t$multi_quest[0]\n" if $debug{single_read};
		return;
	}

	my @map_ends = ();
	my @rd_ends = ();
	my $mq_1 = "";

	for my $j (0 .. $#multi_quest) {

		print {$debug{multi_multiquest}} "$j of $#multi_quest:\t$multi_quest[$j]\n" if $debug{multi_multiquest};

		my @mq = split "\t", $multi_quest[$j];
		if ($j == 0) { $mq_1 = $mq[7] }
		push @map_ends, ($mq[4], $mq[5]);
		push @rd_ends, ($mq[2], $mq[3]);

		if ($j == $#multi_quest) {

			my $map_min = min(@map_ends);
			my $map_max = max(@map_ends);
			my $rd_min  = min(@rd_ends);
			my $rd_max  = max(@rd_ends);
			my $l = abs(($map_max - $map_min) - ($rd_max - $rd_min));  # $l is small for short divided and big for long divided reads
			my $circle = abs($seq_length{$mq[7]} - ($map_max - $map_min)) < 20 ? "Circle"  : "C";

			print {$debug{into_multiquest}} "$mq[0]; L: $l; Mx: $map_max; Mn: $map_min; Rx: $rd_max; Rn: $rd_min\n" if $debug{into_multiquest};

			for my $k (1 .. $#multi_quest) {
				my @mq_check = split "\t", $multi_quest[$k];
				if ($mq_check[7] ne $mq_1) { $circle = "Contig_Join"; }
			}

			my $ll = $opt{min_match_length};  # an arbitrary limit separating short and long divided reads; "Distance" between map positions

			if ($l < $ll) {
				$short_div_read++;
				for my $k (0 .. $#multi_quest) {
					$multi_quest[$k] .= "\t" . $circle;
					my @mq = split "\t", $multi_quest[$k];
					push @{$bt_r_multi_short{$mq[7]}}, $multi_quest[$k];
					print {$debug{short_read}} "\n$mq[0]; Distance $l; mm:$map_min; $map_max; $rd_min; $rd_max; $circle\n" if ($debug{short_read} && $k ==0);

				}
			}
			elsif ($l >= $ll) {
				$long_div_read++;
				for my $k (0 .. $#multi_quest) {
					$multi_quest[$k] .= "\t" . $circle;
					my @mq = split "\t", $multi_quest[$k];
					if (abs($mq[4] - $mq[5]) > $opt{min_match_length}) {
						push @{$bt_r_multi_long{$mq[7]}}, $multi_quest[$k];
						# print $lmr "M>>>:$k:  $multi_quest[$k];\n\tL $l; mm:$map_min; $map_max; $rd_min; $rd_max\n";
						print {$debug{long_read}} "\n$mq[0]; Distance $l; mm:$map_min; $map_max; $rd_min; $rd_max; $circle\n"	if ($debug{long_read} && $k ==0);
						print {$debug{long_read}} "\t\tMap_No: $k\t$multi_quest[$k];\n" if $debug{long_read};
					}
				}
			}
		}
	}
}
# Processing btop reads end

# Postscript printing begin
sub PS_init {
	$map_pages{front} = 0;
	foreach my $key (@seqs, "front") {
		for my $i (0 .. $map_pages{$key}) {
			my $p = sprintf("Ps-page_%s_%03d", $key, $i);
			my $file = "$p.ps";
			open my $fh, '>', $file
			or die "Cannot open $file: $!";
			$PSFH{$p} = $fh;   # store file handle in hash; filename omitting .ps

			my ($a,$b) = ($i*$axis_width*$opt{map_scale},($i+1)*$axis_width*$opt{map_scale});
			my $b1 = min($b,$seq_length{$key});
			my $axis_end = min($axis_width,($b1-$a)/$opt{map_scale});
			my $t = ($b1-$a)/$opt{map_scale};
			while($a  =~ s/(\d+)(\d\d\d)/$1\,$2/){}; # a very elegant way to make 1000 separators!
			while($b  =~ s/(\d+)(\d\d\d)/$1\,$2/){}; # a very elegant way to make 1000 separators!
			while($b1 =~ s/(\d+)(\d\d\d)/$1\,$2/){}; # a very elegant way to make 1000 separators!
			my $ae_mark_no = $b1 ne $b ? "			newpath $axis_end 30 add -6 moveto ($b1) show" : "";
			my ($x,$y) = (10750 * $opt{page_scale},7600 * $opt{page_scale}); 
			my $mid = 0.5*$y;
			my $page_no = $i+1;
			my $page_of = $map_pages{$key} + 1;

			# Arrowheads:
			my $ahx = $opt{l_w}*2.5;
			my $ahy = $ahx*0.5;
			my $print_arrow = "\n	gsave\n	$opt{l_w} setlinewidth
			20 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 1 setlinejoin stroke
			grestore \n" ;

			my $print_arrow_hd = "\n	gsave\n	$opt{l_w} setlinewidth
			2 0 rmoveto $ahx 1 sub 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 1 setlinejoin stroke
			grestore \n" ;

			my $print_arrow_hd_rv = "\n	gsave\n	$opt{l_w} setlinewidth
			2 0 rmoveto -$ahx 1 add 0 rlineto $ahx $ahy rlineto -1 -$ahy rlineto 1 -$ahy rlineto -$ahx $ahy rlineto 1 setlinejoin stroke
			grestore \n" ;

			my $print_arrow_ct = "\n	gsave\n	$opt{l_w} setlinewidth
			20 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 12 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 1 setlinejoin stroke
			grestore \n" ;

			#	% all arrow_inv is obsulute for now:
			my $print_arrow_inv = "\n	gsave\n	$opt{l_w} setlinewidth
			20 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 12 0 rlineto -$ahx $ahy rlineto 1 -$ahy rlineto  -1 -$ahy rlineto $ahx $ahy rlineto 1 setlinejoin stroke
			grestore \n" ;

			my $print_x_axis =  $key ne ("front")? "			 newpath 0 0 moveto $axis_end 0 lineto 1 setlinecap stroke\n			 $ae_mark_no\n" : '';
			my $print_title  =  $key ne ("front")? "			$color{8}\n			newpath 250 $opt{page_scale} mul $y 70 sub moveto (Mapping $sample ONP reads on: $key.      Page $page_no of $page_of pages.      From base $a to base $b1.) show\n  ": '';

			print $fh "%!PS-Adobe-3.0
			<< /PageSize [$x $y] /Orientation 0 >> setpagedevice

			% Frame
			$color{9}
			$opt{l_w} 3 mul setlinewidth
			newpath 0 $y   moveto $x $y lineto 1 setlinecap stroke
			newpath 0 0     moveto $x 0 lineto 1 setlinecap stroke
			newpath 0 0     moveto 0   $y lineto 1 setlinecap stroke
			newpath $x 0   moveto $x   $y lineto 1 setlinecap stroke

			% Font:
			/Courier findfont   % Get the basic font
			24 scalefont            % Scale the font to X points
			setfont                 % Make it the current font

			% Title
			$print_title
			250 $opt{page_scale} mul $mid translate
			$color{7}
			$opt{l_w} setlinewidth


			% Base, x_axis
			$print_x_axis

			$opt{l_w} 2 div setlinewidth

			% Arrows:
			/arrow {$print_arrow} def
			/arrow_ct {$print_arrow_ct} def
			/arrow_hd {$print_arrow_hd} def
			/arrow_hd_rv {$print_arrow_hd_rv} def
			/arrow_inv {$print_arrow_inv} def
			/arrow_circular {$print_arrow
				gsave
				$opt{l_w} 4 mul setlinewidth
				30 0 rmoveto 0 0 rlineto 1 setlinecap stroke
				grestore
				gsave
				$opt{l_w} 2 mul setlinewidth
				$color{9}
				30 0 rmoveto 0 0 rlineto 1 setlinecap stroke
				grestore
			} def

			";

			# Print x-axis - except for the front page:
			unless ($key eq "front") {
				my $sfine = 1000/$opt{map_scale}; # distance betw marks
				my $fine = 10000/$opt{map_scale}; # distance betw marks
				my $coarse = 50000/$opt{map_scale}; # distance betw marks
				my $mark_a = $axis_end/$sfine;
				my $mark_b = $axis_end/$fine;
				my $mark_c = $axis_end/$coarse;

				print $fh "			$opt{l_w} 6 div setlinewidth\n";
				for my $j (1 .. $mark_a ) {
					my $mark = $j*$sfine;
					print $fh "		newpath $mark  0   moveto $mark   3 lineto  1 setlinecap stroke\n";
				}

				print $fh "			$opt{l_w} 2 div  setlinewidth\n";
				for my $j (1 .. $mark_b ) {
					my $mark = $j*$fine;
					print $fh "		newpath $mark  -4   moveto $mark   4 lineto  1 setlinecap stroke\n";
				}
				print $fh "			$opt{l_w} setlinewidth\n";
				for my $j (0 .. $mark_c ) {
					my $mark = $j*$coarse;
					my $mark_no = $i*$opt{map_scale}*$axis_width + $j*$opt{map_scale}*$coarse;
					while($mark_no =~ s/(\d+)(\d\d\d)/$1\,$2/){}; # a very elegant way to make 1000 separators!
					print $fh "newpath $mark  -6   moveto $mark   6 lineto  1 setlinecap stroke\n";
					print $fh "newpath ($mark_no) dup stringwidth pop 2 div $mark exch sub 20 moveto show \n";
				}
			}
			# Notes on center text: dup; dublicate stack, stringwith returns x & y of with; pop deletes the y; 2 div halves (for centering) give the x-base (here $mark);
			# contd: ech: exchange the two values in stack before substracting what is now the half of the length of the text string.
			print $fh "
			% Font:
			/Courier findfont   % Get the basic font
			12 scalefont            % Scale the font to X points
			setfont                 % Make it the current font
			$opt{l_w} 3 div setlinewidth
			";
		}
	}
	return();
}

sub PRINT_parameters {	# Printing PS frontpage report - begin.

	my $page_x = 750;
	my $page_y = 3000;
	my $a_scale = 2;

	my $add1 = 25  * $opt{page_scale};
	my $add2 = $add1 / $a_scale;
	my $add3 = 50  * $opt{page_scale};
	my $add4 = 200 * $opt{page_scale};
	my $add5 = 500 * $opt{page_scale};
	my $add6 = 120 * $opt{page_scale};
	my $add7 = 675 * $opt{page_scale};
	my $add9 = 1000 * $opt{page_scale};
	my $line_sp = -30;
	my $line_no = 0;

	my @p_print;

	# ---------- helper call ----------

	my $ps = PS(line_sp => -30);

	# ---------- header ----------

	$ps->push("%\tPrint parameters & Annotation");
	$ps->push("R\t$page_x $page_y translate");
	$ps->push("R\t3 3 scale");

	$ps->font('Helvetica',72);
	$ps->color($color{2});
	$ps->text(700*$opt{page_scale} ,"ReadMap Analysis of $ARGV[1]");
	$ps->color($color{8});


	# ---------- parameters ----------

	$ps->{line_no} = 10;

	$ps->title("Printing Parameters:");
	$ps->hr($add9);

	for my $p (@PARAM_META) {
	    my $val = $opt{ $p->{key} };

	    my $formatted = $p->{fmt}
	        ? $p->{fmt}->($val)
	        : $val;

		my $default = $opt_default{ $p->{key} };

		$ps->row2(
		    $add1,
		    $add5,
		    $p->{cli},
		    "$formatted   (def: $default)"
		);
	}


	# ---------- Read Lines ----------

	$ps->nl for 1..3;

	$ps->title("Read Lines:");
	$ps->hr($add9);

	$ps->push("V\tgsave");
	$ps->push("V\t$opt{l_w} 2 mul setlinewidth");

	$ps->color($color{1});
	$ps->push("L\t$add1 ".$ps->YL." $add6 ".$ps->YL);
	$ps->text($add4,"Undivided reads");

	$ps->color($color{2});
	$ps->push("L\t$add1 ".$ps->YL." $add6 ".$ps->YL);
	$ps->text($add4,"Divided reads, short distance");

	$ps->color($color{3});
	$ps->push("L\t$add1 ".$ps->YL." $add6 ".$ps->YL);
	$ps->text($add4,"Divided reads, short distance, reverse direction");

	$ps->color($color{4});
	$ps->push("L\t$add1 ".$ps->YL." $add6 ".$ps->YL);
	$ps->text($add4,"Divided reads, long distance");

	$ps->push("V\tgrestore");

	# ---------- Arrows ----------

	$ps->nl for 1..3;

	$ps->title("Arrows:");
	$ps->hr($add9);

	my $line_no_restore = $ps->{line_no};

	$ps->push("V\tgsave");
	$ps->push("V\t$a_scale $a_scale scale");

	# forward arrows
	$ps->push("V\t$color{1} $add2 ".($ps->Y/$a_scale)." moveto $arrow{hd}");
	$ps->push("V\t$color{4} ".($add2+20)." ".($ps->Y/$a_scale)." moveto $arrow{hd}");
	$ps->push("V\t$color{2} ".($add2+40)." ".($ps->Y/$a_scale)." moveto $arrow{hd}");
	$ps->nl;

	# reverse arrows
	$ps->push("V\t$color{1} $add2 ".($ps->Y/$a_scale)." moveto $arrow{hd_rv}");
	$ps->push("V\t$color{4} ".($add2+20)." ".($ps->Y/$a_scale)." moveto $arrow{hd_rv}");
	$ps->push("V\t$color{2} ".($add2+40)." ".($ps->Y/$a_scale)." moveto $arrow{hd_rv}");
	$ps->nl;

	# special arrows
	$ps->push("V\t$color{4} $add2 ".($ps->Y/$a_scale)." moveto $arrow{nm}");
	$ps->nl;

	$ps->push("V\t$color{4} $add2 ".($ps->Y/$a_scale)." moveto $arrow{ct}");
	$ps->nl;

	$ps->push("V\t$color{4} $add2 ".($ps->Y/$a_scale)." moveto $arrow{circular}");
	$ps->nl;

	# restore scale
	$ps->push("V\t1 $a_scale div 1 $a_scale div scale");

	# rotated arrow
	$ps->push("V\tgsave $color{3} newpath $add1 5 add ".$ps->Y." moveto -30 rotate 2 2 scale $arrow{nm} grestore");
	$ps->nl;

	# restore vertical position
	$ps->{line_no} = $line_no_restore;

	$ps->color($color{8});

	$ps->text($add4,"Arrow - this read continues, forward direction");
	$ps->text($add4,"Arrow - this read continues, reverse direction");
	$ps->text($add4,"Arrow - this read joins to an unexpected position");
	$ps->text($add4,"Arrow - this read joins to an end of another contig");
	$ps->text($add4,"Arrow - this read joins to the begin/end of the same contig");
	$ps->text($add4,"[this indicates circularity]");
	$ps->text($add4,"Indicator for inverted sequence");

	# ---------- annotations legend ----------

	$ps->nl for 1..3;
	$ps->title("Annotations:");
	$ps->hr($add9);

	$ps->text($add1,"GenBank annotations are colored according to their");
	$ps->text($add1,"probability of causing genomic instability:");

	for (
	[$color{10},"Gene, other"],
	[$color{11},"* RNA coding gene"],
	[$color{14},"# Phage related gene"],
	[$color{15},"+ Transposase related gene"],
	[$color{12},"^ Pseudogene"],
	[$color{13},"? Hypothetical gene"],
	){
		my ($c,$txt)=@$_;
		$ps->push("V\t$c");
		$ps->text($add3,$txt);
	}

	# ---------- project info ----------

	# Column shift:
	$ps->{line_no} = 0;
	$ps->push("V	1500 $opt{page_scale} mul 0 ".$ps->YL." sub translate ");
	$ps->color($color{10});
	$ps->{line_no} = 10;

	$ps->title("Project and Assembly Information:");
	$ps->hr($add9);

	if ($ref_seq_type eq 'gb') {
		$ps->text($add1,"Organism:                $organism");
		$ps->text($add1,$bio_project);
		$ps->text($add1,$bio_sample);
		$ps->text($add1,$sra);
		$ps->text($add1,$asm_meth);
		$ps->text($add1,$seq_tech);
	}
	else {
		$ps->text($add1,"The reference sequence is in FASTA format, first line reads:");
		$ps->text($add3,$rs_name);
	}

	# ---------- reference lengths ----------

	$ps->nl for 1..3;
	$ps->title("Acc. numbers and lengths of reference sequence:");
	$ps->hr($add9);

	for my $seq (@seqs) {
		my $sl = $seq_length{$seq};
		1 while $sl =~ s/(\d+)(\d\d\d)/$1,$2/;
		$ps->text($add1,sprintf("%-20s","$seq is ").sprintf("%16s",$sl)." bases long");
	}


	# ---------- statistics ----------

	$read_tot =	$single_read_count +	$short_div_read +	$long_div_read;
	$read_tot ||= 1;

	my $si_p = sprintf "%12.2f",$single_read_count*100/$read_tot;
	my $sh_p = sprintf "%12.2f",$short_div_read*100/$read_tot;
	my $lo_p = sprintf "%12.2f",$long_div_read*100/$read_tot;
	my $in_p = sprintf "%12.2f",$arrow_count*100/$read_tot;

	$ps->nl for 1..3;
	$ps->title("Reads mapped statistics:");
	$ps->hr($add9);

	$ps->row2($add1,$add5,"Total reads mapped$cov:", sprintf("%14d",$read_tot));

	for (
	[$color{1},"Undivided reads:",              sprintf("%14d",$single_read_count),$si_p],
	[$color{2},"Divided reads, short distance:",sprintf("%14d",$short_div_read),$sh_p],
	[$color{4},"Divided reads, long distance:", sprintf("%14d",$long_div_read),$lo_p],
	[$color{3},"Reads indicating inversions:",  sprintf("%14d",$arrow_count),$in_p],
	){
		my ($col,$txt,$n,$p)=@$_;

		$ps->color($col);

		$ps->push(
		"3\t$add3\x1F".$ps->Y."\x1F$txt\t".
		"$add5\x1F".$ps->Y."\x1F$n\t".
		"$add7\x1F".$ps->Y."\x1F$p%"
		);

		$ps->nl;
	}

	# ---------- footer ----------

	REPORT_FOOTER_PS($ps);

	# ---------- emit ----------

	my $p = "Ps-page_".$front."_000";

	for my $line (@{ $ps->{p_print} }) {

		my @ppl = split /\t/, ($line // '');
		$#ppl = 3;
		$_ //= '' for @ppl;

		my $tag = $ppl[0];

		if ($tag eq '%') {
			print { $PSFH{$p} } "%\t$ppl[1]\n";
		}
		elsif ($tag eq 'R' or $tag eq 'V') {
			print { $PSFH{$p} } "\t$ppl[1]\n";
		}
		elsif ($tag eq 'L') {
			print { $PSFH{$p} }
			"\tnewpath $ppl[1] moveto $ppl[2] lineto 0 setlinecap stroke\n";
		}
		elsif ($tag =~ /^\d+$/) {

			my $count = $tag;

			for my $j (1..$count) {
				next unless defined $ppl[$j] && length $ppl[$j];
				my ($x,$y,$txt) = split /\x1F/, $ppl[$j], 3;
				$x   //= 0;
				$y   //= 0;
				$txt //= '';
				print { $PSFH{$p} }	"\t$x $y moveto ($txt) show\t";
			}
			print { $PSFH{$p} } "\n";
		}
	}
	print { $PSFH{$p} } "showpage\n";
}	

sub REPORT_FOOTER_PS {

	my ($ps,%args)=@_;

	my $add1 = 25 * $opt{page_scale};

	$ps->color($color{10});
	$ps->nl; $ps->nl;

	$ps->title("Stamp:");
	$ps->hr(1000);

	my $ts = strftime("%Y-%m-%d %H:%M:%S", localtime);

	$ps->text($add1,"Sample: $sample");
	$ps->text($add1,"Analysis time: $ts");

	(my $prog = $0)=~s{\\}{/}g;
	$ps->text($add1,"Program: $prog");

	(my $cmd = $0) =~ s{\\}{/}g;
	$ps->text($add1,"Command line: $cmd @ARGV");

	use Cwd qw(getcwd);
	my $cwd = getcwd();
	$ps->text($add1,"Working directory: $cwd");

	my $host = `hostname`;
	chomp $host;
	$ps->text($add1,"Host: $host");

	$ps->text($add1,"Operating system: $^O");

	$ps->text($add1,"Reference sequence: $opt{seq_file}") if $opt{seq_file};
	$ps->text($add1,"Reads file: $opt{btop_file}")       if $opt{btop_file};

}		#	Printing PS frontpage report - end.


sub Print_annot {	#		Printing annotations (if genbank format ref sequence) - begin

	PS_all ("% Font and Line width:\n	/Courier findfont \n	3 scalefont setfont\n	$opt{l_w} 3 div setlinewidth ");
	my $as = 10; # annotation space
	my ($fr_p,$fr_m) = (0,0); # counter for frame plus or minus
	my ($e_p,$e_m) = (0,0); # counter for frame plus or minus
	my $y; # Y-position
	for my $i (0..$#extract_anno) {

		my ($ct,$strand,$begin,$end,$type,$name,$group)
		= split "\t", $extract_anno[$i], 7;

		# text fields
		$_ //= '' for ($ct,$strand,$type,$name,$group);

		# numeric fields (REAL FIX)
		for ($begin,$end) {
			$_ = 0 unless defined $_ && $_ =~ /^\d+$/;
		}

		$begin /= $opt{map_scale};
		$end   /= $opt{map_scale};

		$name = substr ($name,-5);
		if ($strand) {
			if ($begin-$e_m >	$as){$fr_m = 0 }
			$y = -44 - 6*$fr_m;
			if ($end-$begin < $as) { $fr_m++; $e_m = $end; if ($fr_m == 3) {$fr_m = 0 } }
		}
		else {
			if ($begin-$e_p >	$as){$fr_p = 0 }
			$y = -18 - 6 * $fr_p;
			if ($end-$begin < $as) { $fr_p++; $e_p = $end; if ($fr_p == 3) {$fr_p = 0 } }
		}
		my $sub = ($strand ? -3.5 : 2);
		my $c = $color{10};
		if ($group eq "RNA")    {$c = $color{11};  $name =  $name."*"}
		if ($group eq "pseudo") {$c = $color{12};  $name =  $name."^"}
		if ($group eq "hyp")    {$c = $color{13};  $name =  $name."?"}
		if ($group eq "phage")  {$c = $color{14};  $name =  $name."#"}
		if ($group eq "IS")     {$c = $color{15};  $name =  $name."+"}

		my $pp =  int($begin/$axis_width);
		my $pp_beg = $begin - $pp*$axis_width;
		my $pp_end = $end - $pp*$axis_width;
		my $p = "Ps-page_".$ct."_".sprintf("%03d",$pp);

		unless ($PSFH{$p}) {die "CT:  $ct; P: $p\n$c\nI: $i; E_A: $extract_anno[$i];\n	newpath $pp_beg $y moveto $pp_end $y lineto 0 setlinecap stroke\n\n";}

		print { $PSFH{$p} } "$c\nnewpath $pp_beg $y moveto $pp_end $y lineto 0 setlinecap stroke\n";
		print { $PSFH{$p} } "newpath $pp_beg  $y $sub add moveto ($name) show\n";

	}
	return ();
}	#		Printing annotations (if genbank format ref sequence) - end

#		Printing reads - begin
sub Init_mapstrings {  # Initialize $k_max empty strings for printing reads
	for (keys %seqs) {
		my $key = $_;
		my $l = $seq_length{$key}/$opt{map_scale} + 2*$opt{space};
		my $string = " " x $l;
		for my $i (1 .. $k_max) {
			$mapstring_plus{$key}[$i]  = $string;
			$mapstring_minus{$key}[$i] = $string;
		}
	}
	return;
}

sub Print_reads
{
	my ($begin,$end,$len,$endsp,$k,$ext2);

	PS_all ("\n%Writing long divided reads:\n	$color{4}");
	PS_all ("% Font and Line width:\n	/Courier findfont \n	12 scalefont setfont\n	$opt{l_w} setlinewidth ");

	foreach my $key (@seqs) {
		$cont_cov{$key} ||= 0; # Initialize coverage counting
		my $c = 0;
		while ($bt_r_multi_long{$key}[$c]) {
			my @map = split "\t", $bt_r_multi_long{$key}[$c];
			
			$cont_cov{$key} += abs($map[4]-$map[5]);
			
			if ($map[1]==1 ) {	# Fw strand
				# Note: $map[4] > $map[5]!!
				# find available space
				$begin=int($map[4]/$opt{map_scale});
				$end  =int($map[5]/$opt{map_scale});
				$endsp=int($map[5]/$opt{map_scale})+$opt{space};
				$k=1;
				while (substr($mapstring_plus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					$len = $endsp-$begin;
					my $replace = "-"x$len;
					if ($begin+$len < length($mapstring_plus{$key}[$k])) {
						substr($mapstring_plus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd};
				if ($map[8] eq "Circle") { $ext2 = $arrow{circular};}
				if ($map[8] eq "Contig_Join") { $ext2 = $arrow{ct};}
			}
			if ($map[1]==-1 ) {	# Rev strand
				# Note: $map[5] > $map[4]!!
				# find available space
				$begin=int($map[5]/$opt{map_scale});
				$end  =int($map[4]/$opt{map_scale});
				$endsp=int($map[4]/$opt{map_scale})+$opt{space};
				$k=1;
				while (substr($mapstring_minus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					$len = $endsp-$begin;
					my $replace = "-"x$len;
					if ($begin+$len < length($mapstring_minus{$key}[$k])) {
						substr($mapstring_minus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd_rv};
				if ($map[8] eq "Circle") {$ext2 = $arrow{circular};}
				if ($map[8] eq "Contig_Join") { $ext2 = $arrow{ct};}
			}

			# print postscript now?
			#which page?
			my $pp =  int($begin/$axis_width);
			my $pp_beg = $begin - $pp*$axis_width;
			my $pp_end = $end - $pp*$axis_width;

			my $p0 = "Ps-page_".$key."_".sprintf("%03d",$pp);
			$pp++;
			my $p1 = "Ps-page_".$key."_".sprintf("%03d",$pp);

			my $pp_add = 0;
			my $mp_add = 0;
			#	$pp++;
			#	my $p1 = "Ps-page_".$key."_".sprintf("%03s",$pp);
			my $y = ($opt{space2reads}+$k*$opt{line_space})*$map[1];
			if ( $pp_end > $axis_width ) {
				$pp_add = ($y>0 ? -6 : -3 );
				$mp_add = ($y>0 ? 8 : 8 );
				my $pp_end_nxt = $pp_end-$axis_width;
				$pp_end = $axis_width;
				print { $PSFH{$p0} } "newpath $pp_end    $pp_add add   $y moveto  $ext2\n";
				#	print { $PSFH{$p0} } "newpath $pp_end    $pp_add add   $y moveto  $ext2\n";
				print { $PSFH{$p1} } "newpath -28            $y moveto  $ext2\n";
				print { $PSFH{$p1} } "newpath 0 $y moveto $pp_end_nxt $y lineto 0 setlinecap stroke\n";
				print { $PSFH{$p1} } "newpath $pp_end_nxt 10 add  $pp_add add $y 3 sub moveto ($map[0]) show\n";
			}
			if ($pp_beg < 10) { print { $PSFH{$p0} }  "newpath -80 $y moveto $ext2 \n"; }
			if ($seq_length{$key}/$opt{map_scale} - int($end) < 10) { print { $PSFH{$p0} }  "newpath $pp_end 60 add  $y moveto $ext2\n"	; }
			print { $PSFH{$p0} } "newpath $pp_beg $y moveto $pp_end $y lineto 0 setlinecap stroke\n";
			print { $PSFH{$p0} } "newpath $pp_end 10 add  $mp_add add $y 3 sub moveto ($map[0]) show\n";
			$c++;
		}
	}

	PS_all ("\n%Writing short divided reads:\n	$color{2}");
	# test sort-reads:
	foreach my $key (@seqs)  {
		my (@map,@sr_min,@sr_max,@sr_beg,@sr_end,@sr_name,@sr_mark_a,@sr_strand,@sr_count);
		my $src = -1;		# src: short read count
		my $srrc = 0;		# srr: short read region count
		my $rd_prev = "";
		for my $rec ( @{ $bt_r_multi_short{$key} } ) {
			@map = split "\t", $rec;
			if ($map[0] ne $rd_prev) {
				$src++;
				$sr_name[$src] = $map[0];
				$srrc = 0;
				$sr_min[$src] = min ( $map[4], $map[5] );
				$sr_mark_a[$src] = "$map[4]\t$map[5]\t";
				$rd_prev = $map[0];
			}
			$sr_strand[$src][$srrc] =  $map[1];
			$sr_beg[$src][$srrc] = min($map[4],$map[5]);
			$sr_end[$src][$srrc] = max($map[4],$map[5]);;

			$sr_min[$src] = min ( $map[4], $map[5], $sr_min[$src]);
			$sr_max[$src] = max ( $map[4], $map[5], $sr_max[$src]);

			$sr_mark_a[$src] .= "$map[4]\t$map[5]\t";

			$sr_count[$src] = $srrc;
			$srrc++;
		}

		# Print reads:
		for my $i (0 .. $src)	{
			
			$cont_cov{$key} += ($sr_max[$i] - $sr_min[$i]);
			
			$begin = $sr_min[$i]/$opt{map_scale};
			$end   = $sr_max[$i]/$opt{map_scale};
			$endsp = $end+$opt{space};
			$len   = $endsp-$begin;

			# find available space; Primary string is the string of first==longest segment:
			if ($sr_strand[$i][0] ==1 ) {

				$k=1;
				while (substr($mapstring_plus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					my $replace = "-"x$len;
					if ($begin+$len < length($mapstring_minus{$key}[$k])) {
						substr($mapstring_plus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd};
				if ($map[8] eq "Circle") { $ext2 = $arrow{circular};}
				if ($map[8] eq "Contig_Join") { $ext2 = $arrow{ct};}
			}

			if ($sr_strand[$i][0]==-1 ) {
				$k=1;
				while (substr($mapstring_minus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					my $replace = "-"x$len;

					if ($begin+$len < length($mapstring_minus{$key}[$k])) {
						substr($mapstring_minus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd_rv};				
				if ($map[8] eq "Circle") {$ext2 = $arrow{circular};}
				if ($map[8] eq "Contig_Join") { $ext2 = $arrow{ct};}
			}

			# print postscript now?
			#which page?
			my $pp =  int($begin/$axis_width);
			my $pp_beg = $begin - $pp*$axis_width;
			my $pp_end = $end - $pp*$axis_width;
			my $p0 = "Ps-page_".$key."_".sprintf("%03d",$pp);
			my $pp1 = $pp+1;
			my $p1 = "Ps-page_".$key."_".sprintf("%03d",$pp1);
			my $y = ($opt{space2reads}+$k*$opt{line_space})*$sr_strand[$i][0];
			my $pp_add = 0;
			my $mp_add = 0;

			# Print lines
			for my $j ( 0 .. $sr_count[$i] ) {
				my $line_beg = $sr_beg[$i][$j]/$opt{map_scale} - $pp*$axis_width;
				my $line_end = $sr_end[$i][$j]/$opt{map_scale} - $pp*$axis_width;

				if ( $line_beg > $axis_width && $line_end > $axis_width ) {
					my $l_b = $line_beg- $axis_width;
					my $l_e = $line_end - $axis_width;
					if ($sr_strand[$i][$j] eq $sr_strand[$i][0]) {
						print { $PSFH{$p1} } "newpath $l_b $y moveto $l_e $y lineto 0 setlinecap stroke \n";				
					}
					else {
						print { $PSFH{$p1} } "gsave\n$color{3} \n	$opt{l_w} 2 mul setlinewidth	\n";
						print { $PSFH{$p1} } "newpath $l_b $y moveto $l_e $y lineto 0 setlinecap stroke \n";
						# missing arrow here?						
						print { $PSFH{$p1} } "grestore \n	";
					}
				}

				if ( $line_beg < $axis_width && $line_end < $axis_width ) {
					if ($sr_strand[$i][$j] eq $sr_strand[$i][0]) {
						print { $PSFH{$p0} } "newpath $line_beg $y moveto $line_end $y lineto 0 setlinecap stroke \n";
					}
					else {  #  An Inversion
						print { $PSFH{$p0} } "gsave\n$color{3} \n	$opt{l_w} 6 mul setlinewidth	\n";
						print { $PSFH{$p0} } "newpath $line_beg $y moveto $line_end $y lineto 0 setlinecap stroke\n ";
						print { $PSFH{$p0} } "newpath $line_beg 90 sub $y 55 add moveto -30 rotate  3 3 scale $arrow{nm}\n";
						print { $PSFH{$p0} } "grestore \n	";
						$arrow_count++;
					}
				}

				if ( $line_beg < $axis_width && $line_end >= $axis_width ) {		# missing arrow here? What going on here?
					$pp_add = ($y>0 ? -6 : -3 );
					$mp_add = ($y>0 ? 8 : 8 );
					my $line_end_nxt = $line_end - $axis_width;
					$line_end = $axis_width;
					print { $PSFH{$p0} }  "newpath $line_beg $y moveto $line_end $y lineto 0 setlinecap stroke\n";
					
					print { $PSFH{$p1} } "newpath 0 $y moveto $line_end_nxt $y lineto 0 setlinecap stroke\n";
					print { $PSFH{$p0} } "newpath $line_end  $pp_add add $y moveto $ext2\n";
					print { $PSFH{$p1} } "newpath -28 $y moveto $ext2\n";
					print { $PSFH{$p0} } "newpath $line_end 10 add $mp_add add $y 3 sub moveto ($sr_name[$i]) show\n";
				}

			}

			# Print Marks:
			my @sr_mark = split "\t",$sr_mark_a[$i] ;
			@sr_mark = sort {$a <=> $b} @sr_mark;

			my $mark_height = 2*$opt{l_w};
			print { $PSFH{$p0} } "$opt{l_w} 3 div setlinewidth	\n";
			if ( $pp_end >  $axis_width	) {print { $PSFH{$p1} } "$opt{l_w} 3 div setlinewidth	\n"; }
			for my $i (0 .. $#sr_mark) {
				if ($sr_mark[$i]){
					my $mark = $sr_mark[$i]/$opt{map_scale} - $pp*$axis_width ; ## -0.5;
					if ($mark < $axis_width) {
						print { $PSFH{$p0} } "newpath $mark $y $mark_height sub moveto $mark $y $mark_height add lineto 0 setlinecap stroke\n";
					}
					if ($mark >= $axis_width) {
						$mark = $mark - $axis_width;
						print { $PSFH{$p1} } "newpath $mark $y $mark_height sub moveto $mark $y $mark_height add lineto 0 setlinecap stroke\n";
					}
				}
			}
			print { $PSFH{$p0} } "$opt{l_w} setlinewidth	\n";
			if ( $pp_end >  $axis_width	) {print { $PSFH{$p1} } "$opt{l_w} setlinewidth	\n"; }

			#	Print Read_tag:
			if ($pp_beg < 10) { print { $PSFH{$p0} }  "newpath -80 $y moveto $ext2 \n"; }
			
			if ($seq_length{$key}/$opt{map_scale} - int($end) < 10) { print { $PSFH{$p0} }  "newpath $pp_end 60 add  $y moveto $ext2\n"	; }


			if ( $pp_end <= $axis_width	) { 
				print { $PSFH{$p0} } "newpath $pp_end 10 add $y 3 sub moveto ($sr_name[$i]) show\n"; 				
				# print { $PSFH{$p0} } "newpath $pp_end    $pp_add add   $y moveto  $ext2 \n";
				}
			if ( $pp_end >  $axis_width	) { print { $PSFH{$p1} } "newpath $pp_end 10 add $axis_width sub $y 3 sub moveto ($sr_name[$i]) show\n";  }

		}
	}

	PS_all ("\n%Writing single reads:\n	$color{1}");
	foreach my $key (@seqs)  {
		# push @{$bt_r_single{$key}},"END \n";
		my $c=0;
		while ($bt_r_single{$key}[$c]) {
			my @map = split "\t",$bt_r_single{$key}[$c];
			
			$cont_cov{$key} += abs($map[4]-$map[5]);
			
			if ($map[1]==1 ) {
				# Note: $map[4] > $map[5]!!
				# find available space
				$begin=int($map[4]/$opt{map_scale});
				$end  =int($map[5]/$opt{map_scale});
				$endsp=int($map[5]/$opt{map_scale})+$opt{space};
				$k=1;
				while (substr($mapstring_plus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					$len = $endsp-$begin;
					my $replace = "-"x$len;
					# print tst2 "+ $key; $seq_length{$key}; $k; $begin; $len; \n";
					if ($begin+$len < length($mapstring_plus{$key}[$k])) {
						substr($mapstring_plus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd};
			}
			if ($map[1]==-1 ) {
				# Note: $map[5] > $map[4]!!
				# find available space
				$begin=int($map[5]/$opt{map_scale});
				$end  =int($map[4]/$opt{map_scale});
				$endsp=int($map[4]/$opt{map_scale})+$opt{space};
				$k=1;
				while (substr($mapstring_minus{$key}[$k],$begin,$endsp-$begin) =~ /-/) {
					$k++;
				}
				$k = min ($k, $k_max);
				if ($k < $k_max) {	#reserve space
					$len = $endsp-$begin;
					my $replace = "-"x$len;
					if ($begin+$len < length($mapstring_minus{$key}[$k])) {
						substr($mapstring_minus{$key}[$k],$begin,$len,$replace);
					}
				}
				$ext2 = $arrow{hd_rv};
			}
			# print postscript now?
			#which page?
			my $pp =  int($begin/$axis_width);
			my $pp_beg = $begin - $pp*$axis_width;
			my $pp_end = $end - $pp*$axis_width;
			my $p0 = "Ps-page_".$key."_".sprintf("%03d",$pp);
			$pp++;
			my $p1 = "Ps-page_".$key."_".sprintf("%03d",$pp);
			my $y = ($opt{space2reads}+$k*$opt{line_space})*$map[1];
			my $pp_add = 0;

			if ( $pp_end > $axis_width ) {
				my $pp_end_nxt = $pp_end-$axis_width;
				$pp_end = $axis_width;
				$pp_add = ($y>0 ? -6 : -3 );

				print { $PSFH{$p0} } "newpath $pp_end $pp_add add      $y  moveto $ext2\n";
				print { $PSFH{$p1} } "newpath -28            $y  moveto $ext2\n";
				print { $PSFH{$p1} } "newpath 0 $y moveto $pp_end_nxt $y lineto 1 setlinecap stroke\n";
				print { $PSFH{$p1} } "newpath $pp_end_nxt 10 add $y 3 sub moveto ($map[0]) show\n";

				$pp_add = ($y>0 ? 8 : 8 );
			}
			print { $PSFH{$p0} }  "newpath $pp_beg $y moveto $pp_end $y lineto 1 setlinecap stroke\n";
			print { $PSFH{$p0} }  "newpath $pp_end 10 add $pp_add add $y 3 sub moveto ($map[0]) show\n";
			$single_read_count++;
			$c++;
		}
	}
	print "\nReads mapped$cov:\n\tUndivided reads:               ", sprintf("%10d",$single_read_count)," reads\n";
	print "\tDivided reads, short distance: ", sprintf("%10d",$short_div_read)   ," reads\n";
	print "\tDivided reads, long distance:  ", sprintf("%10d",$long_div_read)    ," reads\n";
	print "\tIndicated inversions:          ", sprintf("%10d",$arrow_count)      ," reads\n\n";

	return();
}
#		Printing reads - end


# --------------------------------------------------
# Postscript helpers
# --------------------------------------------------

sub PS {

	my (%opt) = @_;

	my $self = {
		p_print => [],
		line_no => $opt{line_no} // 0,
		line_sp => $opt{line_sp} // -30,
	};

	bless $self, 'PS';
	return $self;
}

sub PS::Y  { $_[0]->{line_no} * $_[0]->{line_sp} }
sub PS::YL { $_[0]->{line_no} * $_[0]->{line_sp} + 5 }
sub PS::nl { $_[0]->{line_no}++ }
sub PS::push { push @{$_[0]->{p_print}}, $_[1] }

sub PS::font {
	my ($self,$f,$s)=@_;
	$self->push("V\t/$f findfont $s scalefont setfont");
}

sub PS::color {
	my ($self,$c)=@_;
	$self->push("V\t$c");
}

sub PS::title {
	my ($self,$txt)=@_;
	$self->font('Courier',36);
	$self->push("1\t0\x1F".$self->Y."\x1F$txt");
	$self->nl;
	$self->font('Courier',24);
}

sub PS::text {
	my ($self,$x,$txt)=@_;
	$self->push("1\t$x\x1F".$self->Y."\x1F$txt");
	$self->nl;
}

sub PS::row2 {
	my ($self,$x1,$x2,$l,$r)=@_;
	$self->push("2\t$x1\x1F".$self->Y."\x1F$l\t$x2\x1F".$self->Y."\x1F$r");
	$self->nl;
}

sub PS::hr {
	my ($self,$add9)=@_;
	$self->push("L\t-10 ".($self->Y+15)." $add9 ".($self->Y+15));
	$self->nl;
}

sub PS_all {
	my ($txt) = @_;
	foreach my $key (@seqs) {
		for my $i (0 .. $map_pages{$key}) {
			my $p = "Ps-page_" . $key . "_" . sprintf("%03d", $i);
			# print $p "\n$txt\n";  # File handle not accessible; needs refactoring
			print { $PSFH{$p} } "\n$txt\n";
		}
	}
	return;
}


sub PS_finish
{
	PRINT_parameters();

	my ($IN, $OUT);

	foreach my $key (@seqs) {
		for my $i (0 .. $map_pages{$key}) {
			my $p = "Ps-page_".$key."_".sprintf("%03d",$i);
			print { $PSFH{$p} } "showpage\n";
		}
	}

	# IMPORTANT: flush all page files
	for my $p (keys %PSFH) {
		close $PSFH{$p};
	}

	my @files = ('Ps-page_front_000.ps');

	foreach my $key (@seqs) {
		for my $i (0 .. $map_pages{$key}) {
			my $p = "Ps-page_" . $key . "_" . sprintf("%03d",$i);
			push @files, "$p.ps";
		}
	}

	open $OUT, '>', "$opt{out}.ps"
	or die "Cannot write $opt{out}.ps: $!";

	local @ARGV = @files;

	while (<>) {
		print $OUT $_;
	}

	close $OUT;

	if (-e "$opt{out}.pdf") {
		unless (rename "$opt{out}.pdf", "$opt{out}.pdf.bak" ) {
			print "Cannot write to existing $opt{out}.pdf (possibly open in another program): \n$!\n\n";
			return;
		}
	}

	my $rc = system("ps2pdf", "$opt{out}.ps");

	if ($rc == -1) { warn "Failed to launch ps2pdf: $!\n"; 	}
	else {
		my $exit = $rc >> 8;
		if ($exit == 0) { print "\nConverted all panels in $opt{out}.ps to $opt{out}.pdf\n"; }
		else { warn "ps2pdf / Ghostscript failed (exit code $exit)\n"; }
	}

	return;
}

sub define_PS_settings	{
	# Setting colors:
	$color{1}  = "0.0 0.5 0.0 setrgbcolor"; # Undivided reads
	$color{2}  = "0.0 0.0 0.5 setrgbcolor"; # Divided reads, short distance
	$color{3}  = "0.8 0.0 0.0 setrgbcolor"; # Divided reads, short distance, revers direction
	$color{4}  = "0.5 0.0 0.0 setrgbcolor"; # Divided reads, long distance
	$color{7}  = "0.2 0.2 0.2 setrgbcolor"; # Axis, dark grey
	$color{8}  = "0.1 0.0 0.6 setrgbcolor"; # Title
	$color{9}  = "0.9 0.9 0.9 setrgbcolor"; # A pale-grey frame
	# Define colors for annotations:
	$color{10} = "0.0 0.0 0.6 setrgbcolor";	# Normal
	$color{11} = "0.8 0.0 0.0 setrgbcolor";	# RNA
	$color{12} = "0.4 0.0 0.3 setrgbcolor";	# Pseudo
	$color{13} = "0.4 0.2 0.0 setrgbcolor";	# Hypothetical
	$color{14} = "0.7 0.15 0.0 setrgbcolor";	# Phage
	$color{15} = "0.8 0.2 0.2 setrgbcolor";	# IS

	# Define arrows:
	$arrow{nm}     			= "12 0 rmoveto arrow    ";  				# ** arrow
	$arrow{ct}  				= "12 0 rmoveto arrow_ct ";  				# ** arrow_ct
	$arrow{hd}  				= "12 0 rmoveto arrow_hd ";  				# ** arrow_hd
	$arrow{hd_rv}  	  	= "12 0 rmoveto arrow_hd_rv ";	  	# ** arrow_hd_rv
	$arrow{circular}  	= "12 0 rmoveto arrow_circular ";  	# ** arrow_circular
	$arrow{inv}        	= "12 0 rmoveto arrow_inv ";  	    # ** arrow_inv

	return();
}
# Postscript printing end

# --------------------------------------------------
# if extract_reads
# --------------------------------------------------

sub extract_reads {

	my ($reads_file,$contig,$start,$end) = @_;

	my (@btop_list) = get_reads_in_region($contig,$start,$end);
	my $pf = join ("_", "extract-reads-list",$contig,$start,$end) ;
	open my $fh_list, '>', $pf. ".txt" or die "Cannot open $pf: $!";
	$fh_list->autoflush(1);

	my $no_of_reads = $#btop_list +1;
	print $fh_list "# read\tstart\tend\n# $contig\t$start\t$end";

	my %wanted = map { $_->{read} => 1 } @btop_list;
	my $remaining = scalar keys %wanted;
	my $fh = open_maybe_gz($reads_file);

	# detect format
	my $first = <$fh>;
	die "Empty reads file\n" unless defined $first;
	seek($fh, 0, 0);

	my $is_fastq = ($first =~ /^@/);
	my $is_fasta = ($first =~ /^>/);

	die "Unknown format\n" unless $is_fastq || $is_fasta;

	# end detect format


	# extract reads
	my $reads_extracted = 0;
	if ($is_fastq) {

		print "\nExtract $no_of_reads reads in $reads_file contig: $contig from $start to $end)\n\tand save in: $pf.fastq\n";

		open my $fh_seq, '>', $pf.".fastq" or die "Cannot open $pf.fastq: $!";

	while (1) {
	    my $h = <$fh>;
	    last unless defined $h;

	    next unless $h =~ /^@(\S+)/;
	    my $id = $1;

	    my @seq_lines;
	    my @qual_lines;

	    # read sequence lines
	    my $plus;
	    while (my $line = <$fh>) {
	        if ($line =~ /^\+/) {
	            $plus = $line;
	            last;
	        }
	        push @seq_lines, $line;
	    }

	    last unless defined $plus;

	    # compute sequence length
	    my $seq_len = 0;
	    for my $line (@seq_lines) {
	        chomp(my $tmp = $line);
	        $seq_len += length($tmp);
	    }

	    # read quality lines
	    my $qual_len = 0;
	    while ($qual_len < $seq_len) {
	        my $qline = <$fh>;
	        last unless defined $qline;
	        chomp(my $tmp = $qline);
	        $qual_len += length($tmp);
	        push @qual_lines, $qline;
	    }

	    # print if wanted
	    if (exists $wanted{$id}) {
	        print $fh_seq $h;
	        print $fh_seq @seq_lines;
	        print $fh_seq $plus;
	        print $fh_seq @qual_lines;

	        delete $wanted{$id};
	        $remaining--;
	        $reads_extracted++;

	        last if $remaining == 0;
	    }
	}

	close $fh_seq;

	} else {
		print "\nExtract $no_of_reads reads in $reads_file contig: $contig from $start to $end)\n\tand save in: $pf.fasta\n";
		open my $fh_seq, '>', $pf.".fasta" or die "Cannot open $pf.fasta: $!";

		my ($header, $info, $seq) = ("", "", "");

		while (my $line = <$fh>) {
		    chomp $line;
		    if ($line =~ /^>(\S+)/) {
		        # flush previous entry
		        if ($header && exists $wanted{$header}) {
		            print $fh_seq "$info\n$seq\n";
		            $reads_extracted++;
		            delete $wanted{$header};
		            $remaining--;
		            last if $remaining == 0;
		        }
		        $header = $1;
		        $info   = $line;
		        $seq    = "";
		    }
		    else {
		        $seq .= $line;
		    }
		}

		# handle last record (only if still needed)
		if ($header && exists $wanted{$header}) {
		    print $fh_seq "$info\n$seq\n";
		    $reads_extracted++;
		}		
		close $fh_seq;

	}
	return;
}

sub get_reads_in_region {
	my ($contig, $start, $end) = @_;

	my $w1 = int($start / 1000);
	my $w2 = int($end   / 1000);

	my %seen;
	my @hits;

	for my $w ($w1 .. $w2) {
		next unless exists $read_extract{$contig}{$w};

		for my $r (@{ $read_extract{$contig}{$w} }) {

			next if $seen{$r->{read}}++;   # ← removes duplicates early

			if ($r->{end} >= $start && $r->{start} <= $end) {
				push @hits, $r;
			}
		}
	}


	return (@hits);
}

# --------------------------------------------------
# end extract reads
# --------------------------------------------------


# --------------------------------------------------
# Help, Usage and other auxcillaries
# --------------------------------------------------


sub usage {

    my ($msg) = @_;
    print "\nERROR: $msg\n" if $msg;

    print <<"USAGE";

Usage:
  read_print.pl <seq_file> <btop_file> [options]

Required:
  seq_file        Reference sequence (GenBank or FASTA formatted)
  btop_file       Alignment file, long reads mapped with locally installed blast with the 
                  specific -outfmt setting as shown in this example [final btop is optional]:
	                blastn -db blastdb -query all_reads.fa -out file.btop 
	                -outfmt '6 delim=	qseqid sframe qstart qend sstart send qlen sseqid btop'
USAGE

    print "Options:\n";

for my $p (@PARAM_META) {
    my $default = fmt_opt($opt_default{ $p->{key} });

    printf "  %-30s %-30s (default: %s)\n",
        $p->{cli},
        $p->{help},
        $default;
}

printf "  %-30s %-30s \n", "-o  --out", "Output prefix";
print <<'EXTRA';

  -e, --extract-reads FILE,CONTIG,START,END
                              Extract reads from FILE mapped to CONTIG
                              within START - END.
                              May be specified multiple times.
  -a, --report_analysis       Report, suited for batch analysis
  -d, --debug                 Activate debugging
  -h, --help                  Show this help



Examples:
  read_print.pl ref.gb reads.btop
  read_print.pl ref.fa reads.btop -r 4000 -m 1000 -x 20 
      (include reads > 4kbp [fx. Pacbio reads], each match > 1 kbp, reduce coverage for large datasets)
  read_print.pl ref.fa reads.btop -o result -d 
      (set output prefix, start debug) 
  read_print.pl ref.fa reads.btop -e reads.fastq,chr1,50000,100000 
      (extract region specific reads for further analysis) 
  read_print.pl ref.fa reads.btop -e reads.fastq,chr1,50000,100000 -e reads.fastq,chr2,25000,75000
EXTRA

    exit;
}

sub Report_Analysis {
	# Print a summary of this analysis in the $opt{out}.Analysis.txt with this format:
	# First line:
	# bio_project\t bio_sample\t sra acc.\t assembly method \t sequencing technolgy\t organism \t sequence length \t coverage \tNo. of reads\t single reads \t short divided reads\t long divided reads\t inversion indicators
	# Following lines, one for each contig in the reference sequence:
	# 15 x "\t" Acc. no.\t sequence length\t coverage
	
	open my $fh, ">", "$opt{out}_Analysis.txt" or die "Cannot open output file: $!";
	
	my @bp = split ":", text($bio_project );   
	my @bs = split ":", text($bio_sample  );   
	my @sr = split ":", text($sra         );   
	my @am = split ":", text($asm_meth    );   
	my @st = split ":", text($seq_tech    );   
	
	my $cont_cov;
	foreach my $key (@seqs) {$cont_cov += $cont_cov{$key} / $seq_length_total}
	$cont_cov = sprintf "%8.2f",$cont_cov;
	my @cols =  ( $opt{out}, $bp[-1], $bs[-1], $sr[-1], $am[-1], $st[-1], $organism, $seq_length_total,  $cont_cov, $read_tot, $single_read_count, $short_div_read, $long_div_read, $arrow_count  );
	@cols = map { text($_) } @cols;
	chomp @cols;
	print $fh join("\t", @cols), "\n";	
	foreach my $key (@seqs) {
		
		my $cc = sprintf "%8.2f", $cont_cov{$key}/$seq_length{$key}  ;
		@cols = ("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", $key, $seq_length{$key}, $cc );
		print $fh join("\t", @cols), "\n";
	}
	print $fh "\n";
	return();
}


# --------------------------------------------------
# Small helpers
# --------------------------------------------------

sub open_maybe_gz {
    my ($file) = @_;

    if ($file =~ /\.gz$/) {
        open(my $fh, "-|", "gzip -dc $file") or die  $!;
        return $fh;
    } else {
        open(my $fh, "<", $file) or die $file, $!;
        return $fh;
    }
}

sub OPEN_DEBUG {
	for my $name (qw(log single_read short_read long_read bt_read multiquest to_multiquest multi_multiquest into_multiquest )) {
		open $debug{$name}, '>', $name."_debug.txt"
		or die "Cannot open $name _debug.txt: $!";
		$debug{$name}->autoflush(1);
	}
	print {$debug{single_read}} "Debug file reporting reads mapped at a single location without interupts\nSR: only one line in btop file\nM: 0: only one of several btop lines mapped"
}

sub logmsg {
	my ($msg) = @_;
	my $ts = strftime("%Y-%m-%d %H:%M:%S", localtime);
	print"$msg: [$ts] \n" unless $opt{debug};
	return ("$msg: [$ts] \n")
}

sub min {
	my @vals = grep defined, @_;
	return undef unless @vals;
	my $m = shift @vals;
	for (@vals) {
		$m = $_ if $_ < $m;
	}
	return $m;
}

sub max {
	my @vals = grep defined, @_;
	return undef unless @vals;
	my $m = shift @vals;
	for (@vals) {
		$m = $_ if $_ > $m;
	}
	return $m;
}

sub text { defined $_[0] ? $_[0] : '' }

sub num  { defined $_[0] && $_[0] =~ /^-?\d+(?:\.\d+)?$/ ? $_[0] : 0 }

