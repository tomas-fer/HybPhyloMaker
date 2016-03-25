#!/usr/bin/perl
use strict;
use Tie::File;
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned
use File::Copy ;

#### debugged 11th July, 2012 by Patrick Kueck -> Fasta Files could not be handled correctly in subroutine &fas
#### updated 18th July, 2012 by Patrick Kueck -> Inclusion of the ability to determine invariant positions of the whole data set and of single defined clades
#### updated 09th August, 2012 by Patrick Kueck -> Inclusion of aminoacid handling
#### updated 10th August, 2012 by Patrick Kueck -> Inclusion of svg output and partition handling
#### updated 12th September, 2012 by Patrick Kueck -> Inclusion of determination of compositional bias (RCFV Value), Absolute deviation of single states, and skew-calculation (for nucleotide data)
#### debugged 23rd October, 2012 by Patrick Kueck
#### updated 14th November, 2012 by Patrick Kueck -> inclusion of pairwise comparisons (c-value calculations of p-distances and ti/tv ratio; missing data )
#### debugged 26th November, 2012 by Patrick Kueck -> '0' sign now allowed in alignment infile name, error prompt included for alignment infiles without file specific suffix (.aln, .fas, .phy).
#### updated 7th december, 2012 by Patrick Kueck -> Proportion of invariabel sites refers only to the total number of sequences which include at least one informative character state
#### debugged 4th march, 2013 by Patrick Kueck -> bugs fixed from email T.Struck 1.3.13 & 3.3.13
#### debugged 20th march, 2013 by Patrick Kueck & Torsten Struck -> Calculation bug for C-Value determination fixed
#### check 20th march, 2013 by Patrick Kueck & Torsten Struck -> All BaCoCa Calculations
#### updated 23th march, 2013 by Patrick Kueck -> calculation and print out of taxon specific skew-values
#### updated 21st june, 2013 by Patrick Kueck -> pathfile handling enabled; clustalfile handling excluded
#### updated 24th june, 2013 by Patrick Kueck -> enabling the generation of heat maps using R commands
#### updated 1st july, 2013 by Patrick Kueck -> integration of -r option to help menu / generate heat maps for missing data overlap
#### debugged 18th july, 2013 by Patrick Kueck -> Genepartition calculations with less than 2 informative taxon sequences
#### debugged 18th october, 2013 by Patrick Kueck -> Ambiguity character 'B' was missing for nucleotide data analyses
#### debugged 27th august, 2013 by Patrick Kueck -> line 1942 -> N informative sites has to be higher than 1 -> N clae taxa > 2 instead > 1
#### debugged 16th june, 2015 by Patrick Kueck & Torsten Struck -> Substitution of lowercase characters within sequence by uppercase characters

print	"\n\n\t-------------------------------------------",
		"\n\t             BaCoCa.v1.107\n\t",
		"\n\tDeveloped by Patrick Kueck & Torsten Struck,\n\t",
		"written by Patrick Kueck (ZFMK Bonn, 2012/13)\n\t",
		"-------------------------------------------\n";



####################################
## Definition of global variables
my	$single_file		=	""			; # skript abort if input file is undefined via the '-i option'
my	$subclades			=	"None"		; # textfile with predefined clades can be given with the '-c option'
my	$partitions		=	"None"		; # textfile with predefined sequence partitions can be given with the '-p option'
my	$seqtype			=	'nucleotide';

my	(
		%filelist,				# key: filenumber; value: filename
		%taxa_of_file,			# key: filename; value: taxa of file
		%state_of_taxon_of_position,
		%seq_ranges_of_part,
		%freq_of_file_and_char,
		%clade_of_number,
		%file_of_number,
		%filelist,
		%seq_of_tax,
		%taxa_of_file_clade
	) ;
####################################



####################################
## Read IN ARGV
&argv_handling ( \$single_file, \$subclades, \$partitions ) ;
####################################



####################################
## Read IN Infile (specified by -i) and store taxonspecific
## sequences in hash %seq_of_tax (keys: taxon; value: sequence
&data_read_in	( \$single_file, \%seq_of_tax ) ;
####################################



####################################
## store filename as hasvalue, key1: filenumber
## store taxonnames in an extra array (@taxa)
my $counter_files					= 1 ;
$file_of_number{$counter_files}	=  $single_file;
####################################



####################################
## 
&data_check ( \$single_file, \%seq_of_tax, \%state_of_taxon_of_position, \%seq_ranges_of_part, \$seqtype ) ;
####################################



####################################
## If partitionfile is specified via command '-p'
## check format and exclude substring info
## store single partition info in %seq_ranges_of_part
unless ( $partitions =~ /^None$/ ){
	
	open INpart, "<$partitions" or die "\n\t!FILE-ERROR!: Cannot open IN partitionfile : ", $partitions, " !\n" ;
	while ( my $line = <INpart> ){
		
		if	( $line =~ /^\w+,( )?\w+( )?=( )?\d+( )?-( )?\d+(\\\d)?$/g	 ){
			
			my @line_parts	= split ",", $line ;
			$line			= $line_parts[1] ;
		}
		
		chomp $line ;
		$line =~ s/ |\s+//g ;
		
		if 		( $line =~ /^\w+=\d+-\d+\\\d$/g	){ die "\n\t!FILE-ERROR: Unknown partition format in line ", $line, "\n\tCannot take Codon positions into account!\n\tCorrect format of partition file:\n\tPARTITIONNAME1 = STARTPOSITION - ENDPOSITION\n\tPARTITIONNAME2 = STARTPOSITION - ENDPOSITION...\n\n\tSTARTPOSITION and ENDPOSITION = integers" }
		unless	( $line =~ /^\w+=\d+-\d+$/g		){ die "\n\t!FILE-ERROR: Unknown partition format in line ", $line, "\n\tCorrect format of partition file:\n\tPARTITIONNAME1 = STARTPOSITION - ENDPOSITION\n\tPARTITIONNAME2 = STARTPOSITION - ENDPOSITION...\n\n\tSTARTPOSITION and ENDPOSITION = integers" }
		
		$line =~ s/=|-/:/g ;
		
		my @part_data = split ":", $line ;
		if ( $part_data[1]	> $part_data[2] ){ die "\n\t!PARTITION-ERROR!: Endposition of partition ", $part_data[0], " < Startposition\n" }
		
		$counter_files++;
		$file_of_number{$counter_files} = $part_data[0] ;
		$seq_ranges_of_part{$part_data[0]}{start}	= $part_data[1] - 1 ;
		$seq_ranges_of_part{$part_data[0]}{end}  	= $part_data[2] - 1 ;
		
		
		print "\tFound Partition: ", $part_data[0], "\tposition: ", $part_data[1], "-", $part_data[2], "\n" ;
	}
	close INpart ;
	
	if ( $partitions =~ /\/+/ ){
		
		my @part_file_path_parts	= split "/", $partitions ;
		$partitions				= $part_file_path_parts[-1] ;
	}
	#print "\npartition file\t", $partitions,"\n";
}
####################################



####################################
## Generate common result folder for BaCoCa analyses
mkdir "BaCoCa_Results", 0755 ;
####################################



####################################
## in hash %sequence_of_taxon (subroutine &data_read_in)
## Afterwards, check for correct taxon and sequence conditions
FILE: # Anchor point. if incorrect formated filestructures are found in the subroutine &data_check, the script skips directly to the next inputfile
for my $file_number ( sort {$a<=>$b} keys %file_of_number ){
	
	&calc_freq	( \$file_of_number{$file_number}, \%seq_ranges_of_part, \%state_of_taxon_of_position, \$subclades, \%freq_of_file_and_char, \%clade_of_number, \%taxa_of_file_clade, \$seqtype ) ;
}
####################################



my @taxa_all = keys %seq_of_tax ;
&print_out	( \%freq_of_file_and_char, \%clade_of_number, \%file_of_number, \%taxa_of_file_clade, \@taxa_all, \$seqtype, \%seq_ranges_of_part ) ;
&end;




exit;


sub synopsis{
		
		print	"\n\n\tSYNOPSIS\n\t-------------------------------------------",
				"\n\tStart (UNIX):\n\tperl BaCoCa.v1.105.pl -i <I> -c <C> -p <P>\n",
				"\n\tStart (WINDOWS):\n\tBaCoCa.v1.105.pl -i <I> -c <C> -p <P> -r\n\n",
				"\n\t<I>: Inputfile (FASTA, CLUSTAL or PHYLIP)",
				"\n\t<C>: Subcladefile (optional)",
				"\n\t<P>: Partitionfile (optional)",
				"\n\t-r: Generate Heat Maps using R",
				"\n\n\tFor further info about optional files type:",
				"\n\t(perl) BaCoCa.v1.105.pl -help\n\t",
				"\n\n\tFor further info about BaCoCa type:",
				"\n\t(perl) BaCoCa.v1.105.pl -v\n\t",
				"-------------------------------------------\n\n" ;
}


sub help{

	print
<<help

        ----------------BaCoCa HELP-----------------
	BaCoCa is a software to execute base composition
	analyses of nucleotide and aminoacid data. Beside
	estimations of single base, gap and ambiguity freq-
	uencies, BaCoCa additionally determines frequencies 
	of purines and pyrimidines for nucleotide data and 
	frequencies of hydrophilicity, polarity, and atomic 
	charges for aminoacid data. Furthermore, BaCoCa
	calculates the amount of compositional bias by 
	measuring the relative composition frequency var-
	iability (RCFV) for complete data sets and for each
	taxon, the absolute deviations of single nucleotide
	and aminoacid states, the absolute deviations of 
	different state combinations (i.e. AT and Y or R),
	different skew values to detect biases between
	two nucleotide frequencies (not implemented for
	aminoacid data), and performs a homogeneity test
	of base frequencies like included in PAUP.
	
	BaCoCa is able to perform all analyses for com-
	plete data sets, predefined taxon subclades, and
	predefined data partitions in one single run.
	
	Subclade definitions:
	----------------------
	To perform BaCoCa analyses on predefined
	subclades use the '-c' option:
		
		-c 'cladefile.txt'
	
	The 'cladefile.txt' has to be in the following
	format:
		
		Subcladename1,taxon1,taxon2,...,taxonX
		Subcladename2,taxon1,taxon3,...,taxonY
		...
		SubcladenameZ,taxon4,taxon3,...,taxonZ
	
	Be aware that taxa of defined subclades are
	identic to names of the corresponding sequence
	datafile (take care of upper and lower cases).
	Subcladenames are important to assign subclade
	specific BaCoCa output files.
	----------------------
	
	Partition definitions:
	----------------------
	To perform BaCoCa analyses on defined data
	partitions use the '-p' option:
		
		-p 'partitionfile.txt'
	
	The 'partitionfile.txt' should be in the
	following format:
		
		Partitionname1 = Startposition - Endposition
		Partitionname2 = Startposition - Endposition
		...
		PartitionnameZ = Startposition - Endposition
	
	but can be in a RAxML adapted partition format
	as well:
		
		AlphanumericSigns, Part1 = Start - End
		AlphanumericSigns, Part2 = Start - End
		...
		AlphanumericSigns, Part2 = Start - End
		
	Startposition and Endpositions have to be given
	as integers and are not allowed to overlap be-
	tween different partitions. Blanks have not to
	be inserted (negligible).
	----------------------
	
	Generate result heat maps using R (-r)
	--------------------------------------
	Only for BaCoCa.vX.X.r.pl! To perform heat map 
	analyses using R, R as well as additional PERL 
	and R packages have to be installed
	(Linux (Ubuntu 12.04 or higher) / Mac):
	
		sudo apt-get install libstatistics-R-perl
		sudo apt-get install r-base
		sudo apt-get install r-cran-gplots
	
	If additional packages are installed, BaCoCa will
	perform heat map analyses with the -r command.
	--------------------------------------
	
	For further information about the algorithm, usage, 
	and options of BaCoCa read the BaCoCa manual.
	
	For further questions or bug reports, please write to
	patrick_kueck\@web.de
	
help
;

}


sub version{
	
	print
<<preface
	
	--------------------BaCoCa PREFACE---------------------
	
	Version     : 1.105
	Language    : PERL
	Last Update : 18th October, 2013
	Author      : Patrick Kück, ZFMK Bonn, GERMANY
	e-mail      : patrick_kueck\@web.de
	Homepage    : http://www.zfmk.de
	
	This program is free software; you can distribute it 
	and/or modify it under the terms of the GNU General Public 
	License as published by the Free Software Foundation ; 
	either version 2 of the License, or (at your option) any 
	later version.

	This program is distributed in the hope that it will be 
	useful, but WITHOUT ANY WARRANTY; without even the 
	implied warranty of MERCHANTABILITY or FITNESS FOR A 
	PARTICULAR PURPOSE. See the GNU General Public License for 
	more details. 

	You should have received a copy of the GNU General Public 
	License along with this program; if not, write to the Free 
	Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, 
	USA.
	
	For further free downloadable programs visit:
	http://www.zfmk.de/web/ZFMK_Mitarbeiter/KckPatrick/Software/index.de.html
	
	------------------------------------------------------------

preface
; 
}


sub argv_handling{
	
	
	my $sref_single_file	= $_[0]  ;
	my $sref_clades		= $_[1]  ;
	my $sref_partition	= $_[2]  ;
	
	
	####################################
	## ARGV handling: type via commandline 
	## e.g.: 'perl basefrequency_estimator.pl -i inputfile.fas -c cladefile.txt
	my ( $commandline ) = join "", @ARGV ;
	
	if ( $commandline ){ 
		
		$commandline		=~	s/ |\s+//g ;
		my		@commands	=	split "-", $commandline ; 
		shift	@commands; # Remove first empty element due to split "-", which is the first sign in the line
		
		
		for my $single_command ( sort @commands ){
			
			if 		( $single_command =~ /^i/i		){ ( $$sref_single_file	= $single_command	) =~ s/^.// } # Single Data file: -i *.phy or *.fas or *.aln
			elsif 	( $single_command =~ /^c/i		){ ( $$sref_clades		= $single_command	) =~ s/^.// } # Clade file: -c *.txt
			elsif 	( $single_command =~ /^p/i		){ ( $$sref_partition		= $single_command	) =~ s/^.// } # Partition file: -p *.txt
			elsif 	( $single_command =~ /^help$/i	){ &help and exit }
			elsif 	( $single_command =~ /^v$/i		){ &version and exit }
			else 										 { print "\n\t!COMMAND-ERROR!: ", $single_command, " unknown!" and &synopsis }
		}
		
		print	"\n\tBaCoCa Filehandling:\n\t";
		unless ( $$sref_single_file			){ print	"!FILE-ERROR!: No infile (.aln, .phy, .fas) specified\n\t"	; &synopsis; exit	}	else { print "Infile: ",				$$sref_single_file,	"\n\t" }
		if ( $$sref_clades		=~ /^None$/	){ print	"Clades Specified: No\n\t"														}	else { print "Clades specified: ",		$$sref_clades,		"\n\t" }
		if ( $$sref_partition	=~ /^None$/	){ print	"Partition File: No\n\t"														}	else { print "Partitions specified: ",	$$sref_partition,		"\n\t" }
	}
	else{ &synopsis ; exit }
	####################################
}


sub data_read_in{
	
	my $sref_file			=	$_[0] ;
	my $href_seq_of_tax	=	$_[1] ;
	
	print "\n\tReadin data...";
	
	####################################
	## Read IN different possible formats, untie linefeeds, identify single taxa and associated sequences, and store both in the hash
	## $href_seq_of_tax (key:taxon; value: corresponding nucleotide sequence
	if 		( $$sref_file =~ /^.*\.phy$/				){ &tie_linefeeds( \$$sref_file ) ; &phy ( \$$sref_file, \%$href_seq_of_tax ) }
	elsif 	( $$sref_file =~ /^.*\.fas$|^.*\.FASTA$/	){ &tie_linefeeds( \$$sref_file ) ; &fas ( \$$sref_file, \%$href_seq_of_tax ) }
	else 	{ die "\n\t!FILE-ERROR!: Cannot readin ", $$sref_file, " !\n\tMaybe forgotten file suffix (.fas, .phy, .aln) ?\n" }
	
	if 		( $$sref_file =~ /\/+/ ){
		
		my @infile_path_parts	= split "/", $$sref_file ;
		$$sref_file			= $infile_path_parts[-1] ;
	}
	#print "\ninfile\t", $$sref_file,"\n";
	####################################
	
	
	
	######################################################################## subroutines within subroutine data_read_in
	sub phy{
		
		my $sref_phy_file	= $_[0] ;
		my $href_seq_of_tax	= $_[1] ;
		my %seen_phylip_taxon ;
		
		
		########################################################################
		## READ IN phylip interleaved and non-interleaved formated files
		## store single taxa and their associated sequence in has $href_seq_of_tax
		open INphy , $$sref_phy_file or die "\n\t!FILE-ERROR!: Cannot find ", $$sref_phy_file, " !\n" ;
		chomp (my @all_lines_phy = <INphy>) and close INphy ; 
		
		
		
		####################################
		## extract the first line (infoline) and determine the number of taxa ($info_line[0])
		( my $infoline	= shift @all_lines_phy ) =~ s/\s+/ / ;
		my	@infos_line	=  split " ", $infoline ;
		####################################
		
		
		
		####################################
		## phylip files can be in interleaved and non-interleaved format
		## interleaved:	tax1 ACGT...	# first part of all taxa
		## 				tax2 ACGT...
		##								#''space line''
		##				tax1 CCCC...	# Second part of all taxa
		##				tax2 GGGG...
		## to concatenate sequence parts correctly store single lines in separate hashkeys
		## until number of taxa is reached ($c equal $infos_line[0]). afterwards remove the following spaceline and concatenate
		## next lines to their corresponding taxon sequences inferred from the first round and so on...
		## If phylip file is in non-interleaved format, the while lopp stops automatically after the first foreach loop
		my 	%seq_phy = () ;
		while ( @all_lines_phy ){
			
			for ( my $c=1; $c<=$infos_line[0]; $c++ ){ my $seq_line_phy = shift @all_lines_phy ; push ( @{$seq_phy{$c}} , $seq_line_phy ) }
			shift @all_lines_phy ;
		}
		####################################
		
		
		
		####################################
		## join single sequence parts of each taxon (if interleaved there are multiple key values)
		## taxonnames are in the same line as sequenceinformation (separated by one or multiple whitespaces), therefore
		## substitute multiple whitespaces into one whitespace, join all sequence parts to one string (only important if interleaved format)
		## split taxonnames from sequenceinformation and store both in the hashreference %href_seq_of_tax (key: taxon; value: complete sequence)
		for my $line_c ( sort {$a<=>$b} keys %seq_phy ){ 
			
			my	@seq_single_parts				=	exists($seq_phy{$line_c}) ? @{$seq_phy{$line_c}} :( ) ;
				$seq_single_parts[0]			=~	s/\s+/ / ;
			my	$seq_complete					=	join	""		, @seq_single_parts ;
				@seq_single_parts				=	split	" "		, $seq_complete ;
			my	$taxon							=	shift @seq_single_parts ;
			
			unless ( $seen_phylip_taxon{$taxon} ){ $seen_phylip_taxon{$taxon}++ } else{ die "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in file ", $$sref_phy_file, " !\n" }
			
				$seq_complete					=	join	""		, @seq_single_parts ;
				$href_seq_of_tax->{$taxon}	=	$seq_complete ;
		}
		########################################################################
	}
	
	sub fas{
		
		my $sref_fas_file		= $_[0] ;
		my $href_seq_of_tax		= $_[1] ;
		my %seen_fasta_taxon ;
		
		####################################
		## READ IN fasta interleaved and non-interleaved formated files
		## store single taxa and their associated sequence in has $href_seq_of_tax
		my 	%seq_of_tax ; my $taxon ; my $seq_complete ;
		
		open INfas, $$sref_fas_file or die "\n\t!FILE-ERROR!: Cannot find ", $$sref_fas_file, " !\n" ;
		while ( my $line = <INfas> ){
			
			chomp	$line ;
			
			if ( $line =~ /^\>/ )	{ ( $taxon = $line ) =~ s/^\>// ; unless ($seen_fasta_taxon{$taxon}){ $seq_complete = () ; $seen_fasta_taxon{$taxon}++} else{ die "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in file ", $$sref_fas_file, " !\n" } } 
			else 					{ $seq_complete = $seq_complete.$line and $href_seq_of_tax->{$taxon} = $seq_complete} 
		}
		close INfas ;
		####################################
	}
	
	sub tie_linefeeds{
		
		my $sref_tie_file = $_[0] ;
		
		####################################
		## Untie linefeeds of Inout files
		TIE:
		(tie ( my @data, 'Tie::File', $$sref_tie_file )) ;
		
		die  "\n\t!FILE-ERROR!: $$sref_tie_file is empty!\n" if 0 == @data ;
		
		map { s/\r\n/\n/g } @data ;
		map { s/\r/\n/g   } @data ;
		
		untie @data ;
		####################################
	}
	########################################################################
}


sub data_check{
	
	my $sref_file							=	$_[0] ;
	my $href_seq_of_tax						=	$_[1] ;
	my $href_state_of_taxon_of_position		=	$_[2] ;
	my $href_seq_ranges_of_part				=	$_[3] ;
	my $sref_type_seq						=	$_[4] ;
	
	my %seen_taxa;
	my $seq_length; my $Nstate_seq1 ;
	
	print "\n\tCheck sequence states ",  $$sref_file ;
	
	####################################
	## Checking sequence and taxon signs
	for my $taxon ( keys %$href_seq_of_tax ){ 
		
		####################################
		## Checking Taxon Names
		my @taxon_signs = split "", $taxon ;
		for my $tax_sign (@taxon_signs){ unless	( $tax_sign =~ /\w|\d|\s/	){ die "\n\t!FILE-ERROR!: Forbidden Sign in taxon name ", $taxon, " !\n" } }
		####################################
		
		####################################
		## Checking sequences
		if		( $href_seq_of_tax->{$taxon}	=~ /^\n$|^$/	){ die "\n\t!FILE-ERROR!: Sequence missing in ", $$sref_file, " for taxon ", $taxon, " !\n" }
		unless	( $seen_taxa{$taxon} ){ $seen_taxa{$taxon}++ } else { die "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in file ", $$sref_file, " !\n" }
		
		
		$href_seq_of_tax->{$taxon}	=~ s/\Q($href_seq_of_tax->{$taxon})/\U$1/gi ;
		$href_seq_of_tax->{$taxon}	=~ s/u/T/gi ;
		$href_seq_of_tax->{$taxon}	=~ s/(A|C|G|T|N|Y|R|W|S|K|M|D|V|H|I|E|L|Q|F|P|B|X)/\U$1/gi ;
		my @base_states					= split "",	$href_seq_of_tax->{$taxon} ; unless ($Nstate_seq1){ $Nstate_seq1 = @base_states }
		my $N_states					= @base_states ;
		
		
		unless ( $N_states == $Nstate_seq1 ){ 
			
			for my $seentax (sort keys %seen_taxa){ print "\n\t", $seentax, " ", $seen_taxa{$seentax}, "bp" }
			
			if ( $N_states < $Nstate_seq1 ){ my $N_difference_bp = $Nstate_seq1 - $N_states; die "\n\t!FILE-ERROR!: Sequences of infile have unequal sequence lengths !\n\tSequence of ", $taxon, " ", $N_difference_bp, " shorter as sequences before!\n" }
			if ( $N_states > $Nstate_seq1 ){ my $N_difference_bp = $N_states - $Nstate_seq1; die "\n\t!FILE-ERROR!: Sequences of infile have unequal sequence lengths !\n\tSequence of ", $taxon, " ", $N_difference_bp, " larger as sequences before!\n" }
		}
		
		
		my $counter_baseposition		= 0 ;
		
		for my $base_state ( @base_states ){
			
			unless ( $base_state =~ /A|C|G|T|N|Y|R|W|S|K|M|D|V|H|I|E|L|Q|F|P|B|X|-|\?/i	){ die "\n\t!FILE-ERROR!: Forbidden sign in sequence", $taxon,":\n\t", $base_state, " in file ", $$sref_file, " !\n" }
			unless ( $base_state =~ /A|C|G|T|N|Y|R|W|S|K|M|D|V|H|B|-|\?/i					){ $$sref_type_seq = 'aminoacid' }
			
			if 		( $base_state =~ /-/  )	{ $href_state_of_taxon_of_position->{$taxon}{$counter_baseposition} = 'GAP' }
			elsif 	( $base_state =~ /\?/ )	{ $href_state_of_taxon_of_position->{$taxon}{$counter_baseposition} = 'MIS' }
			else 							{ $href_state_of_taxon_of_position->{$taxon}{$counter_baseposition} = $base_state }
			
			#print "\n", $taxon, " :Ngaps: ", $href_state_of_taxon_of_position->{$taxon}{$counter_baseposition};
			$counter_baseposition++
		}
		$seq_length = $counter_baseposition ;
		####################################
	}
	
	print "\n\tType of data: ", $$sref_type_seq, "\n\n";
	
	####################################
	## if different sequence lengths detected, set length to 'variable' -> No determination of invariant positions possible!
	$href_seq_ranges_of_part->{$$sref_file}{start} = 0 ; $href_seq_ranges_of_part->{$$sref_file}{end} = $seq_length-1 ; #print $seq_length; exit;
	####################################
	
	####################################
}


sub calc_freq{
	
	my $sref_file							=	$_[0] ;
	my $href_seq_ranges_of_part				=	$_[1] ;
	my $href_state_of_taxon_of_position		=	$_[2] ;
	my $sref_cladefile						=	$_[3] ;
	my $href_frequeny_of_file_of_character	=	$_[4] ;
	my $href_clade_of_number				=	$_[5] ;
	my $href_taxa_of_file_clade				=	$_[6] ;
	my $sref_seqtype						=	$_[7] ;
	
	
	
	my $cladecounter = 0 ;
	my @taxa = keys %$href_state_of_taxon_of_position ;
	my $N_taxa_complete_data = @taxa ;
	my @taxa_informative;
	
	
	#####################################
	## store infile taxa in a hash %seen_taxa to check
	## if defined taxon clades are all in the infile
	my %seen_taxa;
	for my $taxon (@taxa){$seen_taxa{$taxon}++}
	@{$href_taxa_of_file_clade->{$$sref_file}{All_Taxa}} = @taxa ;
	#####################################
	
	
	
	####################################
	## Calculate and print frequencies of complete data files
	$href_clade_of_number->{$cladecounter}{$$sref_file} = 'All_Taxa';
	&calculate_frequencies	( \$$sref_file, \@taxa, \%$href_frequeny_of_file_of_character, \$cladecounter, \$$sref_seqtype, \$N_taxa_complete_data, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number, \@taxa_informative  ) ;
	####################################
	
	
	
	####################################
	## Calculate and print c-values and standard deviation for nucleotide data
	if ( $$sref_seqtype =~ /nucleotide/ ){
		
		&pairwise_comparisons_nuc		( \$$sref_file, \@taxa, \@taxa_informative, \%$href_frequeny_of_file_of_character, \$cladecounter, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number  )
	}
	else { &pairwise_comparisons_aa	( \$$sref_file, \@taxa, \%$href_frequeny_of_file_of_character, \$cladecounter, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number ) }
	####################################
	
	
	
	########################################################################## Calculate Basefrequencies of predefined taxon clades
	## Calculate and print frequencies for specified clades
	unless ( $$sref_cladefile =~ /^None$/ ){
		
		
		#####################################
		## Read IN specified Clade file .txt, break up comma separated taxon names of each clade (1 clade per line)
		## Calculate and print Clade-frequencies
		open INfile, "<$$sref_cladefile" or die "\n\tFILE-ERROR!: Clade-File ", $$sref_cladefile, " cannot be found !" ;
		while ( my $line = <INfile> ){
			
			$cladecounter++;
			
			chomp $line;
			my @clade_taxa		= split ",", $line ;
			my $subclade_code	= shift @clade_taxa;
			@{$href_taxa_of_file_clade->{$$sref_file}{$subclade_code}} = @clade_taxa;
			
			CLADE:
			for my $taxon (@clade_taxa){
				
				unless ( $seen_taxa{$taxon} ){ warn "\n\t!TAXON-ERROR!: Clade taxon ",$taxon, " not found in inputfile ", $$sref_file, "!\n";  $cladecounter++ and next CLADE}
			}
			
			@taxa_informative = () ;
			
			$href_clade_of_number->{$cladecounter}{$$sref_file} = $subclade_code;
			&calculate_frequencies	( \$$sref_file, \@clade_taxa, \%$href_frequeny_of_file_of_character, \$cladecounter, \$$sref_seqtype, \$N_taxa_complete_data, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number, \@taxa_informative  ) ;
			
			if ( $$sref_seqtype =~ /nucleotide/ ){
				
					&pairwise_comparisons_nuc	( \$$sref_file, \@clade_taxa, \@taxa_informative, \%$href_frequeny_of_file_of_character, \$cladecounter, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number  )
			}
			else {	&pairwise_comparisons_aa	( \$$sref_file, \@clade_taxa, \%$href_frequeny_of_file_of_character, \$cladecounter, \%$href_seq_ranges_of_part, \%$href_state_of_taxon_of_position, \%$href_clade_of_number ) }
		}
		#####################################
	}
	##########################################################################
	
	
	
	
	
	######################################################################## subroutines
	sub calculate_frequencies{
		
		my $sref_infile				= $_[0] ;
		my $aref_clade				= $_[1] ;
		my $href_frequencies			= $_[2] ;
		my $sref_counter				= $_[3] ;
		my $sref_type_seq				= $_[4] ;
		my $sref_N_tax_complete		= $_[5] ;
		my $href_ranges_of_part		= $_[6] ;
		my $href_state_of_tax_pos	= $_[7] ;
		my $href_clade_of_numbers	= $_[8] ;
		my $aref_taxa_informative	= $_[9] ;
		
		
		############################################################################################################################################################################################################################
		################################################################################ DETERMINE SEQUENCE LENGTH #################################################################################################################
		############################################################################################################################################################################################################################
		
		######################
		## Count characters of each sequence position (important for determination of invariant positions)
		my $sequ_length	= $href_ranges_of_part->{$$sref_infile}{end} - $href_ranges_of_part->{$$sref_infile}{start} +1 ; #print "start", $href_ranges_of_part->{$$sref_infile}{start}, "end", $href_ranges_of_part->{$$sref_infile}{end} ,"length", $sequ_length; exit;
		$href_frequencies->{$$sref_infile}{$$sref_counter}{length} = $sequ_length;
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ ASSIGN CHARACTER STATES ###################################################################################################################
		############################################################################################################################################################################################################################
		
		######################
		## Assign single state chracter frequencies of each taxon to the hash %href_frequencies and count single states
		my %seen_state; 
		for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
			
			for my $taxon ( @$aref_clade ){
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$href_state_of_tax_pos->{$taxon}{$position}}++ ; #print $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$href_state_of_tax_pos->{$taxon}{$position}};
				
				$seen_state {$position}{$href_state_of_tax_pos->{$taxon}{$position}}++
			}
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ CHECK INFORMATIVE TAXA ####################################################################################################################
		############################################################################################################################################################################################################################
		print "\n\tCheck for complete uninformative taxon sequences ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile};
		
		######################
		## Check for NUCLEOTIDE taxon sequences which only consists of missing data, ambiguity characters and indel states
		## These taxa are not further applicated ( not included in @$aref_taxa_informative )
		if ( $$sref_type_seq eq 'nucleotide' ){
			
			print "\n\t", $$sref_type_seq, "\n\tsequence length: ", $sequ_length ;
			
			for my $taxon ( @$aref_clade ){
				
				if ( (	$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{A} +
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{C} +
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{G} +
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{T} )
						>= 1
					){ 
						push @$aref_taxa_informative, $taxon;
				}
				else { print "\n\t...uninformative sequence found for ", $taxon }
				
				######################
				## Count total number of uninformative sites of taxon
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} =  $sequ_length - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{A} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{C} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{G} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{T} ;
				unless ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} = 0 }
				######################
				
				
				######################
				## Calculate mean proportion of uninformative sites of taxon
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{uninformative} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} / $sequ_length ) ;
				######################
				
				
				######################
				## Count total number of uninformative sites over all clade taxa
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{uninformative} += ( $sequ_length - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{A} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{C} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{G} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{T} ) ;
				######################
			}
		}
		######################
		
		
		
		######################
		## Check for AMINOACID taxon sequences which only consists of missing data, ambiguity characters and indel states
		## These taxa are not further applicated ( not included in @$aref_taxa_informative )
		elsif ( $$sref_type_seq eq 'aminoacid' ){
			
			print "\n\t", $$sref_type_seq, "\n\tsequence length: ", $sequ_length ;
			
			for my $taxon ( @$aref_clade ){
				
				my $N_uninf = $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{X} + $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{MIS} + $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{GAP} ;
				
				
				if ( $N_uninf < $sequ_length ){
					
					push @$aref_taxa_informative, $taxon;
					my $N_inf_char = $sequ_length - $N_uninf ;
					print "\n\tinf taxon: $taxon inf char: $N_inf_char"
				}
				else { print "\n\t...uninformative sequence found for ", $taxon }
				
				######################
				## Count total number of uninformative sites of taxon
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} =  ($href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{X} + $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{MIS} + $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{GAP} );
				unless ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} = 0 }
				######################
				
				
				######################
				## Calculate mean proportion of uninformative sites of taxon
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{uninformative} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} / $sequ_length ) ;
				######################
				
				
				
				######################
				## Count total number of uninformative sites over all clade taxa
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{uninformative} += $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{uninformative} ;
				######################
			}
		}
		
		######################
		## Calculate mean proportion of uninformative sites over all clade taxa
		my $N_taxa_clade = @$aref_clade ;
		my $N_characters_clade = $N_taxa_clade * $sequ_length ;
		$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{uninformative} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{uninformative} / $N_characters_clade ) ;
		######################
		
		
		
		######################
		## Informative taxa info
		my $N_taxa	= @$aref_taxa_informative ;
		$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Ntaxa}{informative} = $N_taxa ;
		######################
		
		print "\n\tN informative taxa: $N_taxa\n" ;
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ INVARIABEL SITE POSITIONS #################################################################################################################
		############################################################################################################################################################################################################################
		print "\n\tCheck invariable states ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile};
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa > 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
				
				
				######################
				## STATE
				for my $state ( 'A', 'C', 'G', 'T' ){
					
					if ( $seen_state{$position}{$state} == @$aref_taxa_informative ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{state} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{state}++;
					}
					else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{state} = '0'}
				}
				######################
				
				
				
				######################
				## CLASS
				if ( $seen_state{$position}{A} + $seen_state{$position}{G} == @$aref_taxa_informative ){ 
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{class} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{class}++
				}
				elsif  ( $seen_state{$position}{C} + $seen_state{$position}{T} == @$aref_taxa_informative ){ 
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{class} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{class}++
				}
				else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{class} = '0' }
				######################
			}
			
			for my $code ( 'state', 'class' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$code} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$code} / $sequ_length ) }
		}
		######################
		
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa < 2
		elsif ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative < 2 ) ){
			
			for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
				
				for my $state ( 'state', 'class' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{$state} = '0' ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$state} = 'NA'
				}
			}
		}
		######################
		
		
		
		######################
		## AMINOACID DATA if informative taxa > 1
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative > 1 ) ){
				
			for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
				
				
				######################
				## STATE
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){
					
					if ( $seen_state{$position}{$state} == @$aref_taxa_informative ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{state} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{state}++;
					}
					else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{state} = '0' }
				}
				######################
				
				
				
				######################
				## STRUCTURE
				if ( $seen_state{$position}{A} + $seen_state{$position}{W} + $seen_state{$position}{M} + $seen_state{$position}{I} + $seen_state{$position}{L} + $seen_state{$position}{F} + $seen_state{$position}{P} == @$aref_taxa_informative ){ 
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{structure} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{structure}++;
				}
				elsif  ( $seen_state{$position}{C} + $seen_state{$position}{G} + $seen_state{$position}{T} + $seen_state{$position}{N} + $seen_state{$position}{Y} + $seen_state{$position}{R} + $seen_state{$position}{S} + $seen_state{$position}{K} + $seen_state{$position}{D} + $seen_state{$position}{V} + $seen_state{$position}{H} + $seen_state{$position}{E} + $seen_state{$position}{Q} == @$aref_taxa_informative ){ 
				
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{structure} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{structure}++;
				}
				else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{structure} = '0' }
				######################
				
				
				
				######################
				## POLARITY
				if ( $seen_state{$position}{A} + $seen_state{$position}{G} + $seen_state{$position}{W} + $seen_state{$position}{M} + $seen_state{$position}{V} + $seen_state{$position}{I} + $seen_state{$position}{L} + $seen_state{$position}{F} + $seen_state{$position}{P} == @$aref_taxa_informative ){ 
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{polarity} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{polarity}++; 
				}
				elsif  ( $seen_state{$position}{C} + $seen_state{$position}{T} + $seen_state{$position}{N} + $seen_state{$position}{Y} + $seen_state{$position}{R} + $seen_state{$position}{S} + $seen_state{$position}{K} + $seen_state{$position}{D} + $seen_state{$position}{H} + $seen_state{$position}{E} + $seen_state{$position}{Q} == @$aref_taxa_informative ){ 
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{polarity} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{polarity}++;
				}
				else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{structure} = '0' }
				######################
				
				
				
				######################
				## CHARGE
				if ( $seen_state{$position}{R} + $seen_state{$position}{K} == @$aref_taxa_informative ){ 
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{charge} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{charge}++;
				}
				elsif  ( $seen_state{$position}{D} + $seen_state{$position}{E} == @$aref_taxa_informative ){ 
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{charge} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{charge}++;
				}
				elsif  ( $seen_state{$position}{A} + $seen_state{$position}{C} + $seen_state{$position}{G} + $seen_state{$position}{T} + $seen_state{$position}{N} + $seen_state{$position}{Y} + $seen_state{$position}{W} + $seen_state{$position}{S} + $seen_state{$position}{M} + $seen_state{$position}{V} + $seen_state{$position}{H} + $seen_state{$position}{I} + $seen_state{$position}{L} + $seen_state{$position}{Q} + $seen_state{$position}{F} + $seen_state{$position}{P} == @$aref_taxa_informative ){ 
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{charge} = '1'; $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{charge}++;
				}
				else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{charge} = '0' }
				######################
			}
			
			for my $code ( 'state', 'structure', 'polarity', 'charge' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$code} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$code} / $sequ_length ) }
		}
		######################
		
		
		
		######################
		## AMINOACID DATA if informative taxa < 2
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative < 2 ) ){
			
			for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
				
				for my $state ( 'state', 'structure', 'polarity', 'charge' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$position}{$state} = '0' ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{pinv}{$state} = 'NA'
				}
			}
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ CALCULATE STATE FREQUENCIES ###############################################################################################################
		############################################################################################################################################################################################################################
		print "\n\tCalculate frequencies ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile} ;
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa > 0
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative > 0 ) ){
			
			###################### N Sites nucleotide data
			for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B','GAP', 'MIS' ){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					if ( $state =~ /^A$|^C$|^G$|^T$/ ){ 
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax}					+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						
						if ( $state =~ /^A$|^T$/ ){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{AT}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{AT}					+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						
						if ( $state =~ /^A$|^G$/ )	{ 
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{purine}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{purine}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} 
						}
						else {
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{pyrimidine}	+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{pyrimidine}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} 
						}
					}
					
					elsif ( $state =~ /^N$|^Y$|^R$|^W$|^S$|^K$|^M$|^D$|^V$|^H$/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{AMB}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{AMB}					+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
					}
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{Nstate}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax}						+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}						+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
				}
			}
			######################
			
			
			
			###################### Mean Sites nucleotide data
			for my $state ( 'A', 'C', 'G', 'T', 'purine', 'pyrimidine', 'AT', 'AMB', 'MIS', 'GAP'){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					if ( $state =~ /^A$|^C$|^G$|^T$|purine|pyrimidine|AT/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} 	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites} ) ;
						
						if ( $state =~ /^A$|^C$|^G$|^T$/){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} 
						}
					}
					
					else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state}	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{Nstate} ) }
				}
				
				if ( $state =~ /^A$|^C$|^G$|^T$|purine|pyrimidine|AT/ ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state}					+= sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}	/ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} ) ;
					
					if ( $state =~ /^A$|^C$|^G$|^T$/){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state}			 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state} / $N_taxa )
					}
				}
				else{
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state}					+= sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}	/ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} )
				}
			}
			
			for my $state ( 'alltax' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{alltax}	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} / $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} ) }
			######################
		}
		######################
		
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa = 0
		elsif ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative == 0 ) ){
				
				for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB', 'alltax'	){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state} = 0		}
				for my $state ( 'GAP', 'AMB', 'MIS', 'alltax'												){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state} = 'NA'	}
		}
		######################
		
		
		
		######################
		## AMINOACID DATA if informative taxa > 0
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative > 0 ) ){
			
			###################### N Sites aminoacid data
			for my $state( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'X', 'GAP', 'MIS' ){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					unless ( $state =~ /^X$|^GAP$|^MIS$/ ){
						
						###### states
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax}						+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						
						###### structure
						if ( $state =~ /^A$|^W$|^M|^I$|^L$|^F$|^P$/ ){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{hydrophobic}	+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{hydrophobic}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						else{
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{hydrophilic}	+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{hydrophilic}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						
						###### polarity
						if ( $state =~ /^A$|^G$|^W|^M$|^V$|^I$|^L$|^F$|^P$/ ){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{nonpolar}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{nonpolar}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						else{ 	
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{polar}			+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{polar}					+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						
						###### charge
						if ( $state =~ /^R$|^K$/ ){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{positive}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{positive}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						elsif ( $state =~ /^D$|^E$/ ){
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{negative}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{negative}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
						else{
							
							$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{neutral}		+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
							$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{neutral}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
						}
					}
					
					###### ambiguity
					elsif ( $state =~ /^X$/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{AMB}				+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{AMB}						+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
					}
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{Nstate}					+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax}							+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}							+= $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state}
				}
			}
			######################
			
			
			
			###################### Mean Sites aminoacid data
			for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'AMB', 'GAP', 'MIS', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'negative', 'neutral' ){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					if ( $state =~ /^A$|^C$|^G$|^T$|^N$|^Y$|^R$|^W$|^S$|^K$|^M$|^D$|^V$|^H$|^I$|^E$|^L$|^Q$|^F$|^P$/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites} ) ;
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state} += $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} ;
					}
					
					elsif ( $state =~ /hydrophobic|hydrophilic|polar|nonpolar|positive|negative|neutral/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites} )
					}
					
					else{
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$state} = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{Nstate} )
					}
				}
				
				if ( $state =~ /^A$|^C$|^G$|^T$|^N$|^Y$|^R$|^W$|^S$|^K$|^M$|^D$|^V$|^H$|^I$|^E$|^L$|^Q$|^F$|^P$|hydrophobic|hydrophilic|polar|nonpolar|positive|negative|neutral/ ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state}	= sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}	/ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} ) ;
					
					if ( $state =~ /^A$|^C$|^G$|^T$|^N$|^Y$|^R$|^W$|^S$|^K$|^M$|^D$|^V$|^H$|^I$|^E$|^L$|^Q$|^F$|^P$/ ){
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state}	= sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state} / $N_taxa )
					}
				}
				else{
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state}	= sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{$state}	/ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} )
				}
			}
			
			for my $state ( 'alltax' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{alltax}	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} / $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} ) }
			######################
		}
		######################
		
		
		
		######################
		## AMINOACID DATA if informative taxa = 0
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative == 0 ) ){
			
			for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'alltax' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state} = 0 }
			for my $state ( 'GAP', 'AMB', 'MIS', 'alltax' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$state} = 'NA' }
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ CHI SQUARE TEST OF HOMOGENEITY ############################################################################################################
		############################################################################################################################################################################################################################
		print "\n\tChi-square test of homogeneity ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile} ;
		
		######################
		## NUCLEOTIDE DATA if informative taxa > 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			######################
			## Homogeneity of Base Frequencies prt1
			for my $state ( 'A', 'C', 'G', 'T', 'AMB', 'MIS', 'GAP' ){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					#####################################
					## Determine expected number of character states for each taxon (N expected) -> equal to PAUP (N expected)!
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state}	= sprintf "%.2f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state} * $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites} ) ;
					#####################################
					
					
					
					#####################################
					## Chi-Square Test -> equal to PAUP (homogeneity of base frequencies)
					unless ($href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} == 0 ){ 
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} += sprintf "%.5f" , ( ( ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} ) ** 2 ) / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} )
					}
					#####################################
				}
			}
			######################
			
			
			
			#####################################
			## Homogeneity of Base Frequencies prt2
			## Calculate N degrees of freedom
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{alltax}	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} / $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} );
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF} = ( @$aref_taxa_informative-1 ) * 3 ;
			#####################################
			
			
			
			##########################################################################
			## Calculate "upper probability of the chi-square distribution (3 degrees "
			## ."of freedom, chi-squared = 6.25): Q = 1-G = $chisprob\n";
			if (($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF} <= 0) || ((abs($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}) - (abs(int($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF})))) != 0)) {
					
					warn "\n\t!WARNING!: Invalid degree of freedom: $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}\n"; # degree of freedom
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{chisprob}  = 'NA';
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} = 'NA';
			}
			else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{chisprob} = &chisqrprob ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}, $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} ) }
			##########################################################################
		}
		
		
		
		######################
		## AMINOACID DATA if informative taxa > 1
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			######################
			## Homogeneity of Base Frequencies prt1
			for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'AMB', 'GAP', 'MIS' ){
				
				for my $taxon ( @$aref_taxa_informative ){
					
					#####################################
					## Determine expected number of character states for each taxon (N expected) -> equal to PAUP (N expected)!
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state}	= sprintf "%.2f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$state} * $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{allsites} ) ;
					#####################################
					
					
					
					#####################################
					## Chi-Square Test -> equal to PAUP (homogeneity of base frequencies)
					unless ($href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} == 0 ){ 
						
						$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} += sprintf "%.5f" , ( ( ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nsites}{$taxon}{$state} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} ) ** 2 ) / $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{Nexp}{$taxon}{$state} )
					}
					#####################################
				}
			}
			######################
			
			
			
			#####################################
			## Homogeneity of Base Frequencies prt2
			## Calculate N degrees of freedom
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{alltax}	 = sprintf "%.4f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nsites}{alltax} / $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{Nstate}{alltax} );
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF} = ( @$aref_taxa_informative-1 ) * 19 ;
			#####################################
			
			
			
			##########################################################################
			## Calculate "upper probability of the chi-square distribution (3 degrees "
			## ."of freedom, chi-squared = 6.25): Q = 1-G = $chisprob\n";
			if (($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF} <= 0) || ((abs($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}) - (abs(int($href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF})))) != 0)) {
					
					warn "\n\t!WARNING!: Invalid degree of freedom: $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}\n"; # degree of freedom
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{chisprob}  = 'NA';
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} = 'NA';
			}
			else{ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{chisprob} = &chisqrprob ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{DF}, $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{ChiSquare} ) }
			##########################################################################
		}
		
		
		
		######################
		## Informative taxa < 2
		if ( @$aref_taxa_informative < 2 ){
			
			for my $state ( 'ChiSquare', 'DF', 'chisprob' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{ht}{$state} = 'NA' }
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ RCFV VALUE ################################################################################################################################
		############################################################################################################################################################################################################################
		print "\n\tRCFV value ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile} ;
		
		######################
		## NUCLEOTIDE DATA if informative taxa > 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			#####################################
			## RCFV Value and absolute deviation of each clade and taxon -> equal to zhong et al. 2011
			for my $taxon (  @$aref_taxa_informative ){
				
				for my $char ( 'A', 'C', 'G', 'T' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV}	+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char} );
				}
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV}		 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV} / $$sref_N_tax_complete );
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{RCFV}				+=					  $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV};
				
				for my $char ( 'purine', 'pyrimidine', 'AMB', 'GAP', 'AT' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char} );
				}
			}
			
			for my $char ( 'A', 'C', 'G', 'T' ){
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			 	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char} / $$sref_N_tax_complete )
			}
			
			for my $char ( 'purine', 'pyrimidine', 'AMB', 'GAP', 'AT' ){
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}				 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char} / $$sref_N_tax_complete )
			}
			#####################################
		}
		######################
		
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa < 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative < 2 ) ){
			
			for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB', 'MIS', 'alltax' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$state} = 'NA' }
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{RCFV} = 'NA'
		}
		######################
		
		
		
		
		######################
		## AMINOACID DATA if informative taxa > 1
		elsif ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			#####################################
			## RCFV Value and absolute deviation of each clade and taxon -> equal to zhong et al. 2011
			for my $taxon (  @$aref_taxa_informative ){
				
				for my $char ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV}	+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{ht}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char} );
				}
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV}		 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV} / $$sref_N_tax_complete );
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{RCFV}				+=					  $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{RCFV};
				
				for my $char ( 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'neutral', 'positive', 'negative', 'AMB', 'GAP' ){
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			+= abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = abs 				( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$char} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$char} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char}	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{absol_dev}{$taxon}{$char} );
				}
			}
			
			for my $char ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}			 	 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char} / $$sref_N_tax_complete )
			}
			
			for my $char ( 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'neutral', 'positive', 'negative', 'AMB', 'GAP' ){
				
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char}				 = sprintf "%.4f",	( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$char} / $$sref_N_tax_complete )
			}
			#####################################
		}
		######################
		
		
		
		######################
		## AMINOACID DATA if informative taxa < 1
		if ( ( $$sref_type_seq =~ /aminoacid/ )  &&  ( @$aref_taxa_informative < 2 ) ){
			
			for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'GAP', 'AMB', 'MIS', 'alltax', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{absol_dev}{$state} = 'NA' }
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{RCFV} = 'NA'
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		
		
		############################################################################################################################################################################################################################
		################################################################################ SKEW VALUES ###############################################################################################################################
		############################################################################################################################################################################################################################
		
		######################
		## NUCLEOTIDE DATA if informative taxa > 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative > 1 ) ){
			
			print "\n\tSKEW values ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile} ;
			
			for my $skew ( 'AG', 'CT', 'AT', 'GC' ){
				
				my @skew_states = split "", $skew ;
				
				############################# hier anstelle clade taxon und ausgeben in taxon_base_frequencies_all_partitions/RCFV
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{skew}{$skew} = ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$skew_states[0]} - $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$skew_states[1]} ) / ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$skew_states[0]} + $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{MeanSites}{$skew_states[1]} ) ;
				$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{skew}{$skew} = sprintf "%.2f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{skew}{$skew} ) ;
			}
			
			for my $taxon (  @$aref_taxa_informative ){
				
				for my $skew ( 'AG', 'CT', 'AT', 'GC' ){
					
					my @skew_states = split "", $skew ;
					
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{skew}{$taxon}{$skew} = ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$skew_states[0]} - $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$skew_states[1]} ) / ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$skew_states[0]} + $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{MeanSites}{$taxon}{$skew_states[1]} ) ;
					$href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{skew}{$taxon}{$skew} = sprintf "%.2f" , ( $href_frequencies->{$$sref_infile}{$$sref_counter}{taxon}{skew}{$taxon}{$skew} ) ;
				}
			}
		}
		######################
		
		
		
		######################
		## NUCLEOTIDE DATA if informative taxa < 1
		if ( ( $$sref_type_seq =~ /nucleotide/ )  &&  ( @$aref_taxa_informative < 1 ) ){
			
			print "\n\tSKEW values ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile} ;
			
			for my $skew ( 'AG', 'CT', 'AT', 'GC' ){ $href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{skew}{$skew} = 'NA' }
		}
		######################
		
		############################################################################################################################################################################################################################
		############################################################################################################################################################################################################################
		
		sub chisqrprob { # Upper probability   X^2(x^2,n)
			
			my ($n,$x) = @_;
			return &precision_string( &subchisqrprob($n, $x))
		}
		
		sub log10 {
			
			my $n = shift;
			return log($n) / log(10);
		}
		
		sub precision_string {
			
			my ($x) = @_;
			#return sprintf "%." . &precision($x) . "f", $x;
			if ($x) { return sprintf "%." . &precision($x) . "f", $x } else { return "0" }
		}
		
		sub precision {
			
			my ($x) = @_;
			return abs int(log10(abs $x) - SIGNIFICANT);
		}
		
		sub subchisqrprob {
			
			my ($n,$x) = @_;
			my $p;
			
			if 		( $x <= 0  ) { $p = 1 } 
			elsif 	( $n >  100) { $p = &subuprob((($x / $n) ** (1/3) - (1 - 2/9/$n)) / sqrt(2/9/$n)) }
			elsif 	( $x >  400) { $p = 0 }
			else 	{
				
				my ($a, $i, $i1);
				if (($n % 2) != 0) {
					
					$p  = 2 * &subuprob(sqrt($x));
					$a  = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
					$i1 = 1;
				}
				else { $p = $a = exp(-$x/2) and $i1 = 2 }
				
				for  ( $i = $i1; $i <= ($n-2); $i += 2) { $a *= $x / $i and $p += $a }
			}
			
			return $p;
		}
		
		sub subuprob {
			
			my ($x) = @_;
			my $p = 0; # if ($absx > 100)
			my $absx = abs($x);
			
			if ($absx < 1.9) {
				
				$p = (1 +
					$absx * (.049867347
						+ $absx * (.0211410061
							+ $absx * (.0032776263
								+ $absx * (.0000380036
									+ $absx * (.0000488906
										+ $absx * .000005383)))))) ** -16/2;
			}
			elsif ($absx <= 100) {
				
				for (my $i = 18; $i >= 1; $i--) { $p = $i / ($absx + $p) }
				$p = exp(-.5 * $absx * $absx) / sqrt(2 * PI) / ($absx + $p);
			}
			
			$p = 1 - $p if ($x<0);
			return $p;
		}
	}
	
	
	sub pairwise_comparisons_aa{
		
		my $sref_infile				= $_[0] ;
		my $aref_clade				= $_[1] ;
		my $href_frequencies			= $_[2] ;
		my $sref_counter				= $_[3] ;
		my $href_ranges_of_part		= $_[4] ;
		my $href_state_of_tax_pos	= $_[5] ;
		my $href_clade_of_numbers	= $_[6] ;
		
		
		
		
		my @taxa	= sort @$aref_clade ;
		
		print "\n\tSequence comparisons ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile}, "\n";
		
		
		##################
		## Generate resultfolder for sequence comparisons
		mkdir "BaCoCa_Results/missing_data_overlap"		, 0755 ;
		##################
		
		
		##################
		## Generate outputfiles for sequence comparisons (missing data positive overlap (po), and missing data negative overlap (no) separately)
		## Print first line for pairwise matrix
		my $outfile_po_over		= "results_missing_data_positive_overlap_".$$sref_infile."_".$$sref_counter.".txt" ;
		my $outfile_no_over		= "results_missing_data_negative_overlap_".$$sref_infile."_".$$sref_counter.".txt" ;
		open OUTpo,	">BaCoCa_Results/missing_data_overlap/$outfile_po_over"		|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/missing_data_overlap/$outfile_po_over !\n" ;
		open OUTno,	">BaCoCa_Results/missing_data_overlap/$outfile_no_over"		|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/missing_data_overlap/$outfile_no_over !\n" ;
		
		#for my $taxon ( @taxa ){ print OUTp "\t", $taxon ; print OUTtitv "\t", $taxon ; print OUTpo "\t", $taxon ; print OUTno "\t", $taxon }
		for my $taxon ( @taxa ){ print OUTp "\t", $taxon ; print OUTpo "\t", $taxon ; print OUTno "\t", $taxon }
		print OUTpo		"\n" ; print OUTno	"\n" ;
		##################
		
		my %hoh_pvalues_of_first_second_taxon ;
		my $outputspace ;
		my $N_sequence_positions ;
		my %prefix_string_of_taxon ;
		
		while ( @taxa ){
			
			$outputspace = $outputspace."\t" ;
			my $first_taxon = shift @taxa ;
			
			$prefix_string_of_taxon{no}{$first_taxon} = $prefix_string_of_taxon{no}{$first_taxon};
			print OUTpo $first_taxon, $prefix_string_of_taxon{po}{$first_taxon}, "\t" ;
			print OUTno $first_taxon, $prefix_string_of_taxon{no}{$first_taxon}, "\t" ;
			
			for my $second_taxon ( @taxa ){
				
				##################
				## counter of missing data positives and negatives have to be set zero at the beginning of each sequence comparison to avoid empty cells in pairwise output matrices
				$hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon} = 0 ;
				$hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon} = 0 ;
				##################
				
				
				my $counter_position;
				for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
					
					if 		( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^GAP$|^MIS$|^X$/i ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^GAP$|^MIS$|^X$/i ) ){ $hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon}++ }
					unless ( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^GAP$|^MIS$|^X$/i ) || ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^GAP$|^MIS$|^X$/i ) ){ $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					
					$counter_position++
				}
				
				#print	"Taxon 1: ", $first_taxon, "\tTaxon 2: ", $second_taxon, "\nPositive Overlap MIS: ", $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon},
				#"\nNegative Overlap MIS: ", $hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon}, "\n\n" ;
				
				$N_sequence_positions = $counter_position ;
				
				
				##################
				## Calculation of p-values of single pairwise comparisons, missing data positive overlap (po), and missing data negative overlap (no)
				#$hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} / $counter_position ) ;
				$hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} / $counter_position ) ;
				$hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} / $counter_position ) ;
				##################
				
				$prefix_string_of_taxon{po}{$second_taxon} = $prefix_string_of_taxon{po}{$second_taxon}."\t".$hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} ;
				$prefix_string_of_taxon{no}{$second_taxon} = $prefix_string_of_taxon{no}{$second_taxon}."\t".$hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} ;
				
				
				##################
				## print OUT of single comparion results
				print OUTpo	"\t", $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon} ;
				print OUTno	"\t", $hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon} ;
				##################
			}
			
			print OUTpo		"\n" ; print OUTno	"\n" ;
		}
		close OUTpo and close OUTno ;
		
		
		
		my $outfile_po = "BaCoCa_Results/missing_data_overlap/$outfile_po_over" ;
		my $outfile_no = "BaCoCa_Results/missing_data_overlap/$outfile_no_over" ;
		##################
	}
	
	
	sub pairwise_comparisons_nuc{
		
		my $sref_infile				= $_[0] ;
		my $aref_clade				= $_[1] ;
		my $aref_clade_inf			= $_[2] ;
		my $href_frequencies		= $_[3] ;
		my $sref_counter			= $_[4] ;
		my $href_ranges_of_part		= $_[5] ;
		my $href_state_of_tax_pos	= $_[6] ;
		my $href_clade_of_numbers	= $_[7] ;
		
		print "\n\tMissing data overlap ",  $$sref_infile, ", ", $href_clade_of_numbers->{$$sref_counter}{$$sref_infile}, "\n";
		
		
		my @taxa	= sort @$aref_clade ;
		my $N_taxa	= @$aref_clade ;
		
		my %seen_inf;
		for my $taxon ( @$aref_clade_inf ){ $seen_inf{$taxon}++ }
		
		##################
		## Generate resultfolder for sequence comparisons
		mkdir "BaCoCa_Results/c_value_calculations"	, 0755 ;
		mkdir "BaCoCa_Results/missing_data_overlap"	, 0755 ;
		##################
		
		
		##################
		## Generate outputfiles for sequence comparisons (p_value, ti/tv, missing data positive overlap (po), and missing data negative overlap (no) separately)
		## Print first line for pairwise matrix
		my $outfile_p_dist		= "results_pairwise_comparisons_pvalues".$$sref_infile."_".$$sref_counter.".txt" ;
		my $outfile_titv_dist	= "results_pairwise_comparisons_titv".$$sref_infile."_".$$sref_counter.".txt" ;
		my $outfile_po_over		= "results_missing_data_positive_overlap_".$$sref_infile."_".$$sref_counter.".txt" ;
		my $outfile_no_over		= "results_missing_data_negative_overlap_".$$sref_infile."_".$$sref_counter.".txt" ;
		
		if ( @$aref_clade_inf > 1 ){ open OUTp,	">BaCoCa_Results/c_value_calculations/$outfile_p_dist"		|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/c_value_calculations/$outfile_p_dist !\n" }
		if ( @$aref_clade_inf > 1 ){ open OUTtitv,	">BaCoCa_Results/c_value_calculations/$outfile_titv_dist"	|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/c_value_calculations/$outfile_titv_dist !\n" }
		open OUTpo,	">BaCoCa_Results/missing_data_overlap/$outfile_po_over"		|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/missing_data_overlap/$outfile_po_over !\n" ;
		open OUTno,	">BaCoCa_Results/missing_data_overlap/$outfile_no_over"		|| die "\n\t!OUTPUT-ERROR!: Can not open file BaCoCa_Results/missing_data_overlap/$outfile_no_over !\n" ;
		
		for my $taxon ( @taxa ){ print OUTpo "\t", $taxon ; print OUTno "\t", $taxon }
		if ( @$aref_clade_inf > 1 ){ for my $taxon ( sort @$aref_clade_inf ){ print OUTp "\t", $taxon ; print OUTtitv "\t", $taxon } }
		print OUTp		"\n" ; print OUTpo		"\n" ; if ( @$aref_clade_inf > 1 ){ print OUTno	"\n" ; print OUTtitv	"\n" }
		##################
		
		
		##################
		## Pairwise sequence comparisons to determine p-values, ti/tv ratios, missing data positive overlap (po), and missing data negative overlap (no) of single comparisons
		## total p-values and ti/tv ratio summarized over all comparisons
		my %hoh_pvalues_of_first_second_taxon ;
		my $N_comparisons		= 0 ;
		my $N_comparisons_inf	= 0 ;
		my ( $total_p, $total_ti_tv_ratio ) ;
		my $outputspace ;
		my $outputspace_inf ;
		my $N_sequence_positions ;
		my $N_comparisons_total;
		my %prefix_string_of_taxon ;
		my %seen_values_po;
		my %seen_values_no;
		
		while ( @taxa ){
			
			$outputspace			= $outputspace."\t" ;
			
			
			my $first_taxon = shift @taxa ;
			
			$prefix_string_of_taxon{no}{$first_taxon} = $prefix_string_of_taxon{no}{$first_taxon};
			print OUTpo $first_taxon, $prefix_string_of_taxon{po}{$first_taxon}, "\t" ;
			print OUTno $first_taxon, $prefix_string_of_taxon{no}{$first_taxon}, "\t" ;
			
			if ( ( $seen_inf{$first_taxon} ) ){
				
				$outputspace_inf = $outputspace_inf."\t" ;
				print OUTp $first_taxon, $outputspace_inf ; print OUTtitv $first_taxon, $outputspace_inf
			}
			
			for my $second_taxon ( @taxa ){
				
				$N_comparisons++ ;
				if (( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $N_comparisons_inf++ }
				
				##################
				## counter of p-values, missing data positives and negatives have to be set zero at the beginning of each sequence comparison to avoid empty cells in pairwise output matrices
				## counter of ti and tv are extra handled (see below)
				$hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} = 0 ;
				$hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} = 0 ;
				$hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} = 0 ;
				##################
				
				
				my $counter_position;
				for my $position ( $href_ranges_of_part->{$$sref_infile}{start} .. $href_ranges_of_part->{$$sref_infile}{end} ){
					
					if ( ( $href_state_of_tax_pos->{$first_taxon}{$position} eq $href_state_of_tax_pos->{$second_taxon}{$position} ) && ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon}++ ; $total_p++ }
					
					if 		( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^A$|^G$/i ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^C$|^T$/i ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{tv}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^C$|^T$/i ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^A$|^G$/i ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{tv}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^C$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^T$/i     ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^T$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^C$/i     ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^A$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^G$/i     ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^G$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^A$/i     ) ){ if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon}++ } $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^G$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^G$/i     ) ){ $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^A$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^A$/i     ) ){ $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^C$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^C$/i     ) ){ $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^T$/i     ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^T$/i     ) ){ $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}++ }
					
					
					elsif 	( ( $href_state_of_tax_pos->{$first_taxon}{$position} =~ /^GAP$|^MIS$|^N$/i ) && ( $href_state_of_tax_pos->{$second_taxon}{$position} =~ /^GAP$|^MIS$|^N$/i ) ){ $hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon}++ }
					
					$counter_position++
				}
				
				
				##################
				## Calculation of p-values of single pairwise comparisons, missing data positive overlap (po), and missing data negative overlap (no)
				if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){ $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} / $counter_position ) }
				
				$hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} / $counter_position ) ;
				$hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} = ( $hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} / $counter_position ) ;
				##################
				
				
				#print	"Taxon 1: ", $first_taxon, "\tTaxon 2: ", $second_taxon, "\nP-value: ",  $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon}, "\nTi-value: ",
				#		$hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon}, "\nTv-value: ", $hoh_pvalues_of_first_second_taxon{tv}{$first_taxon}{$second_taxon}, "\n\n" ;
				
				
				if ( ( $seen_inf{$first_taxon} )  && ( $seen_inf{$second_taxon} ) ){
					
					##################
					## Calculation of ti/tv for single pairwise comparisons
					## if value of transversions equal zero, the number of transitions is devided by 0.001 to avoid devition error
					## if the number of transitions equal zero, the ti/tv ratio is also zero
					if ( ( $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon} ) && ( $hoh_pvalues_of_first_second_taxon{tv}{$first_taxon}{$second_taxon} ) ){
						
						$hoh_pvalues_of_first_second_taxon{titv}{$first_taxon}{$second_taxon} = $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon} / $hoh_pvalues_of_first_second_taxon{tv}{$first_taxon}{$second_taxon}
					}
					elsif ( $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon} ){
						
						$hoh_pvalues_of_first_second_taxon{titv}{$first_taxon}{$second_taxon} = $hoh_pvalues_of_first_second_taxon{ti}{$first_taxon}{$second_taxon} / 0.001
					}
					else{ $hoh_pvalues_of_first_second_taxon{titv}{$first_taxon}{$second_taxon} = 0 }
					##################
					
					
					##################
					## total value of ti/tv summarized over all comparisons
					$total_ti_tv_ratio += $hoh_pvalues_of_first_second_taxon{titv}{$first_taxon}{$second_taxon} ;
					##################
					
					
					##################
					## print OUT of single comparion results
					if ( @$aref_clade_inf > 1 ){
						
						print OUTp		"\t", $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} ;
						print OUTtitv	"\t", $hoh_pvalues_of_first_second_taxon{titv   }{$first_taxon}{$second_taxon}
					}
					
					$N_comparisons_total += $counter_position ;
				}
				
				$prefix_string_of_taxon{po}{$second_taxon} = $prefix_string_of_taxon{po}{$second_taxon}."\t".$hoh_pvalues_of_first_second_taxon{po     }{$first_taxon}{$second_taxon} ;
				$prefix_string_of_taxon{no}{$second_taxon} = $prefix_string_of_taxon{no}{$second_taxon}."\t".$hoh_pvalues_of_first_second_taxon{no     }{$first_taxon}{$second_taxon} ;
				
				$seen_values_po{$hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon}}++ ;
				$seen_values_no{$hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon}}++ ;
				
				print OUTpo	"\t", $hoh_pvalues_of_first_second_taxon{po}{$first_taxon}{$second_taxon} ;
				print OUTno	"\t", $hoh_pvalues_of_first_second_taxon{no}{$first_taxon}{$second_taxon} ;
				##################
				
				$N_sequence_positions = $counter_position ;
			}
			
			if ( ( @$aref_clade_inf > 1 ) && ( $seen_inf{$first_taxon} ) ){	print OUTp		"\n" ; print OUTtitv	"\n" }
																				print OUTpo	"\n" ; print OUTno		"\n" ;
		}
		close OUTpo and close OUTno ;
		
		
		my $outfile_po = "BaCoCa_Results/missing_data_overlap/$outfile_po_over" ;
		my $outfile_no = "BaCoCa_Results/missing_data_overlap/$outfile_no_over" ;
		
		my $N_diff_po = keys %seen_values_po ;
		my $N_diff_no = keys %seen_values_no ;
		##################
		
		
		
		if ( @$aref_clade_inf > 2 ){
			
			##################
			## Calculation of mean p-values and mean ti/tv ratio
			my $mean_p						= $total_p  / $N_comparisons_total;
			my $mean_ratio_total_ti_tv	= $total_ti_tv_ratio / $N_comparisons_inf ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{meanvalue}{titv  } = $mean_ratio_total_ti_tv ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{meanvalue}{pvalue} = $mean_p ;
			##################
			
			
			##################
			## 
			my @taxa_inf = sort @$aref_clade_inf ;
			my ( $sq_total_p , $sq_total_ti_tv ) ;
			while ( @taxa_inf ){
				
				my $first_taxon = shift @taxa_inf ; for my $second_taxon ( @taxa_inf ){
					
					$sq_total_p		+= ( $hoh_pvalues_of_first_second_taxon{pvalues}{$first_taxon}{$second_taxon} - $mean_p )**2 ;
					$sq_total_ti_tv	+= ( $hoh_pvalues_of_first_second_taxon{titv   }{$first_taxon}{$second_taxon} - $mean_ratio_total_ti_tv )**2 ;
				}
			}
			##################
			
			
			##################
			## Calculation of standard derivation for p-values and ti/tv ratios
			my $standard_derivation_p		= sqrt ( $sq_total_p / ( $N_comparisons_inf - 1 ) ) ;
			my $standard_derivation_ti_tv	= sqrt ( $sq_total_ti_tv / ( $N_comparisons_inf - 1 ) ) ;
			
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{c_value      }{total } = ( $standard_derivation_ti_tv / $standard_derivation_p ) ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{stdderivation}{pvalue} = $standard_derivation_p ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{stdderivation}{titv  } = $standard_derivation_ti_tv ;
			##################
			
			
			##################
			## Print SQvalue and Standard Derivation
			#print OUTp		"\n\nSQ Total Mean:\t", $sq_total_p,		"\nStandard Derivation:\t", $standard_derivation_p ;
			#print OUTtitv	"\n\nSQ Total Mean:\t", $sq_total_ti_tv,	"\nStandard Derivation:\t", $standard_derivation_ti_tv ;
			close OUTp and close OUTtitv ;
			##################
		}
		
		elsif ( @$aref_clade_inf < 3 ){
			
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{c_value      }{total } = 'NA' ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{stdderivation}{pvalue} = 'NA' ;
			$href_frequencies->{$$sref_infile}{$$sref_counter}{clade}{stdderivation}{titv  } = 'NA' ;
		}
	}
}


sub print_out{
	
	my $href_freq_of_file_char	=	$_[0] ;
	my $href_clade_of_number		=	$_[1] ;
	my $href_file_of_number		=	$_[2] ;
	my $href_taxa_of_file_clade	=	$_[3] ;
	my $aref_all_taxa				=	$_[4] ;
	my $sref_seqtype				=	$_[5] ;
	my $href_seq_ranges_of_part	=	$_[6] ;
	
	#my $heat_map_infile_suffix = ".txt" ;
	#if ( $$sref_r_heatmaps eq 'yes' ){ $heat_map_infile_suffix = ".csv" }
	
	#print "\n\n\tGenerate outputfolder..." ;
	####################################
	## Generate OUTputfolder
	mkdir "BaCoCa_Results/chisquare_test_homogeneity_taxa"			, 0755 ;
	mkdir "BaCoCa_Results/taxon_basefrequencies_single_partions"	, 0755 ;
	mkdir "BaCoCa_Results/taxon_basefrequencies_all_partions"		, 0755 ;
	mkdir "BaCoCa_Results/compositional_bias"						, 0755 ;
	mkdir "BaCoCa_Results/invariant_alignment_positions"			, 0755 ;
	mkdir "BaCoCa_Results/invariant_alignment_positions/fasta"		, 0755 ;
	mkdir "BaCoCa_Results/invariant_alignment_positions/txt"		, 0755 ;
	mkdir "BaCoCa_Results/invariant_alignment_positions/svg"		, 0755 ;
	
	####################################
	
	
	print "\n\tPrint results..." ;
	######################################################################################################################################################################################################################## Output Nucleotid data
	if ( $$sref_seqtype =~ /nucleotide/ ){
		
		open	OUTall,	">BaCoCa_Results/summarized_frequencies.txt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/summarized_frequencies.txt !\n" ;
		print	OUTall	"Basefrequencies\nFile\tClade\tN taxa (informative)\tp(A)\tp(C)\tp(G)\tp(T|U)\tp(AT)\tp(Purine)\tp(Pyrimidine)\tcsRCFV (A)\tcsRCFV (C)\tcsRCFV (G)\tcsRCFV (T|U)\tcsRCFV (AT)\tcsRCFV (Purine)\tcsRCFV (Pyrimidine)\tSkew (A-G)\tSkew (C-T)\tSkew (A-T)\tSkew (G-C)\tChi-square value\tDegrees of freedom\tp-value\tRCFV Value\tc-value\tMean p-distances\tMean ti/tv\tStd. Deviation p-distances\tStd. Deviation ti/tv\tp(Sites)\tp(Gap)\tp(Ambiguity)\tcsRCFV (Gap)\tcsRCFV (Ambiguity)\tp_uninformative sites\tp_inv(State)\tp_inv(Class)\n";
		
		for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ 
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){
				
				mkdir "BaCoCa_Results/skew_values"					, 0755 ;
	
				print "\n\t", $href_file_of_number->{$file_number}, ", ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}} ;
				############################################################################################################ begin summarized all for nc
				## Print OUT summarized results for nucleotide data
				print OUTall $href_file_of_number->{$file_number}, "\t", $clade_number, " ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}, "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{Ntaxa}{informative} ;
				
				for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} }
				
				for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{absol_dev}{$state} }
				
				for my $skew ( 'AG', 'CT', 'AT', 'GC' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{skew}{$skew} }
				
				for my $state ( 'ChiSquare', 'DF', 'chisprob' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{$state} }
				
				for my $state ( 'RCFV' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} }
				
				for my $state ( 'total' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{c_value}{$state} }
				
				for my $state ( 'pvalue', 'titv' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{meanvalue}{$state} }
				
				for my $state ( 'pvalue', 'titv' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{stdderivation}{$state} }
				
				for my $state ( 'alltax', 'GAP', 'AMB' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} }
				
				for my $state ( 'GAP', 'AMB' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{absol_dev}{$state} }
				
				for my $state ( 'uninformative' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} }
				
				for my $state ( 'state', 'class' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{pinv}{$state} }
				
				print OUTall "\n";
				############################################################################################################ end summarized all for nc
				
				
				
				############################################################################################################ begin homogeneity among taxa for nc
				## Print OUT Results of chis-quare tests of homogeneity among taxa
				my $outfile_ht = $href_file_of_number->{$file_number}."_clade".$clade_number."_".$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}.".txt" ;
				open OUTht,  ">BaCoCa_Results/chisquare_test_homogeneity_taxa/$outfile_ht" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/chisquare_test_homogeneity_taxa/$outfile_ht !\n" ;
				
				print OUTht "File: ", $href_file_of_number->{$file_number}, "\tClade: ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}, "\nBase frequencies of taxa:\nTaxon\tp(A)\tp(C)\tp(G)\tp(T|U)\tN sites\n";
				
				
				my $aref_taxa = exists($href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}}) ? \@{$href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}}} :( ) ;
				for my $taxon ( @$aref_taxa ){
					
					print OUTht	$taxon, "\t",	$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{MeanSites}{$taxon}{A}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{MeanSites}{$taxon}{C}, "\t",
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{MeanSites}{$taxon}{G}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{MeanSites}{$taxon}{T}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{allsites}, "\n"
				}
				
				print OUTht		"Mean\t",		$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{A}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{C}, "\t",
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{G}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{T}, "\t", 
													$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{alltax}, "\n\n",
													"Chi-Square test of homogeneity of base frequencies across taxa:\n\nTaxon\tp(A)\tp(C)\tp(G)\tp(T|U)\n" ;
				
				for my $taxon ( @$aref_taxa ){
					
					print OUTht	$taxon, "Observed\t",	$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{A}, "\t", 
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{C}, "\t",
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{G}, "\t", 
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{T}, "\n",
								$taxon, "Expected\t",		$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nexp}{$taxon}{A}, "\t",
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nexp}{$taxon}{C}, "\t",
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nexp}{$taxon}{G}, "\t",
															$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nexp}{$taxon}{T}, "\n"
				}
				
				print OUTht	"\nChi-Square = ",	$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{ChiSquare}, 
								"(df = ", 				$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{DF},
								"), P = ",				$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{chisprob} ;
								if ($href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{chisprob}<0.01) {print OUTht "\nNOTE: Given a significance level of 0.01 the composition is significantly diviating from homogeneity!\n"} else {print OUTht "\nNOTE: Given a significance level of 0.01 the composition is NOT significantly diviating from homogeneity!\n"}
				print OUTht	"\nWarning: This test ignores correlation due to phylogenetic structure.\n";
				############################################################################################################ end homogeneity among taxa for nc
			}
			
			
			
			
			############################################################################################################ begin taxon frequencies of single partitions for nc
			if ( $file_number == 1 ){
				
				my $outfile_tf_prt = $href_file_of_number->{$file_number}."_summarized_taxon_basefrequencies.txt" ;
				open	OUTprt,	">BaCoCa_Results/taxon_basefrequencies_single_partions/$outfile_tf_prt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/taxon_basefrequencies_single_partions/$outfile_tf_prt !\n" ;
				print	OUTprt	"Basefrequencies\nFile\tTaxon\tp(A)\tp(C)\tp(G)\tp(T|U)\tp(AT)\tp(Purine)\tp(Pyrimidine)\tp(Gap)\tp(Ambiguity)\tRCFV Value (taxon specific)\tAbsolute Deviation (A)\tAbsolute Deviation (C)\tAbsolute Deviation (G)\tAbsolute Deviation (T|U)\tAbsolute Deviation (AT)\tAbsolute Deviation (Purine)\tAbsolute Deviation (Pyrimidine)\tAbsolute Deviation (Gap)\tAbsolute Deviation (Ambiguity)\tp_uninformative sites\n" ;
				
				my $aref_taxa = exists($href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{0}{$href_file_of_number->{$file_number}}}) ? \@{$href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{0}{$href_file_of_number->{$file_number}}}} :( ) ;
				for my $taxon ( @$aref_taxa ){
					
					print OUTprt	$href_file_of_number->{$file_number}, "\t", $taxon ;
					for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB', 'RCFV' ){ print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state} }
					for my $state ( 'A', 'C', 'G', 'T', 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB' ){ print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{absol_dev}{$taxon}{$state} }
					for my $state ( 'uninformative' ){ print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state} }
					print OUTprt "\n"
				}
				close OUTprt
			}
			############################################################################################################ end taxon frequencies of single partitions for nc
		}
		
		
		
		############################################################################################################ begin taxon frequencies (specific) of all partitions for nc
		for my $state ( 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB', 'RCFV', 'MIS', 'uninformative' ){
			
			my $outfile	= "BaCoCa_Results/taxon_basefrequencies_all_partions/".$state."_frequencies_all_partitions.txt";
			open	OUTprt,	">$outfile" or die "\n\t!OUTPUT-ERROR!: Cannot open file $outfile !\n" ;
			
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ print OUTprt "\t", $href_file_of_number->{$file_number} } print OUTprt "\n" ;
			
			for my $taxon ( @$aref_all_taxa ){ print OUTprt $taxon ;
				
				for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
					
					print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state} ;
				}
				print OUTprt "\n";
			}
			close OUTprt;
		}
		
		for my $state ( 'AT', 'purine', 'pyrimidine', 'GAP', 'AMB' ){
			
			my $outfile = "BaCoCa_Results/compositional_bias/".$state."_absolute_deviation_frequencies_all_partitions.txt";
			open	OUTprt,	">$outfile" or die "\n\t!OUTPUT-ERROR!: Cannot open file $outfile !\n" ;
			
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ print OUTprt "\t", $href_file_of_number->{$file_number} } print OUTprt "\n" ;
			
			for my $taxon ( @$aref_all_taxa ){ print OUTprt $taxon ;
				
				for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
					
					print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{absol_dev}{$taxon}{$state} ;
				}
				print OUTprt "\n"
			}
			close OUTprt ;
		}
		
		for my $skew ( 'AG', 'CT', 'AT', 'GC' ){
			
			my $outfile = "BaCoCa_Results/skew_values/".$skew."_skew_values_all_partitions.txt";
			open	OUTprt,	">$outfile" or die "\n\t!OUTPUT-ERROR!: Cannot open file $outfile !\n" ;
			
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ print OUTprt "\t", $href_file_of_number->{$file_number} } print OUTprt "\n" ;
			
			for my $taxon ( @$aref_all_taxa ){ print OUTprt $taxon ;
				
				for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ 
					
					print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{skew}{$taxon}{$skew};
				}
				print OUTprt "\n"
			}
			close OUTprt ;
		}
		############################################################################################################ end taxon frequencies (specific) of all partitions for nc
		
		
		
		############################################################################################################ begin invariant alignment partitions for nc
		for my $state ( 'state', 'class' ){
			
			my $outfile = $state."_invariant_msa_positions";
			open	OUTfas,	">BaCoCa_Results/invariant_alignment_positions/fasta/${outfile}.fas" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/invariant_alignment_positions/${outfile}.fas !\n" ;
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){
				
				print OUTfas ">", $href_clade_of_number->{$clade_number}{$href_file_of_number->{1}}, "\n";
				
				for my $pos ( $href_seq_ranges_of_part->{$href_file_of_number->{1}}{start} .. $href_seq_ranges_of_part->{$href_file_of_number->{1}}{end} ){ #print $pos,"\n";####################################################
					
					if ( $href_freq_of_file_char->{$href_file_of_number->{1}}{$clade_number}{pinv}{$pos}{$state} == 1 ){ print OUTfas 'T' }else { print OUTfas 'A' }
				} print OUTfas "\n"
			}
			close OUTfas;
			
			open	OUTinv,	">BaCoCa_Results/invariant_alignment_positions/txt/${outfile}.txt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/invariant_alignment_positions/${outfile}.txt !\n" ;
			print	OUTinv	"Invariant alignment positions of ", $state, "\n\n" ;
			
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){ print OUTinv "\tClade ", $clade_number, ": ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{1}}} print OUTinv "\n";
			
			for my $pos ( $href_seq_ranges_of_part->{$href_file_of_number->{1}}{start} .. $href_seq_ranges_of_part->{$href_file_of_number->{1}}{end} ){ print OUTinv $pos+1;
				
				for my $clade_number ( sort {$a<=>$b} keys %$href_clade_of_number ){
					
					print OUTinv "\t", $href_freq_of_file_char->{$href_file_of_number->{1}}{$clade_number}{pinv}{$pos}{$state}
				} print OUTinv "\n" ;
			}
			close OUTinv;
			
			################################################
			## svg output
			for my $file_number (sort {$a<=>$b} keys %$href_file_of_number ){
				
				my $svg_filename = "invariant_alignment_positions/svg/".$state."_". $href_file_of_number->{$file_number}."_invariant_msa_positions.svg" ;
				my @position_numbers ;
				
				for my $pos ( $href_seq_ranges_of_part->{$href_file_of_number->{$file_number}}{start} .. $href_seq_ranges_of_part->{$href_file_of_number->{$file_number}}{end} ){ $pos += 1; push @position_numbers, $pos }
				
				&svg_plot ( \%$href_clade_of_number, \%$href_freq_of_file_char, \%$href_file_of_number, \$state, \$svg_filename, \$file_number, \@position_numbers )
			}
		}
		############################################################################################################ end invariant alignment partitions for nc
	}
	
	
	
	
	
	
	######################################################################################################################################################################################################################## Output Aminoacid data
	else{
		open	OUTall,	">BaCoCa_Results/summarized_frequencies.txt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/summarized_frequencies.txt !\n" ;
		print	OUTall	"Basefrequencies\nFile\tClade\tN taxa (informative)\tp(A)\tp(V)\tp(L)\tp(I)\tp(P)\tp(M)\tp(F)\tp(W)\tp(G)\tp(S)\tp(T)\tp(C)\tp(N)\tp(Q)\tp(Y)\tp(D)\tp(E)\tp(K)\tp(R)\tp(H)\tp(Hydrophobic)\tp(Hydrophilic)\tp(Polar)\tp(Nonpolar)\tp(Positive)\tp(Neutral)\tp(Negative)\tRCFV Value\tcsRCFV (A)\tcsRCFV (V)\tcsRCFV (L)\tcsRCFV (I)\tcsRCFV (P)\tcsRCFV (M)\tcsRCFV (F)\tcsRCFV (W)\tcsRCFV (G)\tcsRCFV (S)\tcsRCFV (T)\tcsRCFV (C)\tcsRCFV (N)\tcsRCFV (Q)\tcsRCFV (Y)\tcsRCFV (D)\tcsRCFV (E)\tcsRCFV (K)\tcsRCFV (R)\tcsRCFV (H)\tcsRCFV (Hydrophobic)\tcsRCFV (Hydrophilic)\tcsRCFV (Polar)\tcsRCFV (Nonpolar)\tcsRCFV (Positive)\tcsRCFV (Neutral)\tcsRCFV (Negative)\tChi-square value\tDegrees of freedom\tp-value\tp(sites)\tp(Gap) [-]\tp(Ambiguity) [X]\tp(Missing Data) [?]\tcsRCFV (Gap)\tcsRCFV (Ambiguity)\tp_uninformative sites\tp_inv(state)\tp_inv(structure)\tp_inv(polarity)\tp_inv(charge)\n";
		
		for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
			
			############################################################################################################ begin summarized all for aa
			## Print OUT summarized results for aminoacid data
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){
				
				print "\n\t", $href_file_of_number->{$file_number}, ", ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}} ;
				print OUTall $href_file_of_number->{$file_number}, "\t", $clade_number, " ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}, "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{Ntaxa}{informative} ;
				
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'RCFV' ){ 
					
					print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state}
				}
				
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative' ){ 
					
					print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{absol_dev}{$state}
				}
				
				for my $state ( 'ChiSquare', 'DF', 'chisprob' ){
					
					print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{$state}
				}
				
				for my $state ( 'alltax', 'GAP', 'AMB', 'MIS' ){ 
					
					print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state}
				}
				
				for my $state ( 'GAP', 'AMB' ){ 
					
					print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{absol_dev}{$state}
				}
				
				for my $state ( 'uninformative' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} }
				
				for my $state ( 'state', 'structure', 'polarity', 'charge' ){ print OUTall "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{pinv}{$state} }
				print OUTall "\n";
				############################################################################################################ end summarized all for aa
				
				
				
				############################################################################################################ begin homogeneity among taxa for aa
				## Print OUT Results of chis-quare tests of homogeneity among taxa
				my $outfile_ht = $href_file_of_number->{$file_number}."_clade".$clade_number."_".$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}.".txt" ;
				open OUTht,  ">BaCoCa_Results/chisquare_test_homogeneity_taxa/$outfile_ht" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/chisquare_test_homogeneity_taxa/$outfile_ht !\n" ;
				
				print OUTht "File: ", $href_file_of_number->{$file_number}, "\tClade: ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}, "\nBase frequencies of taxa:\nTaxon";
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht "\tp(", $state, ")"}
				print OUTht "\tN sites\n";
				
				
				my $aref_taxa = exists($href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}}) ? \@{$href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{$clade_number}{$href_file_of_number->{$file_number}}}} :( ) ;
				for my $taxon ( @$aref_taxa ){
					
					print OUTht $taxon;
					for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht	"\t",	$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{MeanSites}{$taxon}{$state} }
					print OUTht	"\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{allsites}, "\n"
				}
				
				print OUTht "Mean";
				
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{$state} } 
				print OUTht "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{MeanSites}{alltax}, "\n\n",
							"Chi-Square test of homogeneity of base frequencies across taxa:\n\nTaxon";
				
				for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht "\tp(", $state, ")" }
				print OUTht "\n" ;
				
				for my $taxon ( @$aref_taxa ){
					
					print OUTht $taxon, "Observed";
					
					for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nsites}{$taxon}{$state} }
					print OUTht "\n" , $taxon, "Expected";
					for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P' ){ print OUTht "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{taxon}{Nexp}{$taxon}{$state} }
					print OUTht "\n" ;
				}
				
				print OUTht	"\nChi-Square = ",	$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{ChiSquare}, 
								"(df = ", 				$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{DF},
								"), P = ",				$href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{chisprob} ;
								if ($href_freq_of_file_char->{$href_file_of_number->{$file_number}}{$clade_number}{clade}{ht}{chisprob}<0.01) {print OUTht "\nNOTE: Given a significance level of 0.01 the composition is significantly diviating from homogeneity!\n"} else {print OUTht "\nNOTE: Given a significance level of 0.01 the composition is NOT significantly diviating from homogeneity!\n"}
				print OUTht	"\nWarning: This test ignores correlation due to phylogenetic structure.\n";
				############################################################################################################ end homogeneity among taxa for aa
			}
			
			
			
			############################################################################################################ begin taxon frequencies of partitions for aa
			if ( $file_number == 1 ){
				
				my $outfile_tf_prt = $href_file_of_number->{$file_number}."_summarized_taxon_basefrequencies.txt" ;
				open	OUTprt,	">BaCoCa_Results/taxon_basefrequencies_single_partions/$outfile_tf_prt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/taxon_basefrequencies_single_partions/$outfile_tf_prt !\n" ;
				print	OUTprt	"Basefrequencies\nFile\tTaxon\tp(A)\tp(V)\tp(L)\tp(I)\tp(P)\tp(M)\tp(F)\tp(W)\tp(G)\tp(S)\tp(T)\tp(C)\tp(N)\tp(Q)\tp(Y)\tp(D)\tp(E)\tp(K)\tp(R)\tp(H)\tp(Hydrophobic)\tp(Hydrophilic)\tp(Polar)\tp(Nonpolar)\tp(Positive)\tp(Neutral)\tp(Negative)\tp(gap)\tp(Ambiguity) [X]\tRCFV Value (taxon specific)\tAbsolute Deviation (A)\tAbsolute Deviation (V)\tAbsolute Deviation (L)\tAbsolute Deviation (I)\tAbsolute Deviation (P)\tAbsolute Deviation (M)\tAbsolute Deviation (F)\tAbsolute Deviation (W)\tAbsolute Deviation (G)\tAbsolute Deviation (S)\tAbsolute Deviation (T)\tAbsolute Deviation (C)\tAbsolute Deviation (N)\tAbsolute Deviation (Q)\tAbsolute Deviation (Y)\tAbsolute Deviation (D)\tAbsolute Deviation (E)\tAbsolute Deviation (K)\tAbsolute Deviation (R)\tAbsolute Deviation (H)\tAbsolute Deviation (Hydrophobic)\tAbsolute Deviation (Hydrophilic)\tAbsolute Deviation (Polar)\tAbsolute Deviation (Nonpolar)\tAbsolute Deviation (Positive)\tAbsolute Deviation (Neutral)\tAbsolute Deviation (Negative)\tAbsolute Deviation (Gap)\tAbsolute Deviation (Ambiguity) [X]\tp_uninformative sites\n" ;
				
				my $aref_taxa = exists($href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{0}{$href_file_of_number->{$file_number}}}) ? \@{$href_taxa_of_file_clade->{$href_file_of_number->{$file_number}}{$href_clade_of_number->{0}{$href_file_of_number->{$file_number}}}} :( ) ;
				for my $taxon ( @$aref_taxa ){
					
					print OUTprt	$href_file_of_number->{$file_number}, "\t", $taxon ;
					for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'GAP', 'AMB', 'RCFV' ){
						
						print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state}
					}
					
					for my $state ( 'A', 'C', 'G', 'T', 'N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'I', 'E', 'L', 'Q', 'F', 'P', 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'GAP', 'AMB' ){
						
						print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{absol_dev}{$taxon}{$state}
					}
					
					for my $state ( 'uninformative' ){ print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state} }
					
					print OUTprt "\n"
				}
				close OUTprt
			}
			############################################################################################################ end taxon frequencies of partitions for aa
		}
		
		
		
		############################################################################################################ begin taxon frequencies (specific) of all partitions for aa
		for my $state ( 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'GAP', 'AMB', 'RCFV', 'MIS', 'uninformative' ){
			
			my $outfile = "BaCoCa_Results/taxon_basefrequencies_all_partions/".$state."_frequencies_all_partitions.txt";
			open	OUTprt,	">$outfile" or die "\n\t!OUTPUT-ERROR!: Cannot open file $outfile !\n" ;
			
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ print OUTprt "\t", $href_file_of_number->{$file_number} } print OUTprt "\n" ;
			
			for my $taxon ( @$aref_all_taxa ){ print OUTprt $taxon ;
				
				for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
					
					print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{MeanSites}{$taxon}{$state} ;
				}
				print OUTprt "\n"
			}
			close OUTprt ;
		}
		
		for my $state ( 'hydrophobic', 'hydrophilic', 'polar', 'nonpolar', 'positive', 'neutral', 'negative', 'GAP', 'AMB' ){
			
			my $outfile = "BaCoCa_Results/compositional_bias/".$state."_absolute_deviation_frequencies_all_partitions.txt";
			open	OUTprt,	">$outfile" or die "\n\t!OUTPUT-ERROR!: Cannot open file $outfile !\n" ;
			
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){ print OUTprt "\t", $href_file_of_number->{$file_number} } print OUTprt "\n" ;
			
			for my $taxon ( @$aref_all_taxa ){ print OUTprt $taxon ;
				
				for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
					
					print OUTprt "\t", $href_freq_of_file_char->{$href_file_of_number->{$file_number}}{0}{taxon}{absol_dev}{$taxon}{$state} ;
				}
				print OUTprt "\n"
			}
			close OUTprt ;
		}
		############################################################################################################ end taxon frequencies (specific) of all partitions for aa
		
		
		
		############################################################################################################ begin invariant alignment partitions for aa
		for my $state ( 'state', 'structure', 'polarity', 'charge' ){
			
			my $outfile = $state."_invariant_msa_positions";
			open	OUTfas,	">BaCoCa_Results/invariant_alignment_positions/fasta/${outfile}.fas" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/invariant_alignment_positions/${outfile}.fas !\n" ;
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){
				
				print OUTfas ">", $href_clade_of_number->{$clade_number}{$href_file_of_number->{1}}, "\n";
				
				for (my $pos = 1; $pos <= $href_freq_of_file_char->{$href_file_of_number->{1}}{0}{length}; $pos++){
					
					if ( $href_freq_of_file_char->{$href_file_of_number->{1}}{$clade_number}{pinv}{$pos}{$state} == 1 ){ print OUTfas 'T' }else { print OUTfas 'A' }
				} print OUTfas "\n"
			}
			close OUTfas;
			
			open	OUTinv,	">BaCoCa_Results/invariant_alignment_positions/txt/${outfile}.txt" or die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/invariant_alignment_positions/${outfile}.txt !\n" ;
			print	OUTinv	"Invariant alignment positions of ", $state, "\n\n" ;
			
			
			for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){ print OUTinv "\tClade ", $clade_number, ": ", $href_clade_of_number->{$clade_number}{$href_file_of_number->{1}}} print OUTinv "\n";
			
			for (my $pos = 1; $pos <= $href_freq_of_file_char->{$href_file_of_number->{1}}{0}{length}; $pos++){ print OUTinv $pos;
				
				for my $clade_number (sort {$a<=>$b} keys %$href_clade_of_number ){
					
					print OUTinv "\t", $href_freq_of_file_char->{$href_file_of_number->{1}}{$clade_number}{pinv}{$pos}{$state}
				} print OUTinv "\n" ;
			}
			close OUTinv;
			
			################################################
			## svg output
			for my $file_number ( sort {$a<=>$b} keys %$href_file_of_number ){
				
				my $svg_filename = "invariant_alignment_positions/svg/".$state."_". $href_file_of_number->{$file_number}."_invariant_msa_positions.svg" ;
				
				my @position_numbers;
				for my $pos ( $href_seq_ranges_of_part->{$href_file_of_number->{$file_number}}{start} .. $href_seq_ranges_of_part->{$href_file_of_number->{$file_number}}{end} ){ $pos += 1; push @position_numbers, $pos }
				
				&svg_plot ( \%$href_clade_of_number, \%$href_freq_of_file_char, \%$href_file_of_number, \$state, \$svg_filename, \$file_number, \@position_numbers )
			}
		}
		############################################################################################################ end invariant alignment partitions for aa
	}
	
	open OUTinfo, ">BaCoCa_Results/explanation_guide.txt" || die "\n\t!OUTPUT-ERROR!: Cannot open file BaCoCa_Results/explanation_guide.txt !\n" ;
	
	print OUTinfo	"summarized_frequencies.txt\n",
					"In the summary_frequency.txt file the columns provide the following data in this order:\n",
					"The first column is the partition or complete data set name followed by...\n",
					"...the name of the user defined clade,\n",
					"...the number of informative taxa in that partition and clade,\n",
					"...the frequencies of individual nucleotides or amino acids (e.g. p(A)),\n",
					"...the frequencies of classes of nucleotides or amino acids (e.g. p(AT) or p(polar)),\n",
					"...the character-state specific RCFV (csRCFV) values for individual nucleotides or amino acids or for corresponding classes in the same order as the frequencies  (see manual and manuscript),\n",
					"...for nucleotides only skew-values as described in the manual and manuscript,\n",
					"...results of the chi-square tests with the actual chi-square value and the degree of freedom as well as the corresponding p-value,\n",
					"      NOTE: For example, for a significance level of 0.01 each p-value below 0.01 means that the composition is significantly diviating from homogeneity!\n",
					"...the RCFV value for the complete data set or partition  (see manual and manuscript),\n",
					"...for nucleotides only the c-value of saturation including the corresponding mean p and ti/tv ratio as well as their standard deviations (see manual and manuscript),\n",
					"...the frequencies of defined sites for all informative taxa (p(sites)) as well as of GAPs, ambiguities, and missing data for all taxa,\n",
					"...the character-state specific RCFV (csRCFV) values for defined sites for all informative taxa as well as of GAPs, ambiguities, and missing data for all taxa in the same order as the frequencies  (see manual and manuscript),\n",
					"...the frequencies of all uninformative characters and of invariant sites for states and classes (see manual).\n\n",
					"taxon_basefrequencies_single_partions/*.summarized_taxon_base_frequencies.txt\n",
					"In this file the columns provide the following data in this order:\n",
					"The first column is the complete data set name followed by...\n",
					"...the name of the taxon,\n",
					"...the frequencies of individual nucleotides or amino acids (e.g. p(A)), of classes of nucleotides or amino acids (e.g. p(AT) or p(polar)) and of GAPs and ambiguities\n",
					"...the RCFV value for the complete data set or partition  (see manual and manuscript),\n",
					"...the absolute deviation from the mean for individual nucleotides or amino acids, for corresponding classes or for GAPs and ambiguities in the same order as the frequencies,\n",
					"...the frequencies of all uninformative characters.\n\n",
					"In the folders taxon_basefrequencies_all_partions, skew_values, and compositional bias taxon versus partition matrices are deposited for different parameters such as AT content, uninformative sites, or RCFV values. For a detailed discussion refer to the manual and the manuscript\n\n",
					"In the folders missing_data_overlap and c_value_calculations taxon versus taxon matrices are deposited for different partitions and user-defined taxon subsets for negative and positive shared data overlap or uncorrected pairwise evolutionary distances and ti/tv ratio, respectively. For a detailed discussion refer to the manual and the manuscript\n\n",
					"In the chisquare_test_homogeneity_taxa folder, the chi-square analysis and the results are given as plain text files in the same format as used in PAUP.\n";
	close OUTinfo ;
	
	
	
	sub svg_plot{
		
		my $href_clade_number	= $_[0] ;
		my $href_freq_file		= $_[1] ;
		my $href_number_file	= $_[2] ;
		my $sref_state			= $_[3] ;
		my $sref_filename		= $_[4] ;
		my $sref_file_number	= $_[5] ;
		my $aref_position_numb	= $_[6] ;
		
		my %hoh_color_of_state = (
									'state'		=> 'darkred',
									'class'		=> 'royalblue',
									'structure'	=> 'orange' ,
									'polarity'	=> 'green' ,
									'charge'	=> 'violet' ,
								);
		
		my $ncolumns = $href_freq_file->{$href_file_of_number->{$$sref_file_number}}{0}{length}	; $ncolumns--	; $ncolumns	*= 10  ;
		my $nrows    = keys	%$href_clade_number													; $nrows--		; $nrows		*= 10  ;
		
		my $init_line = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>' ;
		my $gen_line  = '<!-- created by BaCoCa.v0.41beta.pl -->'                ;
		
		open my $fh_matrix_out , ">BaCoCa_Results/$$sref_filename" ; 
		
		my $width = $ncolumns + 11 ;
		my $height= $nrows    + 11 ;
		my @clade = () ;
		 
		
		
		
		
		print $fh_matrix_out <<FRAME;
$init_line
$gen_line

<svg transform="translate(10,70)"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   version="1.0"
   width="$width"
   height="$height"
   id="svg2">

  <defs
     id="defs4" />

<rect
     width="$width"
     height="$height"
     x="0"
     y="0"
     style="opacity:1;fill:white;fill-opacity:0;stroke:black;stroke-width:1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect001" />


FRAME
		
		
		my $y = 1 ;
		
		for my $cladenr ( sort {$a<=>$b} keys %$href_clade_number ) {
			
			my @characs;
			
			push @clade , $href_clade_number->{$cladenr}{$href_number_file->{$$sref_file_number}} ;
			for (my $pos = 0; $pos < $ncolumns; $pos++){ push @characs, $href_freq_file->{$href_number_file->{$$sref_file_number}}{$cladenr}{pinv}{$pos}{$$sref_state} }
			my $x = 1 ;
			
			for my $char ( @characs ) {
				
				my $id = int rand 100000 ;
				
				if ( $char =~ /^1$/ ) {
				
print $fh_matrix_out <<RECT;

<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:$hoh_color_of_state{$$sref_state};fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
				$x +=10
				}
				
				elsif ( $char =~ /^0$/ ) {
				
print $fh_matrix_out <<RECT;

<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:white;fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
				$x +=10
				}
				
			}
			$y += 10
		}
		
		
		$y = 10 ;
		
		#write clades to matrix
		
print $fh_matrix_out <<TAXA;	
 <text
     x="-5"
     y="10"
     style="font-size:8px;font-style:normal;font-weight:normal;text-align:end;text-anchor:end;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
    id="text1892"
    
TAXA
		
		print $fh_matrix_out    "xml:space=\"preserve\">" ;
		
		
		for my $clade ( @clade ) {
			
			$clade =~ s/_/ /g ;
			
			my $id = int rand 100000 ;
			
			print $fh_matrix_out "<tspan\n       x=\"-5\"\n       y=\"$y\"\n       style=\"text-align:end;text-anchor:end\"\n       id=\"tspan$id\">$clade</tspan>" ;
			
			$y += 10 ;
		}
		
		# legende
		$y += 5;
		my $id = int rand 100000 ;
		print $fh_matrix_out "<tspan\n       x=\"105\"\n       y=\"$y\"\n       style=\"text-align:end;text-anchor:end\"\n       id=\"tspan$id\">Invariant Positions ($$sref_state)</tspan>" ;
		$id = int rand 100000 ;
		
		$y -= 2.5;
		
print $fh_matrix_out <<TAXA;	
</text>

TAXA

my $x = 10     ;
for my $bin ( @$aref_position_numb ) {  
	
	my $id = int rand 100000 ;
	
print $fh_matrix_out <<BINS;
<text
     x="$x"
     y="-5"
     transform="rotate(-90,$x,-5)"     	
     style="font-size:8px;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
     id="text8747"
     xml:space="preserve"><tspan
       x="$x"
       y="-5"
       id="tspan$id">$bin</tspan></text>
BINS
	
	$x += 10 ;
	
}




		
print $fh_matrix_out <<CIRCLE;

<circle
     r="4"
     cx="112"
     cy="$y"
     style="fill:$hoh_color_of_state{$$sref_state};fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
CIRCLE
		
		
print $fh_matrix_out <<FINISH;

</svg>

	
FINISH
	}
	
}


sub end{
		
		TIMER:
		# set timer
		my ( $user, $system, $cuser, $csystem ) = times ;
		
		print "\n\n\t-------------------------------------------\n\tBye Bye... BaCoCa\n\t-------------------------------------------\n\t";
		
		print <<TIME;
		
		***  time used: $user sec  ***
		
TIME
		

		
		
		
}

































































