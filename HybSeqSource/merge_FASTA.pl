#!/usr/bin/env perl

############     merge_FASTA.pl     ############
#
#	This Perl script concatenates several FASTA sequence
#	files to one large FASTA sequence file.
#
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de



############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed and
#	the name of the output file as arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	The output file name has to be passed to the script via -o:
#	-o output_file.txt
#
#	If no output file name is passed, the output is written into:
#	merged.fas
#
#	For example you can type the following command:
#	perl merge_FASTA.pl -i file1.fas -i file2.fas -o output_merged.fas
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	Multiple files and combinations of all kinds of arguments are
#	allowed:
#	perl merge_FASTA.pl -i file1.fas -I list_of_files.txt -o output_merged.fas




@input_files=();
$output_file_name="merged.fas";
$|=1;
###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl merge_fasta.pl -argument1 -argument2 ...\n\n";
	exit;
	}

$argv="";
foreach(@ARGV)
	{
	$argv.=$_;
	}
@arguments=split('-',$argv);

foreach(@arguments)
	{
	if($_=~/^ *i/)
		{
		$_=~s/^ *i//;
		$_=~s/ //g;
		push(@input_files,$_);
		}
	elsif($_=~/^ *I/)
		{
		$_=~s/^ *I//;
		$_=~s/ //g;
		open(FILES_IN,"$_");
		while(<FILES_IN>)
			{
			unless($_=~/^\s*$/)
				{
				$_=~s/\s//sg;
				push(@input_files,$_);
				}
			}
		}
	elsif($_=~/^ *o/)
		{
		$_=~s/^ *o//;
		$_=~s/ //g;
		$output_file_name=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
	exit;
	}

###   PRINT ARGUMENTS   ###
print"The following files will be processed:\n";
foreach(@input_files)
	{
	if(-e $_)
		{
		print"$_\n";
		push(@input_files_ok,$_);
		}
	else
		{
		print"could not find file: $_. It will be ignored.\n";
		}
	}
print"Output file name: $output_file_name\n\n";

###   START   ###
$total_seqs=0;
open(OUT,">$output_file_name");
foreach$file(@input_files_ok)
	{
	$seqs_in_file=0;
	print"processing $file";
	open(IN,$file);
	while(<IN>)
		{
		$line=$_;
		if($_=/^>/)
			{
			$seqs_in_file++;
			print OUT $line;
			}
		elsif($_=~/^\s*$/)
			{
			print OUT $line;
			}
		}
	print"\tSequences: $seqs_in_file\n";
	$total_seqs+=$seqs_in_file;
	}
print"Total sequences: $total_seqs\n\n";
exit;