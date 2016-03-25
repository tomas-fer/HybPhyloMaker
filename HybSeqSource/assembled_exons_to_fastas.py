#!/usr/bin/env python
from os import path
from os import makedirs
from re import sub   #This imports regular expression usage
from optparse import OptionParser    #Imports the option parser module

###### OPTIONS and USAGE ######
parser = OptionParser(usage = """assembled_exons_to_fastas.py -l PLSX_list -f FASTA_file -d OUT_directory [-g] [-r] [-n reference_NAME]

assembled_exons_to_fastas.py -- Processes the BLAT output from several samples 
    comparing reference-guided contig assemblies with the reference (probe) 
    contigs. The output is a directory containing a fasta file for each
    reference contig holding the sequence for each sample.

Copyright (c) 2014 Kevin Weitemier
Version 1.0

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details. You
should have received a copy of the GNU General Public License along with this
program.  If not, see <http://www.gnu.org/licenses/>.

If you find this program useful, please cite:

K. Weitemier, S.C.K. Straub, R. Cronn, M. Fishbein, R. Schmickl, A. McDonnell,
and A. Liston. 2014. Hyb-Seq: Combining target enrichment and genome skimming
for plant phylogenomics. Applications in Plant Sciences 2(9): 1400042.

Input - A file containing a list of every .pslx file to be processed (one file
    per line), and a fasta file of the targeted contigs or exons. The .pslx
    files to be processed need to have been created with the targeted exons
    input as the 'database' or 'targets.'""")
parser.add_option("-l", action="store", type="string", dest="conname",
    help="File containing a list of the .pslx files to process; one name per line", default="")
parser.add_option("-f", action="store", type="string", dest="faname",
    help="Fasta file of targeted contigs", default="")
parser.add_option("-d", action="store", type="string", dest="dirname",
    help="Name of directory for output files", default="Exons_to_be_aligned")
parser.add_option("-g", action="store_true", default=False, dest="filler",
    help="""Insert gaps? If the sample sequence matches the reference across two or more
blocks, a string of Ns will be inserted equal to the size of the gap(s). In
certain cases this may match the size of introns. Without this separate blocks will be concatenated. Default=False""")
parser.add_option("-r", action="store_true", default=False, dest="ref_out",
    help="""Output the reference sequence as an entry in each sequence alignment? default=False""")
parser.add_option("-n", action="store", type="string", dest="ref_name", default="Reference",
    help="""If the reference sequence is output (-r), what name will be given to this sample
in the alignment files? Default=Reference""")
(options, args) = parser.parse_args()

# Makes sure all filenames are given
if options.conname == "":
    parser.error("Please include a .pslx list file using -l.")
if options.faname == "":
    parser.error("Please include a fasta file using -f.")
if path.exists(options.dirname):
    parser.error("The output directory already exists. Please remove it or provide a new name using the -d option.")

###### OPENING INPUT/OUTPUT FILES ######
ConfigFile = open(options.conname, 'r')
FastaFile = open(options.faname, 'r')
makedirs(options.dirname)

# Opening the contigs file and makeing a dictionary of each sequence
ReferenceContigs = {}
FaLine = FastaFile.readline()
while FaLine:
    FaLine = FaLine.strip()
    if not FaLine.startswith('>'):
        sys.exit("""The fasta file is not formatted correctly. Please be sure that each entry is given on exactly two lines: the ID line and the sequence line.""")
    ID = FaLine.lstrip('>')
    FaLine = FastaFile.readline().strip()
    Seq = FaLine
    ReferenceContigs[ID] = Seq
    FaLine = FastaFile.readline()
FastaFile.close()

Contigs = {}
for exon in ReferenceContigs:
    Contigs[exon] = {}

# Opening the pslx list and grabbing the entries
pslxFiles = []
ConLine = ConfigFile.readline()
while ConLine:
    ConLine = ConLine.strip()
    pslxFiles.append(ConLine)
    ConLine = ConfigFile.readline()
ConfigFile.close()

Names = []
# Processing each pslx file
for Filename in pslxFiles:
    Name = sub('Final_Assembly_',r'',Filename)
    Name = sub('.pslx',r'',Name)
    Names.append(Name)
    File = open(Filename, 'r')
    pslxLine = File.readline()
    while pslxLine:
        pslxLine = pslxLine.strip()
        if not ',' in pslxLine:
            pslxLine = File.readline()
            continue
        Fields = pslxLine.split("\t")
        Length = int(Fields[0]) + int(Fields[1])
        ThisExon = Fields[13]
        MySeq = ''
        if int(Fields[17]) == 1 or not options.filler:
            MySeq = Fields[21].replace(",", "")
        else:
            blockSizes = Fields[18].split(",")
            blockStarts = Fields[19].split(",")
            Blocks = Fields[21].split(',')
            MySeq = Blocks[0]
            for i in range(int(Fields[17])-1):
                GapSize = (int(blockStarts[i+1]) - (int(blockStarts[i]) + int(blockSizes[i])))
                GapFiller = 'N' * GapSize
                MySeq = MySeq + GapFiller + Blocks[i+1]
        if not ThisExon in Contigs:
            sys.exit("The contigs in the fasta file don't match those in the .pslx files. Be sure the names of the contigs don't contain spaces.")
        if Name in Contigs[ThisExon]:
            if Length > Contigs[ThisExon][Name][0]:
                Contigs[ThisExon][Name] = [Length, MySeq]
        else:
            Contigs[ThisExon][Name] = [Length, MySeq]
        pslxLine = File.readline()
    File.close()

for exon in Contigs:
    OutExon = sub(',', r'...', exon)
    OutName = options.dirname + '/' + "To_align_" + OutExon + ".fasta"
    OutFile = open(OutName, 'w')
    Filler = 'n' * len(ReferenceContigs[exon])
    for Name in Names:
        if Name in Contigs[exon]:
            OutFile.write(">%s\n%s\n" % (Name, Contigs[exon][Name][1]))
        else:
            OutFile.write(">%s\n%s\n" % (Name, Filler))
    if options.ref_out:
        OutFile.write(">%s\n%s\n" % (options.ref_name, ReferenceContigs[exon]))
    OutFile.close()
#EOF
