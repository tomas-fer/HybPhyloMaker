# HybPipe
Set of bash scripts for analysis of HybSeq raw data. Consists of several steps:
0. Download FASTQ files from Illumina BaseSpace storage
1. Processing raw reads (PhiX removal, adaptor removal, quality filtering, summary statistics)
2. Read mapping to 'pseudoreference' in Geneious, export consensus file
3. Generate PSLX files using BLAT
