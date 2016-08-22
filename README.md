# HybPhyloMaker
Set of bash scripts for analysis of HybSeq raw data. Consists of several steps:  
0. Download FASTQ files from Illumina BaseSpace storage  
1. Processing raw reads (PhiX removal, adaptor removal, quality filtering, summary statistics)  
_Intermediate manual step - read mapping to 'pseudoreference' in Geneious, export consensus file_  
2. Recognize sequences matching probes (generate PSLX files using BLAT)  
3. Create alignments for all genes  
4. Treat missing data, select best genes  
5. Generate FastTree or RAxML gene trees + trees-alignment properties  
6. Root gene trees with outgroup, combine gene trees into a single file  
7. Estimate species tree (ASTRAL, ASTRID, MRL, concatenation) 
8. Subselect suitable genes and repeat steps 6+7 
  
Uses many additional software that must be installed and put in the PATH prior tu run scripts.  
Also utilizes many scripts developed by others (located in HybSeqSource folder). PLEASE CITE APPROPRIATELY THOSE SCRIPTS WHEN USING HybPhyloMaker!  
Read manual located in docs folder before running HybPhyloMaker.  

