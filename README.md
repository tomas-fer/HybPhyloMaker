# HybPhyloMaker
   
> Fér T. & Schmickl R. (2018): [HybPhyloMaker: target enrichment data analysis from raw reads to species trees.](http://journals.sagepub.com/doi/10.1177/1176934317742613) *Evolutionary Bioinformatics* 14: 1-9. doi: 10.1177/1176934317742613  
   
**Set of bash scripts for analysis of HybSeq data (from raw reads to species trees).** Consists of several steps:   
  

0:  Prepare FASTQ files to folder (optionally download files from Illumina BaseSpace storage)  
1:  Processing raw reads (PhiX removal, adaptor removal, quality filtering, summary statistics)  
2:  Mapping reads to reference (using Bowtie2/BWA), create consensus sequence  
3:  Recognize sequences matching probes (generate PSLX files using BLAT)  
4:  Create alignments for all genes (+ optionally correct reading frame and/or select low heterozygosity loci)  
5:  Treat missing data, select best genes  
6:  Generate FastTree or RAxML gene trees + calculate/plot trees-alignment properties  
7:  Root gene trees with outgroup, combine gene trees into a single file  
8:  Estimate species tree (ASTRAL, ASTRID, MRL, BUCKy, concatenation)  
9:  Subselect suitable genes and repeat steps 7+8  
10:Subselect trees based on samples presence, collapse unsupported branches  
  
Uses many additional software that must be installed and put in the PATH prior to run scripts (see [Table](docs/HybPhyloMaker_software.pdf) located in docs folder and consider to run [install_software.sh](install_software.sh)).  
Also utilizes many scripts developed by others (located in [HybSeqSource folder](HybSeqSource)). PLEASE CITE APPROPRIATELY THOSE SCRIPTS WHEN USING HybPhyloMaker!  

Read [manual](docs/HybPhyloMaker_manual1.6.4.pdf) located in docs folder before running HybPhyloMaker.  

# HybPhyloMaker workflow
![HybPhyloMaker workflow](https://github.com/tomas-fer/HybPhyloMaker/blob/master/docs/HybPhyloMaker_workflow.png)

# Papers citing HybPhyloMaker
1. Carlsen MM, Fér T, Schmickl R, Leong-Škorničková J, Newman M, Kress WJ. 2018. [Resolving the rapid plant radiation of early diverging lineages in the tropical Zingiberales: Pushing the limits of genomic data](https://www.sciencedirect.com/science/article/pii/S1055790317309296). Molecular Phylogenetics and Evolution 128:55-68. doi: 10.1016/j.ympev.2018.07.020  
2. Herrando-Moraira S, Cardueae Radiations Group. 2018. [Exploring data processing strategies in NGS target enrichment to disentangle radiations in the tribe Cardueae (Compositae)](https://www.sciencedirect.com/science/article/pii/S1055790318302501). Molecular Phylogenetics and Evolution 128:69-87. doi: 10.1016/j.ympev.2018.07.012  
3. Villaverde T, Pokorny L, Olsson S, Rincón-Barrado M, Johnson MG, Gardner EM, Wickett NJ, Molero J, Riina R, Sanmartín I. 2018. [Bridging the micro- and macroevolutionary levels in phylogenomics: Hyb-Seq solves relationships from populations to species and above](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15312). New Phytologist 220:636–650. doi: 10.1111/nph.15312  
4. Jones KE, Fér T, Schmickl RE, Dikow RB, Funk VA, Herrando-Moraira S, Siniscalchi CM, Susanna A, Slovák M, Thapa R, Watson LE, Mandel JR. 2019. [An empirical assessment of a single family-wide hybrid capture locus set at multiple evolutionary timescales in Asteraceae](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/aps3.11295). Applications in Plant Sciences 7(10):e11295. doi: 10.1002/aps3.11295  
5. Karimi N, Grover CE, Galagher JP, Wendel JF, Ané C, Baum DA. 2020. [Reticulate evolution helps explain apparent homoplasy in floral biology and pollination in baobabs (_Adansonia_; Bombacoideae; Malvaceae)](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syz073/5613901). Systematic Biology 69:462–478. doi: 10.1093/sysbio/syz073  
6. Mao Y, Hou S, Shi J, Economo EP. 2020. [TREEasy: an automated workflow to infer gene trees, species trees, and phylogenetic networks from multilocus data](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13149). Molecular Ecology Resources 20:832-840. doi: 10.1111/1755‐0998.13149  
7. Tomasello S, Karbstein K, Hodač L, Paetzold C, Hörandl E. 2020. [Phylogenomics unravels Quaternary vicariance and allopatric speciation patterns in temperate‐montane plant species: a case study on the _Ranunculus auricomus_ species complex](https://doi.org/10.1111/mec.15458). Molecular Ecology 29:2031-2049. doi: 10.1111/mec.15458  
8. Rejlová L, Böhmová A, Chumová Z, Hořčicová Š, Josefiová J, Schmidt PA, Trávníček P, Urfus T, Vít P, Chrtek J. 2020. [Disparity between morphology and genetics in _Urtica dioica_ (Urticaceae)](https://academic.oup.com/botlinnean/advance-article/doi/10.1093/botlinnean/boaa076/5903206). Botanical Journal of the Linnean Society. doi: 10.1093/botlinnean/boaa076
9. Karbstein K, Tomasello S, Hodač L, Dunkel FG, Daubert M, Hörandl E. 2020. [Phylogenomics supported by geometric morphometrics reveals delimitation of sexual species within the polyploid apomictic _Ranunculus auricomus_ complex (Ranunculaceae)](https://www.biorxiv.org/content/10.1101/2020.01.07.896902v1.full). bioRxiv. doi: 10.1101/2020.01.07.896902  
10. Knyshov A, Gordon ERL, Christiane Weirauch C. 2020. [New alignment-based sequence extraction software (ALiBaSeq) and its utility for deep level phylogenetics](https://www.biorxiv.org/content/10.1101/2020.04.27.064790v1.full). bioRxiv. doi: 10.1101/2020.04.27.064790  
11. Chumová Z,  Záveská E,  Ponert J,  Schmidt PA,  Trávníček P. 2020. [Partial endoreplication stimulates diversification in the species-richest lineage of orchids](https://www.biorxiv.org/content/10.1101/2020.05.12.091074v1). bioRxiv. doi: 10.1101/2020.05.12.091074  
12. Ojeda DI, Koenen E, Cervantes S, de la Estrella M, Banguera-Hinestroza E, Janssens SB, Migliore J, Demenou B, Bruneau A, Forest F, Hardy OJ. 2020 [Phylogenomics within the _Anthonotha_ clade (Detarioideae, Leguminosae) reveals a high diversity in floral trait shifts and a general trend towards organ number reduction](https://www.biorxiv.org/content/10.1101/511949v1.full). bioRxiv. doi: 10.1101/511949
