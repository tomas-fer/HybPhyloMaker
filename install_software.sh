#INSTALL SOFTWARE NECESSARY FOR HybPhyloMaker
#and clone HybPhyloMaker GitHub repository (incl. test dataset)
#this should work on major Linux distribution (tested on Debian)

#Install software using apt-get
apt-get install -y python
apt-get install -y python3
apt-get install -y perl
apt-get install -y parallel
apt-get install -y bowtie2
apt-get install -y samtools
apt-get install -y fastx-toolkit
apt-get install -y openjdk-7-jre
apt-get install -y mafft
apt-get install -y fasttree
apt-get install -y r-base
apt-get install -y git

#Install R packages
R -q -e "install.packages('ape', repos='http://cran.rstudio.com/')"
R -q -e "install.packages('seqinr', repos='http://cran.rstudio.com/')"
R -q -e "install.packages('data.table', repos='http://cran.rstudio.com/')"

#Install other software
mkdir install
cd install

#RAxML
git clone https://github.com/stamatak/standard-RAxML
cd standard-RAxML
make -f Makefile.gcc
rm *.o
make -f Makefile.SSE3.gcc
rm *.o
make -f Makefile.AVX.gcc
rm *.o
make -f Makefile.PTHREADS.gcc
rm *.o
make -f Makefile.SSE3.PTHREADS.gcc
rm *.o
make -f Makefile.AVX.PTHREADS.gcc
rm *.o
cp raxmlHPC* /usr/local/bin
cd ..

#MstatX
git clone https://github.com/gcollet/MstatX
cd MstatX
make
cp mstatx /usr/local/bin
cd ..

#TrimAl
git clone https://github.com/scapella/trimal
cd trimal/source
make
cp trimal /usr/local/bin
cd ../..

#Newick Utilities
wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
tar xfz newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
rm newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
cp newick-utils-1.6/src/nw_* /usr/local/bin

#bam2fastq
wget https://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz
tar xfz bam2fastq-1.1.0.tgz
rm bam2fastq-1.1.0.tgz
cd bam2fastq-1.1.0
make
cp bam2fastq /usr/local/bin
cd ..

#BLAT
#see, e.g., http://nix-bio.blogspot.cz/2013/10/installing-blat-and-blast.html for installing tips
wget https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
unzip blatSrc35.zip
rm blatSrc35.zip
cd blatSrc
echo MACHTYPE #if you get x86_64-pc-linux-gnu continue with next line, if something else write appropriate short name
MACHTYPE=x86_64
export MACHTYPE
mkdir lib/$MACHTYPE
mkdir -p ~/bin/$MACHTYPE
make
cp ~/bin/$MACHTYPE/blat /usr/local/bin
cd ..

#p4
#see http://p4.nhm.ac.uk/installation.html
apt-get install -y python-numpy
apt-get install -y python-scipy
apt-get install -y libgsl0-dev
apt-get install -y python-dev
git clone https://github.com/pgfoster/p4-phylogenetics
cd p4-phylogenetics
python setup.py build
python setup.py install
python setup.py build_ext -i
cd ..

#Leave 'install' directory
cd ..

#Check if everything is installed correctly
echo -e "\nSoftware installed...checking for binaries in PATH"
for i in parallel bowtie2 samtools bam2fastq java fastx_collapser perl blat mafft python python3 trimal mstatx fasttree nw_reroot nw_topology raxmlHPC raxmlHPC-PTHREADS R p4; do
	command -v $i >/dev/null 2>&1 || { echo -n $i; echo >&2 "...not found"; }
done
echo "Necessary software installed, possible errors indicated above."
echo -e "If you don't see any *not found* your system is now ready to run HybPhyloMaker!\n"

#Clone HybPhyloMaker GitHub repository
git clone https://github.com/tomas-fer/HybPhyloMaker
cd HybPhyloMaker
chmod +x *.sh
chmod +x HybSeqSource/ASTRID

