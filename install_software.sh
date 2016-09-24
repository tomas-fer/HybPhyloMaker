#INSTALL SOFTWARE NECESSARY FOR HybPhyloMaker
#and clone HybPhyloMaker GitHub repository (incl. test dataset)
#This should work on major Linux distribution (tested on Debian), but carefully set appropriate installer and names of some libraries!!!
#Be sure that you have installed gcc, gcc-c++, make before running this script

#Change name of your default package management tool (apt-get on Debian/Ubuntu, zypper on OpenSUSE, yum on Fedora/CentOS/RHEL/Scientific)
installer=apt-get

#Install software using
$installer install -y python
$installer install -y python3 #Does not work on CentOS???
$installer install -y perl
#$installer install -y parallel #better to install from source, see below
#$installer install -y bowtie2 #better to install from source, see below
#$installer install -y samtools #better to install from source, see below
#$installer install -y fastx-toolkit #better to install from source, see below
$installer install -y openjdk-7-jre #java-1.7.0-openjdk.x86_64 in Fedora; java-1_7_0-openjdk in OpenSUSE
#$installer install -y mafft #better to install from source, see below
#$installer install -y fasttree #better to install from source, see below
$installer install -y r-base #R in Fedora
$installer install -y git
$installer install -y libpng-dev #libpng-devel on Fedora/CentOS/OpenSUSE
$installer install -y zlib1g-dev #zlib-devel on Fedora/CentOS/OpenSUSE
$installer install -y wget
$installer install -y tar

#Install R packages
R -q -e "install.packages('ape', repos='http://cran.rstudio.com/')"
R -q -e "install.packages('seqinr', repos='http://cran.rstudio.com/')"
R -q -e "install.packages('data.table', repos='http://cran.rstudio.com/')"

#Install other software
mkdir install
cd install

#MAFFT
wget http://mafft.cbrc.jp/alignment/software/mafft-7.305-without-extensions-src.tgz
tar -xvf mafft-7.305-without-extensions-src.tgz
rm mafft-7.305-without-extensions-src.tgz
cd mafft-7.305-without-extensions/core
make
make install
cd ../..

#GNU parallel
wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar xjvf parallel-latest.tar.bz2
rm parallel-latest.tar.bz2
cd parallel*
./configure
make
make install
cd ..

#Samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xjvf samtools-1.3.1.tar.bz2
rm samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --without-curses
make
make install
cd ..

#FastTree
mkdir FastTree
cd FastTree
wget http://www.microbesonline.org/fasttree/FastTree
chmod +x FastTree
cp FastTree /usr/local/bin
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm
chmod +x FastTreeMP
cp FastTreeMP /usr/local/bin
cd ..

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
echo $MACHTYPE #if you get x86_64-pc-linux-gnu continue with next line, if something else write appropriate short name
MACHTYPE=x86_64
export MACHTYPE
mkdir lib/$MACHTYPE
mkdir -p ~/bin/$MACHTYPE
make
cp ~/bin/$MACHTYPE/blat /usr/local/bin
cd ..

#Fastx-toolkit
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
tar -xvf libgtextutils-0.7.tar.gz
rm libgtextutils-0.7.tar.gz
cd libgtextutils-0.7
./configure
make
make install
cd ..
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
tar xjvf fastx_toolkit-0.0.14.tar.bz2
rm fastx_toolkit-0.0.14.tar.bz2
cd fastx_toolkit-0.0.14
./configure
make
make install
cd ..

#Bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-source.zip
unzip bowtie2-2.2.9-source.zip
rm bowtie2-2.2.9-source.zip
cd bowtie2-2.2.9
make
cp bowtie2* /usr/local/bin
cd ..

#p4
#see http://p4.nhm.ac.uk/installation.html
$installer install -y python-numpy
$installer install -y python-scipy
$installer install -y libgsl0-dev #gsl-devel in Fedora and OpenSUSE
$installer install -y python-dev #python-devel in OpenSUSE
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
for i in parallel bowtie2 samtools bam2fastq java fastx_collapser perl blat mafft python python3 trimal mstatx FastTree nw_reroot nw_topology raxmlHPC raxmlHPC-PTHREADS R p4; do
	command -v $i >/dev/null 2>&1 || { echo -n $i; echo >&2 "...not found"; }
done
echo "Necessary software installed, possible errors indicated above."
echo -e "If you don't see any *not found* your system is now ready to run HybPhyloMaker!\n"

#Clone HybPhyloMaker GitHub repository
git clone https://github.com/tomas-fer/HybPhyloMaker
cd HybPhyloMaker
chmod +x *.sh
chmod +x HybSeqSource/ASTRID

