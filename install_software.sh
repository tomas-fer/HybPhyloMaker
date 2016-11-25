#INSTALL SOFTWARE NECESSARY FOR HybPhyloMaker
#and clone HybPhyloMaker GitHub repository (incl. test dataset)
#This should work on major Linux distribution (tested on Debian, Ubuntu, OpenSUSE, Fedora, and CentOS)
#Without changes works only on 64-bit platforms (x86_64)

#Carefully set your distribution
distribution=Debian #one of: Debian (also for Ubuntu), OpenSUSE, Fedora, CentOS)
#Change name of your default package management tool (apt-get on Debian/Ubuntu, zypper on OpenSUSE, yum on Fedora/CentOS/RHEL/Scientific)
installer=apt-get #one of apt-get, zypper, yum

#------INSTALL SOFTWARE FROM REPOSITORIES------
#Install software using default package manager specified above

#Compilation utilities, i.e., gcc, g++, make
if [[ $distribution =~ "Debian" ]]; then
	for i in gcc g++ make; do
		if ! [ -x "$(command -v $i)" ]; then
			$installer install -y $i
		fi
	done
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	for i in gcc-c++ g++ make; do
		if ! [ -x "$(command -v $i)" ]; then
			$installer install -y $i
		fi
	done
fi

#Perl
if ! [ -x "$(command -v perl)" ]; then
	$installer install -y perl
fi

#Python
if ! [ -x "$(command -v python)" ]; then
	$installer install -y python
fi

#Python3
if ! [ -x "$(command -v python3)" ]; then
	if [[ $distribution =~ "CentOS" ]]; then
		$installer install -y epel-release #Only for CentOS
		$installer install -y python34 #Only for CentOS
	else
		$installer install -y python3 #Does not work on CentOS
	fi
fi

#Java
if [[ $distribution =~ "Debian" ]]; then
	$installer install -y openjdk-7-jre #Debian/Ubuntu
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	$installer install -y java-1.7.0-openjdk.x86_64 #Fedora/CentOS
elif [[ $distribution =~ "OpenSUSE" ]]; then
	$installer install -y java-1_7_0-openjdk #OpenSUSE
fi


#zlib library
if [[ $distribution =~ "Debian" ]]; then
	$installer install -y libpng-dev #Debian/Ubuntu
	$installer install -y zlib1g-dev #Debian/Ubuntu
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	$installer install -y libpng-devel #Fedora/CentOS/OpenSUSE
	$installer install -y zlib-devel #Fedora/CentOS/OpenSUSE
fi

#pkg-config
if [[ $distribution =~ "Debian" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	$installer install -y pkg-config #Debian/Ubuntu/OpenSUSE
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	$installer install -y pkgconfig #Fedora/CentOS
fi

#R
if [[ $distribution =~ "Debian" ]]; then
	$installer install -y r-base-dev #Debian/Ubuntu
	$installer install -y gfortran #Debian/Ubuntu
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	$installer install -y R #Fedora/CentOS
	$installer install -y gcc-gfortran #Fedora/CentOS
elif [[ $distribution =~ "OpenSUSE" ]]; then
	$installer install -y R-base-devel #OpenSUSE
	$installer install -y gcc-fortran #necessary in OpenSUSE for compilation of some R packages
fi

#R packages
for Rpackage in ape seqinr data.table; do
	R -q -e "is.element('$Rpackage', installed.packages()[,1])" > testpackage
	if grep -Fxq "[1] FALSE" testpackage; then
		R -q -e "install.packages('$Rpackage', repos='http://cran.rstudio.com/')"
	else
		echo -e "R package $Rpackage already installed"
	fi
done
rm testpackage

#wget, tar, bzip2, bc, git
for i in wget tar bzip2 bc git; do
	if ! [ -x "$(command -v $i)" ]; then
		$installer install -y $i
	fi
done

##------UNCOMMENT NEXT LINES (the whole block) IF YOU WISH TO INSTALL THIS SOFTWARE FROM REPOSITORIES------
#If commented then newest versions will be installed later from source, see below
# if [[ $distribution =~ "Debian" ]]; then
	# $installer install -y parallel #Debian/Ubuntu/Fedora
	# $installer install -y bowtie2 #Debian/Ubuntu
	# $installer install -y samtools #Debian/Ubuntu/CentOS/Fedora
	# $installer install -y fastx-toolkit #Debian/Ubuntu
	# $installer install -y mafft #Debian/Ubuntu
	# $installer install -y fasttree #Debian/Ubuntu
# elif [[ $distribution =~ "OpenSUSE" ]]; then
	# $installer install -y gnu_parallel #OpenSUSE
# elif [[ $distribution =~ "Fedora" ]]; then
	# $installer install -y parallel #Debian/Ubuntu/Fedora
	# $installer install -y samtools #Debian/Ubuntu/CentOS/Fedora
	# $installer install -y fastx_toolkit #Fedora; if commented newest version will be installed later from source, see below
# elif [[ $distribution =~ "CentOS" ]]; then
	# $installer install -y samtools #Debian/Ubuntu/CentOS/Fedora
# fi

#------INSTALL OTHER SOFTWARE------
mkdir install
cd install

#MAFFT
if ! [ -x "$(command -v mafft)" ]; then
	wget http://mafft.cbrc.jp/alignment/software/mafft-7.305-without-extensions-src.tgz
	tar -xvf mafft-7.305-without-extensions-src.tgz
	rm mafft-7.305-without-extensions-src.tgz
	cd mafft-7.305-without-extensions/core
	make
	make install
	cd ../..
fi

#GNU parallel
if ! [ -x "$(command -v parallel)" ]; then
	wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
	tar xjvf parallel-latest.tar.bz2
	rm parallel-latest.tar.bz2
	cd parallel*
	./configure
	make
	make install
	cd ..
fi

#Samtools
if ! [ -x "$(command -v samtools)" ]; then
	wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
	tar xjvf samtools-1.3.1.tar.bz2
	rm samtools-1.3.1.tar.bz2
	cd samtools-1.3.1
	./configure --without-curses
	make
	make install
	cd ..
fi

#FastTree
if ! [ -x "$(command -v FastTree)" ]; then
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
fi

#RAxML
git clone https://github.com/stamatak/standard-RAxML
cd standard-RAxML
if ! [ -x "$(command -v raxmlHPC)" ]; then
	make -f Makefile.gcc
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-SSE3)" ]; then
	make -f Makefile.SSE3.gcc
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC­AVX)" ]; then
	make -f Makefile.AVX.gcc
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-PTHREADS)" ]; then
	make -f Makefile.PTHREADS.gcc
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-PTHREADS-SSE3)" ]; then
	make -f Makefile.SSE3.PTHREADS.gcc
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC­PTHREADS-AVX)" ]; then
	make -f Makefile.AVX.PTHREADS.gcc
	rm *.o
fi
cp raxmlHPC* /usr/local/bin
cd ..

#MstatX
if ! [ -x "$(command -v mstatx)" ]; then
	git clone https://github.com/gcollet/MstatX
	cd MstatX
	make
	cp mstatx /usr/local/bin
	cd ..
fi

#TrimAl
if ! [ -x "$(command -v trimal)" ]; then
	git clone https://github.com/scapella/trimal
	cd trimal/source
	make
	cp trimal /usr/local/bin
	cd ../..
fi

#Newick Utilities
if ! [ -x "$(command -v nw_reroot)" ]; then
	wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	tar xfz newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	rm newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	cp newick-utils-1.6/src/nw_* /usr/local/bin
fi

#bam2fastq
if ! [ -x "$(command -v bam2fastq)" ]; then
	wget https://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz
	tar xfz bam2fastq-1.1.0.tgz
	rm bam2fastq-1.1.0.tgz
	cd bam2fastq-1.1.0
	make
	cp bam2fastq /usr/local/bin
	cd ..
fi

#BLAT
#see, e.g., http://nix-bio.blogspot.cz/2013/10/installing-blat-and-blast.html for installing tips
if ! [ -x "$(command -v blat)" ]; then
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
fi

#FastUniq
if ! [ -x "$(command -v fastuniq)" ]; then
	wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz
	tar xfz FastUniq-1.1.tar.gz
	rm FastUniq-1.1.tar.gz
	cd FastUniq/source
	make
	cp fastuniq /usr/local/bin
	cd ../..
fi

#Bowtie2
if ! [ -x "$(command -v bowtie2)" ]; then
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-source.zip
	unzip bowtie2-2.2.9-source.zip
	rm bowtie2-2.2.9-source.zip
	cd bowtie2-2.2.9
	make
	cp bowtie2* /usr/local/bin
	cd ..
fi

#ococo
if ! [ -x "$(command -v ococo)" ]; then
	git clone --recursive https://github.com/karel-brinda/ococo
	cd ococo
	make -j
	cp ococo /usr/local/bin
	cd ..
fi

#p4 (only necessary for combining bootstrap support in Astral and Astrid trees)
#see http://p4.nhm.ac.uk/installation.html
#compilation probably fails on Fedora/CentOS/OpenSUSE - unable to find gsl!!! Solution???
if ! [ -x "$(command -v p4)" ]; then
	$installer install -y python-numpy #Debian/Ubuntu/OpenSUSE
	#$installer install -y numpy #CentOS, Fedora
	$installer install -y python-scipy #Debian, OpenSUSE
	#$installer install -y scipy #CentOS, Fedora
	$installer install -y libgsl0-dev #Debian
	#$installer install -y gsl-devel #CentOS, Fedora and OpenSUSE
	$installer install -y python-dev #Debian
	#$installer install -y python-devel #CentOS, Fedora and OpenSUSE
	git clone https://github.com/pgfoster/p4-phylogenetics
	cd p4-phylogenetics
	python setup.py build
	python setup.py install
	python setup.py build_ext -i
	cd ..
fi

#ococo (necessary for majority rule consensus building from mapped reads in BAM file)
#see https://github.com/karel-brinda/ococo
if ! [ -x "$(command -v ococo)" ]; then
	git clone --recursive https://github.com/karel-brinda/ococo
	cd ococo
	make -j
	cp ococo /usr/local/bin
	cd ..
fi

#Leave 'install' directory
cd ..

#Check if everything is installed correctly
echo -e "\nSoftware installed...checking for binaries in PATH"
for i in parallel bowtie2 ococo samtools bam2fastq java fastuniq perl blat mafft python python3 trimal mstatx FastTree nw_reroot nw_topology raxmlHPC raxmlHPC-PTHREADS R p4; do
	command -v $i >/dev/null 2>&1 || { echo -n $i; echo >&2 "...not found"; }
done
#Check R packages
for Rpackage in ape seqinr data.table; do
	R -q -e "aa <- file('Rtest', open='wt'); sink(aa, type='message'); require($Rpackage); sink(type='message'); close(aa)" > /dev/null
	if grep -Fq "no package called" Rtest; then
		echo -e "R package $Rpackage...not found"
	elif grep -Fq "Error" Rtest; then
		echo -e "R package $Rpackage...unable to load"
	fi
done
rm Rtest
echo "Necessary software installed, possible errors indicated above."
echo -e "If you don't see any *not found* your system is now ready to run HybPhyloMaker!\n"

#Clone HybPhyloMaker GitHub repository
git clone https://github.com/tomas-fer/HybPhyloMaker
cd HybPhyloMaker
chmod +x *.sh
chmod +x HybSeqSource/ASTRID
