##########################################################################################################################
# INSTALL SOFTWARE NECESSARY FOR running HybPhyloMaker                                                                   #
# and clone HybPhyloMaker GitHub repository (incl. test dataset)                                                         #
# This should work on major Linux distribution (tested on Debian, Ubuntu, OpenSUSE, Fedora, CentOS and Scientific Linux) #
# Without changes works only on 64-bit platforms (x86_64)                                                                #
# This script MUST be run with root privileges!                                                                          #
#                                                                                                                        #
# Tomas Fer, 2017, 2018                                                                                                  #
# tomas.fer@natur.cuni.cz                                                                                                #
# https://github.com/tomas-fer/HybPhyloMaker                                                                             #
# v.1.6.2                                                                                                                #
##########################################################################################################################

#Carefully set your distribution
distribution=Debian #one of: Debian (also for Ubuntu), OpenSUSE, Fedora, CentOS (also for Scientific Linux)
#Change name of your default package management tool (apt-get on Debian/Ubuntu, zypper on OpenSUSE, yum on Fedora/CentOS/Scientific/RHEL)
installer=apt-get #one of apt-get, zypper, yum

echo -e "\n************************************************"
echo -e "Installation script for HybPhyloMaker is running"
echo -e "All necessary software will be installed."
echo -e "See installation logs in 'install' folder.\n"
echo -e "Your distribution: $distribution"
echo -e "Your package manager: $installer"
echo -e "************************************************\n"

mkdir -p install
cd install

#------INSTALL SOFTWARE FROM REPOSITORIES------
#Install software using default package manager specified above

#Compilation utilities, i.e., gcc, g++, make, autoconf
if [[ $distribution =~ "Debian" ]]; then
	for i in gcc g++ make automake; do
		if ! [ -x "$(command -v $i)" ]; then
			echo -e "Installing '$i'"
			$installer install -y $i &> ${i}_install.log
		fi
	done
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	for i in make automake; do
		if ! [ -x "$(command -v $i)" ]; then
			echo -e "Installing '$i'"
			$installer install -y $i &> ${i}_install.log
		fi
	done
	if ! [ -x "$(command -v gcc)" ]; then
		echo -e "Installing 'gcc'"
		$installer install -y gcc &> gcc_install.log
	fi
	if ! [ -x "$(command -v g++)" ]; then
		echo -e "Installing 'g++'"
		$installer install -y gcc-c++ &> g++_install.log
	fi
fi
if ! [ -x "$(command -v autoconf)" ]; then
	echo -e "Installing 'autoconf'"
	$installer install -y autoconf &> autoconf_install.log
fi
if ! [ -x "$(command -v autoreconf)" ]; then
	echo -e "Installing 'dh-autoreconf'"
	$installer install -y dh-autoreconf &> dh-autoreconf_install.log
fi

#Perl
if ! [ -x "$(command -v perl)" ]; then
	echo -e "Installing 'perl'"
	$installer install -y perl &> perl_install.log
fi

#Python
if ! [ -x "$(command -v python)" ]; then
	echo -e "Installing 'python'"
	$installer install -y python &> python_install.log
fi

#Python3
if ! [ -x "$(command -v python3)" ]; then
	if [[ $distribution =~ "CentOS" ]]; then
		echo -e "Installing 'python3'"
		$installer install -y epel-release &> python3_install.log #Only for CentOS
		$installer install -y python34 &>> python3_install.log #Only for CentOS
	else
		echo -e "Installing 'python3'"
		$installer install -y python3 &> python3_install.log #Does not work on CentOS
		
	fi
fi

#Pip
if [[ $distribution =~ "Debian" ]]; then
	if ! [ -x "$(command -v pip)" ]; then
		echo -e "Installing 'pip'"
		$installer install -y python-pip &> pip_install.log
	fi
	pip install --upgrade pip &>> pip_install.log
else
	if ! [ -x "$(command -v pip2.7)" ]; then
		echo -e "Installing 'pip'"
		$installer install -y python-pip &> pip_install.log
	fi
	pip2.7 install --upgrade pip &>> pip_install.log
fi

#Pip3 (also required for 'biopython' and 'kindel' installation, see below)
if ! [ -x "$(command -v pip3)" ]; then
	if [[ $distribution =~ "CentOS" ]]; then
		echo -e "Installing 'pip3'"
		$installer install -y python34-devel &>> python3_install.log #Only for CentOS
		$installer install -y python34-pip &> pip3_install.log #Only for CentOS
	elif [[ $distribution =~ "Debian" ]]; then
		echo -e "Installing 'pip3'"
		$installer install -y python3-dev &>> python3_install.log #Does not work on CentOS
		$installer install -y python3-pip &> pip3_install.log #Does not work on CentOS
	elif [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "Fedora" ]]; then
		echo -e "Installing 'pip3'"
		$installer install -y python3-devel &>> python3_install.log #Does not work on CentOS
		$installer install -y python3-pip &> pip3_install.log #Does not work on CentOS
	fi
fi
if [[ $distribution =~ "Fedora" ]]; then
	python3 -m pip install --upgrade pip &>> pip3_install.log
else
	pip3 install --upgrade pip &>> pip3_install.log
fi

#Biopython
if [ ! `python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null` ]; then
	echo -e "Installing 'biopython for python3'"
	pip3 install biopython &> biopython_install.log
fi

#Java
if ! [ -x "$(command -v java)" ]; then
	if [[ $distribution =~ "Debian" ]]; then
		echo -e "Installing 'java'"
		distrib=$(cat /etc/*release | grep ^ID= | cut -d'=' -f2)
		if [[ $distrib =~ "debian" ]]; then
			debver=$(cat /etc/debian_version | cut -d"." -f1)
			if [ "$debver" -eq "9" ]; then
				$installer install -y openjdk-8-jre &> java_install.log #Debian9/Ubuntu
			else
				$installer install -y openjdk-7-jre &> java_install.log #Debian9/Ubuntu
			fi
		elif [[ $distrib =~ "ubuntu" ]]; then
			ubuver=$(lsb_release -r -s | cut -d"." -f1)
			if [ "$ubuver" -eq "16" ]; then
				$installer install -y openjdk-8-jre &> java_install.log #Debian9/Ubuntu
			else
				$installer install -y openjdk-7-jre &> java_install.log #Debian9/Ubuntu
			fi
		fi
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		echo -e "Installing 'java'"
		$installer install -y java-1.7.0-openjdk.x86_64 &> java_install.log #Fedora/CentOS
	elif [[ $distribution =~ "OpenSUSE" ]]; then
		echo -e "Installing 'java'"
		$installer install -y java-1_7_0-openjdk &> java_install.log #OpenSUSE
	fi
fi

#zlib library
if [ ! "$(whereis zlib | grep /)" ]; then
	if [[ $distribution =~ "Debian" ]]; then
		echo -e "Installing 'zlib'"
		$installer install -y libpng-dev &> zlib_install.log #Debian/Ubuntu
		$installer install -y zlib1g-dev &>> zlib_install.log #Debian/Ubuntu
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
		echo -e "Installing 'zlib'"
		$installer install -y libpng-devel &> zlib_install.log #Fedora/CentOS/OpenSUSE
		$installer install -y zlib-devel &>> zlib_install.log #Fedora/CentOS/OpenSUSE
	fi
fi

#pkg-config
if [ ! "$(whereis pkg-config | grep /)" ]; then
	if [[ $distribution =~ "Debian" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
		echo -e "Installing 'pkg-config'"
		$installer install -y pkg-config &> pkg-config_install.log #Debian/Ubuntu/OpenSUSE
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		echo -e "Installing 'pkg-config'"
		$installer install -y pkgconfig &> pkg-config_install.log #Fedora/CentOS
	fi
fi

#libtool
if [ ! "$(whereis libtool | grep /)" ]; then
	echo -e "Installing 'libtool'"
	$installer install -y libtool &> libtool_install.log
fi

#R
#Comment for Ubuntu/Debian: you should install the newest version of R by adding CRAN mirror to /etc/apt/sources.list
#Look at the end of this script for an advice how to do that...
#Newer version (i.e., at least v3.3) should be installed before running this script!
if ! [ -x "$(command -v R)" ]; then
	if [[ $distribution =~ "Debian" ]]; then
		echo -e "Installing 'R'"
		$installer install -y r-base-dev &> R_install.log #Debian/Ubuntu
		$installer install -y gfortran &> gfortran_install.log #Debian/Ubuntu
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		echo -e "Installing 'R'"
		$installer install -y R &> R_install.log #Fedora/CentOS
		$installer install -y gcc-gfortran &> gfortran_install.log #Fedora/CentOS
	elif [[ $distribution =~ "OpenSUSE" ]]; then
		echo -e "Installing 'R'"
		$installer install -y R-base-devel &> R_install.log #OpenSUSE
		$installer install -y gcc-fortran &> gfortran_install.log #necessary in OpenSUSE for compilation of some R packages
	fi
fi

#R packages
for Rpackage in ape seqinr data.table openxlsx; do
	R -q -e "is.element('$Rpackage', installed.packages()[,1])" > testpackage
	if grep -Fxq "[1] FALSE" testpackage; then
		echo -e "Installing '$Rpackage for R'"
		R -q -e "install.packages('$Rpackage', repos='http://cran.rstudio.com/')" &> R_${Rpackage}_install.log
	else
		echo -e "R package $Rpackage already installed"
	fi
done
rm testpackage

#wget, tar, bzip2, bc, git
for i in wget tar bzip2 bc git; do
	if ! [ -x "$(command -v $i)" ]; then
		echo -e "Installing '$i'"
		$installer install -y $i &> ${i}_install.log
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

#MAFFT
if ! [ -x "$(command -v mafft)" ]; then
	echo -e "Installing 'mafft'"
	wget http://mafft.cbrc.jp/alignment/software/mafft-7.305-without-extensions-src.tgz &> mafft_install.log
	tar -xvf mafft-7.305-without-extensions-src.tgz 1>/dev/null
	rm mafft-7.305-without-extensions-src.tgz
	cd mafft-7.305-without-extensions/core
	make &>> ../../mafft_install.log
	make install &>> ../../mafft_install.log
	cd ../..
fi

#GNU parallel
if ! [ -x "$(command -v parallel)" ]; then
	echo -e "Installing 'GNU parallel'"
	wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 &> GNUparallel_install.log
	tar xjvf parallel-latest.tar.bz2 1>/dev/null
	rm parallel-latest.tar.bz2
	cd parallel*
	./configure &>> ../GNUparallel_install.log
	make &>> ../GNUparallel_install.log
	make install &>> ../GNUparallel_install.log
	cd ..
fi

#Samtools
if ! [ -x "$(command -v samtools)" ]; then
	echo -e "Installing 'samtools'"
	wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 &> samtools_install.log
	tar xjvf samtools-1.3.1.tar.bz2 1>/dev/null
	rm samtools-1.3.1.tar.bz2
	cd samtools-1.3.1
	./configure --without-curses &>> ../samtools_install.log
	make &>> ../samtools_install.log
	make install &>> ../samtools_install.log
	cd ..
fi

#FastTree
if ! [ -x "$(command -v FastTree)" ]; then
	echo -e "Installing 'FastTree'"
	mkdir FastTree
	cd FastTree
	wget http://www.microbesonline.org/fasttree/FastTree &> ../fasttree_install.log
	chmod +x FastTree
	cp FastTree /usr/local/bin
	wget http://www.microbesonline.org/fasttree/FastTree.c &>> ../fasttree_install.log
	gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm &>> ../fasttree_install.log
	chmod +x FastTreeMP
	cp FastTreeMP /usr/local/bin
	cd ..
fi

#EMBOSS
if ! [ -x "$(command -v transeq)" ]; then
	echo -e "Installing 'EMBOSS'"
	wget ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz &> EMBOSS_install.log
	tar xfz emboss-latest.tar.gz 1>/dev/null
	rm emboss-latest.tar.gz
	cd EMBOSS*
	./configure --without-x &>> ../EMBOSS_install.log
	make &>> ../EMBOSS_install.log
	ldconfig
	make install &>> ../EMBOSS_install.log
	ldconfig
	make install &>> ../EMBOSS_install.log
	cd ..
fi

#RAxML
git clone https://github.com/stamatak/standard-RAxML &> raxml_install.log
cd standard-RAxML
if ! [ -x "$(command -v raxmlHPC)" ]; then
	echo -e "Installing 'raxmlHPC'"
	make -f Makefile.gcc &>> ../raxml_install.log
	rm *.o
fi
if [[ `grep sse3 /proc/cpuinfo` ]]; then
	if ! [ -x "$(command -v raxmlHPC-SSE3)" ]; then
		echo -e "Installing 'raxmlHPC-SSE3'"
		make -f Makefile.SSE3.gcc &>> ../raxml_install.log
		rm *.o
	fi
fi
if [[ `grep avx /proc/cpuinfo` ]]; then
	if ! [ -x "$(command -v raxmlHPC­AVX)" ]; then
		echo -e "Installing 'raxmlHPC­AVX'"
		make -f Makefile.AVX.gcc &>> ../raxml_install.log
		rm *.o
	fi
fi
	if ! [ -x "$(command -v raxmlHPC-PTHREADS)" ]; then
	echo -e "Installing 'raxmlHPC-PTHREADS'"
	make -f Makefile.PTHREADS.gcc &>> ../raxml_install.log
	rm *.o
fi
if [[ `grep sse3 /proc/cpuinfo` ]]; then
	if ! [ -x "$(command -v raxmlHPC-PTHREADS-SSE3)" ]; then
		echo -e "Installing 'raxmlHPC-PTHREADS-SSE3'"
		make -f Makefile.SSE3.PTHREADS.gcc &>> ../raxml_install.log
		rm *.o
	fi
fi
if [[ `grep avx /proc/cpuinfo` ]]; then
	if ! [ -x "$(command -v raxmlHPC­PTHREADS-AVX)" ]; then
		echo -e "Installing 'raxmlHPC­PTHREADS-AVX'"
		make -f Makefile.AVX.PTHREADS.gcc &>> ../raxml_install.log
		rm *.o
	fi
fi
cp raxmlHPC* /usr/local/bin
cd ..

#MstatX
if ! [ -x "$(command -v mstatx)" ]; then
	echo -e "Installing 'MstatX'"
	git clone https://github.com/gcollet/MstatX &> mstatx_install.log
	cd MstatX
	make &>> ../mstatx_install.log
	cp mstatx /usr/local/bin
	cd ..
fi

#TrimAl
if ! [ -x "$(command -v trimal)" ]; then
	echo -e "Installing 'TrimAl'"
	git clone https://github.com/scapella/trimal &> trimal_install.log
	cd trimal/source
	make &>> ../../trimal_install.log
	cp trimal /usr/local/bin
	cd ../..
fi

#Newick Utilities
if ! [ -x "$(command -v nw_reroot)" ]; then
	echo -e "Installing 'Newick Utilities'"
	# wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz &> newickutil_install.log
	# tar xfz newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz 1>/dev/null
	# rm newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	# cp newick-utils-1.6/src/nw_* /usr/local/bin
	if ! [ -x "$(command -v bison)" ]; then
		$installer install -y bison &> bison_install.log
	fi
	if ! [ -x "$(command -v flex)" ]; then
		$installer install -y flex &> flex_install.log
	fi
	if ! [ -x "$(command -v autoreconf)" ]; then
		$installer install -y dh-autoreconf &> autoreconf_install.log
	fi
	git clone https://github.com/tjunier/newick_utils  &> newickutil_install.log
	cd newick_utils/ &>> newickutil_install.log
	libtoolize &>> newickutil_install.log
	aclocal &>> newickutil_install.log
	autoheader &>> newickutil_install.log
	autoreconf -fi &>> newickutil_install.log
	./configure &>> newickutil_install.log
	make &>> newickutil_install.log
	make check &>> newickutil_install.log
	make install &>> newickutil_install.log
	ldconfig &>> newickutil_install.log
	cd ..
fi

#bam2fastq
if ! [ -x "$(command -v bam2fastq)" ]; then
	echo -e "Installing 'bam2fastq'"
	wget https://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz &> bam2fastq_install.log
	tar xfz bam2fastq-1.1.0.tgz 1>/dev/null
	rm bam2fastq-1.1.0.tgz
	cd bam2fastq-1.1.0
	make &>> ../bam2fastq_install.log
	cp bam2fastq /usr/local/bin
	cd ..
fi

#BLAT
#see, e.g., http://nix-bio.blogspot.cz/2013/10/installing-blat-and-blast.html for installing tips
if ! [ -x "$(command -v blat)" ]; then
	echo -e "Installing 'BLAT'"
	wget https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip &> blat_install.log
	unzip blatSrc35.zip 1>/dev/null
	rm blatSrc35.zip
	cd blatSrc
	#echo $MACHTYPE #if you get x86_64-pc-linux-gnu continue with next line, if something else write appropriate short name
	MACHTYPE=x86_64
	export MACHTYPE
	mkdir -p lib/$MACHTYPE
	mkdir -p ~/bin/$MACHTYPE
	make &>> ../blat_install.log
	cp ~/bin/$MACHTYPE/blat /usr/local/bin
	cd ..
fi

#FastUniq
if ! [ -x "$(command -v fastuniq)" ]; then
	echo -e "Installing 'FastUniq'"
	wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz &> fastuniq_install.log
	tar xfz FastUniq-1.1.tar.gz 1>/dev/null
	rm FastUniq-1.1.tar.gz
	cd FastUniq/source
	make &>> ../fastuniq_install.log
	cp fastuniq /usr/local/bin
	cd ../..
fi

#Bowtie2
if ! [ -x "$(command -v bowtie2)" ]; then
	echo -e "Installing 'Bowtie2'"
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-source.zip &> bowtie2_install.log
	unzip bowtie2-2.2.9-source.zip 1>/dev/null
	rm bowtie2-2.2.9-source.zip
	cd bowtie2-2.2.9
	make &>> ../bowtie2_install.log
	cp bowtie2* /usr/local/bin
	cd ..
fi

#bwa
if ! [ -x "$(command -v bwa)" ]; then
	echo -e "Installing 'bwa'"
	wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.16a.tar.bz2 &> bwa_install.log
	tar jxf bwa-0.7.16a.tar.bz2 1>/dev/null
	rm bwa-0.7.16a.tar.bz2
	cd bwa-0.7.16a
	make &>> ../bwa_install.log
	cp bwa /usr/local/bin
	cd ..
fi

#ococo (necessary for majority rule consensus building from mapped reads in BAM file)
#see https://github.com/karel-brinda/ococo
if ! [ -x "$(command -v ococo)" ]; then
	echo -e "Installing 'OCOCO'"
	git clone --recursive https://github.com/karel-brinda/ococo &>> ococo_install.log
	cd ococo
	make -j &>> ../ococo_install.log
	make install  &>> ../ococo_install.log
	cd ..
fi

#p4 (only necessary for combining bootstrap support in Astral and Astrid trees)
#see http://p4.nhm.ac.uk/installation.html
#For compilation on Fedora/CentOS/OpenSUSE you need to specify where 'gsl' is installed (in setup.py) - modification of 'setup.py' is included below
if ! [ -x "$(command -v p4)" ]; then
	echo -e "Installing 'p4'"
	if [[ $distribution =~ "Debian" ]]; then
		if ! [[ `pip show numpy | grep Version` ]]; then
			pip install numpy &> numpy_install.log
			#$installer install -y python-numpy &> numpy_install.log #Debian/Ubuntu/OpenSUSE
		fi
		if  ! [[ `pip show scipy | grep Version` ]]; then
			pip install scipy &> scipy_install.log
			$installer install -y python-scipy &> scipy_install.log #Debian, OpenSUSE
		fi
	elif [[ $distribution =~ "OpenSUSE" ]]; then
		if ! [[ `pip2.7 show numpy | grep Version` ]]; then
			#pip2.7 install numpy &> numpy_install.log
			$installer install -y python-numpy-devel &> numpy_install.log #Debian/Ubuntu/OpenSUSE
		fi
		if  ! [[ `pip2.7 show scipy | grep Version` ]]; then
			pip2.7 install scipy &> scipy_install.log
			#$installer install -y python-scipy &> scipy_install.log #Debian, OpenSUSE
		fi
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		if ! [[ `pip2.7 show numpy | grep Version` ]]; then
			#pip2.7 install numpy &> numpy_install.log
			$installer install -y numpy &> numpy_install.log #CentOS, Fedora
		fi
		if  ! [[ `pip2.7 show scipy | grep Version` ]]; then
			#pip2.7 install scipy &> scipy_install.log
			$installer install -y scipy &> scipy_install.log #CentOS, Fedora
		fi
	fi
	
	if [[ $distribution =~ "Debian" ]]; then
		$installer install -y libgsl0-dev &> libgsl_install.log #Debian
		$installer install -y python-dev &> python-dev_install.log #Debian
	elif [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		$installer install -y gsl-devel &> libgsl_install.log #CentOS, Fedora and OpenSUSE
		$installer install -y python-devel &> python-dev_install.log #CentOS, Fedora and OpenSUSE
	fi
	
	if [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		$installer install redhat-rpm-config &> rpm-config_install.log
	fi
	#install python module 'future'
	if [[ $distribution =~ "Debian" ]]; then
		if ! [[ `pip show future | grep Version` ]]; then
			pip install future &> python-future_install.log
		fi
	else
		if ! [[ `pip2.7 show future | grep Version` ]]; then
			pip2.7 install future &> python-future_install.log
		fi
	fi
	git clone https://github.com/pgfoster/p4-phylogenetics &> p4_install.log
	cd p4-phylogenetics
	#Modify setup.py to be able to find gsl
	if [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "Fedora" ]]; then
		replace="my_include_dirs = [\'\/usr\/include\/\']"
		#sed -i.bak "46s/.*/$replace/" setup.py
		#sed -i.bak2 "47s/# //" setup.py
	fi
	python setup.py build &>> ../p4_install.log
	python setup.py install &>> ../p4_install.log
	python setup.py build_ext -i &>> ../p4_install.log
	cd ..
fi

#kindel (necessary for majority rule consensus building from mapped reads in BAM file)
#requires v0.1.4 (higher version will not work as they are missing option for threshold!)
#see https://pypi.python.org/pypi/kindel
if ! [ -x "$(command -v kindel)" ]; then
	echo -e "Installing 'kindel'"
	# if [[ $distribution =~ "Debian" ]]; then
		# $installer install -y python3-dev &> python3-dev_install.log #Debian
	# elif [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		# $installer install -y python3-devel &> python3-dev_install.log #CentOS, Fedora and OpenSUSE
	# fi
	if [[ $distribution =~ "Fedora" ]]; then
		python3 -m pip install 'kindel==0.1.4' &>> kindel_install.log
	else
		pip3 install 'kindel==0.1.4' &>> kindel_install.log
	fi
fi

#other python modules (mainly for PartitionFinder)
if [[ $distribution =~ "Debian" ]]; then
	if ! [[ `pip show pandas | grep Version` ]]; then
		echo -e "Installing 'pandas for python'"
		pip install pandas &> python-pandas_install.log
		#$installer install -y python-pandas &> python-pandas_install.log #Debian
	fi
	if ! [[ `pip show scikit-learn | grep Version` ]]; then
		echo -e "Installing 'scikit-learn for python'"
		pip install scikit-learn &> python-sklearn_install.log
		#$installer install -y python-sklearn &> python-sklearn_install.log #Debian
	fi
	if ! [[ `pip show tables | grep Version` ]]; then
		echo -e "Installing 'tables for python'"
		pip install tables &> python-tables.log
	fi
	if ! [[ `pip show parsing | grep Version` ]]; then
		echo -e "Installing 'parsing for python'"
		pip install parsing &> python-parsing.log
	fi
elif [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	if ! [[ `pip2.7 show pandas | grep Version` ]]; then
		echo -e "Installing 'pandas for python'"
		pip2.7 install pandas &> python-pandas_install.log #CentOS, Fedora and OpenSUSE
	fi
	if ! [[ `pip2.7 show scikit-learn | grep Version` ]]; then
		echo -e "Installing 'scikit-learn for python'"
		pip2.7 install scikit-learn &> python-sklearn_install.log #CentOS, Fedora and OpenSUSE
	fi
	if ! [[ `pip show tables | grep Version` ]]; then
		echo -e "Installing 'tables for python'"
		pip2.7 install tables &> python-tables.log
	fi
	if ! [[ `pip show parsing | grep Version` ]]; then
		echo -e "Installing 'parsing for python'"
		pip2.7 install parsing &> python-parsing.log
	fi
fi

#OpenMPI
if [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	module load mpi/openmpi-x86_64 2>/dev/null
	if ! [ -x "$(command -v mpicc)" ]; then
		echo -e "\nDo you want to install ExaML [y/n]:"
		read -s -n 1 instexaml
		if [[ $instexaml =~ "y" ]]; then
			echo -e "Installing 'OpenMPI'"
			$installer -y install openmpi openmpi-devel &> openmpi_install.log
			$installer -y install environment-modules &>> openmpi_install.log
			echo -e "\nPlease logout and login back and run again 'install_software.sh' in order to complete ExaML installation" && exit 3
		fi
	else
		instexaml=yes
	fi
else
	if ! [ -x "$(command -v mpicc)" ]; then
		echo -e "Installing 'OpenMPI'"
		if [[ $distribution =~ "Debian" ]]; then
			instexaml=yes
			$installer -y install mpi-default-dev &> openmpi_install.log
		elif [[ $distribution =~ "OpenSUSE" ]]; then
			instexaml=yes
			wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz &> openmpi_install.log
			tar -xvf openmpi-2.1.1.tar.gz 1>/dev/null
			cd openmpi-2.1.1/
			./configure &>> openmpi_install.log
			make &>> openmpi_install.log
			make install &>> openmpi_install.log
			ldconfig &>> openmpi_install.log
			cd ..
		fi
	fi
fi

#ExaML
if [[ $instexaml =~ "y" ]]; then
	if ! [ -x "$(command -v examl)" ]; then
		echo -e "Installing 'ExaML'"
		if [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
			module load mpi/openmpi-x86_64
		fi
		git clone https://github.com/stamatak/ExaML &> examl_install.log
		cd ExaML/parser
		make -f Makefile.SSE3.gcc &> ../../examl_install.log
		cp parse-examl /usr/local/bin
		cd ../examl
		make -f Makefile.SSE3.gcc &> ../../examl_install.log
		cp examl /usr/local/bin
		if [[ `grep avx /proc/cpuinfo` ]]; then
			make -f Makefile.AVX.gcc &> ../../examl_install.log
			cp examl-AVX /usr/local/bin
		fi
		cd ../..
	fi
fi

#BUCKy
if ! [ -x "$(command -v bucky)" ]; then
	echo -e "Installing 'BUCKy'"
	wget http://dstats.net/download/http://www.stat.wisc.edu/~ane/bucky/v1.4/bucky-1.4.4.tgz &> bucky_install.log
	tar -xzvf bucky-1.4.4.tgz 1>/dev/null
	rm bucky-1.4.4.tgz
	cd bucky-1.4.4/src/
	make &>> ../../bucky_install.log
	cp mbsum /usr/local/bin
	cp bucky /usr/local/bin
	cd ../..
fi

#Leave 'install' directory
cd ..

#Check if everything is installed correctly
echo -e "\n**************************************************************"
echo -e "Software installed...checking for binaries in PATH"
rm not_installed.txt 2>/dev/null
for i in parallel bowtie2 bwa ococo kindel samtools transeq bam2fastq java fastuniq perl blat mafft python python3 trimal mstatx FastTree nw_reroot nw_topology raxmlHPC raxmlHPC-PTHREADS examl R p4 bucky; do
	#command -v $i >/dev/null 2>&1 || { echo -n $i; echo >&2 "...not found"; }
	command -v $i >/dev/null 2>&1 && echo ${i}...OK || { echo -n $i; echo >&2 "...not found"; echo $i >> not_installed.txt; }
done
sed -i.bak 's/transeq/EMBOSS/' not_installed.txt 2>/dev/null
sed -i.bak2 's/nw_reroot//' not_installed.txt 2>/dev/null
sed -i.bak4 's/nw_topology/NewickUtilities/' not_installed.txt 2>/dev/null
sed -i.bak3 's/parallel/GNUparallel/' not_installed.txt 2>/dev/null
sed -i.bak4 '/^$/d' not_installed.txt 2>/dev/null
rm -f *.bak* 2>/dev/null
echo -e "\n**************************************************************"
echo -e "List of software which is not installed:"
if [ ! -f not_installed.txt ]; then
	echo -e "ALL NECESSARY SOFTWARE INSTALLED"
else
	cat not_installed.txt 2>/dev/null
	echo -e "Check particular log file(s) to look for possible solution..."
fi
echo -e "**************************************************************"

#Check R packages
echo -e "\n**************************************************************"
echo -e "Checking R packages"
for Rpackage in ape seqinr data.table; do
	R -q -e "aa <- file('Rtest', open='wt'); sink(aa, type='message'); require($Rpackage); sink(type='message'); close(aa)" > /dev/null
	if grep -Fq "no package called" Rtest; then
		echo -e "R package $Rpackage...not found"
	elif grep -Fq "Error" Rtest; then
		echo -e "R package $Rpackage...unable to load"
	else
		echo -e "R package $Rpackage...OK"
	fi
done
rm Rtest
echo -e "**************************************************************"

echo -e "\n**************************************************************"
echo -e "Necessary software installed, possible errors indicated above."
echo -e "**************************************************************"
echo -e "\nIf you don't see any *not found* your system is now ready to run HybPhyloMaker!"
echo -e "Consult appropriate '_install.log' (in case *not found* was reported) to solve the installation problem."
echo -e "If there is a problems with R packages, installation of newer R version might solve the problem."

#Clone HybPhyloMaker GitHub repository
if [[ ! -d HybPhyloMaker ]]; then
	echo -e "\nCloning HybPhyloMaker GitHub repository..."
	git clone https://github.com/tomas-fer/HybPhyloMaker &> HybPhyloMaker_gitcloning.log
	cd HybPhyloMaker
	chmod +x *.sh
	chmod +x HybSeqSource/ASTRID
fi

echo -e "\nInstalation script finished.\n"

# Successfully tested on
# - Ubuntu 16.04
# - Ubuntu 14.04 LTS #Newer R version necessary, see below!!!
# - Debian 9.0
# - Debian 8.6 #Newer R version necessary, see below!!!
# - OpenSUSE 42
# - CentOS 7.2
# - Fedora 24
# - Scientific Linux 7.2


#Tips for older versions - run these command before running this 'install_software.sh'

#Ubuntu 14.04
#**********************************************************************************
#too old R version (see https://cran.r-project.org/bin/linux/ubuntu/ how to update)
# codename=$(lsb_release -c -s)
# echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu $codename/" | sudo tee -a /etc/apt/sources.list > /dev/null
# apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
# add-apt-repository ppa:marutter/rdev
# apt-get update
# apt-get upgrade

#CentOS 6
#**********************************************************************************
#too old g++ (4.7 and later required for proper building of 'OCOCO')
#see http://ask.xmodulo.com/upgrade-gcc-centos.html how to upgrade on CentOS
# wget http://people.centos.org/tru/devtools-1.1/devtools-1.1.repo -P /etc/yum.repos.d
# sh -c 'echo "enabled=1" >> /etc/yum.repos.d/devtools-1.1.repo'
# yum install -y devtoolset-1.1
# scl enable devtoolset-1.1 bash
#
#and too old glibc. GLIBC_2.14 is required - and it is not easy to upgrade, better to upgrade to CentOS 7
#NewickUtilities will not work!

#**********************************************************************************
#Debian 7, 8
#too old R version (see https://cran.r-project.org/bin/linux/debian/ how to update)
# codename=$(lsb_release -c -s)
# echo "deb http://cran.cnr.berkeley.edu/bin/linux/debian ${codename}-cran3/" | sudo tee -a /etc/apt/sources.list > /dev/null
# apt-key adv --keyserver keys.gnupg.net --recv-key 6212B7B7931C4BB16280BA1306F90DE5381BA480
# apt-get update
# apt-get upgrade
#
#only Debian 7 - too old glibc. GLIBC_2.14 is required - and it is not easy to upgrade, better to upgrade to Debian 8
#NewickUtilities will not work!
