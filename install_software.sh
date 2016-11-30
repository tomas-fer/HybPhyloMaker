#########################################################################################################
# INSTALL SOFTWARE NECESSARY FOR running HybPhyloMaker                                                  #
# and clone HybPhyloMaker GitHub repository (incl. test dataset)                                        #
# This should work on major Linux distribution (tested on Debian, Ubuntu, OpenSUSE, Fedora, and CentOS) #
# Without changes works only on 64-bit platforms (x86_64)                                               #
# This script MUST be run with root privileges!                                                         #
#                                                                                                       #
# Tomas Fer, 2016                                                                                       #
# tomas.fer@natur.cuni.cz                                                                               #
# https://github.com/tomas-fer/HybPhyloMaker                                                            #
# v.1.3.2                                                                                               #
#########################################################################################################

#Carefully set your distribution
distribution=Debian #one of: Debian (also for Ubuntu), OpenSUSE, Fedora, CentOS)
#Change name of your default package management tool (apt-get on Debian/Ubuntu, zypper on OpenSUSE, yum on Fedora/CentOS/RHEL/Scientific)
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

#Compilation utilities, i.e., gcc, g++, make
if [[ $distribution =~ "Debian" ]]; then
	for i in gcc g++ make; do
		if ! [ -x "$(command -v $i)" ]; then
			echo -e "Installing '$i'"
			$installer install -y $i &> ${i}_install.log
		fi
	done
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	for i in gcc-c++ g++ make; do
		if ! [ -x "$(command -v $i)" ]; then
			echo -e "Installing '$i'"
			$installer install -y $i &> ${i}_install.log
		fi
	done
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

#Java
if [[ $distribution =~ "Debian" ]]; then
	echo -e "Installing 'java'"
	$installer install -y openjdk-7-jre &> java_install.log #Debian/Ubuntu
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	echo -e "Installing 'java'"
	$installer install -y java-1.7.0-openjdk.x86_64 &> java_install.log #Fedora/CentOS
elif [[ $distribution =~ "OpenSUSE" ]]; then
	echo -e "Installing 'java'"
	$installer install -y java-1_7_0-openjdk &> java_install.log #OpenSUSE
fi

#zlib library
if [[ $distribution =~ "Debian" ]]; then
	echo -e "Installing 'zlib'"
	$installer install -y libpng-dev &> zlib_install.log #Debian/Ubuntu
	$installer install -y zlib1g-dev &>> zlib_install.log #Debian/Ubuntu
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	echo -e "Installing 'zlib'"
	$installer install -y libpng-devel &> zlib_install.log #Fedora/CentOS/OpenSUSE
	$installer install -y zlib-devel &>> zlib_install.log #Fedora/CentOS/OpenSUSE
fi

#pkg-config
if [[ $distribution =~ "Debian" ]] || [[ $distribution =~ "OpenSUSE" ]]; then
	echo -e "Installing 'pkg-config'"
	$installer install -y pkg-config &> pkg-config_install.log #Debian/Ubuntu/OpenSUSE
elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
	echo -e "Installing 'pkg-config'"
	$installer install -y pkgconfig &> pkg-config_install.log #Fedora/CentOS
fi

#R
#Comment for Ubuntu: you should install the newest version of R by adding CRAN mirror to /etc/apt/sources.list
#(see, e.g., https://cran.r-project.org/bin/linux/ubuntu/README.html)
# codename=$(lsb_release -c -s)
# echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu $codename/" | sudo tee -a /etc/apt/sources.list > /dev/null
# apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
# add-apt-repository ppa:marutter/rdev
# apt-get update
# apt-get upgrade
# apt-get install -y r-base r-base-dev
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

#R packages
for Rpackage in ape seqinr data.table; do
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

#RAxML
git clone https://github.com/stamatak/standard-RAxML &> raxml_install.log
cd standard-RAxML
if ! [ -x "$(command -v raxmlHPC)" ]; then
	echo -e "Installing 'raxmlHPC'"
	make -f Makefile.gcc &>> ../raxml_install.log
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-SSE3)" ]; then
	echo -e "Installing 'raxmlHPC-SSE3'"
	make -f Makefile.SSE3.gcc &>> ../raxml_install.log
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC­AVX)" ]; then
	echo -e "Installing 'raxmlHPC­AVX'"
	make -f Makefile.AVX.gcc &>> ../raxml_install.log
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-PTHREADS)" ]; then
	echo -e "Installing 'raxmlHPC-PTHREADS'"
	make -f Makefile.PTHREADS.gcc &>> ../raxml_install.log
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC-PTHREADS-SSE3)" ]; then
	echo -e "Installing 'raxmlHPC-PTHREADS-SSE3'"
	make -f Makefile.SSE3.PTHREADS.gcc &>> ../raxml_install.log
	rm *.o
fi
if ! [ -x "$(command -v raxmlHPC­PTHREADS-AVX)" ]; then
	echo -e "Installing 'raxmlHPC­PTHREADS-AVX'"
	make -f Makefile.AVX.PTHREADS.gcc &>> ../raxml_install.log
	rm *.o
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
	wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz &> newickutil_install.log
	tar xfz newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz 1>/dev/null
	rm newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	cp newick-utils-1.6/src/nw_* /usr/local/bin
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

#ococo (necessary for majority rule consensus building from mapped reads in BAM file)
#see https://github.com/karel-brinda/ococo
if ! [ -x "$(command -v ococo)" ]; then
	echo -e "Installing 'OCOCO'"
	wget https://github.com/karel-brinda/ococo/archive/0.1.2.4.tar.gz &> ococo_install.log
	tar xfz 0.1.2.4.tar.gz 1>/dev/null
	cd ococo-0.1.2.4/ext
	git clone https://github.com/samtools/htslib &>> ococo_install.log
	cd ..
	make -j &>> ../ococo_install.log
	cp ococo /usr/local/bin
	cd ..
fi

#p4 (only necessary for combining bootstrap support in Astral and Astrid trees)
#see http://p4.nhm.ac.uk/installation.html
#For compilation on Fedora/CentOS/OpenSUSE you need to specify where 'gsl' is installed (in setup.py) - modification of 'setup.py' is included below
if ! [ -x "$(command -v p4)" ]; then
	echo -e "Installing 'p4'"
	if [[ $distribution =~ "Debian" ]]; then
		$installer install -y python-numpy &> numpy_install.log #Debian/Ubuntu/OpenSUSE
		$installer install -y python-scipy &> scipy_install.log #Debian, OpenSUSE
	elif [[ $distribution =~ "OpenSUSE" ]]; then
		$installer install -y python-numpy-devel &> numpy_install.log #Debian/Ubuntu/OpenSUSE
		$installer install -y python-scipy &> scipy_install.log #Debian, OpenSUSE
	elif [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		$installer install -y numpy &> numpy_install.log #CentOS, Fedora
		$installer install -y scipy &> scipy_install.log #CentOS, Fedora
	fi
	
	if [[ $distribution =~ "Debian" ]]; then
		$installer install -y libgsl0-dev &> libgsl_install.log #Debian
		$installer install -y python-dev &> python-dev_install.log #Debian
	elif [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "Fedora" ]] || [[ $distribution =~ "CentOS" ]]; then
		$installer install -y gsl-devel &> libgsl_install.log #CentOS, Fedora and OpenSUSE
		$installer install -y python-devel &> python-dev_install.log #CentOS, Fedora and OpenSUSE
	fi
	
	if [[ $distribution =~ "Fedora" ]]; then
		$installer install redhat-rpm-config &> rpm-config_install.log
	fi
	
	git clone https://github.com/pgfoster/p4-phylogenetics &> p4_install.log
	cd p4-phylogenetics
	#Modify setup.py to be able to find gsl
	if [[ $distribution =~ "OpenSUSE" ]] || [[ $distribution =~ "CentOS" ]] || [[ $distribution =~ "Fedora" ]]; then
		replace="my_include_dirs = [\'\/usr\/include\/\']"
		sed -i.bak "45s/.*/$replace/" setup.py
		sed -i.bak2 "46s/# //" setup.py
	fi
	python setup.py build &>> ../p4_install.log
	python setup.py install &>> ../p4_install.log
	python setup.py build_ext -i &>> ../p4_install.log
	cd ..
fi

#Leave 'install' directory
cd ..

#Check if everything is installed correctly
echo -e "\n**************************************************************"
echo -e "Software installed...checking for binaries in PATH"
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
echo -e "**************************************************************"
echo -e "\nIf you don't see any *not found* your system is now ready to run HybPhyloMaker!"
echo -e "Consult appropriate '_install.log' (in case *not found* was reported) to solve the installation problem."
echo -e "If there is a problems with R packages, installation of newer R version might solve the problem."

#Clone HybPhyloMaker GitHub repository
echo -e "\nCloning HybPhyloMaker GitHub repository..."
git clone https://github.com/tomas-fer/HybPhyloMaker &> HybPhyloMaker_gitcloning.log
cd HybPhyloMaker
chmod +x *.sh
chmod +x HybSeqSource/ASTRID

echo -e "\nInstalation script finished.\n"

# Successfully tested on
# - Ubuntu 14.04 #Newer R version necessary!!!
# - Debian 8.6
# - OpenSUSE 42
# - CentOS 7.2
# - Fedora 24


# Tips for older versions

#CentOS 6
#too old g++ (4.7 and later required for proper building of 'OCOCO')
#see http://ask.xmodulo.com/upgrade-gcc-centos.html how to upgrade on CentOS
#wget http://people.centos.org/tru/devtools-1.1/devtools-1.1.repo -P /etc/yum.repos.d
#sh -c 'echo "enabled=1" >> /etc/yum.repos.d/devtools-1.1.repo'
#yum install devtoolset-1.1
#scl enable devtoolset-1.1 bash

#Ubuntu
#too old R version
