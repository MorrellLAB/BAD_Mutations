#!/bin/bash
#	A shell script to download dependencies
#	for LRT_Predict
#	Author: Paul Hoffman

#	This script requires Git, Wget, and CMake v3.0
#	or higher to run.

#	Make directory for dependencies

set -e
set -u
set -o pipefail

#	Set our variables from input
#	The first variable is the dependencies diretory, and everything else
#	are the dependencies that need to be downloaded
DEPSDIR=$1
MISSING=${@:1}

#	Make the dependencies directory and cd into it
mkdir -p $DEPSDIR
cd $DEPSDIR

for x in $MISSING
do
	case $x in
		"tBLASTx" )
			cd $DEPSDIR
			#	Install BLAST
				#	Pull BLAST from NCBI with wget
				#	Designed to always pull latest version 
				#	unless NCBI changes their ftp heirarchy
			wget --no-directories --progress=bar -r -A.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
				#	Get rid of all OS-specifc installers
			find . -maxdepth 1 -not -name "*src.tar.gz" -name "*.tar.gz" -delete
				#	Extract the source file
				#	Designed to always install the latest verion of BLAST
				#	unless NCBI changes their file naming scheme
			tar -xvzf ncbi*
				#	Install BLAST
			cd `find . -maxdepth 1 -type d -name "ncbi*"`/c++
			./configure
			make
				#	Cleanup tarball
			cd $DEPSDIR
			rm ncbi*tar.gz
			;;
		"PRANK" )
			cd $DEPSDIR
			#	Install PRANK
				#	Pull PRANK from Wasabi
				#	Designed to always pull latest version
				#	unless Wasabi changes their downloads page
			wget --no-directories --progress=bar -r -np -l 1 -A.tgz http://wasabiapp.org/download/prank
				#	Get rid of all OS-specific installers
				find . -maxdepth 1 -not -name "prank.source*" -name "*.tgz" -delete
				#	Extract the source file
				#	Designed to always install the lastest version of PRANK
				#	unless Wasabi changes their file naming scheme
			tar -xvzf prank*
				#	Install PRANK
			cd prank-msa/src
			make
				#	Cleanup tarball and robots.txt
			cd $DEPSDIR
			rm prank*.tgz robots.txt
			;;
		"requests" )
			cd $DEPSDIR
			#	Install Requests
				#	Pull Requests from GitHub
				#	Requires Git to be installed
				#	Does not require a GitHub account
			git clone git://github.com/kennethreitz/requests.git
				#	Install Requests
			cd requests
			REQ=`pwd`
			mkdir $REQ/modules
			MOD=$REQ/modules
			export PYTHONPATH="$PYTHONPATH:$REQ:$MOD"
			python setup.py build
			python setup.py install --install-base="." --install-lib='$base/modules' --install-scripts='$base/bin' --install-data='$base/data'/ --install-headers='$base/include/'
			cd $DEPSDIR
			;;
		"HyPhy" )
			cd $DEPSDIR
			#	Install HyPhy
				#	Pull HyPhy from GitHub
				#	Requires Git to be installed
				#	Does not require a GitHub account
			git clone https://github.com/veg/hyphy.git
				#	Installing HyPhy requires CMake 3.0 or higher to be installed
			cd hyphy
			cmake -DINSTALL_PREFIX=./ ./
			make install
			cd $DEPSDIR
			;;
		"Bio" )
			cd $DEPSDIR
			#	Install BioPython
				#	Pull BioPython with wget
				#	Designed to always pull latest version 
				#	unless BioPython changes their downloads page
			wget -O - http://biopython.org/DIST > biopython.txt
			BIO=`grep -E -o 'biopython-[0-9]\.[0-9a-b\.]*tar\.gz' biopython.txt | tail -1`
			wget --no-directories --progress=bar --timestamping -r http://biopython.org/DIST/$BIO
				#	Extract the source file
				#	Designed to always install the latest version of BioPython
				#	unless BioPython changes their file naming scheme
			tar -xvzf $BIO
				#	Install BioPython
			cd biopython*
			python setup.py build
			python setup.py test
			python setup.py install --home=./
				#	Cleanup tarball and biopython.txt
			cd $DEPSDIR
			rm biopython*.tar.gz biopython.txt
			;;
		* )
			echo "Nothing missing"
			break
			;;
	esac
done
