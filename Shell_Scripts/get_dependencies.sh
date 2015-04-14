#!/bin/bash
#	Written by Paul Hoffman

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
BLAST=`pwd`
cd ../../
#export PATH=$PATH:<your additions seperated with colons>:

#http://ssh-commands.blogspot.com/2013/01/the-ultimate-wget-download-guide-with.html

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
PRANK=`pwd`
cd ../../

#	Install HyPhy
	#	Pull HyPhy from GitHub
	#	Requires Git to be installed
	#	Does not require a GitHub account
git clone https://github.com/veg/hyphy.git
	#	Installing HyPhy requires CMake 3.0 or higher to be installed
cd hyphy
cmake -DINSTALL_PREFIX=./ ./
make install
HYPHY=`pwd`
cd ../

#	Install BioPython
	#	Pull BioPython with wget
	#	Designed to always pull latest version 
	#	unless BioPython changes their downloads page
mkdir BioPython
cd BioPython
wget --no-directories --progress=bar --timestamping -r -np -l 1 -A.tar.gz http://biopython.org/DIST
	#	Get rid of all files that aren't latest
find . -maxdepth 1 -not -name "biopython*" -name "*.tar.gz" -delete
mv `ls -t *.tar.gz | head -n 1` ../
cd ../
rm -rf BioPython
	#	Extract the source file
	#	Designed to always install the latest version of BioPython
	#	unless BioPython changes their file naming scheme
tar -xvzf biopython*
	#	Install BioPython
	#	Requires superuser access
cd biopython*
python setup.py build
python setup.py test
python setup.py install --home=./
BIO=`pwd`
cd ../

export PATH=$PATH:$BLAST:$PRANK:$HYPHY:$BIO
