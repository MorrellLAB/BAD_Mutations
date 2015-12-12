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
DEPSDIR="$1"
MISSING="${@:1}"

#	Make the dependencies directory and cd into it
mkdir -p "${DEPSDIR}"
cd "${DEPSDIR}"

#   A function to check for, download, and isntall Pip
function getPip() {
    if ! `command -v pip > /dev/null 2> /dev/null`
    then
        wget https://bootstrap.pypa.io/get-pip.py
        python get-pip.py --user
        echo export PATH='$PATH':"${HOME}"/.local/bin >> "${HOME}"/.bash_profile
        source "${HOME}"/.bash_profile
        rm get-pip.py
    fi
}

for x in "${MISSING}"
do
	case "${x}" in
		"tBLASTx" )
			cd "${DEPSDIR}"
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
			cd "${DEPSDIR}"
			rm ncbi*tar.gz
			;;
		"PRANK" )
			cd "${DEPSDIR}"
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
			cd "${DEPSDIR}"
			rm prank*.tgz robots.txt
			;;
        "pasta" )
            #   Install Pasta
                #   Pull Pasta from GitHub
                #	Requires Git to be installed
				#	Does not require a GitHub account
            git clone https://github.com/smirarab/pasta.git
                #   Pull sate-tools from GitHub
                #	Requires Git to be installed
				#	Does not require a GitHub account
            if [[ $( uname ) -eq "Linux" ]] # If we are on Linux
            then # Get the linux version of sate-tools
                git clone https://github.com/smirarab/sate-tools-linux.git
            elif [[ $( uname ) -eq "Darwin" ]] # If we're on Mac OS X
            then # Get the Mac OS X version of sate-tools
                https://github.com/smirarab/sate-tools-mac.git sate-tools
            else # If we aren't on either
                echo "Please use Mac OS X or Linux"
                exit 1
            fi
            cd pasta
            python setup.py develop --user
            cd "${DEPSDIR}"
		"requests" )
			cd "${DEPSDIR}"
			#	Install Requests
				#	Check for Pip
				#	Install Requests using Pip
			getPip
            pip install --user requests
			;;
		"HyPhy" )
			cd "${DEPSDIR}"
			#	Install HyPhy
				#	Pull HyPhy from GitHub
				#	Requires Git to be installed
				#	Does not require a GitHub account
			git clone https://github.com/veg/hyphy.git
				#	Installing HyPhy requires CMake 3.0 or higher to be installed
			cd hyphy
			cmake -DINSTALL_PREFIX=./ ./
			make install
			cd "${DEPSDIR}"
			;;
		"Bio" )
			cd "${DEPSDIR}"
			#	Install BioPython
				#	Check for Pip
				#	Install BioPython using Pip
			getPip
            pip install biopython
			;;
		* )
			echo "Nothing missing"
			break
			;;
	esac
done
