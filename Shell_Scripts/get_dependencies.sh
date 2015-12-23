#!/bin/bash
#	A shell script to download dependencies
#	for LRT_Predict
#	Author: Paul Hoffman

#	This script requires Git, Wget, and CMake v3.0
#	or higher to run.

#	Make directory for dependencies

set -e
set -o pipefail

#	Set our variables from input
#	The first variable is the dependencies diretory, and everything else
#	are the dependencies that need to be downloaded
DEPSDIR="$1"
declare -a MISSING=("${@:2}")

#	Make the dependencies directory and cd into it
echo ${DEPSDIR}
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
				#	Use OS-specific installer
			if [[ $(uname) == "Linux" ]]
			then
				find . -maxdepth 1 -not -name "*linux*" -name "*.tar.gz" -delete
			elif [[ $(uname) == "Darwin" ]]
			then
				find . -maxdepth 1 -not -name "*macosx*" -name "*.tar.gz" -delete
			else
				echo "Not on a supported operating system!"
				exit 1
			fi
				#	Extract the source file
				#	Designed to always install the latest verion of BLAST
				#	unless NCBI changes their file naming scheme
			tar -xvzf ncbi*
				#	Move NCBI BLAST+ to a standard folder
			mv ncbi* ncbi_blast+
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
		"PASTA" )
			#   Install Pasta
				#	Download Pasta from GitHub using wget
			wget https://github.com/smirarab/pasta/archive/master.zip
			unzip master.zip
				#   Download sate-tools from GitHub
			if [[ $( uname ) == "Linux" ]] # If we are on Linux
			then # Get the linux version of sate-tools
				wget -O sate-tools.zip https://github.com/smirarab/sate-tools-linux/archive/master.zip
			elif [[ $( uname ) == "Darwin" ]] # If we're on Mac OS X
			then # Get the Mac OS X version of sate-tools
				wget -O sate-tools.zip https://github.com/smirarab/sate-tools-mac/archive/master.zip
			else # If we aren't on either
				echo "Please use Mac OS X or Linux"
				exit 1
			fi
				#	Extract sate-tools
			unzip sate-tools.zip
				#	Cleanup the zip files
			rm -rf *.zip
				#	Move sate-tools to a properly named directory
			mv sate-tools* $(echo sate-tools* | cut -f 1,2,3 -d '-' )
			cd pasta-master
			python setup.py develop --user
			cd "${DEPSDIR}"
			;;
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
			wget https://github.com/veg/hyphy/archive/master.zip
			unzip master.zip
				#	Installing HyPhy requires CMake 3.0 or higher to be installed
			cd hyphy-master
				#	Actually install HyPhy
			sed 's/EXCLUDE_FROM_ALL//g' CMakeLists.txt > CMAKE.txt
			mv CMAKE.txt CMakeLists.txt
				#	Allow non-root users to install HyPhy
			cmake -DINSTALL_PREFIX=$(pwd)
				#	Install HyPhy
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
