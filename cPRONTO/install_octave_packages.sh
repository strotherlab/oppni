#!/bin/sh
# script for installing octave support packages packages required for OPPNI
# L. Mark Prati  mprati@research.baycrest.org
#
echo -e "\nChecking for Octave support packages. Missing packages will be installed from Octave forge to ~/octave"
echo -e "Be patient, this process may take up to 20 minutes to complete...\n"
sleep 1 
cd ~
#attemp module load command first
module "load nixpkgs/16.09 gcc/7.3.0 octave/4.4.1" 1> /dev/null 2>&1;
status_octave=$?
#if it fails check if an octave version is on the the PATH
if [ "$status_octave" -ne 0 ]; then
	which octave 1> /dev/null 2>&1;
	status_octave=$?
	if [ "$status_octave" -ne 0 ]; then
	    echo -e "\n\e[31mOctave package installation failed to complete ...." 	
    	echo "Unable to load or locate octave-cli"    	
        echo -e 'Install Octave 4.4.1 then rerun "install_octave_packages.sh\e[0m"\n'
        exit
	else
		echo -e "\nA local instal of Octave has been found. Pakage installation is continuing ....\n"
	fi
fi 
#test verson installed and isuue warning NOTE version does not affect package script install.
octave_version="$(octave-cli --no-init-file --eval 'disp(version)')"
#version_invalid="$octave_version<4.4.1" | bc -l
#echo "valid = $version_invalid"
#if [ $(echo "$octave_version < 4.4.1"|bc) ]; then   	
if [ $octave_version != '4.4.1' ]; then
	echo "WARNING - OPPNI requires Octave version 4.4.1 or higher, version detected is: $octave_version)'"	
fi	 
echo "Searching Octave forge for packages ...')"
if ls ~/octave/io* 1> /dev/null 2>&1; then
	echo "Skipping io - local packages found"
	status_io=0
else	
	echo "Installing io package! You can ignore depreciated warning messages during installation ..."
	sleep 1
	octave-cli --no-init-file --eval "pkg install -forge io"
	status_io=$?
fi
if ls ~/octave/control* 1> /dev/null 2>&1; then
	echo "Skipping control - local packages found"
	status_ctrl=0
else	
	echo "Installing control package ..."
	octave-cli --no-init-file --eval "pkg install -forge control"
	status_ctrl=$?
fi	
if ls ~/octave/struct* 1> /dev/null 2>&1; then
	echo "Skipping struct - local packages found"
	status_struc=0
else	
	echo "Installing struct package ..."
	octave-cli --no-init-file --eval "pkg install -forge struct"
	status_struc=$?
fi
if ls ~/octave/statistics* 1> /dev/null 2>&1; then
	echo "Skipping statistics - local packages found"
	status_stats=0
else	
	echo "Installing statistics package ..."
	octave-cli --no-init-file --eval "pkg install -forge statistics"
	status_stats=$?
fi	
if ls ~/octave/signal* 1> /dev/null 2>&1; then
	echo "Skipping signal - local packages found"
	status_sig=0
else	
	echo "Installing signal package ..."
	octave-cli --no-init-file --eval "pkg install -forge signal"
	status_sig=$?
fi
if ls ~/octave/optim* 1> /dev/null 2>&1; then
	echo "Skipping optim - local packages found"
	status_opt=0
else	
	echo "Installing optim package ..."
	octave-cli --no-init-file --eval "pkg install -forge optim"
	status_opt=$?
fi
if [ "$status_io" -ne 0 ] || [ "$status_ctrl" -ne 0 ] || [ "$status_struc" -ne 0 ] || [ "$status_stats" -ne 0 ] || [ "$status_sig" -ne 0 ] || [ "$status_opt" -ne 0 ] 
then
    echo -e "\e[31mOctave package installation failed to complete ...."
    echo "The following package(s) need to be re-installed:"
    if [ "$status_io" -ne 0 ]; then
    	echo "Pacakge: io"
    fi
    if [ "$status_ctrl" -ne 0 ]; then
    	echo "Pacakge: control"
    fi
	if [ "$status_struc" -ne 0 ]; then
    	echo "Pacakge: struct"
	fi
    if [ "$status_stats" -ne 0 ]; then
    	echo "Pacakge: statistics"
	fi
    if [ "$status_sig" -ne 0 ]; then
    	echo "Pacakge: signal"
	fi
    if [ "$status_opt" -ne 0 ]; then
    	echo "Pacakge: optim"
	fi
	echo -e 'Try re-running "install_octave_packages.sh\e[0m"'
else
	echo "Octave package installation check has completed ...."
	echo "Any required missing packages have been installed"
	ls ~/octave 
fi    
exit