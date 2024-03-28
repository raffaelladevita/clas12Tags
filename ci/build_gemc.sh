#!/usr/bin/env zsh

# Purpose: compiles gemc and installs it in gemc

# Container run:
# docker run -it --rm --platform linux/amd64 jeffersonlab/gemc:dev-g4v10.7.4-fedora36-cvmfs sh
# git clone http://github.com/gemc/clas12Tags /root/clas12Tags && cd /root/clas12Tags
# ./ci/build_gemc.sh


# if we are in the docker container, we need to load the modules
if [[ -z "${DISTTAG}" ]]; then
    echo "\nNot in container"
else
    echo "\nIn container: ${DISTTAG}"
    source  /app/localSetup.sh
fi

function compileGEMC {
	# getting number of available CPUS
	copt=" -j"`getconf _NPROCESSORS_ONLN`" OPT=1"
	echo
	echo Compiling GEMC with options: "$copt"
	scons SHOWENV=1 SHOWBUILD=1 $copt
	# checking existence of executable
	ls gemc
	if [ $? -ne 0 ]; then
		echo gemc executable not found
		exit 1
	fi
}

cd source
compileGEMC
echo "copying gemc to "$GEMC for experiment tests
cp gemc $GEMC
echo
echo "content of "$GEMC":"
ls -lrt $GEMC
