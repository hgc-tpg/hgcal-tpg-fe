#!/bin/bash
#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
# #Array: Size=$#, an element=$1, all element = $@

printf "Start Running job at: ";/bin/date
printf "Worker node hostname: ";/bin/hostname


###########################################################
# Set CMSSW enviroment for accessible root,boost and others
###########################################################
HOME=/home/hep/idas
currdir=$PWD
echo HOME: $HOME
echo currdir: $currdir
echo _CONDOR_SCRATCH_DIR: ${_CONDOR_SCRATCH_DIR}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $HOME/CMSSW/CMSSW_14_0_0_pre1/src
eval `scramv1 runtime -sh`
ls
cd ${currdir}
echo "CMSSW_BASE : $CMSSW_BASE"
printf "SCRAM :  ${SCRAM_ARCH}\n"
which root
echo $ROOTSYS
##########################################################


###########################################################
# Get the input arguments
###########################################################
echo "All arguements: "$@
echo "Number of arguements: "$#
# mass=$1
# year=$2
# channel=$3
# cattype=$4
# rundir=${_CONDOR_SCRATCH_DIR}
# random=${RANDOM}
# clusproc=$5
# iloop=$6

infile=$1
ofindex=$2
nevents=$3
sidelength=$4
ofextn=$5
sampletype=$6
iloop=$7
clusproc=$8
rundir=${_CONDOR_SCRATCH_DIR}
random=${RANDOM}
##########################################################


###########################################################
# Set the output and running/execution directory
###########################################################
outdir=$sampletype
outpath=Result_iter$iloop/$outdir
condorOutDir=$HOME/stage2_emulation_results/$outpath
testdir=$rundir/$outpath

if [ ! -d $testdir ] ; then
    mkdir -p $testdir
else
    rm $testdir/*.*
fi

if [ ! -d $condorOutDir ] ; then
    mkdir -p $condorOutDir
fi
##########################################################


##########################################################
# Go inside run/exec directory and download dependencies
##########################################################
#mv local.tar.gz $testdir/
cd $testdir
rsync -avP  $HOME/EmulatorChain/hgcal-tpg-fe/condor/tmpLog_s2emu_iter$iloop/local.tar.gz .
tar --strip-components=1 -zxf local.tar.gz
ls -lah
##########################################################


##########################################################
# run
##########################################################
make clean
make
echo "stage2SemiEmulator.exe $infile $ofindex $nevents $sidelength $ofextn $sampletype \> out_${$sampletype}_${ofextn}_${clusproc}.log "
./stage2SemiEmulator.exe $infile $ofindex $nevents $sidelength $ofextn $sampletype > out_${sampletype}_${ofextn}_${clusproc}.log 2>&1
##########################################################


##########################################################
# transfer output file
##########################################################
rsync -avP  TPGS2Emu_tree_${ofextn}_${ofindex}.root $condorOutDir
rsync -avP  stage2SemiEmulator_${ofextn}_${ofindex}.root $condorOutDir
rsync -avP  local.tar.gz $condorOutDir/../
rsync -avP  out_${sampletype}_${ofextn}_${clusproc}.log $HOME/EmulatorChain/hgcal-tpg-fe/condor/tmpLog_s2emu_iter$iloop/logs/
printf "Done transfer: ";/bin/date
cd $testdir
printf "Job Done at: ";/bin/date
##########################################################
