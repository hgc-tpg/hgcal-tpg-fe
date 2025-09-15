#!/bin/bash

CURRDIR="$PWD"
islxplus=`hostname | grep lxplus | wc -l`
ishgcbeamtestpc=`hostname | grep hgcbeamtestpc | wc -l`

if [ $islxplus -eq 1 ] ; then
    downloaddir=/tmp/$USER
else
    downloaddir=/tmp
fi

if [ ! -d $downloadir/hgcal10glinkreceiver ] ; then
    cd $downloaddir
    git clone https://gitlab.cern.ch/pdauncey/hgcal10glinkreceiver.git
else
    cd $downloaddir/hgcal10glinkreceiver
    git pull
fi

# if [ ! -d $downloadir/slow_control_configuration ] ; then
#     cd $downloadir
#     git clone https://gitlab.cern.ch/hgcal-daq-sw/serenity-daq/slow_control_configuration.git
# else
#     cd $downloadir/slow_control_configuration
#     git pull
# fi


cd $CURRDIR
rsync -avP $downloaddir/hgcal10glinkreceiver/common .
#rsync -avP $downloadir/slow_control_configuration cfgmap/

if [ -L dat ] ; then
    unlink dat
fi

if [ $islxplus -eq 1 ] ; then
    ln -s /eos/cms/store/group/dpg_hgcal/tb_hgcal/2024/BeamTestSep/HgcalBeamtestSep2024 dat
fi
if [ $ishgcbeamtestpc -eq 1 ] ;  then
    ln -s /till/HgcalBeamtest2024_TEST dat
fi
