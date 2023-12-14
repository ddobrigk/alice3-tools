#!/bin/bash
#
# This macro runs all steps in a strangeness tracking effort
#
# Step 1: create simulation via o2-sim
#
# Basic configurations
NEVENTS=50000
NCPUS=4

DIR=`pwd`

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size 10000000000"

echo "CPU #${1} starting"
mkdir job_${1}
cp config_custom_Pythia.cfg job_${1}
cp configPythia.ini job_${1}
cp generator_pythia8_PbPb.C job_${1}
cp pythia8_hi.cmnd job_${1}
cp pythia8_pp.cmnd job_${1}
cd job_${1}

##run the simulation
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
time o2-sim -j ${NCPUS} --field +15U -e TGeant3 -n ${NEVENTS} -g external --configFile configPythia.ini -m ZDC --configKeyValues "Diamond.width[2]=6.;Diamond.width[0]=0.005;Diamond.width[1]=0.005;" -o o2sim
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"

echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
echo "Doing master cleanup"
#rm o2sim_Kine.root
rm o2sim.root
rm o2sim_par.root
#rm o2sim_geometry.root
#rm o2sim_grp.root
rm MCStepLoggerVolMap.dat
rm o2sim_proc-cut.dat
rm MCStepLoggerSenVol.dat
rm gphysi.dat
echo "Done! Enjoy!"
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
