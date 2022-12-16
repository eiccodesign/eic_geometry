#!/bin/bash

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -part <particle> -p <momentum> "
  echo "  OPTIONS: "
  echo "    -n,--nevents     Number of events"
  echo "    -part,--particle particle type"
  echo "    -p, --momentum   particle momentum (GeV)"
  exit
}

# Default simulation parameters
particle="pi-"
beam_energy=50
num_events=1000
theta_min=2.83 # in degrees
theta_max=2.83 # in degrees
phi_min=0. # in degrees
phi_max=360. # in degrees 
distribution=3 # Momentum distribution: 0=fixed, 1=uniform, 2=Gaussian, 3=log uniform
# Note: Log uniform is energy independent and will select energies in {2, 4, 8, 16, 32, 64, 128} GeV
physics_list="FTFP_BERT_HP"

while [ True ]; do
if [ "$1" = "--help" -o "$1" = "-h" ]; then
   print_the_help
   shift
elif [ "$1" = "-part" -o "$1" = "--particle" ]; then
   particle=$2
   shift 2 # past argument
elif [ "$1" = "-n" -o "$1" = "--nevents" ]; then
   num_events=$2
   shift 2 # past argument
elif [ "$1" = "-p" -o "$1" = "--momentum" ]; then
   beam_energy=$2
   shift 2 # past argument
else
   break
fi
done

# Output file names
info_string="${particle}_${beam_energy}GeV_theta_${theta_min}-${theta_max}deg"
hepmcfile="gen_${info_string}.hepmc"
simfile="sim_${info_string}.edm4hep.root"
recofile="reco_${info_string}.edm4hep.root"

# Generating hepmc file
root -l -b -q "${DETECTOR_PATH}/hepmc_generation/gen_particles.cxx(\
${num_events},\
\"${hepmcfile}\",\
\"${particle}\",\
${theta_min},\
${theta_max},\
${phi_min},\
${phi_max},\
${beam_energy},\
${distribution})"

# Running simulation
npsim \
   --compactFile ${DETECTOR_PATH}/hadron_endcap.xml \
   --numberOfEvents ${num_events} \
   --physicsList ${physics_list} \
   --inputFiles ${hepmcfile} \
   --outputFile ${simfile}  || exit
   
# Deleting hepmcfile  
rm ${hepmcfile}

# Running reconstruction
export JUGGLER_SIM_FILE=${simfile} JUGGLER_REC_FILE=${recofile} JUGGLER_N_EVENTS=${num_events}
gaudirun.py ${DETECTOR_PATH}/scripts/hadron_endcap_reco.py
