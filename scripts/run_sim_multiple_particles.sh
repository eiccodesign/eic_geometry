#!/bin/bash

function print_the_help {
  echo "USAGE: ${0} -n <nevents> -part <"\"particle\""> -p <momentum> -thmin <theta_min> -thmax <theta_max> -dist <"\"dist\""> "
  echo "  OPTIONS: "
  echo "    -n,--nevents            Number of events"
  echo "    -part,--particle        particle type"
  echo "    -pmin, --energy_min     particle min. momentum (GeV)"
  echo "    -pmax, --energy_max     particle min. momentum (GeV)"
  echo "    -dist, --distribution   energy distribution (fixed, uniform, gaussian, log10continuous, log10discrete)"
  echo "                            for fixed distribution, make sure pmin and pmax are equal"
  echo "                            Set to log10continuous by default"
  echo "    -thmin, --theta_min     Minimum theta (degrees)"
  echo "    -thmax, --theta_max     Maximum theta (degrees)" 
  echo "                            theta (polar angle) will generate uniformly between theta_min and theta_max"
  echo "    -jid, --jid             job ID"  
  echo "    -nmin, --nmin           minimum number of particles"
  echo "    -nmax, --nmax           maximum number of particles"  
  exit
}

# Input simulation parameters
particle="neutron"
energy_min=100
energy_max=100
num_events=10
theta_min=0.0 # in degrees
theta_max=0.3 # in degrees
phi_min=0. # in degrees
phi_max=360. # in degrees 
distribution="fixed"
physics_list="FTFP_BERT"
job_id=000
min_num_particles=1
max_num_particles=10


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
elif [ "$1" = "-pmin" -o "$1" = "--energy_min" ]; then
   energy_min=$2
   shift 2 # past argument
elif [ "$1" = "-pmax" -o "$1" = "--energy_max" ]; then
   energy_max=$2
   shift 2 # past argument
elif [ "$1" = "-thmin" -o "$1" = "--theta_min" ]; then
   theta_min=$2
   shift 2 # past argument
elif [ "$1" = "-thmax" -o "$1" = "--theta_max" ]; then
   theta_max=$2
   shift 2 # past argument
elif [ "$1" = "-dist" -o "$1" = "--distribution" ]; then
   distribution=$2
   shift 2 # past argument
elif [ "$1" = "-jid" -o "$1" = "--jid" ]; then
   job_id=$2
   shift 2 # past argument
elif [ "$1" = "-nmin" -o "$1" = "--nmin" ]; then
   min_num_particles=$2
   shift 2 # past argument
elif [ "$1" = "-nmax" -o "$1" = "--nmax" ]; then
   max_num_particles=$2
   shift 2 # past argument
else
   break
fi
done

# Output file names
info_string="${particle}_${distribution}_${energy_min}GeV-${energy_max}GeV_theta_${theta_min}deg-${theta_max}deg_${min_num_particles}-${max_num_particles}particles_${job_id}"
hepmcfile="gen_${info_string}.hepmc"
simfile="sim_${info_string}.edm4hep.root"
recofile="reco_${info_string}.edm4hep.root"

# Generating hepmc file
root -l -b -q "${DETECTOR_PATH}/hepmc_generation/gen_multiple_particles.cxx(\
${num_events},\
\"${hepmcfile}\",\
\"${particle}\",\
${theta_min},\
${theta_max},\
${phi_min},\
${phi_max},\
${energy_min},\
${energy_max},\
\"${distribution}\",\
${min_num_particles},\
${max_num_particles})"

# Running simulation
npsim \
   --compactFile ${DETECTOR_PATH}/${DETECTOR}.xml \
   --numberOfEvents ${num_events} \
   --physicsList ${physics_list} \
   --inputFiles ${hepmcfile} \
   --outputFile ${simfile}  || exit

# Running reconstruction
export JUGGLER_SIM_FILE=${simfile} JUGGLER_REC_FILE=${recofile} JUGGLER_N_EVENTS=${num_events}
gaudirun.py ${DETECTOR_PATH}/scripts/${DETECTOR}_reco.py