#!/bin/sh

# Use either run_sim_hepmc or run_sim_gps
simulation="${INSERT_PATH}/scripts/run_sim_hepmc.sh"
energies=(1 2 5 10 20 30 40 50 60 70 80 90 100)

for energy in "${energies[@]}"
do
	sed -i "s/beam_energy=.*/beam_energy=${energy}/" ${simulation}
	echo "Starting energy ${energy} GeV"
	time source ${simulation}
done