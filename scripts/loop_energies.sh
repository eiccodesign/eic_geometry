#!/bin/bash

simulation="${DETECTOR_PATH}/scripts/run_sim_hepmc.sh"
energies=(10 20 30 40 50 60 70 80 90 100)

for energy in "${energies[@]}"
do
	sed -i "s/beam_energy=.*/beam_energy=${energy}/" ${simulation}
	echo "Starting energy ${energy} GeV"
	time ${simulation}
done