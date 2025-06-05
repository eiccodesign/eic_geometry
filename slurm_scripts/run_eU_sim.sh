#!/bin/bash

output_dir="/p/lustre1/milton3/eU_generation/eU_data"
EIC_DIR="/p/lustre1/milton3/eU_generation"
DETECTOR_DIR="/p/lustre1/milton3/eU_generation/local/share/epic/epic_farforward_sansLYSO.xml"
num_events_per_job=5000
num_jobs=20
mkdir -p "$output_dir"

# Loop through all .root files in the input directory
for ((i=0; i<num_jobs; i++))
do
    start_event=$((i*$num_events_per_job))
    end_event=$((start_event+num_events_per_job))
    info="eU_events_${start_event}-${end_event}"

    # Set output file name
    simfile="sim_${info}.edm4hep.root"
    recofile="reco_${info}.edm4hep.root"

    # Create a SLURM job script
    job_script=$(mktemp)

    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --job-name=${info}
#SBATCH --output=${output_dir}/${info}.out
#SBATCH --error=${output_dir}/${info}.err
#SBATCH --time=40:00:00  # Set the job time limit to 40 hours

export GENSCRIPTNAME="${info}.sh"
echo \${GENSCRIPTNAME}
if [ -f "\${GENSCRIPTNAME}" ]; then
    echo "\${GENSCRIPTNAME} exists and will be removed..."
    rm "\${GENSCRIPTNAME}"
fi
cd ${EIC_DIR}
#Write to Script
echo "#!/bin/bash" > \${GENSCRIPTNAME}
echo -en "\n" >> \${GENSCRIPTNAME}
echo "source ${EIC_DIR}/setup_env.sh" >> \${GENSCRIPTNAME}
echo "npsim --compactFile ${DETECTOR_DIR} --numberOfEvents ${num_events_per_job} --skipNEvents ${start_event} --inputFiles ab_filtered_eU_0.hepmc --outputFile ${simfile}"  >> \${GENSCRIPTNAME}
echo "eicrecon -Pdd4hep:xml_files=${DETECTOR_DIR} -Ppodio:output_file=${recofile} -Ppodio:output_collections=HcalFarForwardZDCRecHits,MCParticles ${simfile}" >> \${GENSCRIPTNAME}
echo "mv "${simfile}" "${output_dir}/""  >> \${GENSCRIPTNAME}
echo "mv "${recofile}" "${output_dir}/""  >> \${GENSCRIPTNAME}
chmod 700 \${GENSCRIPTNAME}
bash ${EIC_DIR}/eic-shell -- ./\${GENSCRIPTNAME}
EOT

    # Submit the SLURM job
    sbatch "$job_script"


done

echo "All jobs submitted."
