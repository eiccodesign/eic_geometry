#!/bin/bash

input_dir="/p/lustre2/milton3/epic_zdc_lambda/zdc_lambdadecay_h4_50GeV-300GeV_log10continuous_0.0-0.3mrad_fullacceptance_12_05_2024/gensim/"
output_dir="/p/lustre2/milton3/epic_zdc_lambda/zdc_lambdadecay_h4_50GeV-300GeV_log10continuous_0.0-0.3mrad_fullacceptance_12_05_2024/sebouh_eicrecon_recosim/"
setup_script="/p/lustre2/milton3/epic_zdc_lambda/thisepic.sh"
eicrecon_script="/p/lustre2/milton3/epic_zdc_lambda/local/bin/eicrecon-this.sh"
source ${setup_script}
source ${eicrecon_script}
EIC_DIR="/p/lustre2/milton3/epic_zdc_lambda/"

mkdir -p "$output_dir"

# Loop through all .root files in the input directory
for simfile in "$input_dir"/*.root; do
    # echo "Starting ${simfile}"
    # Extract the INFO part from the simfile name
    base_name=$(basename "$simfile" .root)

    info=$(echo "$base_name" | sed 's/^sim_//')

    # Set output file name
    recofile="reco_${info}.root"

    job_script=$(mktemp)

    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --job-name=reco_${info}
#SBATCH --output=${output_dir}/reco_${info}.out
#SBATCH --error=${output_dir}/reco_${info}.err
#SBATCH -t 03:00:00
export GENSCRIPTNAME="${base_name}.sh"
echo \${GENSCRIPTNAME}
if [ -f "\${GENSCRIPTNAME}" ]; then
    echo "\${GENSCRIPTNAME} exists and will be removed..."
    rm "\${GENSCRIPTNAME}"
fi
cd ${EIC_DIR}
#Write to Script
echo "#!/bin/bash" > \${GENSCRIPTNAME}
echo -en "\n" >> \${GENSCRIPTNAME}
echo "source ${setup_script}" >> \${GENSCRIPTNAME}
echo "source ${eicrecon_script}" >> \${GENSCRIPTNAME}
echo "eicrecon -Pdd4hep:xml_files=${DETECTOR_PATH}/${DETECTOR}.xml -Ppodio:output_file=${recofile} -Ppodio:output_collections=HcalFarForwardZDCRecHits,HcalFarForwardZDCClusters,MCParticles ${simfile}" >> \${GENSCRIPTNAME}
echo "mv "${recofile}" "${output_dir}/""  >> \${GENSCRIPTNAME}
chmod 700 \${GENSCRIPTNAME}
bash ${EIC_DIR}/eic-shell -- ./\${GENSCRIPTNAME}
EOT

    # Submit the SLURM job
    sbatch "$job_script"

    # Remove the temporary job script
    # rm "$job_script"

done

echo "All jobs submitted."
