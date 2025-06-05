#!/bin/bash

input_dir="/p/lustre1/milton3/lambda_zdc_production/generate_data/eic/zdc+20cmPbWO4_300GeV_dynamicrange_pi0_h4_10GeV-300GeV_log10discrete_0.0-0.23deg_09_11_2024/gensim"
output_dir="/p/lustre1/milton3/lambda_zdc_production/generate_data/eic/zdc+20cmPbWO4_300GeV_dynamicrange_pi0_h4_10GeV-300GeV_log10discrete_0.0-0.23deg_09_11_2024/recosim_sigma1_thresh4_DR250GeV_16bit"
setup_script="/p/lustre1/milton3/lambda_zdc_production/generate_data/eic/setup_env.sh"
EIC_DIR="/p/lustre1/milton3/lambda_zdc_production/generate_data/eic"
num_events=10000
mkdir -p "$output_dir"

# Loop through all .root files in the input directory
for simfile in "$input_dir"/*.root; do
    # Extract the INFO part from the simfile name
    base_name=$(basename "$simfile" .root)
    info=$(echo "$base_name" | sed 's/^sim_//')

    # Set output file name
    recofile="reco_${info}.root"


    # Create a SLURM job script
    job_script=$(mktemp)

    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --job-name=reco_${info}
#SBATCH --output=${output_dir}/reco_${info}.out
#SBATCH --error=${output_dir}/reco_${info}.err

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
echo "source ${EIC_DIR}/setup_env.sh" >> \${GENSCRIPTNAME}
echo "export JUGGLER_SIM_FILE=${simfile}"  >> \${GENSCRIPTNAME}
echo "export JUGGLER_REC_FILE=${recofile}"  >> \${GENSCRIPTNAME}
echo "export JUGGLER_N_EVENTS=${num_events}"  >> \${GENSCRIPTNAME}
echo "gaudirun.py ${EIC_DIR}/local/share/eic_geometry/scripts/zdc_reco.py"  >> \${GENSCRIPTNAME}
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
