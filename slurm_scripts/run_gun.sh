#!/bin/bash

# NOTE: THE COMMANDS NEED TO BE MODIFIED FOR YOUR USE CASE
# THIS USES A NEUTRON EVENT WITH EICRECON
#===========================
# User Specified Directories
#===========================
export EIC_DIR=/p/lustre2/milton3/epic_zdc_lambda
export OUTPUT_DIR=/p/lustre2/milton3/epic_zdc_lambda
export DETECTOR_PATH=${EIC_DIR}/local/share/eic_geometry
export num_events="50"
export theta_min=0
export theta_max=3
export PARTICLE="neutron"
export ENERGY_MIN="50"
export ENERGY_MAX="300"
export DISTRIBUTION="discrete" # //Momentum distribution options: fixed, uniform, Gaussian, discrete, log10continuous

JOB_ARRAY_ID=${SLURM_ARRAY_TASK_ID}                           
SLURM_JOB_ID=${SLURM_JOB_ID}

function print_the_help {
  echo"USAGE: ${0} -n <nevents> -d <output_dir> -p <particle> "
  echo "  OPTIONS: "
  echo "    -n,--nevents     Number of events to generate"
  echo "    -d,--directory   Directory for storing output root files, as well as temporary job scripts"
  echo "    -j,--jobname     Name of Job used for labelling large batches of submissions"
  echo "    -pmin           Minimum particle momentum (GeV)"
  echo "    -pmax           Maximum particle momentum (GeV)"
  echo "    -th_max           Theta Maximum (deg)"
  echo "    -th_min          Theta Minimum(deg)"
  echo "    -part,--particle    Particle Species"
  echo "                     Allowed Particles: pi0, pi+, pi-, ka0, ka+, ka-, proton, neutron, e-, e+,photon"
  exit
}
POSITIONAL=()
while [[ $# -gt 0 ]] 
do
  key="$1"
  case ${key} in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
     -d|--directory)
      export OUTPUT_DIR="$2"
      shift #past argument
      shift #past value
      ;;
    -j|--jobname)
      export JOB_NAME="$2"
      shift #past argument
      shift #past value
      ;;
    -part|--particle)
      export PARTICLE="$2"
      shift #past argument
      shift #past value
      ;;
    -n|--nevents)
      export num_events="$2"
      shift #past argument
      shift #past value
      ;;
    -th_max)
      export theta_max="$2"
      shift #past argument
      shift #past value
      ;;
    -th_min)
      export theta_min="$2"
      shift #past argument
      shift #past value
      ;;
    -pmin)
      export ENERGY_MIN="$2"
      shift #past argument
      shift #past value
      ;;
    -pmax)
      export ENERGY_MAX="$2"
      shift #past argument
      shift #past value
      ;;
     -jid)
      export jid="$2"
      shift #past argument
      shift #past value                                                                                                                         
      ;;

    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $2"
      print_the_help
      shift # past argument
      break
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters
#Check the Directories
printf "\n EIC_DIR set to ${EIC_DIR} \n" 
printf "\n Output Files will be moved to ${OUTPUT_DIR} \n" 

#=====================================================
# Crete Directories for storing Root Files and Scripts
#=====================================================
#Directory for run scripts
export TEMPSCRIPT_DIR=${OUTPUT_DIR}/tempscripts
if [ ! -d "${TEMPSCRIPT_DIR}" ]; then
    mkdir -p ${TEMPSCRIPT_DIR}
fi

#Reco Root Directory 
export RECO_DIR=${OUTPUT_DIR}/recosim
if [ ! -d "${RECO_DIR}" ]; then
    mkdir -p ${RECO_DIR}
fi

#Gen Root Directory
export GEN_DIR=${OUTPUT_DIR}/gensim
if [ ! -d "${GEN_DIR}" ]; then
    mkdir -p ${GEN_DIR}
fi

#Gen Hepmc Directory
export HEPMC_DIR=${GEN_DIR}/hepmc
if [ ! -d "${HEPMC_DIR}" ]; then
    mkdir -p ${HEPMC_DIR}
fi

NAME_TAG="${SLURM_JOB_ID}_${JOB_ARRAY_ID}"
echo "________ ${NAME_TAG}"
#================================================
# Create file to execute upon entering eic-shell
#================================================
#Make the script that runs inside the container
export GENSCRIPTNAME="${NAME_TAG}.sh"
if [ -f "${GENSCRIPTNAME}" ]; then
    echo "${GENSCRIPTNAME} exists and will be removed..."
    rm "${GENSCRIPTNAME}"
fi
cd ${EIC_DIR}
#Write to Script
echo "#!/bin/bash" > ${GENSCRIPTNAME}
echo -en "\n" >> ${GENSCRIPTNAME}
echo "source ${EIC_DIR}/thisepic.sh" >> ${GENSCRIPTNAME}
echo "source ${EIC_DIR}/local/bin/eicrecon-this.sh" >> ${GENSCRIPTNAME}
echo "/usr/bin/which eicrecon" >> ${GENSCRIPTNAME}
echo "bash  ${EIC_DIR}/run_sim_epic_nolambda.sh -part \"${PARTICLE}\" -n ${num_events} -thmin ${theta_min} -thmax ${theta_max} -pmin ${ENERGY_MIN} -pmax ${ENERGY_MAX} -dist \"${DISTRIBUTION}\" -jid ${NAME_TAG}"  >> ${GENSCRIPTNAME}
chmod 700 ${GENSCRIPTNAME}
bash ${EIC_DIR}/eic-shell -- ./${GENSCRIPTNAME}
echo "bash  ${EIC_DIR}/run_sim_epic_nolambda.sh -part \"${PARTICLE}\" -n ${num_events} -thmin ${theta_min} -thmax ${theta_max} -pmin ${ENERGY_MIN} -pmax ${ENERGY_MAX} -dist \"${DISTRIBUTION}\" -jid ${NAME_TAG}"
 
#=========================================
# Move Output Files to Correct Directories
#=========================================
#Make sure that similar file does not exist to avoid complications
printf "\n Moving files to ${OUTPUT_DIR} \n"
info_string="${PARTICLE}_${DISTRIBUTION}_${ENERGY_MIN}GeV-${ENERGY_MAX}GeV_theta_${theta_min}mrad-${theta_max}mrad_${NAME_TAG}"
mv ${GENSCRIPTNAME} ${TEMPSCRIPT_DIR}
mv ${EIC_DIR}/reco_${info_string}.edm4hep.root ${RECO_DIR}
mv ${EIC_DIR}/sim_${info_string}*.root ${GEN_DIR}
mv ${EIC_DIR}/gen_${info_string}*.hepmc ${HEPMC_DIR}

#Make the stored files read/writeable
chmod -R 777 ${RECO_DIR}
chmod -R 777 ${GEN_DIR}
chmod -R 777 ${HEPMC_DIR}
