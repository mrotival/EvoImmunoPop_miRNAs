#!/bin/sh

# for FEATURE in Gene miR miRGene
#     do
#       sbatch -p geh --qos=geh --mem-per-cpu=100000 -J ${FEATURE}QTL -o "/pasteur/homes/mrotival/JobOutput2/JobOutput_${FEATURE}QTL_FDR.log" /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/cis_xQTL_FDR_tars_runner.sh ${FEATURE}
#     done

source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module load gcc
module load R/3.5.0

# PERMDIR="/pasteur/scratch/mrotival/sQTL"
# mkdir $PERMDIR
LOGDIR="/pasteur/homes/mrotival/JobOutput2/"

FEATURE=$1

FILENAME=${LOGDIR}/scripts/run_FDR_${FEATURE}QTLs.R
echo "QTL_type='${FEATURE}'" > ${FILENAME}
head -n 2  /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/Load${FEATURE}Feature.R >> ${FILENAME}
cat /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/FDR_compute.R >> ${FILENAME}
Rscript ${FILENAME} ||exit 1
echo "R finished Running on Node "
srun hostname || exit 2
exit 0

