#!/bin/sh

# for FEATURE in Gene miR miRGene
#     do
#         for CHR in `seq 1 22` ;
#         do
#             echo "$FEATURE $CHR" ;
#             sbatch --array=0-1 -p dedicated --qos=fast --mem-per-cpu=40000 -J ${FEATURE}QTL -o "/pasteur/homes/mrotival/JobOutput2/JobOutput_${FEATURE}QTL%a_${CHR}.log" /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/cis_xQTL_tars_runner.sh $CHR ${FEATURE}
#         done ;
#     done

# for FEATURE in Gene miR miRGene
#     do
#         for CHR in `seq 1 22` ;
#         do
#             echo "$FEATURE $CHR" ;
#             sbatch --array=2-100 -p dedicated --qos=fast --mem-per-cpu=40000 -J ${FEATURE}QTL -o "/pasteur/homes/mrotival/JobOutput2/JobOutput_${FEATURE}QTL%a_${CHR}.log" /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/cis_xQTL_tars_runner.sh $CHR ${FEATURE}
#         done ;
#     done

source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module load gcc
module load R/3.5.0

# PERMDIR="/pasteur/scratch/mrotival/sQTL"
# mkdir $PERMDIR
LOGDIR="/pasteur/homes/mrotival/JobOutput2/"

PERM=${SLURM_ARRAY_TASK_ID}
CHR=$1
FEATURE=$2

FILENAME=${LOGDIR}/scripts/run_${FEATURE}QTLs_chr${CHR}_perm${PERM}_ALL_allCond.R
echo "CHR=${CHR}" > ${FILENAME}
echo "perm=${PERM}" >> ${FILENAME}
cat /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/Load${FEATURE}Feature.R >> ${FILENAME}
cat /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/15_eQTL_comparison/eQTL_mapping/QTL_mapping_ALL-cond-chr_perm_tars.R >> ${FILENAME}

Rscript ${FILENAME} ||exit 1
echo "R finished Running on Node "
srun hostname || exit 2
exit 0

