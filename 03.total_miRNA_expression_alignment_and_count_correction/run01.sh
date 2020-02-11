#!/bin/zsh

######################################################################################################
##WARNING : some of them will crash because of the cluster (10/1000), just relaunch and it will pass##
######################################################################################################

source /local/gensoft2/adm/etc/profile.d/modules.sh
module load samtools/1.9
module load bowtie/1.1.1
module load Python/2.7.8
module load bedtools/2.25.0
source /pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/python_virtual_environment/mirna_python/bin/activate
python3 /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/03.total_miRNA_expression_alignment_and_count_correction/01.alignement_on_hg19.py ${SLURM_ARRAY_TASK_ID}
