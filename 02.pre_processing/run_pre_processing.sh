#!/bin/zsh
source /local/gensoft2/adm/etc/profile.d/modules.sh
module load fastx_toolkit
source /pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/python_virtual_environment/mirna_python/bin/activate
python3 /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/02.pre_processing/01.trimming_reads.py ${SLURM_ARRAY_TASK_ID}
