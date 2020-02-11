#!/bin/zsh
source /local/gensoft2/adm/etc/profile.d/modules.sh
module load fastqc/0.11.5
source /pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/python_virtual_environment/mirna_python/bin/activate
python3 /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/01.description_of_raw_data/05.launch_fastqc.py ${SLURM_ARRAY_TASK_ID}
