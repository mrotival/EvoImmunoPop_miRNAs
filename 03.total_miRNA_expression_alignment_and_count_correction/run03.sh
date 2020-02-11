#!/bin/zsh

source /pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/python_virtual_environment/mirna_python/bin/activate
python3 /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/03.total_miRNA_expression_alignment_and_count_correction/03.get_count_after_alignement.py ${SLURM_ARRAY_TASK_ID}
