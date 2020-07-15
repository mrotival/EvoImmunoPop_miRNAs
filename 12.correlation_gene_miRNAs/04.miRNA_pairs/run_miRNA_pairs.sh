#!/bin/zsh
SCRIPTDIR="/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts_to_share/12.correlation_gene_miRNAs/04.miRNA_pairs/"
LOGDIR="${SCRIPTDIR}/Logs/"


for condition in 1 2 3 4 5 
do
  {
    echo '#!/bin/zsh'
    echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
    echo "module load R/3.5.0"
    echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts_to_share/12.correlation_gene_miRNAs/04.miRNA_pairs/Correlate_miRNA_pair_with_gene_expression.R $condition "'$SLURM_ARRAY_TASK_ID'
  }>tempScript.sh
  chmod +x tempScript.sh
  sbatch -o ${LOGDIR}/miR_pair_${condition}_%a.log -J miR_pairs --mem=20G --array=1-390 tempScript.sh
  rm tempScript.sh
done
