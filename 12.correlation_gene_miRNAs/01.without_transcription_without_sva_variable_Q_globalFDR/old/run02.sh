#!/bin/zsh

for condition in 1 2 3 4 5
do
  {
    echo "#!/bin/zsh"
    echo "#SBATCH -o glmnet_condition_${condition}_%a"
    echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
    echo "module load R/3.5.0"
    echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/12.correlation_gene_miRNAs/without_transcription_without_sva/02.correlation_allmiRNAs_geneFPKM_withoutTranscrtiption_withoutSVA.R $condition"' $SLURM_ARRAY_TASK_ID'
  }>tempScript.sh
  chmod +x tempScript.sh
  sbatch --mem=20G --array=0-19 tempScript.sh
  rm tempScript.sh
done
