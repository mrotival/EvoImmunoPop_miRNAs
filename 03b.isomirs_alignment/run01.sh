#!/bin/zsh

for condition in 1 2 3 4 5
do
  {
    echo "#!/bin/zsh"
    echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
    echo "module load R/3.5.0"
    echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/20.Gene_regulation_by_miRNA/scripts/01.correlation_miRNA_geneFPKM_transcriptionadjusted.R $condition "'$SLURM_ARRAY_TASK_ID'
  }>tempScript.sh
  chmod +x tempScript.sh
  sbatch --array=1-658 tempScript.sh
  rm tempScript.sh
done
