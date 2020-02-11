#!/bin/zsh

for condition in {1..5}
do
  for chromosome in {1..22}
  do
    {
      echo "#!/bin/zsh"
      echo "#SBATCH -J mirQTL_"$condition"_"$chromosome"_%a"
      echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
      echo "module load R/3.5.0"
      echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/09.mirQTL/02.compute_mirQTL.R $chromosome $condition "'$SLURM_ARRAY_TASK_ID'
    } > eqtl.sh
    chmod +x eqtl.sh
    sbatch --qos=fast -p dedicated --array=0-1000 eqtl.sh
  done
done
