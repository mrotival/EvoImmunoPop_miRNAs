#!/bin/zsh

for condition in {1..5}
do
  for chromosome in {1..22}
  do
    {
      echo "#!/bin/zsh"
      echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
      echo "module load R/3.1.2"
      echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/10.isomiRQTL/02.compute_isomiRQTL.R $chromosome $condition "'$SLURM_ARRAY_TASK_ID'
    } > eqtl.sh
    chmod +x eqtl.sh
    sbatch --qos=fast -p dedicated --array=0-1000 eqtl.sh
  done
done
