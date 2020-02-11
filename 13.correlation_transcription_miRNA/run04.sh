#!/bin/zsh

for condition in 1 2 3 4 5
do
  {
    echo "#!/bin/zsh"
    echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
    echo "module load R/3.5.0"
    echo "Rscript /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/13.correlation_transcription_miRNA/04.compute_variance_of_gene_transcription_explained_by_miRNAs_CARSCORE.R $condition"
  }>carscore.sh
  chmod +x carscore.sh
  sbatch carscore.sh
  rm carscore.sh
done
