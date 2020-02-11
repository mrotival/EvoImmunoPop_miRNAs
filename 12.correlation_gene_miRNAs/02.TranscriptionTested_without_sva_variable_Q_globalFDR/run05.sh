#!/bin/zsh
SCRIPTDIR="/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/12.correlation_gene_miRNAs/TranscriptionTested_without_sva_variable_Q_globalFDR/"
LOGDIR="${SCRIPTDIR}/Logs/"

for condition in 1 2 3 4 5
do
  {
    echo '#!/bin/zsh'
    echo "#SBATCH -o ${LOGDIR}/CARscore_TRtest_${condition}.log"
    echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
    echo "module load R/3.5.0"
    echo "Rscript ${SCRIPTDIR}/05.compute_variance_of_gene_explained_by_miRNAs_CARSCORE.R $condition"
  }>carscore.sh
  chmod +x carscore.sh
  sbatch -p geh --qos=geh -J CAR_TRtest carscore.sh
  rm carscore.sh
done
