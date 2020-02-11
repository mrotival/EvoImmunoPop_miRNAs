#!/bin/zsh

SCRIPTDIR="/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/12.correlation_gene_miRNAs/without_transcription_without_sva_variable_Q_globalFDR/"
LOGDIR="${SCRIPTDIR}/Logs/"


for perm in 0 1
do
	for condition in 1 2 3 4 5
	do
	  {
		echo '#!/bin/zsh'
		echo "#SBATCH -o ${LOGDIR}/glmnet_${condition}_${perm}_%a.log"
		echo "source /local/gensoft2/adm/etc/profile.d/modules.sh"
		echo "module load R/3.5.0"
		echo "Rscript ${SCRIPTDIR}/02.correlation_allmiRNAs_geneFPKM_withoutTranscription_withoutSVA.R $condition $perm"' $SLURM_ARRAY_TASK_ID'
	  }>tempScript.sh
	  chmod +x tempScript.sh
	  sbatch -J noTR_miRGene --mem=20G --array=0-99 tempScript.sh 
	  rm tempScript.sh
	done
done