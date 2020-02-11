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
      for i in `seq 0 99`
      do
        lastLine=`tail -n 1 ${LOGDIR}/glmnet_${condition}_${perm}_${i}.log`
        [[ $lastLine != '[1] "done"' ]] && echo "$perm $condition $i $lastLine"
        [[ $lastLine != '[1] "done"' ]] && sbatch  -J noTR_miRGene -p geh --qos geh --mem=40G --array=$i tempScript.sh
        done
    rm tempScript.sh
    done
done


# cd /Volumes/evo_immuno_pop/Maxime/miRNA_V2/scripts/12.correlation_gene_miRNAs/without_sva_adjustment_variable_Q_globalFDR/Logs
# for i in `ls` 
# do
#   lastLine=`tail -n 1 $i`
#   [[ $lastLine != '[1] "done"' ]] && echo "$i $lastLine"
# done
# 
# 
# 
# rht='[1] "done"'
# 
# [ $lastLine = $rht ] && echo yes
# [ $lastLine == $rht ] && echo yes
# [[ $lastLine = $rht ]] && echo yes
# [[ $lastLine == $rht ]] && echo yes
# [[ $lastLine == "$rht" ]] && echo yes

        echo "#!/bin/zsh"
