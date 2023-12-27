#!/bin/sh

#PBS -l mem=50gb,nodes=1:ppn=4,walltime=24:00:00 
#PBS -m abe 
#PBS -M lileichinaus@gmail.com
#PBS -q small
#PBS -e /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/script/predict_list01
#PBS -o /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/script/predict_list01

module load bad_mutations/1.0

FASTAFILES=($(cat /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_fasta_tree.list01))
CURRENT_FASTAFILES=${FASTAFILES[${PBS_ARRAYID}]}
GENES=($(awk -F"/" '{print $NF}' /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_fasta_tree.list01|awk -F"." '{print $1}'))
CURRENT_GENES=${GENES[${PBS_ARRAYID}]}
MSAFILES=($(cat /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_MSA.list01))
CURRENT_MSAFILES=${MSAFILES[${PBS_ARRAYID}]}
TREEFILES=($(cat /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_tree.list01))
CURRENT_TREEFILES=${TREEFILES[${PBS_ARRAYID}]}
SUBSFILES=($(cat /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/file_list/sorted_sub.list01))
CURRENT_SUBFILES=${SUBSFILES[${PBS_ARRAYID}]}

python /home/morrellp/llei/BAD_Mutations/BAD_Mutations.py  -v DEBUG\
 predict\
 -c /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Config/BAD_Mutations_Config_Ath.txt\
 -f ${CURRENT_FASTAFILES}\
 -a ${CURRENT_MSAFILES}\
 -r ${CURRENT_TREEFILES}\
 -s ${CURRENT_SUBFILES}\
 -o /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list01
 2> /panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/log_predict/list01/${CURRENT_GENES}_Predictions.log
