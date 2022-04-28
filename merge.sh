#Script to merge results files of CPH analysis
#To run: sh merge.sh 

awk 'NR == 1 || FNR > 1'  ${COHORT}_DEMENTIA_cxphres_*.txt > ${COHORT}_DEMENTIA_cxphres_all.txt 
