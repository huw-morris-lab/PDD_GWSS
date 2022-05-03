#Written by Maryam Shoai (UCL)
#Edit END with number of raw files generated

END=
for n in $(seq 0 $END); do
	echo "R CMD BATCH --no-save Dementia_${COHORT}_$n.r"
done > 3_submitcxph_allscripts.sh
