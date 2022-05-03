#Written by Maryam Shoai (UCL)

END=
for n in $(seq 0 $END); do

	echo "#!/usr/bin/Rscript
	library(data.table)
	library(survival)
	library(dplyr)

	#Load genetic data
	subset<- fread(\"$DIR/$COHORT/RESULTS/${COHORT}_file${n}.raw\")
	#Keep IID and genetic data columns
	subset<- subset[,-c(3,4,5,6)]

	#Load clinical data
	clinical <- fread(\"$DIR/${COHORT}_clinical.txt\")

	#Load Principal Components
	PCs <- fread(\"DIR/${COHORT}_PCA.eigenvec\")

	#Select first 10 principal components
	PCs <- PCs %>%
		select(V2:V12) %>%
		rename(IID = V2,
			PC1 = V3,
			PC2 = V4,
			PC3 = V5,
			PC4 = V6,
			PC5 = V7,
			PC6 = V8,
			PC7 = V9,
			PC8 = V10,
			PC9 = V11,
			PC10 = V12)
	
	#Inner join all datasets 
	TABLE <- clinical %>%
		inner_join(PCs, by = \"IID\") %>%
		inner_join(subset, by = \"IID\")

	TABLE<- as.data.frame(TABLE)
	
	coefficients<-as.data.frame(matrix(ncol= 9))
 
	names(coefficients) <- c(\"SNP\",\"Coeff\", \"se\", \"Pvalue\", \"Cox.zphPVal\", \"N\", \"ov.lik.ratio\",\"logrank\", \"r2\" )
 
	for (i in 18:ncol(TABLE)) {
	print(colnames(TABLE)[i])
  	snp <- TABLE[,c(i,1:17)]
	j= i-17

#Test Cox model - if there is an error then put NoConverge into results
  	testCox = try(coxph(Surv(snp\$timeToEvent_dementia, snp\$event_dementia) ~ snp[,1] + snp\$age_onset + snp\$gender + snp\$PC1+ snp\$PC2 + snp\$PC3 + snp\$PC4 + snp\$PC5, data=snp), 
  	              silent = T)
  	
  	if(class(testCox)[1]==\"try-error\"){
  	  coefficients[j,1]<- paste(colnames(TABLE)[i])
  	  coefficients[j,2:9] <- \"NoConverge\"
  	} else {
  	
  	#If Cox model is not an error, run and get the coefficients
  	model.cox<- coxph(Surv(snp\$timeToEvent_dementia, snp\$event_dementia) ~ snp[,1] + snp\$age_onset + snp\$gender + snp\$PC1+ snp\$PC2 + snp\$PC3 + snp\$PC4 + snp\$PC5, data=snp)
  	
  	testkmz <- try(kmz<- cox.zph(model.cox, transform = \"km\"), silent = T)
  	
  	coefficients[j,1]<- paste(colnames(TABLE)[i])
  	coefficients[j,2]<- summary(model.cox)\$coefficients[1,1]
  	coefficients[j,3]<- summary(model.cox)\$coefficients[1,3]
  	coefficients[j,4]<- summary(model.cox)\$coefficients[1,5]
  	coefficients[j,6]<- model.cox\$n
  	coefficients[j,7]<- summary(model.cox)\$logtest[[1]]
  	coefficients[j,8]<- summary(model.cox)\$sctest[[1]]
  	coefficients[j,9]<- summary(model.cox)\$rsq[[1]] # nagelkerke r square
  	
  	#However sometimes the kmz transformation results in an error
  	#If it is an error, put solveError into the results table, otherwise get the Cox.zphPVal
  	if(class(testkmz)[1] == \"try-error\"){
  	coefficients[j,5]<- \"solveError\"
  	} else {
  	coefficients[j,5]<- kmz\$table[1,3]
  	}
  	}
  	}

	fwrite(coefficients, \"${COHORT}_DEMENTIA_cxphres_${n}.txt\", row.names=FALSE, sep=\"\t\", quote= FALSE)" > Dementia_${COHORT}_${n}.r
done

