COHORT=cohort_name
DIR=path/to/imputed/files
bcftools concat $DIR/chr1.dose.vcf.gz $DIR2/chr2.dose.vcf.gz $DIR/chr3.dose.vcf.gz $DIR/chr4.dose.vcf.gz $DIR/chr5.dose.vcf.gz $DIR/chr6.dose.vcf.gz $DIR/chr7.dose.vcf.gz $DIR/chr8.dose.vcf.gz $DIR/chr9.dose.vcf.gz $DIR/chr10.dose.vcf.gz $DIR/chr11.dose.vcf.gz $DIR/chr12.dose.vcf.gz $DIR/chr13.dose.vcf.gz $DIR/chr14.dose.vcf.gz $DIR/chr15.dose.vcf.gz $DIR/chr16.dose.vcf.gz $DIR/chr17.dose.vcf.gz $DIR/chr18.dose.vcf.gz $DIR/chr19.dose.vcf.gz $DIR/chr20.dose.vcf.gz $DIR/chr21.dose.vcf.gz $DIR/chr22.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools norm -Ou -m -any |
bcftools norm -Ou -f human_g1k_v37.fasta |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o $COHORT.HRC.allchromosomes.R2_0.3.vcf.gz
