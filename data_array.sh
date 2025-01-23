#!/bin/bash
#prepare:target data: vcf(correct header;sorted)  base data:gwas(localization)  .cov(in advance)  .pheno(in advance)
source /home/bioinfo/software/environment.vne
#without imputation:  bash /home/wangfan819/PGD/BRI002/gdmPRS/PRSbeta.sh /home/wangfan819/PGD/BRI002/GDMmerge.vcf GDM 4041 0.005  0.1  0.1  1e-6 10 0.25 0.05 TRUE TRUE
data=$1 #vcf data
prefix=$2 #prefix of target data
pheno=$3 #prefix of base data
maf=$4
geno=$5
mind=$6
hwe=$7
pcas=$8 #pca number
r2=$9 #LD r2 for Heterozygosity,Relatedness,PCA later
RS=${10} #R2 threshold for clumping
binary=${11} #whether binary traits, True/False
manh=${12} #whether plot predictors, True/False
blank=`printf "\n"`
info=0.8 #imputation r2 cutoff
pre="$prefix.$info"
seed=`echo $RANDOM`
echo "random seed  $seed" > input.record
time1=`date +"%Y-%m-%d %H:%M:%S"`
echo "Start at: $time1"
startTime=`date +%s`
echo "start time: $time1" >> input.record
echo "base data name: $pheno" >> input.record
echo "target data name: $prefix" >> input.record
echo "minor allele frequency cutoff: > $maf" >> input.record
genos=`echo 100-$geno/0.01|bc`
echo "genotyping rate: > $genos%" >> input.record
minds=`echo $mind/0.01|bc`
echo "sample missingness: < $minds%" >> input.record
echo "Hardy-Weinberg Equilibrium Pvalue: > $hwe" >> input.record
echo "principal component number:  $pcas" >> input.record
echo "filter out highly correlated SNPs with LD r2: > $r2" >> input.record
echo "R2 threshold for clumping: $RS" >> input.record
echo "$blank" >> input.record
sleep 3
function rand(){   
    low=$1   
    up=$(($2-$low+1))   
    num=$(($RANDOM+1000000000))
    echo $(($num%$up+$low))   
}
rd=$(rand 1000 10000)
line=`echo $rd`
init=`echo ${data##*/}`
initN=`echo ${init%.vcf*}`
num=`head -n $line $data | sed -n '/#CHROM/='`
sed -n "1,${num}p" $data > ${initN}head.txt
sed "1,${num}d" $data > ${initN}main.vcf
sed 's/chr//g' ${initN}main.vcf | awk '{print $1,$2,$4,$5}' > ${initN}main1.vcf
awk '{print $1}' ${initN}main1.vcf > ${initN}chr.txt
awk '{print $2}' ${initN}main1.vcf > ${initN}pos.txt
awk '{print $3}' ${initN}main1.vcf > ${initN}ref.txt
awk '{print $4}' ${initN}main1.vcf > ${initN}alt.txt
paste -d ":" ${initN}chr.txt ${initN}pos.txt ${initN}ref.txt ${initN}alt.txt > ${initN}SNP.txt
awk '{print $1}' ${initN}SNP.txt | paste - ${initN}main.vcf | awk '{$(3+1) = $1; print $0}' | sed 's/[^ ]* //' > ${initN}main2.vcf
sed 's/ /\t/g' ${initN}main2.vcf  > ${initN}main3.vcf
sort -n -k 2,3 -u ${initN}main3.vcf > ${initN}main4.vcf
cat ${initN}head.txt ${initN}main4.vcf > ${initN}NEW.vcf
bcftools sort ${initN}NEW.vcf -o ${initN}sortNEW.vcf
rm {${initN}chr.txt,${initN}pos.txt,${initN}ref.txt,${initN}alt.txt,${initN}SNP.txt,${initN}head.txt,${initN}main.vcf,${initN}main1.vcf,${initN}main2.vcf,${initN}main3.vcf,${initN}main4.vcf,${initN}NEW.vcf}
plink --keep-allele-order --vcf ${initN}sortNEW.vcf --double-id --make-bed --out $prefix
echo "standard vcf for PRS is getted !"

file1=$prefix".pheno"
file2=$prefix".cov"
if [ -f "$file1" -a -f "$file2" ]
then
#Standard QC
plink \
    --keep-allele-order \
    --bfile $prefix \
    --maf $maf \
    --hwe $hwe \
    --geno $geno \
    --mind $mind \
    --write-snplist \
    --make-just-fam \
    --out ${pre}.QC
echo "STQC done"
else
  echo "please provide ${prefix}.pheno or ${prefix}.cov file in advance!" >> input.record
  exit
fi
#high corrlated
plink --keep-allele-order --bfile $prefix --keep ${pre}.QC.fam  --extract ${pre}.QC.snplist --indep-pairwise 200 50 $r2 --out ${pre}.QC
echo "high corrlated done"
#Heterozygosity
plink --keep-allele-order --bfile $prefix --extract ${pre}.QC.prune.in --keep ${pre}.QC.fam --het --out ${pre}.QC
Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -c 1
echo "Heterozygosity check done"

#mismatch
Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -o $pheno -c 2 
echo "mismatch done"
#sex check
#plink --keep-allele-order --bfile $prefix --extract ${pre}.QC.prune.in --keep ${pre}.valid.sample --check-sex --out ${pre}.QC
#Rscript PRSall.R -p $prefix -i $info -c 3

sed '1d' ${pre}.valid.sample > ${pre}.QC.valid  # perform if no sex check
#Relatedness
plink --keep-allele-order --bfile $prefix --extract ${pre}.QC.prune.in --keep ${pre}.QC.valid --rel-cutoff 0.125 --out ${pre}.QC
echo "Relatedness done"
#final
plink \
   --keep-allele-order \
   --bfile $prefix \
   --make-bed \
   --keep ${pre}.QC.rel.id \
   --out ${pre}.QC \
   --extract ${pre}.QC.snplist \
   --exclude ${pre}.mismatch \
   --a1-allele ${pre}.a1
echo "final QC done"
qcN1=`cat ${pre}.QC.fam|wc -l`
echo "Samples number after QC: $qcN1" >> input.record
qcN2=`cat ${pre}.QC.bim|wc -l`
echo "SNPs number after QC: $qcN2" >> input.record

##########################################
#population, get .eigenvec and .eigenval
plink --keep-allele-order --bfile ${pre}.QC  --indep-pairwise 200 50 $r2  --out $pre  #r2
plink --keep-allele-order --bfile ${pre}.QC  --extract ${pre}.prune.in  --pca $pcas --out $pre   #PCA num variable
echo "population and PCA done"


if [ $binary == T -o $binary == TRUE -o $binary == true -o $binary == True ] 
then
  echo "You select: Binary Traits" >> input.record
  plink \
    --keep-allele-order \
    --bfile ${pre}.QC \
    --clump-p1 1 \
    --clump-r2 $RS \
    --clump-kb 250 \
    --clump ${pheno}.QC \
    --clump-snp-field SNP \
    --clump-field P \
    --out $pre
    echo "plink clumping done"
    awk 'NR!=1{print $3}' ${pre}.clumped >  ${pre}.valid.snp
    echo "get index SNPs done"
    awk '{print $3,$8}' ${pheno}.QC > SNP.pvalue
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -c 5
  plink \
    --keep-allele-order \
    --bfile ${pre}.QC \
    --score ${pheno}.QC 3 5 9 header \
    --q-score-range range_list SNP.pvalue \
    --extract ${pre}.valid.snp \
    --out $pre
    echo "plink PRS done"
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -n $pcas -c 6
    echo "plink best cutoff done"
    rm ${pre}.*.profile
else
    echo "You select: Quantitative Traits" >> input.record
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -o $pheno -c 4  #OR to BETA
  plink \
    --keep-allele-order \
    --bfile ${pre}.QC \
    --clump-p1 1 \
    --clump-r2 $RS \
    --clump-kb 250 \
    --clump ${pheno}.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out $pre
    echo "plink clumping done"
    awk 'NR!=1{print $3}' ${pre}.clumped >  ${pre}.valid.snp
    echo "get index SNPs done"
    awk '{print $3,$8}' ${pheno}.QC.Transformed > SNP.pvalue  #SNPid and pvalue
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -c 5
  plink \
    --keep-allele-order \
    --bfile ${pre}.QC \
    --score ${pheno}.QC.Transformed 3 5 9 header \
    --q-score-range range_list SNP.pvalue \
    --extract ${pre}.valid.snp \
    --out $pre
    #col 3(SNPid) col 5(effect Allel) col 9(BETA) .Transformed has a header
    echo "plink PRS done"
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -n $pcas  -c 7
    echo "plink best cutoff done"
    rm ${pre}.*.profile
fi

plinkN=`cat ${pre}.clumped|wc -l`
echo "predictors number of plink: $plinkN" >> input.record


##########################################
#download PRSice
#wget https://github.com/choishingwan/PRSice/releases/download/2.3.3/PRSice_linux.zip
#unzip PRSice_linux.zip

#get covariate: cov(sex) + pca
Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -n $pcas -c 8
echo "get .covariate done"


#calculate PRSscore
if [ $binary == T -o $binary == TRUE -o $binary == true -o $binary == True ]
then
  Rscript /home/wangfan819/software/PRSice.R \
    --prsice /home/wangfan819/software/PRSice_linux \
    --base ${pheno}.QC.gz \
    --target ${pre}.QC \
    --binary-target T \
    --ld LDref \
    --ld-type bed \
    --pheno ${prefix}.pheno \
    --cov ${pre}.covariate \
    --base-maf MAF:${maf} \
    --stat OR \
    --or \
    --clump-p  1 \
    --clump-kb 250 \
    --clump-r2 $RS \
    --interval 0.00005 \
    --all-score \
    --print-snp \
    --quantile 100 \
    --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
    --quant-ref 60 \
    --out $pre
    echo "PRSice PRS done"
    #Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -c 9
    Rscript /home/wangfan819/PGD/BRI002/gdmPRS/PRSall.R -p $prefix -i $info -n $pcas -c 10
    echo "PRSice best cutoff done"
else
   Rscript /home/wangfan819/software/PRSice.R \
    --prsice /home/wangfan819/software/PRSice_linux \
    --base ${pheno}.QC.gz \
    --target ${pre}.QC \
    --ld LDref \
    --ld-type bed \
    --pheno ${prefix}.pheno \
    --cov ${pre}.covariate \
    --base-maf MAF:${maf} \
    --stat BETA \
    --beta \
    --clump-p  1 \
    --clump-kb 250 \
    --clump-r2 $RS \
    --interval 0.00005 \
    --all-score \
    --print-snp \
    --quantile 100 \
    --quant-break 1,5,10,20,40,60,80,90,95,99,100 \
    --quant-ref 60 \
    --out $pre
    echo "PRSice PRS done"
fi

prsiceN=`cat ${pre}.snp|wc -l`
((prsiceN-=1))
echo "predictors number of PRSice: $prsiceN" >> input.record
echo "$blank" >> input.record

if [ $manh == T -o $manh == TRUE -o $manh == true -o $manh == True ]
then
  awk '{print $3,$1,$2,$8}' ${pheno}.QC > base.plot
  Rscript PRSall.R -p $prefix -i $info -c 11
  rm base.plot
  echo "plot after clumping done"
fi

endTime=`date +%s`
time2=`date +"%Y-%m-%d %H:%M:%S"`
echo "Finish at: $time2"

echo "ALL Score by plink in: plinkprs.all_score" >> input.record
echo "ALL Score by PRSice in: ${pre}.all_score " >> input.record
echo "People risk score in best plink model in: plinkprs.best" >> input.record
echo "People risk score in best prsice model in: prsice2PRS.best" >> input.record
echo "Pictures were saved in: .pdf" >> input.record

echo "finish time: $time2" >> input.record
sumTime=$(( $endTime - $startTime ))
echo "Total run $sumTime seconds"

