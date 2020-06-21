

wDir=`pwd`			# working directory
scriptDir="${wDir}/scripts/"	# directory with custom scripts
bamDir=`pwd`			# directory with mapped bam files

TISSUE="Liver BM Cortex"
CONDITIONS="Liver_WT Liver_ADAR1 BM_WT BM_ADAR1 Cortex_WT1 Cortex_ADAR1 Cortex_WT2 Cortex_ADAR2"
ESbed6="${wDir}/ES.bed6"	# path to bed6 files with all editingsites considered
refGenomFa="${wDir}/mm10.fa";

mkdir -p ${wDir}/SS
cd ${wDir}/SS
for i in $(ls /${bamDir}/*.bam);do
  bn=$(basename $i .bam); 
  d=$(date);
  echo -e "start processing\t$bn\t$d";
  samtools view $i | perl ${scriptDir}/getSS_from_bam.pl > SSraw_${bn}.csv;
done

# remove SS with not more than 10 suporting reads in at least 2 samples (out of the 3 samples) for each condition/tissue and being suported by more than 0 reads in all samples of each condition/tissue
for cond in (CONDITIONS)
do
  inputFiles=`ls SSraw_*${cond}*.csv | perl -e 'while(<>){chomp;push @F, $_;}print join("::",@F)."\n"'`
  perl ${scriptDir}/filter_SS_by_readsupport.pl 10 2 $inputFiles  | sort -k1,1 -k2,2n > SSfiltered_${cond}.bed3
done
cat SSfiltered_*.bed3 | sort -k1,1 -k2,2n | uniq > SSfiltered.bed3

# find editing site ES which are closer than 50 nt away from an splice junction SJ.
cat ${ESbed6} | perl -lane 'print join("\t", $F[0], $F[1], $F[2])' | sort -k1,1 -k2,2n > editingSites.bed
cat ${ESbed6} | sort -k1,1 -k2,2n > editingSites.bed6

# get splice site and editins site paires of less than 50 nt distance
bedtools closest -t all -k 15 -d -a SSfiltered.bed3 -b editingSites.bed | perl -F"\t" -lane 'BEGIN{open SS, "> SS_closeTo_ES.bed3"; open ES, "> ES_closeTo_SS.bed3";open BB, "> SS_ES_pairs.csv";} if (abs($F[-1]) <= 50 && $F[-2] != -1){print SS join("\t", @F[0..2]);print ES join("\t", @F[3..5]); print BB join("\t", @F);}'
cat SS_closeTo_ES.bed3 | sort -k1,1 -k2,2n | uniq > tmp; mv -f tmp SS_closeTo_ES.bed3
cat ES_closeTo_SS.bed3 | sort -k1,1 -k2,2n | uniq > tmp; mv -f tmp ES_closeTo_SS.bed3

# extract all reads from all bam files which overlap at least one ES and one SS  = read of interests ROI
mkdir -p ${wDir}/bam
cd ${wDir}/bam
for i in $(ls ${bamDir}/*.bam)
do
  bn=$(basename $i .bam)
  d=$(date)
  echo -e "processing $bn\t$d"
  bedtools intersect -abam $i -b ${wDir}/SS/SS_closeTo_ES.bed3 > SStmp_${bn}.bam
  bedtools intersect -abam SStmp_${bn}.bam -b ${wDir}/SS/ES_closeTo_SS.bed3 > ${bn}_ROI.bam
  rm -f SStmp_${bn}.bam
done

# split ROI.bam into split and read-through reads
mkdir -p ${wDir}/bam/split_subsets/
for i in $(ls *_ROI.bam)
do
  echo $i
  perl ${scriptDir}/seperate_splices_unspliced_reads_for_SS_ES_paires.pl ${wDir}/SS/SS_ES_pairs.csv $i
done

# use bam-readcount to count editing levels in both subsets
cd ${wDir}/bam/split_subsets/
mkdir -p ${wDir}/bamreadcountData/
ln -s ${$refGenomFa}
ln -s ${$refGenomFa}.fai

for i in `find . -iname *.bam`
do
  echo $i
  samtools index $i
done

for i in `find . -iname *.bam`
do
  bn=$(basename $i .bam)
  dir=$(dirname $i)
  echo $bn | perl -F"_" -lane '$F[2]=~s/ES//;print join("\t", $F[0], $F[2], $F[2])' > loi.txt
  zeit=$(date)
  echo -e "$dir\t$bn\t@ $zeit"
  mkdir -p ../bamreadcountData/${dir}
  bam-readcount --site-list loi.txt -q 20 -b 20 -w 1 -f ${$refGenomFa} ${i} >> ${wDir}/bamreadcountData/${dir}/${bn}_bamreadcount.csv 2>/dev/null
  rm -f loi.txt >/dev/null
done

for i in `find ${wDir}/bamreadcountData/ -iname *_bamreadcount.csv`
do
  gzip $i
done


# make Fisher's exact test on 2x2 table( A. split - no editing; B. split - editing; C. no split - no editing; D. no split - editing)
mkdir -p ${wDir}/table
cd ${wDir}/table

for i in $(ls ${wDir}/bam/bamreadcountData/) 
do
  d=$(date)
  echo -e "processing $i  @ $d"
  perl ${scriptDir}/make2x2table_fromDir_andTest_3.pl $i > ${i}.csv 2> ${i}.log
done

# merge p-values for each tissue separately
for tissue in (TISSUE)
do
  inputFiles=`ls ${tissue}*_ROI.csv | perl -e 'while(<>){chomp;push @F, $_;}print join("::",@F)."\n"'`
  perl ${scriptDir}/FisherMerge_sampleTest.pl $inputFiles  > differentially_editing_splicing_Patterns_${tissue}.bed 2> differentially_editing_splicing_Patterns_${tissue}.csv
done

# concatenate p-value merges for each tissue
cat differentially_editing_splicing_Patterns_*.bed > differentially_editing_splicing_Patterns_concatenated.bed


