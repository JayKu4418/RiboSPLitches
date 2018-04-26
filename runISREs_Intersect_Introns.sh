chroms='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'

for chrom in $chroms
 do
  echo $chrom
  intersectBed -a ../processed_data/Introns_JustmRNAs.bed -b ../temp/Chromosomedata_ISREs/${chrom}_ISRE_locations.bed -wa -wb > ../temp/Chromosomedata_ISREs/Introns_Intersect_${chrom}_ISREs.txt
done
