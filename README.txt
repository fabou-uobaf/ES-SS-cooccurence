# editing splicing co-occurence

## Synopsis

Input a list of editing sites in bed6 format and a list of NGS reads mapped with STAR to the reference genome. From each bam file splice sites are extracted and filtered if they have a editing site in close proximity (+/- 50 nt). Reads overlapping a SS-ES pair are extracted and split in two bam files holding spliced and unspliced reads, respectively. In each of this reads the number of editing events at the annotated editing site is deduced using the tool bam-readcounts. For each sample and each ES-SS pair, independance of ES and SS events are tested using a Fisher Exact Test. Obtainted p-values from each sample is merged in a meta p-value using Fisher's Method.

## Tools pre-installed and available in $PATH

* samtools

* bedtools

* bam-readcounts

* gzip
