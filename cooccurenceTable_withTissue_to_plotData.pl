#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $qval_th  = shift;
my $bed_fh   = shift;
my $files_fh = shift;
my @files_fh = split":", $files_fh;

my %sigs = ();
my %file2tissue = ();
$file2tissue{'RNAseq_2864-WT-Primer2_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_2864-WT-Primer2_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_2868-Hom-Primer4_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_2868-Hom-Primer4_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_2960-Hom-Primer5_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_2960-Hom-Primer5_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_2961-WT-Primer6_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_2961-WT-Primer6_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_298-WT-Primer2_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_298-WT-Primer2_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_300-ADAR2-Primer4_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_300-ADAR2-Primer4_ROI.csv'}->{condition}="dADAR2";
$file2tissue{'RNAseq_3113-Hom-Primer7_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_3113-Hom-Primer7_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_3118-WT-Primer12_ROI.csv'}->{tissue}="Cortex1";
$file2tissue{'RNAseq_3118-WT-Primer12_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_321-WT-Primer5_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_321-WT-Primer5_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_322-ADAR2-Primer6_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_322-ADAR2-Primer6_ROI.csv'}->{condition}="dADAR2";
$file2tissue{'RNAseq_3829-WT-Primer7_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_3829-WT-Primer7_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_3832-ADAR2-Primer12_ROI.csv'}->{tissue}="Cortex2";
$file2tissue{'RNAseq_3832-ADAR2-Primer12_ROI.csv'}->{condition}="dADAR2";
$file2tissue{'RNAseq_BM_ADAR1_littermate21_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_ADAR1_littermate21_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_BM_ADAR1_littermate22_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_ADAR1_littermate22_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_BM_ADAR1_littermate23_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_ADAR1_littermate23_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_BM_WT_littermate21_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_WT_littermate21_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_BM_WT_littermate22_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_WT_littermate22_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_BM_WT_littermate23_ROI.csv'}->{tissue}="Bonemarrow";
$file2tissue{'RNAseq_BM_WT_littermate23_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_Liver_ADAR1_littermat13_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_ADAR1_littermat13_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_Liver_ADAR1_littermate11_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_ADAR1_littermate11_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_Liver_ADAR1_littermate12_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_ADAR1_littermate12_ROI.csv'}->{condition}="dADAR1";
$file2tissue{'RNAseq_Liver_WT_littermate11_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_WT_littermate11_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_Liver_WT_littermate12_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_WT_littermate12_ROI.csv'}->{condition}="Wildtype";
$file2tissue{'RNAseq_Liver_WT_littermate13_ROI.csv'}->{tissue}="Liver";
$file2tissue{'RNAseq_Liver_WT_littermate13_ROI.csv'}->{condition}="Wildtype";


open B, "< $bed_fh" or die "can t open $bed_fh\n";
while(<B>){
  chomp;
  my @F=split"\t",$_;
  my ($chr, $es, $ss) = split":", $F[3];
  if ($F[-2] <= $qval_th){
    if ( defined($sigs{$chr}->{$es}->{$ss}) && $sigs{$chr}->{$es}->{$ss} > $F[-2]){
      $sigs{$chr}->{$es}->{$ss} = $F[-2]
    }
    else{
      $sigs{$chr}->{$es}->{$ss} = $F[-2];
    }
  }
}
close B;

print join("\t", qw/sample condition tissue es ss pval min_pval total_read_count editingRate splitRate/)."\n";
foreach my $file (@files_fh){
  my $sample = $file;
  $sample=~s/RNAseq_//g;
  $sample=~s/_ROI.csv//g;
  my $condition = "unkown";
  if ( defined($file2tissue{$file}->{condition}) ){
    $condition = $file2tissue{$file}->{condition};
  }
  my $tissue = "unkown";
  if ( defined($file2tissue{$file}->{tissue}) ){
    $tissue = $file2tissue{$file}->{tissue};
  }

  open F, "<$file" or die "can t open $file\n";
  while(<F>){
    chomp;
    my @F=split"\t",$_;
    next if ($F[2] eq "ES");
    next unless ( defined($sigs{$F[1]}->{$F[2]}->{$F[3]}) );
    my $es = "$F[1]:$F[2]";
    my $ss = "$F[1]:$F[3]";
    my $pval = $F[9];
    my $total = $F[4]+$F[5]+$F[6]+$F[7]+1;
    my $editing = $F[5]+$F[7];
    my $split   = $F[4]+$F[5];
    my $editingRate = ($editing)/($total);
    my $splitRate   = ($split)/($total);
    print join("\t", $sample, $condition, $tissue, "ES $es", "SS $ss", $pval, $sigs{$F[1]}->{$F[2]}->{$F[3]}, $total, $editingRate, $splitRate)."\n";
  }
}

