#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

my $roi_fh = shift;
my $bam_fh = shift;
my $basename = fileparse($bam_fh,qw/ .bam/);

my $samtools = "samtools";
my $bedtools = "bedtools";
my $tmpESbed_fh           = "ES.bed";
my $tmpSSbed_fh           = "SS.bed";
my $tpmIntermediateBAM_fh = "intermediate.bam";
my $tpmIntersectedBAM_fh  = "intersected.bam";
my $resDir_fh             = "split_subsets/${basename}/";
system("mkdir -p $resDir_fh");

my %data     = ();

open R, "< $roi_fh" or die " can t open $roi_fh\n";
while(<R>){
  chomp;
  my @F  = split"\t", $_;
  my $id = join("_", $F[0], "SS$F[2]", "ES$F[5]");
  $data{$id}->{SS}  = $F[2];
  $data{$id}->{ES}  = $F[5];
  $data{$id}->{chr} = $F[0];
  print STDERR "Warnings:\tES and SS not on same chr: ($F[0] != $F[3])\n" if ($F[0] ne $F[3]);
}

foreach my $id (keys %data){
  # print SS.bed
  open SS, "> $tmpSSbed_fh" or die "can t write to $tmpSSbed_fh\n";
  print SS join("\t", $data{$id}->{chr}, $data{$id}->{SS}-1, $data{$id}->{SS})."\n";
  close SS;
  # print ES.bed
  open ES, "> $tmpESbed_fh" or die "can t write to $tmpESbed_fh\n";
  print ES join("\t", $data{$id}->{chr}, $data{$id}->{ES}-1, $data{$id}->{ES})."\n";
  close ES;

  # make sequential bedtools intersect
  # remove reads which do not have a match at ES (-split option in second intersect call
  system("$bedtools intersect -abam $bam_fh -b $tmpSSbed_fh > $tpmIntermediateBAM_fh");
  system("$bedtools intersect -split -abam $tpmIntermediateBAM_fh -b $tmpESbed_fh > $tpmIntersectedBAM_fh");

  # read in an split intersected bam file for split at position SS (and match at position ES)
  my $fileName_split   = "${id}_${basename}_split.bam";
  my $fileName_unsplit = "${id}_${basename}_unsplit.bam";
  open BS, " | samtools view -bS - > ${resDir_fh}/${fileName_split}" or die "can t write bam file via <.. | samtools view -bS - > ${resDir_fh}/${fileName_split}>\n";
  open BU, " | samtools view -bS - > ${resDir_fh}/${fileName_unsplit}" or die "can t write bam file via <.. | samtools view -bS - > ${resDir_fh}/${fileName_unsplit}>\n";
  open IB, "samtools view -h $tpmIntersectedBAM_fh | " or die "can t read bam file via <samtools view $tpmIntersectedBAM_fh | ..>\n";
  while(<IB>){
    if (m/^@/){
      print BS "$_";
      print BU "$_";
      next;
    }
    chomp;
    my @F=split"\t",$_;
    my ($chr, $bitflag, $left, $cigar) = ($F[2], $F[1], $F[3], $F[5]);

    my %ss_pos = &getSS($cigar, $left);
    
    if ( defined( $ss_pos{$data{$id}->{SS}}) ){
      print BS "$_\n";
    }
    else{
      print BU "$_\n";
    }
  }

  # remove tmp files
  system("rm -f $tmpESbed_fh $tmpSSbed_fh $tpmIntermediateBAM_fh $tpmIntersectedBAM_fh 2>/dev/null");
}

#print Dumper(%data);

## subroutines

sub getSS {
  my $cigar_string = shift;
  my $cumLength    = shift;
  my %ss_pos       = ();

  while($cigar_string=~m/(\d+)([MDX=N])/g){
    my $dist = $1;
    my $type = $2;
    if ($type eq "N"){
      $ss_pos{$cumLength+1}++;
      $ss_pos{$cumLength+$dist}++;
    }
    $cumLength += $dist;
  }
  return(%ss_pos);
}

