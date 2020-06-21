#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %data = ();

while(<>){
  chomp;
  my @F=split"\t",$_;
  my ($chr, $bitflag, $left, $cigar) = ($F[2], $F[1], $F[3], $F[5]);

  next unless ($cigar=~m/\dN/);

  #my $strand = &getStrand($bitflag);
  my @ss_pos = &getSS($cigar, $left);
  foreach my $ss_pos (@ss_pos){
    #$data{$chr}->{$ss_pos}->{$strand}++;
    $data{$chr}->{$ss_pos}++;
  }
}

foreach my $chr (sort keys %data){
  foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}){
    #foreach my $strand (sort keys %{$data{$chr}->{$pos}}){
    #  print join("\t", $chr, $pos, $strand, $data{$chr}->{$pos}->{$strand})."\n";
    #}
    print join("\t", $chr, $pos, $data{$chr}->{$pos})."\n";
  }
}

sub getStrand{
  my $flag = shift;
  my $strand = "+";
  if( ($flag+0) & 16){
    $strand = "-";
  }
  return($strand);
}

sub getSS {
  my $cigar_string = shift;
  my $cumLength    = shift;
  my @ss_pos       = ();

  while($cigar_string=~m/(\d+)([MDX=N])/g){
    my $dist = $1;
    my $type = $2;
    if ($type eq "N"){
      push @ss_pos, $cumLength+1;
      push @ss_pos, $cumLength+$dist;
    }
    $cumLength += $dist;
  }
  return(@ss_pos);
}

