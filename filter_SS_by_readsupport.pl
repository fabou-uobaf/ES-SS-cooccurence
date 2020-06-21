#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $readCount_th    = shift;
my $sampleCount_th  = shift;
my $fhs             = shift;
my @fhs             = split"::", $fhs;

my %data = ();

foreach my $fh (@fhs){
  open I, "< $fh" or die "can t open $fh\n";
  while(<I>){
    chomp;
    my @F=split"\t",$_;
    $data{$F[0]}->{$F[1]}->{$fh} = $F[2];
  }
  close I;
}

foreach my $chr (keys %data){
  foreach my $pos (keys %{$data{$chr}}){
    foreach my $fh (@fhs){
      $data{$chr}->{$pos}->{$fh} = 0 unless ( defined($data{$chr}->{$pos}->{$fh}) )
    }
  }
}

foreach my $chr (keys %data){
  foreach my $pos (keys %{$data{$chr}}){
    my $check = 0;
    foreach my $fh (keys %{$data{$chr}->{$pos}}){
      if ($data{$chr}->{$pos}->{$fh} >= $readCount_th){
        $check++;
      }
      elsif ($data{$chr}->{$pos}->{$fh} <= 0){
        $check--;
      }
    }
    if ($check >= $sampleCount_th){
      print join("\t", $chr, $pos-1, $pos)."\n";
    }
  }
}


#print Dumper(%data);
