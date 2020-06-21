#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dirbase = "../bam/bamreadcountData/";
my $basename = shift;
my $dir = "${dirbase}/${basename}/";

my $rData_fh = "R.data";
my $rResu_fh = "R.res";
my $rScri_fh = "R.script";

my %data = ();

opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @files = grep { !/^\./ && -f "$dir/$_" } readdir($dh);

foreach my $file (@files){
  my $modFile = $file;
  $modFile =~s/$basename/0/;
  my @parts = split"_", $modFile;
  my ($ss, $es, $sample, $split, $sufix) = ($parts[-5], $parts[-4], $parts[-3], $parts[-2], $parts[-1]);
  $sample = $basename;
  my $chr = join("_", @parts[0 .. ($#parts - 5)]);
  $es=~s/ES//;
  $ss=~s/SS//;
  
  my $fh = "${dir}/${file}";
  if (-s $fh){
    if ($file =~ m/.gz$/){
      open FH, "zcat $fh | " or die "can t open $fh via <zcat $fh | ..>\n";
    }
    else{
      open FH, "< $fh" or die "can t open $fh\n";
    }
    while(<FH>){
      chomp;
      my @F=split"\t", $_;
      print STDERR "ERROR:\tError in $file: $_\n" if ($es != $F[1] || $chr ne $F[0]);
      my $refBase = uc($F[2]);
      print STDERR "WARNING:\tRef prob in $file: $_\n" if ($refBase=~!/^[at]$/);
      my $totalCoverage = $F[3];
      my ($at, $gc) = (0,0);
      foreach my $col (4..$#F){
        my @G = split":", $F[$col];
        $at += $G[1] if ( $G[0]=~m/^[at]$/i);
        $gc += $G[1] if ( $G[0]=~m/^[gc]$/i);
      }
      my $acgt = ($at+$gc);
      print STDERR "WARNING:\ttotal cov and at+gc to not match up ($totalCoverage != ($at+$gc) = $acgt) \n#\t$_\n" if ($totalCoverage != ($at+$gc));
      
      $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{refBase} = $refBase ;
      $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{totalCoverage} = $totalCoverage ;
      $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{AT} = $at ;
      $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{GC} = $gc ;
    }
  }
  else{

    $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{totalCoverage} = 0 ;
    $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{AT} = 0 ;
    $data{$basename}->{$chr}->{$es}->{$ss}->{$split}->{GC} = 0 ;

  }
}
closedir $dh;

open RD, "> $rData_fh" or die "can t write to $rData_fh\n";
foreach my $bn (sort keys %data){
  foreach my $chr (sort keys %{$data{$bn}}){
    foreach my $es (sort {$a <=> $b} keys %{$data{$bn}->{$chr}}){
      foreach my $ss (sort {$a <=> $b} keys %{$data{$bn}->{$chr}->{$es}}){
        my $ref_split      = ( defined($data{$bn}->{$chr}->{$es}->{$ss}->{split}) )?($data{$bn}->{$chr}->{$es}->{$ss}->{split}->{AT}):(0);
        my $edited_split   = ( defined($data{$bn}->{$chr}->{$es}->{$ss}->{split}) )?($data{$bn}->{$chr}->{$es}->{$ss}->{split}->{GC}):(0);
        my $ref_unsplit    = ( defined($data{$bn}->{$chr}->{$es}->{$ss}->{unsplit}) )?($data{$bn}->{$chr}->{$es}->{$ss}->{unsplit}->{AT}):(0);
        my $edited_unsplit = ( defined($data{$bn}->{$chr}->{$es}->{$ss}->{unsplit}) )?($data{$bn}->{$chr}->{$es}->{$ss}->{unsplit}->{GC}):(0);

        print RD join("\t", $bn, $chr, $es, $ss, $ref_split, $edited_split, $ref_unsplit, $edited_unsplit)."\n";
      }
    }
  }
}
close RD;

open RS, "> $rScri_fh" or die "can t write to $rScri_fh\n";
print RS '# read in data table ';
print RS "\n";
print RS "filename <- '$rData_fh'  ";
print RS "\n";
print RS 'data <- read.table(file=filename, sep="\t", header=F) ';
print RS "\n";
print RS "\n\n";
print RS '# define functions ';
print RS "\n";
print RS 'FisherVector <- function(X) {fisher.test(matrix(c(as.numeric(X[5]), as.numeric(X[6]), as.numeric(X[7]), as.numeric(X[8])), ncol=2, byrow=F), alternative = "two.sided")$p.value} ';
print RS "\n";
print RS 'Enrichment   <- function(X) {((as.numeric(X[5])+1)/(as.numeric(X[6])+1))/((as.numeric(X[7])+1)/(as.numeric(X[8])+1))} ';
print RS "\n\n";
print RS '# calculate scores ';
print RS "\n";
print RS 'data$enrichmentFactor <- apply(data, 1, Enrichment) ';
print RS "\n";
print RS 'data$p.value <- apply(data, 1, FisherVector) ';
print RS "\n";
print RS 'data$q.value <- p.adjust(data$p.value, method ="fdr") ';
print RS "\n\n";
print RS '# print output ';
print RS "\n";
print RS "resname <- '$rResu_fh' ";
print RS "\n";
print RS 'write.table(as.data.frame(data),file=resname, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE) ';
print RS "\n";
close RS;

system("R CMD BATCH --no-save --no-restore $rScri_fh");

print join("\t", qw/sample chr ES SS unedited_split edited_split unedited_unsplit edited_unsplit enrichment pval qval/)."\n";
open RR, "< $rResu_fh" or die "can t open $rResu_fh\n";
while(<RR>){
  chomp;
  my @F=split"\t",$_;
  my ($bn, $chr, $es, $ss, $ref_split, $edited_split, $ref_unsplit, $edited_unsplit, $oddsratio, $pVal, $qVal) = @F;
  $data{$bn}->{$chr}->{$es}->{$ss}->{oddsratio} = $oddsratio;
  $data{$bn}->{$chr}->{$es}->{$ss}->{pVal} = $pVal;
  $data{$bn}->{$chr}->{$es}->{$ss}->{qVal} = $qVal;
  print join("\t", @F)."\n";
}
close RR;

system("rm -f $rScri_fh $rResu_fh $rData_fh ${rScri_fh}.Rout");
#print Dumper(%data);
