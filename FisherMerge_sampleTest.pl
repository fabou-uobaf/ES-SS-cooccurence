#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $fhs = shift;
my @fhs = split"::",$fhs;
my %data = ();
my %samples = ();

foreach my $fh (@fhs){
  open F, "< $fh" or die "can t open $fh\n";
  while(<F>){
    chomp;
    my @F=split"\t",$_;
    if ($F[0] eq "sample" && $F[4] eq "unedited_split"){
      next;
    }
    my ($sample,$chr,$ESpos, $SSpos,$unedited_split,$edited_split,$unedited_unsplit,$edited_unsplit,$enrichment,$pval,$qval) = @F;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{unedited_split} = $unedited_split;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{edited_split} = $edited_split;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{unedited_unsplit} = $unedited_unsplit;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{edited_unsplit} = $edited_unsplit;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{enrichment} = $enrichment;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval} = $pval;
    $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{qval} = $qval;

    $samples{$sample}++;
  }
  close F;
}

open RD, "> R.data" or die "can t write to R.data\n";
foreach my $chr (keys %data){
  foreach my $ESpos (keys %{$data{$chr}}){
    foreach my $SSpos (keys %{$data{$chr}->{$ESpos}}){
      my @pvalues = ();
      foreach my $sample (sort keys %samples){
        if ( defined($data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval}) ){
          push @pvalues, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval};
        }
        else{
          push @pvalues, 1;
        }
      }
      print RD join("\t", $chr, $ESpos, $SSpos, @pvalues)."\n";
   }
  }
}
close RD;

open RS, "> R.script" or die "can t write to R.script\n";
print RS 'FisherMethodVector <- function(X) {'."\n";
print RS '  line <- as.numeric(X[4:length(X)]);'."\n";
print RS '  line <- replace(line,which(line<.Machine$double.eps),.Machine$double.eps);'."\n";
print RS '  pval <- pchisq(-2*sum(log(line)),df = 2*length(line),lower.tail=FALSE);'."\n";
print RS '  return(pval)'."\n";
print RS '}'."\n";
print RS 'data <- read.table(file="R.data", sep="\t", header = FALSE)'."\n";
print RS 'data$p.value <- apply(data, 1, FisherMethodVector)'."\n";
print RS 'data$q.value <- p.adjust(data$p.value, method = "fdr")'."\n";
print RS 'write.table(as.data.frame(data),file="R.res", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)'."\n";
close RS;

system("R CMD BATCH --no-save --no-restore R.script");

open RR, "< R.res" or die "can t open R.res\n";
while(<RR>){
  chomp;
  my @F=split"\t",$_;
  $data{$F[0]}->{$F[1]}->{$F[2]}->{pval} = $F[-2];
  $data{$F[0]}->{$F[1]}->{$F[2]}->{qval} = $F[-1];
}
close RR;
system("rm -f R.res R.script R.script.Rout R.data");

print STDERR "#".join("\t", qw/chr ESpos SSpos/, map ({"pval_$_"} sort keys %samples), "merged_p.value", "merged_q.value" )."\n";
print STDERR "#".join("\t", "", "->", qw/sampleID unedited_split edited_split unedited_unsplit edited_unsplit enrichment p.value q.value/)."\n"; 
foreach my $chr (sort keys %data){
  foreach my $ESpos (sort {$a <=> $b} keys %{$data{$chr}}){
    foreach my $SSpos (sort {$a <=> $b} keys %{$data{$chr}->{$ESpos}}){
      my @pvalues = ();
      foreach my $sample (sort keys %samples){
        if ( defined($data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval}) ){
          push @pvalues, sprintf "%.3g", $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval};
        }
        else{
          push @pvalues, "NA=1";
        }
      }

      print join("\t", $chr, $ESpos-1, $ESpos, join(":", $chr, $ESpos, $SSpos), sprintf("%.3g", $data{$chr}->{$ESpos}->{$SSpos}->{qval}), "." )."\n";
      print STDERR join("\t", $chr, $ESpos, $SSpos, @pvalues, sprintf("%.3g",$data{$chr}->{$ESpos}->{$SSpos}->{pval}), sprintf("%.3g", $data{$chr}->{$ESpos}->{$SSpos}->{qval}) )."\n"; 

      foreach my $sample (sort keys %samples){
        if ( defined($data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval}) ){
          print STDERR join("\t", "", "->", $sample, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{unedited_split}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{edited_split}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{unedited_unsplit}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{edited_unsplit}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{enrichment}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{pval}, $data{$chr}->{$ESpos}->{$SSpos}->{$sample}->{qval})."\n";
        }
        else{
          print STDERR join("\t", "", "->", $sample, "NA", "NA", "NA", "NA", "NA", "NA", "NA")."\n";
        }
      }
    }
  }
}

#print Dumper(%data);

