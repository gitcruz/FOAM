#!/usr/bin/env perl
 
#use Pod::Usage; # commented this on 31-01-2023
use Getopt::Long;
use File::Basename;
 
$this_program = basename($0);
my $minmatch = 800;
$minqual = 12;
my $target = 0;
my $result = GetOptions(
			'b|buffer=i' => \$printlines,
			'm|minq=i'   => \$minqual,
			'u|uniq=s'   => \$uniq,
			'ref:s' => \$target
		       );
die "Must supply reference for filtering using -ref ref_id\n" if !$target;
my %refseqlen;
while (<>) {
  if (m/SQ.*SN:(\S+).*LN:(\d+)/) {
    $refseqlen{$1}=$2;
  }
  next if ($_=~/^\@/);
  $t++;
  my @samfields = split /\t/;
  my $strand = strand($samfields[1]);
  my $supplemental = supplemental($samfields[1]);
  next if $supplemental;
  my $contig = $samfields[2];
  my $fastqstring = '@'.$samfields[0]."\n".$samfields[9]."\n+\n".$samfields[10]."\n";
  if ($contig eq $target) {
    my $read   = $samfields[0];
    my $from   = $samfields[3];
    my $mapq = $samfields[4];
    next if $mapq < $minqual;
    my $refalnlen    = 0;
    my $matches = 0;
    my $insert = 0;
    my $del = 0;
    my $soft = 0;
    my $hard = 0;
    my $intron = 0;
    if ($samfields[4] < $minq) {
      $lowqual++;
    }
    while ($samfields[5]=~/(\d+)M/g) {
      $matches += $1;
    }
    while ($samfields[5]=~/(\d+)I/g) {
      $insert += $1;
    }
    while ($samfields[5]=~/(\d+)D/g) {
      $del += $1;
    }
    while ($samfields[5]=~/(\d+)S/g) {
      $soft += $1;
    }
    while ($samfields[5]=~/(\d+)H/g) {
      $hard += $1;
    }
    while ($samfields[5]=~/(\d+)N/g) {
      $intron += $1;
    }
    $refalnlen = $matches + $del + $intron;
    $to = $from + $refalnlen -1;
    my $qalnlen = $matches + $insert;
    my $readlen = $matches + $insert + $soft +$hard;
    # M+I+S=ALNQLEN
    # M
    my $lefthard = 0;
    my $righthard = 0;
    my $leftsoft = 0;
    my $rightsoft = 0;
    my $qstart = 0;
    my $qend = $readlen;
    if ($samfields[5]=~s/^(\d+)H//) {
      $lefthard=$1;
    }
    if ($samfields[5]=~/^(\d+)S/) {
      $leftsoft=$1;
    }
    if ($samfields[5]=~s/(\d+)H$//) {
      $righthard=$1;
    }
    if ($samfields[5]=~/(\d+)S$/) {
      $rightsoft=$1;
    }
    if ($strand eq "+") {
      $qstart = 1 + $lefthard +$leftsoft;
      $qend = $readlen - $righthard -$rightsoft;
    } else {
      $qstart = 1 + $righthard + $rightsoft;
      $qend = $readlen - $lefthard - $leftsoft;
    }
    if ($refalnlen >= $minmatch) {
      print STDERR $fastqstring;
    }else{
      print STDOUT $fastqstring;
    }
  }else{
    print STDOUT $fastqstring;
  }
  # print("$read\t$readlen\t$qstart\t$qend\t+\t$contig\t$refseqlen{$contig}\t$from\t$to\t$strand\t$mapq\t");printf("%0.2f\n",(1-(($insert+$del)/(($qalnlen+$refalnlen)/2))));
}
sub supplemental{
  my $flag = shift;
  my $supplemental = 0;
  if ($flag & 2048) {
    $supplemental = 1;
  }
  return $supplemental;		
}

sub strand{
  my $flag = shift;
  my $strand = '+';
  if ($flag & 0x10) {
    $strand = '-';
  }
  return $strand;		
}
