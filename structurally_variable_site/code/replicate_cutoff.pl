#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

#my $Path1 = $ARGV[0];
my $outfile = $ARGV[0]; 
my $outfile2 = $ARGV[1]; 

my $usage = "This script is to calculate the threshold value of L1 distance of single nucleotide between replicates.
usage: $0 <outfile> <outfile2> 

example: perl replicate_cutoff.pl dyn_str_cutoff01.txt dyn_str_cutoff02.txt
";
die $usage if $#ARGV<1;

my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $file1 = ""; my $file2 = ""; my $file3 = ""; my $file4 = "";
my %inf = (); my %minf = (); my %prot = (); my %ginf = (); my %func = (); my %clipdb = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $flag = 0;
my $sta = 0; my $end = 0; my $cell_line = "";

my @cell = ("H9","HEK293","HEK293T","Hela","HepG2","K562");
my $path1 = "/150T/zhangqf2/huangwz/total_smart_icshape/new_smartSHAPE_replicate/";


opendir(TEMPDIR, $path1) or die "can't open it:$path1";
#H9_smartSHAPE_N1.out
#H9_smartSHAPE_N2.out

my @Dir2 = readdir TEMPDIR;
my @fa2 = grep /_smartSHAPE_N1\.out$/, @Dir2;
closedir TEMPDIR;

#int(rand(10000))
#my $cutoff_score = 1.5626; my $len = 20;

my @dist1 = (); my @dist2 = ();
foreach $key (@fa2){
	my @sent1 = split(/_/, $key);
	my $infile1 = $path1.$key;
	my $infile2 = $path1.$sent1[0]."_smartSHAPE_N2.out";
	my $Sinf1 = &get_shape_value($infile1);
	my $Sinf2 = &get_shape_value($infile2);
	my %sinf1 = %{$Sinf1}; my %sinf2 = %{$Sinf2};
	foreach $sid (keys %sinf1){
		if(exists $sinf2{$sid}){
			@sen1 = @{$sinf1{$sid}}; 
			@sen2 = @{$sinf2{$sid}}; 
			if($#sen1 == $#sen2){
				for($i=0; $i<=$#sen1; $i++){
					if(($sen1[$i] ne "NULL")&&($sen2[$i] ne "NULL")){
						if(rand(100000)<10){
							push(@dist1, abs($sen1[$i]-$sen2[$i]));
						}
					}
				}
			}
		}
	}
}
@dist1 = sort {$a<=>$b} @dist1;

open(OUT1, ">", $outfile);
print OUT1 "total number: ",$#dist1+1,"\n";
my $dist1_cutoff = $dist1[int(0.99*($#dist1 + 1))];
print OUT1 "L1 distant 0.99 : ",$dist1_cutoff,"\n";
$dist1_cutoff = $dist1[int(0.975*($#dist1 + 1))];
print OUT1 "L1 distant 0.975 : ",$dist1_cutoff,"\n";
$dist1_cutoff = $dist1[int(0.95*($#dist1 + 1))];
print OUT1 "L1 distant 0.95 : ",$dist1_cutoff,"\n";
$dist1_cutoff = $dist1[int(0.90*($#dist1 + 1))];
print OUT1 "L1 distant 0.90 : ",$dist1_cutoff,"\n";
$dist1_cutoff = $dist1[int(0.85*($#dist1 + 1))];
print OUT1 "L1 distant 0.85 : ",$dist1_cutoff,"\n";
$dist1_cutoff = $dist1[int(0.80*($#dist1 + 1))];
print OUT1 "L1 distant 0.80 : ",$dist1_cutoff,"\n";
close OUT1;

open(OUT1, ">", $outfile2);
print OUT1 join("\n", @dist1),"\n";
close OUT1;

sub get_shape_value{
	my $trx_file = shift;
	my $sen = ""; my $seq = ""; my $sna = ""; 
	my %inf = (); 
	my @sen = (); my @sen1 = (); my @sen2 = ();
	open(FILE3, $trx_file)||die("open $trx_file error!\n");
	#ENST00000416931.1	372	*	NULL	NULL	NULL
	while($sen = <FILE3>){
		chomp($sen);
		@sen1 = split(/\t/,$sen);
		@sen2 = split(/\./,$sen1[0]);
		$inf{$sen2[0]} = [@sen1[3..($#sen1)]];
	}
	close FILE3;
	return \%inf;
}

exit;