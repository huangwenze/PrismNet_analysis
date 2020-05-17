#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;


my $outfile = $ARGV[0];
my $usage = "This script is to merge the overlapped structure change window.
usage: $0 <outfile> 

example: perl merge_peak.pl dyn_str_number.txt
";
die $usage if $#ARGV<0;

my $inpath = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/dynamic_structure05/";
my $outpath1 = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/dynamic_str_region05/filter_win/"; 
my $outpath2 = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/dynamic_str_region05/combine_region/";

my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $file1 = ""; my $file2 = ""; my $file3 = ""; my $file4 = "";
my %inf = (); my %minf = (); my %prot = (); my %ginf = (); my %func = (); my %clipdb = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $flag = 0;
my $sta = 0; my $end = 0; my $cell_line = "";

my @cell = ("H9","HEK293","HEK293T","Hela","HepG2","K562");

opendir(TEMPDIR, $inpath) or die "can't open it:$inpath";
#HEK293_H9_win_icshape_test.txt

my @Dir2 = readdir TEMPDIR;
my @fa2 = grep /_icshape_test\.txt$/, @Dir2;
closedir TEMPDIR;

my $pvalue = 3; my $cutoff_score = 0.2279; my $len = 20;

open(OUT, ">", $outfile);

foreach $key (@fa2){
	my $infile = $inpath.$key;
	@sen1 = split(/_/, $key);
	my $outfile1 = $outpath1.$sen1[0]."_".$sen1[1]."_dystrwin.txt";
	my $outfile2 = $outpath1.$sen1[1]."_".$sen1[0]."_dystrwin.txt";
	my $outfile01 = $outpath2.$sen1[0]."_".$sen1[1]."_dystrsite.txt";
	my $outfile02 = $outpath2.$sen1[1]."_".$sen1[0]."_dystrsite.txt";
	my $num1 = 0; my $num2 = 0; my $num3 = 0;
	open(FILE1, $infile)||die("open $infile error!\n");
	open(OUT1, ">", $outfile1);
	open(OUT2, ">", $outfile2);
	open(OUT3, ">", $outfile01);
	open(OUT4, ">", $outfile02);
	#ENST00000480106	0	0	15	0.1082	0.50736
	#ENST00000480106	1	0	16	0.10144	0.50736
	$sid = ""; $sta = -1; $end = -1; $flag = 0;
	while($sen = <FILE1>){
		chomp($sen);
		$num1 = $num1 + 1;
		@sen2 = split(/\t/, $sen);
		if(($sen2[2] >= $pvalue)&&($sen2[4] > $cutoff_score)){
			$num2 = $num2 + 1;
			print OUT1 $sen2[0],"\t",$sen2[1]+1,"\t",$sen2[1]+20,"\t",$sen2[4],"\n";
			print OUT2 $sen2[0],"\t",$sen2[1]+1,"\t",$sen2[1]+20,"\t",$sen2[4],"\n";
			if($sen2[0] ne $sid){
				if($sid ne ""){
					print OUT3 $sid,"\t",$sta+1,"\t",$end+1,"\t",$flag,"\n";
					print OUT4 $sid,"\t",$sta+1,"\t",$end+1,"\t",$flag,"\n";
					$num3 = $num3 + 1;
				}
				$sid = $sen2[0]; $sta = $sen2[1]; $end = $sen2[1] + 19; $flag = $sen2[4];
			}else{
				if($sen2[1] <= $end){
					$end = $sen2[1] + 19;
					$flag = max($flag, $sen2[4]);
				}else{
					print OUT3 $sid,"\t",$sta+1,"\t",$end+1,"\t",$flag,"\n";
					print OUT4 $sid,"\t",$sta+1,"\t",$end+1,"\t",$flag,"\n";
					$num3 = $num3 + 1;
					$sid = $sen2[0]; $sta = $sen2[1]; $end = $sen2[1] + 19; $flag = $sen2[4];
				}
			}
		}
	}
	print OUT $sen1[0]."\t".$sen1[1],"\t",$num1,"\t",$num2,"\t",$num3,"\n";
	print OUT $sen1[1]."\t".$sen1[0],"\t",$num1,"\t",$num2,"\t",$num3,"\n";
	close FILE1;
	close OUT1;
	close OUT2;	
	close OUT3;
	close OUT4;	
}

close OUT;


exit;