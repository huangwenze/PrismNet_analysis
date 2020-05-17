#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

my $path = $ARGV[0];
my $outfile1 = $ARGV[1];
my $outfile2 = $ARGV[2];

#perl motif_compare2.pl /150T/zhangqf2/huangwz/structure_motif/human_dataset6/motif_kmer/motif_data/total_motif2/ motif_yeast_summary.txt motif_yeast_prob.txt
my $usage = "This script is to get the motif from deeplearn model output file.
usage: $0 <path> <outfile1> <outfile2>
";
die $usage if $#ARGV<2;

opendir(TEMPDIR, $path) or die "can't open it:$path";
#ALKBH5_HEK293_combine_summary.txt
#ALKBH5_HEK293_topmotif.txt

my @Dir = readdir TEMPDIR;
#my $prex = "_seq";
#my @fa1 = grep /${prex}$/, @Dir;
my @fa1 = grep /_combine_summary\.txt$/, @Dir;
#my @fa2 = grep /_HepG2_K562_ics$/, @Dir;
#my @fa = (@fa1, @fa2);
closedir TEMPDIR;

my @cells = ("H9","HEK293","HEK293T","Hela","HepG2","K562");
#my @cells = ("yeast");

open(OUT1, ">", $outfile1);
open(OUT2, ">", $outfile2);

for(my $j=0; $j<=$#fa1; $j++){
	my @sent1 = split(/\_/, $fa1[$j]);
	if(grep { $_ eq $sent1[1] } @cells){
		my $file1 = $path.$fa1[$j];
		my $file2 = $path.$sent1[0]."_".$sent1[1]."_topmotif.txt";
		open(FILE1, $file1)||die("open $file1 error!\n");
		open(FILE2, $file2)||die("open $file2 error!\n");
		#2
		#UGCAUG|UUUUUU	157
		#UGCAUG|PPPPPP	119
		my $sen = <FILE1>;
		my @sent2 = split(/\t/, $sen);
		if($sent2[0] <= 40){
			next;
		}
		my $r = 0;
		while($sen = <FILE1>){
			chomp($sen);
			$r = $r + 1;
			print OUT1 $sent1[0]."_".$sent1[1]."_".$r,"\t",$sen,"\n";
			print OUT2 $sent1[0]."_".$sent1[1]."_".$r;
			my $sen1 = <FILE2>; my $sen2 = <FILE2>; my $sen3 = <FILE2>; my $sen4 = <FILE2>; my $sen5 = <FILE2>; my $sen6 = <FILE2>;
			chomp($sen1); chomp($sen2); chomp($sen3); chomp($sen4); chomp($sen5); chomp($sen6);
			my @sen1 = split(/\t/, $sen1); my @sen2 = split(/\t/, $sen2); my @sen3 = split(/\t/, $sen3); my @sen4 = split(/\t/, $sen4); my @sen5 = split(/\t/, $sen5); my @sen6 = split(/\t/, $sen6);
			for(my $i=0; $i<=9; $i++){
				my $Sum1 = $sen1[$i] + $sen2[$i] + $sen3[$i] + $sen4[$i]; my $Sum2 = $sen5[$i] + $sen6[$i];
				my $new_var = sprintf("\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", $sen1[$i]/$Sum1, $sen2[$i]/$Sum1, $sen3[$i]/$Sum1, $sen4[$i]/$Sum1, $sen5[$i]/$Sum2, $sen6[$i]/$Sum2);
				print OUT2 $new_var;
			}
			print OUT2 "\n";
			if($r >= 5){
				last;
			}
		}
		close FILE1;
	}
}

close OUT1;
close OUT2;

exit;
