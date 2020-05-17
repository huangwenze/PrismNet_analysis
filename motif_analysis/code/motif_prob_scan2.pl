#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;
use experimental qw(smartmatch);
use Parallel::ForkManager;


my ($cell1) = $ARGV[0];
my ($human_motif_name) = $ARGV[1];
my ($motif_file) = $ARGV[2];
my ($trx_motif_file) = $ARGV[3];
my ($clip_file) = $ARGV[4];

my $usage = "This script is to scan the transcriptome by using the PrismNet motif.
usage: $0 <cell1> <human_motif> <bind_peak_file> <trx_motif_file> <clip_file>

example: perl motif_prob_scan2.pl HepG2 SRSF1_combine_1 SRSF1_combine_1_motif.meme SRSF1_combine_1_motif_HepG2_trx_pu2.out SRSF1_HepG2_clipdata.txt
";
die $usage if $#ARGV<4;


my $trx_file = "/150T/zhangqf2/huangwz/Gencode/hg38_transcriptome_std_simple.fa";

#my $human_icbind_total_motif_file = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/motif_in_clippeak/icbind_motif/motif_human_prob_top5.txt";
my $human_icbind_total_motif_file = "motif_human_prob_top5.txt";

my $trx_shape_path = "/150T/zhangqf2/huangwz/total_smart_icshape/new_smartSHAPE_0412/";

my $shape_file1 = $trx_shape_path.$cell1."_smartSHAPE.out";

my $p_cutoff = 0.001;

my $Str_inf = &get_trx_ics($shape_file1); my $Seq_inf = &get_trx_seq($trx_file);
my %seq_inf = %{$Seq_inf};
my $Scode_inf = &get_trx_strcode($Str_inf);
my %scode_inf = %{$Scode_inf};

my %index_trx = ();
my $r = 0;
foreach my $ke1 (keys %scode_inf){
	if((exists $seq_inf{$ke1}) && (length($seq_inf{$ke1}) == length($scode_inf{$ke1}))){
		$r = $r + 1;
		my $k = int($r/1000);
		${$index_trx{$k}}{$ke1} = $r;
	}
}

my $Motif_INF1 = &pick_motif($human_icbind_total_motif_file, $human_motif_name, $motif_file);

open(FILE1, $clip_file)||die("open $clip_file error!\n");
#>ENST00000356638|1680|1779|1
my $num1 = 0; my %clip_inf1 = ();
while(my $sen = <FILE1>){
	chomp($sen);
	if($sen=~m/>/){
		$sen=~s/>//g;
		my @sen1 = split(/\|/, $sen);
		my $sid = $sen1[0];
		if(exists $clip_inf1{$sid}){
			push(@{$clip_inf1{$sid}}, $sen1[1], $sen1[2]);
		}else{
			$clip_inf1{$sid} = [$sen1[1], $sen1[2]];
		}
		$num1 = $num1 + 1;
	}
}
close FILE1;

my ($hscore, $lscore) = &cal_rank_prob($Motif_INF1, 6);

my %scan_inf = %{$hscore};
my $n1 = 0; my $n2 = 0; my $n3 = 0;
open(OUT, ">", $trx_motif_file);
foreach my $ke1 (keys %scode_inf){
	#if($n3 >= 1000){
	#	last;
	#}
	if((exists $seq_inf{$ke1}) && (length($seq_inf{$ke1}) == length($scode_inf{$ke1}))){
		my $len1 = length($seq_inf{$ke1});
		$seq_inf{$ke1} =~s/T/U/g;
		for(my $i=0; $i<=$len1-6; $i++){
			my $subseq01 = substr($seq_inf{$ke1}, $i, 6); my $substr01 = substr($scode_inf{$ke1}, $i, 6);
			my $count = $substr01 =~ tr/N/N/;
			if($count == 0){
				my $sna = $subseq01."|".$substr01;
				if(exists $scan_inf{$sna}){
					my $prob01 = $scan_inf{$sna}[0]; my $pvalue01 = $scan_inf{$sna}[1];
					#my $prob01 = &cal_site_prob($Motif_INF1, $subseq01, $substr01);
					#my $pvalue01 = &cal_pvalue($prob01, $lscore);
					if( $pvalue01 < $p_cutoff ){
						print OUT $ke1,"\t",$i+1,"\t",$i+6,"\t",$subseq01,"\t",$substr01,"\t",$prob01,"\t",$pvalue01,"\t";
						my $flag = 0; my $sta = $i+1; my $end = $i+6;
						if(exists $clip_inf1{$ke1}){
							my @sent1 = @{$clip_inf1{$ke1}};
							for(my $j=0; $j<=$#sent1; $j=$j+2){
								if( abs(($sent1[$j] + $sent1[$j+1])/2 - ($sta + $end)/2) < ($sent1[$j+1] - $sent1[$j])/2 + ($end - $sta)/2 ){
									$flag = 1;
									last;
								}
							}
						}
						print OUT $flag,"\n";
						if($flag > 0){
							$n2 = $n2 + 1;
						}
						$n1 = $n1 + 1;
					}
				}
			}
		}
		$n3 = $n3 + 1;
	}
}
close OUT;

print "Result:\t",$motif_file,"\t",$n1,"\t",$n2,"\n";



sub cal_pvalue{
	my $Score = shift; my $Plist = shift;
	my @plist = @{$Plist};
	my $i = 0;
	for($i=0; $i<=$#plist; $i++){
		if($Score > $plist[$i]){
			last;
		}
	}
	my $Rank = $i/($#plist + 1);
	return($Rank);
}


sub cal_rank_prob{
	my $Motif_inf = shift; my $len1 = shift;
	my @rseq = ("A", "C", "G", "U"); my @rstr = ("P", "U", "N");
	my %inf1 = (); 	my %inf2 = ();
	for(my $i=0; $i<=$#rseq; $i++){
		$inf1{$rseq[$i]} = 1;
	}
	for(my $i=0; $i<=$#rstr; $i++){
		$inf2{$rstr[$i]} = 1;
	}
	my $INF01 = &traversal(\%inf1, \@rseq, $len1-1);
	my $INF02 = &traversal(\%inf2, \@rstr, $len1-1);
	my %inf01 = %{$INF01}; my %inf02 = %{$INF02};
	my %pscore = ();
	foreach my $ke1 (keys %inf01){
		foreach my $ke2 (keys %inf02){
			my $sna = $ke1."|".$ke2;
			$pscore{$sna} = [&cal_site_prob($Motif_inf, $ke1, $ke2)];
		}
	}
	my @plist = (sort {$b <=> $a} values %pscore);
	my $tnum = $#plist + 1;
	my $r = 0;
	foreach my $ke (sort {${$pscore{$b}}[0] <=> ${$pscore{$a}}[0]} keys %pscore){
		$r = $r + 1;
		#$pscore{$ke} = $r/$tnum;
		push(@{$pscore{$ke}}, $r/$tnum);
	}
	return(\%pscore, \@plist);
}

sub traversal{
	my $INF = shift; my $Sent = shift; my $num = shift;
	my %inf1 = %{$INF}; my @sent = @{$Sent};
	my %inf2 = ();
	if($num == 1){
		foreach my $ke (keys %inf1){
			for(my $i=0; $i<=$#sent; $i++){
				my $sna = $ke.$sent[$i];
				$inf2{$sna} = 1;
			}
		}
		return(\%inf2);
	}elsif($num >= 2){
		foreach my $ke (keys %inf1){
			for(my $i=0; $i<=$#sent; $i++){
				my $sna = $ke.$sent[$i];
				$inf2{$sna} = 1;
			}
		}
		my $num2 = $num - 1;
		my $INF2 = &traversal(\%inf2, $Sent, $num2);
		return($INF2);
	}
}

sub cal_site_prob{
	my $Motif_inf = shift; my $Seq1 = shift; my $Shape = shift;
	my %motif_inf = %{$Motif_inf}; 
	my @sequ = split(//, $Seq1); my @shape = split(//, $Shape);
	my $j = 0;
	my $prob = 0;
	for(my $i=0; $i<=$#sequ; $i++){
		$j = $i + 3;
		if($sequ[$i] eq "A"){
			$prob = $prob + log(${$motif_inf{$j}}[0])/log(10);
		}elsif($sequ[$i] eq "C"){
			$prob = $prob + log(${$motif_inf{$j}}[1])/log(10);
		}elsif($sequ[$i] eq "G"){
			$prob = $prob + log(${$motif_inf{$j}}[2])/log(10);
		}else{
			$prob = $prob + log(${$motif_inf{$j}}[3])/log(10);
		}
		if($shape[$i] eq "P"){
			$prob = $prob + log(${$motif_inf{$j}}[4])/log(10);
		}elsif($shape[$i] eq "U"){
			$prob = $prob + log(${$motif_inf{$j}}[5])/log(10);
		}else{
			$prob = $prob + log(0.5)/log(10);
		}
	}
	return($prob);
}

sub pick_motif{
	my $total_icbind_motif_file = shift; my $motif_name = shift; my $motif_meme_file = shift;
	my %inf1 = ();
	my $header = "MEME version 4.10.1 (Release date: Wed Mar 25 11:40:43 2015 +1000)\nstrands: +\n\nMOTIF\t${motif_name}\nletter-probability matrix: alength= 4 w= 10\n";
	open(FILE1, $total_icbind_motif_file)||die("open $total_icbind_motif_file error!\n");
	open(OUT,">",$motif_meme_file);
	print OUT $header;
	#ATXN2_HEK293T_1	0.255682	0.255682	0.244318	0.244318	0.500000
	while(my $sen = <FILE1>){
		chomp($sen);
		my @sent = split(/\t/, $sen);
		if($sent[0] eq $motif_name){
			my $j = 1;
			for(my $i=1; $i<=$#sent; $i=$i+6){
				$inf1{$j} = [max($sent[$i], 0.000001), max($sent[$i+1], 0.000001), max($sent[$i+2], 0.000001), max($sent[$i+3], 0.000001), max($sent[$i+4], 0.000001), max($sent[$i+5], 0.000001)]; #ACGT,PU
				$j++;
				#my $new_var = sprintf("\s%.6f\s\s%.6f\s\s%.6f\s\s%.6f\s\n", $sent[$i], $sent[$i+1], $sent[$i+2], $sent[$i+3]);
				my $new_var = sprintf(" %.6f  %.6f  %.6f  %.6f \n", $sent[$i], $sent[$i+1], $sent[$i+2], $sent[$i+3]);
				print OUT $new_var;
			}
		}
	}
	close FILE1;
	close OUT;
	return(\%inf1);
}

sub icSHAPE_str2{
	my $Shape_value = shift;
	#my $Shape_list = shift; my $sta = shift; my $len = shift;
	my $str_seq = ""; my $str_pro = "";
	my @dvalue = (0.0, 0.088, 0.233, 0.484, 1.0);
	if($Shape_value < $dvalue[0]){
		$str_seq = "N";
		$str_pro =	0;	
	}elsif($Shape_value <= $dvalue[1]){
		$str_seq = "P";
		$str_pro = 1 - ($Shape_value - $dvalue[0])/($dvalue[1] - $dvalue[0])*0.5;
	}elsif($Shape_value <= $dvalue[2]){
		$str_seq = "P";
		$str_pro = 0.5 - ($Shape_value - $dvalue[1])/($dvalue[2] - $dvalue[1])*0.5;
	}elsif($Shape_value <= $dvalue[3]){
		$str_seq = "U";
		$str_pro = ($Shape_value - $dvalue[2])/($dvalue[3] - $dvalue[2])*0.5;
	}else{
		$str_seq = "U";
		$str_pro = 0.5 + ($Shape_value - $dvalue[3])/($dvalue[4] - $dvalue[3])*0.5;
	}
	return ($str_seq, $str_pro);
}

sub icSHAPE_str{
	my $Shape_value = shift;
	#my $Shape_list = shift; my $sta = shift; my $len = shift;
	my $str_seq = ""; my $str_pro = "";
	my @dvalue = (0.0, 0.088, 0.233, 0.484, 1.0);
	if($Shape_value < $dvalue[0]){
		$str_seq = "N";
	}elsif($Shape_value <= $dvalue[2]){
		$str_seq = "P";
	}else{
		$str_seq = "U";
	}
	return ($str_seq);
}


sub get_trx_strcode{
	my $INF1 = shift; 
	my %inf1 = %{$INF1}; my %inf2 = ();
	foreach my $ke (keys %inf1){
		my @sent1 = @{$inf1{$ke}};
		my $scode = "";
		for(my $i=0; $i<=$#sent1; $i++){
			my $code = &icSHAPE_str($sent1[$i]);
			$scode = $scode.$code;
		}
		$inf2{$ke} = $scode;
	}
	return(\%inf2);
}


sub get_trx_ics{
	my $trx_file = shift;
	my $sen = ""; my $seq = ""; my $sna = ""; 
	my %inf = (); 
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $sid = "";
	my $i = 0; my $num = 0;
	open(FILE3, $trx_file)||die("open $trx_file error!\n");
	#ENST00000531760	502	12.4653677062886	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	
	while($sen = <FILE3>){
		chomp($sen);
		$sen =~s/NULL/-1/g;
		@sen1 = split(/\t/,$sen);
		@sen2 = split(/\./,$sen1[0]);
		$inf{$sen2[0]} = [@sen1[3..$#sen1]];
		$i = $i + 1;
	}
	close FILE3;
	return \%inf;
}

sub get_trx_seq{
	my $trx_file = shift;
	my $sen = ""; my $seq = ""; my $sna = ""; 
	my %inf = (); 
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $sid = "";
	my $i = 0; my $num = 0;
	open(FILE3, $trx_file)||die("open $trx_file error!\n");
	#>ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-002|DDX11L1|1657|processed_transcript|
	while($sen = <FILE3>){
		chomp($sen);
		if($sen =~m/^>/){
			$sen =~s/>//g;
			@sen1 = split(/\|/,$sen);
			@sen2 = split(/\./,$sen1[0]);
			$sen = <FILE3>;
			chomp($sen);
			$inf{$sen2[0]} = $sen;
		}
	}
	close FILE3;
	return \%inf;
}

exit;