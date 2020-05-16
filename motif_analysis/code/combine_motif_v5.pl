#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

#my $prot_cell = $ARGV[0];
my $protein_list = $ARGV[0];

my $usage = "This script is to build the motif from PrismNet model output attention file. For the RBP with two or more models, we combined their HARs and build the combined motifs.
usage: $0 <protein_list>
";
die $usage if $#ARGV<0;

my $path1 = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/motif_refine/RBP_motif2/";

my %prot_inf = ();
open(FILE1, $protein_list)||die("open $protein_list error!\n");
#ALKBH5_HEK293
#ATXN2_HEK293T
while(my $sen = <FILE1>){
	chomp($sen);
	my @sent1 = split(/_/, $sen);
	if(exists $prot_inf{$sent1[0]}){
		push(@{$prot_inf{$sent1[0]}}, $sen);
	}else{
		$prot_inf{$sent1[0]} = [$sen];
	}
}
close FILE1;

foreach my $key (keys %prot_inf){
	my @sent1 = @{$prot_inf{$key}};
	if($#sent1 == 0){
		next;
	}
	my @file_list = ();
	for(my $i=0; $i<=$#sent1; $i++){
		my $file1 = $path1.$sent1[$i]."_seq_6kmer_seq.txt";
		push(@file_list, $file1);
	}
	my $Kmer1 = &read_kmer_file(\@file_list);
	my ($Motif_matrix1, $Motif_matrix2) = &combine_kmer($Kmer1, "-ACGT", $key."_combine");
	&motif_print($Motif_matrix1, $key."_combine", 10, $key."_combine_motif_10_seq.meme", "-ACGT");
	my %mmat1 = %{$Motif_matrix1};
	my %mmat2 = %{$Motif_matrix2};
	open(OUT1, ">", $key."_combine_motif_10_str.meme");
	foreach my $k1 (sort {$a<=>$b} keys %mmat2){
		my %tmp = %{$mmat2{$k1}};
		print OUT1 $k1,"\n";
		foreach my $k1 (sort {$a<=>$b} keys %tmp){
			#print OUT1 ${$tmp{$k1}}[0],"|",${$tmp{$k1}}[1],"\t";
			print OUT1 sprintf("%.4f", ${$tmp{$k1}}[0]),"|",sprintf("%.4f", ${$tmp{$k1}}[1]),"\t";
		}
		print OUT1 "\n";
	}
	close OUT1;
	
	open(OUT1, ">", $key."_combine_motif_10_seq.txt");
	foreach my $k1 (sort {$a<=>$b} keys %mmat1){
		my %tmp = %{$mmat1{$k1}};
		#print OUT1 $k1,"\n";
		for(my $i=0; $i<=3; $i++){
			foreach my $k2 (sort {$a<=>$b} keys %tmp){
				#print OUT1 ${$tmp{$k1}}[0],"|",${$tmp{$k1}}[1],"\t";
				print OUT1 sprintf("%.4f", ${$tmp{$k2}}[$i]),"\t";
			}
			print OUT1 "\n";
		}	
	}
	close OUT1;
	open(OUT1, ">", $key."_combine_motif_10_str.txt");
	foreach my $k1 (sort {$a<=>$b} keys %mmat2){
		my %tmp = %{$mmat2{$k1}};
		#print OUT1 $k1,"\n";
		for(my $i=0; $i<=1; $i++){
			foreach my $k2 (sort {$a<=>$b} keys %tmp){
				#print OUT1 ${$tmp{$k1}}[0],"|",${$tmp{$k1}}[1],"\t";
				print OUT1 sprintf("%.4f", ${$tmp{$k2}}[$i]),"\t";
			}
			print OUT1 "\n";
		}	
	}
	close OUT1;
	#last;
}


#id，int 
#label，int 
#Predictscore，float 
#Sequence，str,101 
#Icshape, float,101 
#Saliency,101x5(5v), 101x6(6v,7v)


sub subshape{
	my $Shape_list = shift; my $sta = shift; my $len = shift;
	#my @dvalue = (0.0, 0.088, 0.233, 0.484, 1.0);
	my @sent1 = split(/\|/, $Shape_list);
	my @sent2 = ();
	#shift(@sent1);
	my $i = 0; my $str_seq = "";
	for($i=$sta; $i<=$sta+$len-1; $i++){
		push(@sent2, $sent1[$i]);
	}
	$str_seq = join("|", @sent2);
	return ($str_seq);
}

sub subshape2{
	my $Shape_list = shift; my $sta = shift; my $len = shift;
	my @dvalue = (0.0, 0.088, 0.233, 0.484, 1.0);
	my @sent1 = split(/\|/, $Shape_list);
	shift(@sent1);
	my $i = 0; my $str_seq = "";
	for($i=$sta; $i<=$sta+$len-1; $i++){
		if($sent1[$i] <= $dvalue[1]){
			$str_seq = $str_seq."P";
		}elsif($sent1[$i] <= $dvalue[2]){
			$str_seq = $str_seq."Q";
		}elsif($sent1[$i] <= $dvalue[3]){
			$str_seq = $str_seq."S";
		}else{
			$str_seq = $str_seq."Z";
		}
	}
	return ($str_seq);
}

sub shape2str{
	my $Shape_list = shift; my $Str_list = shift;
	my @dvalue = (0.0, 0.088, 0.233, 0.484, 1.0);
	my @sent1 = split(/\|/, $Shape_list);
	my @sent2 = split(//, $Str_list);
	#shift(@sent1);
	my $i = 0; my $str_seq = "";
	for($i=0; $i<=$#sent1; $i++){
		if($sent1[$i] <= 0){
			$str_seq = $str_seq.$sent2[$i];
		}elsif($sent1[$i] <= $dvalue[2]){
			$str_seq = $str_seq."P";
		}else{
			$str_seq = $str_seq."U";
		}
	}
	return ($str_seq);
}

sub max_index{
	my $list = shift; my $len = shift;
	my @sen1 = @{$list};
	my $i = 0; my $j = 0; my $r = 0; my $index = 0; my $maxn = 0; my $sum = 0;
	for($i=0; $i<=$#sen1-$len+1; $i++){
		$sum = 0;
		for($j=0; $j<$len; $j++){
			$sum = $sum + $sen1[$i+$j];
		}
		if($sum > $maxn){
			$index = $i;
			$maxn = $sum;
		}
	}
	return ($index, $maxn);
}

sub max_per_seq{
	my $Inf = shift; my $len = shift; my $per = shift; my $bind_score = shift;
	#my $list = shift; my $len = shift;
	my %inf = %{$Inf};
	my $key = ""; my @total = (); my $i = 0;
	foreach $key (sort {${$inf{$b}}[0] <=> ${$inf{$a}}[0]} keys %inf){
		if(${$inf{$key}}[0] < $bind_score){
			last;
		}
		my @sent1 = split(/\|/, ${$inf{$key}}[2]);
		shift(@sent1);
		for($i=0; $i<=$#sent1-$len+1; $i++){
			push(@total, sum(@sent1[$i..($i+$len-1)]));
		}
	}
	@total = sort {$b <=> $a} @total;
	my $boun = $total[int(($#total + 1)*$per)-1];
	return ($boun);
}

sub max_per_seq_str{
	my $Inf = shift; my $len = shift; my $per = shift; my $bind_score = shift;
	#my $list = shift; my $len = shift;
	my %inf = %{$Inf};
	my $key = ""; my @total1 = (); my @total2 = (); my $i = 0;
	foreach $key (sort {${$inf{$b}}[0] <=> ${$inf{$a}}[0]} keys %inf){
		if(${$inf{$key}}[0] < $bind_score){
			last;
		}
		my @sent1 = split(/\|/, ${$inf{$key}}[2]);
		shift(@sent1);
		my @sent2 = split(/\|/, ${$inf{$key}}[4]);
		shift(@sent2);
		for($i=0; $i<=$#sent1-$len+1; $i++){
			push(@total1, sum(@sent1[$i..($i+$len-1)]));
		}
		for($i=0; $i<=$#sent2-$len+1; $i++){
			push(@total2, sum(@sent2[$i..($i+$len-1)]));
		}
	}
	@total1 = sort {$b <=> $a} @total1;
	@total2 = sort {$b <=> $a} @total2;
	my $boun1 = $total1[int(($#total1 + 1)*$per)-1];
	my $boun2 = $total2[int(($#total2 + 1)*$per*2)-1];
	return ($boun1, $boun2);
}

sub motif_print{
	#my $ref_ics = "/150T/zhangqf2/huangwz/icshape_value/293T.invivo.out";
	#my $ref_trx = "/150T/zhangqf2/huangwz/CLIP/gencode.v26.transcripts.ics.std.fa";
	my $Motif = shift; my $out_pref = shift; my $motif_len = shift; my $pmotif_file = shift; my $ALPHA = shift;
	my %mot_inf = %{$Motif};
	my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $ics = "";
	my @sen = (); my @sen1 = (); my @sen2 = (); my @sent1 = (); my @sent2 = ();
	my @alpha = split(//, $ALPHA);
	my $sid = ""; my $key;
	my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tsum = 0; my $k1 = 0;
	
	#my $head = "MEME version 4.10.1 (Release date: Wed Mar 25 11:40:43 2015 +1000)\nstrands: +\n\nMOTIF ";
	#my $head_str = "MEME version 4.10.1 (Release date: Wed Mar 25 11:40:43 2015 +1000)\n\nALPHABET\nP\nQ\nS\nZ\nEND ALPHABET\n\nstrands: +\n\nMOTIF ";
	my $head = "MEME version 4.10.1 (Release date: Wed Mar 25 11:40:43 2015 +1000)\n\nALPHABET\n".$alpha[1]."\n".$alpha[2]."\n".$alpha[3]."\n".$alpha[4]."\nEND ALPHABET\n\nstrands: +\n\nMOTIF ";
	
	#open(FILE1, $pf_file)||die("open $pf_file error!\n");
	open(OUT1, ">", $pmotif_file);
	$j = 0;
	foreach $key (sort{$a<=>$b} keys %mot_inf){
		print OUT1 $head;
		print OUT1 $out_pref,$key,"\n";
		my %motifi = %{$mot_inf{$key}};
		@sen1 = @{$motifi{1}};
		$tsum = $sen1[0] + $sen1[1] + $sen1[2] + $sen1[3];
		print OUT1 "letter-probability matrix: alength= 4 w= $motif_len nsites = $tsum\n";
		foreach $i (sort{$a<=>$b} keys %motifi){
			@sen1 = @{$motifi{$i}};
			$tsum = $sen1[0] + $sen1[1] + $sen1[2] + $sen1[3];
			if($tsum > 0){
				my $new_var = sprintf(" %.6f  %.6f  %.6f  %.6f \n", $sen1[0]/$tsum, $sen1[1]/$tsum, $sen1[2]/$tsum, $sen1[3]/$tsum);
				print OUT1 $new_var;
			}else{
				my $new_var = sprintf(" %.6f  %.6f  %.6f  %.6f \n", 0.25, 0.25, 0.25, 0.25);
				print OUT1 $new_var;
			}
		}
		print OUT1 "\n";
	}
	close OUT1;
}

sub kmer_cal{
	my $Inf = shift; my $len = shift;
	my %inf = %{$Inf}; my %kmer_seq = (); my %kmer_loc = ();
	my @sent = ();
	my $i = 0; my $key = ""; my $seq = ""; my $subseq = ""; my $sta = 0; my $end = 0; my $sna = "";
	foreach $key (sort {${$inf{$b}}[0] <=> ${$inf{$a}}[0]} keys %inf){
		@sent = split(/\_/, $key);
		$seq = ${$inf{$key}}[2];
		for($i=0; $i<=length($seq)-$len; $i++){
			$subseq = substr($seq, $i, $len);
			$sta = $sent[1] + $i; $end = $sta + $len - 1; $sna = $sent[0]."_".$sta."_".$end;
			if(exists $kmer_seq{$subseq}){
				$kmer_seq{$subseq} = $kmer_seq{$subseq} + 1;
				push(@{$kmer_loc{$subseq}}, $sna);
			}else{
				$kmer_seq{$subseq} = 1;
				$kmer_loc{$subseq} = [$sna];
			}
		}
	}
	return (\%kmer_seq, \%kmer_loc);
}

sub kmer_cal2{
	my $Inf = shift; my $len = shift;
	my %inf = %{$Inf}; my %kmer_seq = (); my %kmer_loc = ();
	my @sent = ();
	my $i = 0; my $key = ""; my $seq = ""; my $subseq = ""; my $sta = 0; my $end = 0; my $sna = "";
	my $stru = ""; my $substru = ""; my $kmer_name = "";
	foreach $key (sort {${$inf{$b}}[0] <=> ${$inf{$a}}[0]} keys %inf){
		@sent = split(/\_/, $key);
		$seq = ${$inf{$key}}[2];
		#$stru = ${$inf{$key}}[4];
		$stru = &shape2str(${$inf{$key}}[3], ${$inf{$key}}[4]);
		for($i=0; $i<=length($seq)-$len; $i++){
			$subseq = substr($seq, $i, $len);
			$substru = substr($stru, $i, $len);
			$sta = $sent[1] + $i; $end = $sta + $len - 1; $sna = $sent[0]."_".$sta."_".$end;
			$kmer_name = $subseq."|".$substru;
			if(exists $kmer_seq{$kmer_name}){
				$kmer_seq{$kmer_name} = $kmer_seq{$kmer_name} + 1;
				push(@{$kmer_loc{$kmer_name}}, $sna);
			}else{
				$kmer_seq{$kmer_name} = 1;
				$kmer_loc{$kmer_name} = [$sna];
			}
		}
	}
	return (\%kmer_seq, \%kmer_loc);
}

sub read_kmer_file{
	my $file_list = shift;
	my %inf = ();
	my @sen1 = @{$file_list};
	foreach my $file1 (@sen1){
		open(FILE1, $file1)||die("open $file1 error!\n");
		#AUUUGU|UUUUUU	7
		#UUUUGU|UUUUUU	6
		while(my $sen = <FILE1>){
			chomp($sen);
			my @sent1 = split(/\t/, $sen);
			if(exists $inf{$sent1[0]}){
				$inf{$sent1[0]} = $inf{$sent1[0]} + $sent1[1];
			}else{
				$inf{$sent1[0]} = $sent1[1];
			}
		}
		close FILE1;
	}
	return(\%inf);
}

sub combine_kmer{
	my $Inf = shift; my $ALPHA = shift; my $protein_name = shift;
	my %inf = %{$Inf}; my %cinf = (); my %mot_inf = (); my %mot_str_inf = (); my %cinf_con = ();
	#my %kmer_loc = %{$Kmer_loc}; my %data_inf = %{$Data_Inf};
	my $key = ""; my $k1 = ""; my $r = 0; my $flag = 0; my $exkey = "";
	my $kmer_sht = sum(values %inf)*0.2;
	my $kmer_sh = 0; my $kmer_sum = 0;
	foreach $key ( sort{$inf{$b} <=> $inf{$a}} keys %inf){
		$kmer_sum = $kmer_sum + $inf{$key};
		if($kmer_sum > $kmer_sht){
			$kmer_sh = $inf{$key};
			last;
		}
	}
	if(max(values %inf) <= 5){
		$kmer_sh = 0;
	}
	print $kmer_sh,"\n";
	#my $kmer_sh = 0;
	open(OUT, ">", $protein_name."_summary.txt");
	open(OUT2, ">", $protein_name."_summary2.txt");
	print OUT $kmer_sh,"\n";
	foreach $key ( sort{$inf{$b} <=> $inf{$a}} keys %inf){
		if($inf{$key} <= $kmer_sh){
			last;
		}
		$flag = 0;
		foreach $k1 (sort{${$cinf{$b}}[0] <=> ${$cinf{$a}}[0]} keys %cinf){
			($exkey, $flag) = &tcluster($k1, $key);
			if($flag == 1){
				${$cinf{$k1}}[0] = ${$cinf{$k1}}[0] + $inf{$key};
				push(@{$cinf{$k1}}, $exkey, $inf{$key});
				push(@{$cinf_con{$k1}}, $key);
				last;
			}
		}
		if($flag == 0){
			$cinf{$key} = [$inf{$key}, $inf{$key}];
			$cinf_con{$key} = [$key];
		}
	}
	$r = 0;
	foreach $key ( sort {${$cinf{$b}}[0] <=> ${$cinf{$a}}[0]} keys %cinf){
		$r = $r + 1;
		my ($mot1, $mot2, $num1) = &build_motif2($key, $cinf{$key}, $ALPHA, "-PU");
		#if($r < 10){
			print OUT $key,"\t",$num1,"\n";
			print OUT2 $key,"\t",$num1,"\n";
			for(my $i=0; $i<=$#{$cinf{$key}}; $i++){
				print OUT2 ${$cinf{$key}}[$i],"\t";
			}
			print OUT2 "\n";
		#}
		$mot_inf{$r} = $mot1; 
		$mot_str_inf{$r} = $mot2;
		#last;
	}
	close OUT;
	close OUT2;
	return (\%mot_inf, \%mot_str_inf);
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
		$str_pro = 1 - ($Shape_value - $dvalue[0])/($dvalue[1] - $dvalue[0])*0.8;
	}elsif($Shape_value <= $dvalue[2]){
		$str_seq = "P";
		$str_pro = 0.2 - ($Shape_value - $dvalue[1])/($dvalue[2] - $dvalue[1])*0.2;
	}elsif($Shape_value <= $dvalue[3]){
		$str_seq = "U";
		$str_pro = ($Shape_value - $dvalue[2])/($dvalue[3] - $dvalue[2])*0.2;
	}else{
		$str_seq = "U";
		$str_pro = 0.2 + ($Shape_value - $dvalue[3])/($dvalue[4] - $dvalue[3])*0.8;
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

sub fivechar{
	my $char1 = shift; my $char2 = shift;
	my @sen1 = split(//, $char1); my @sen2 = split(//, $char2); my $i = 0; my $r = 0; my $mismatch = 0; my $flag = 0;
	if(substr($char1, 0, 4) eq substr($char2, 1, 4)){
		return ("-".$char2."---", 1);
	}elsif(substr($char1, 1, 4) eq substr($char2, 0, 4)){
		return ("---".$char2."-", 1);
	}
	for($i=0; $i<=$#sen1; $i++){
		if($sen1[$i] eq $sen2[$i]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - $r;
	if($mismatch == 1){
		return ("--".$char2."--", 1);
	}
	if(substr($char1, 0, 3) eq substr($char2, 2, 3)){
		return ($char2."----", 1);
	}elsif(substr($char1, 2, 3) eq substr($char2, 0, 3)){
		return ("----".$char2, 1);
	}
	for($i=0; $i<=$#sen1-1; $i++){
		if($sen1[$i] eq $sen2[$i+1]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - 1 - $r;
	if($mismatch == 1){
		return ("-".$char2."---", 1);
	}
	for($i=1; $i<=$#sen1; $i++){
		if($sen1[$i] eq $sen2[$i-1]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - 1 - $r;
	if($mismatch == 1){
		return ("---".$char2."-", 1);
	}
	return($char2, 0);
}

sub tcluster{
	my $char1 = shift; my $char2 = shift;
	my @sen1 = split(/\|/, $char1); my @sen2 = split(/\|/, $char2);
	my @sent1 = split(//, $sen1[1]);
	my $r = 0; my $i = 0;
	for($i=0; $i<=$#sent1; $i++){
		if($sent1[$i] eq "P"){
			$r = $r + 1;
		}
	}
	#$r = $r/($#sent1 + 1);
	my ($ch1, $flag1) = &clusterchar1($sen1[0], $sen2[0]);
	my ($ch2, $flag2) = &clusterchar1($sen1[1], $sen2[1]);
	if(($flag1 == 1)&&($flag2 == 1)){
		return($ch1."|".$ch2, 1);
	}
	if(($flag1 == 0)&&($flag2 == 1)&&($r > 4)){
		my ($ch01, $flag01) = &clusterchar2($sen1[0], $sen2[0]);
		if($flag01 == 1){
			return($ch01."|".$ch2, 1);
		}
		if(&mismatch($sen1[0], $sen2[0]) <= $r - 3){
			return("--".$sen2[0]."--|".$ch2, 1);
		}
	}
	return($char2, 0);
}

sub tcluster2{
	my $char1 = shift; my $char2 = shift;
	my @sen1 = split(/\|/, $char1); my @sen2 = split(/\|/, $char2);
	my @sent1 = split(//, $sen1[1]);
	my $r = 0; my $i = 0;
	for($i=0; $i<=$#sent1; $i++){
		if($sent1[$i] eq "P"){
			$r = $r + 1;
		}
	}
	$r = $r/($#sent1 + 1);
	my ($ch1, $flag1) = &clusterchar1($sen1[0], $sen2[0]);
	my ($ch2, $flag2) = &clusterchar1($sen1[1], $sen2[1]);
	if(($flag1 == 1)&&($flag2 == 1)){
		return($ch1."|".$ch2, 1);
	}elsif(($flag1 == 1)&&($flag2 == 0)){
		my ($ch02, $flag02) = &clusterchar2($sen1[1], $sen2[1]);
		if($flag02 == 1){
			return($ch1."|".$ch02, 1);
		}
	}elsif(($flag1 == 0)&&($flag2 == 1)){
		my ($ch01, $flag01) = &clusterchar2($sen1[0], $sen2[0]);
		if($flag01 == 1){
			return($ch01."|".$ch2, 1);
		}
	}
	return($char2, 0);
}

sub clusterchar1{
	my $char1 = shift; my $char2 = shift;
	my $tnum = length($char1);
	my @sen1 = split(//, $char1); my @sen2 = split(//, $char2); my $i = 0; my $r = 0; my $mismatch = 0; my $flag = 0;
	if($char1 eq $char2){
		return ("--".$char2."--", 1);
	}elsif(substr($char1, 0, $tnum-1) eq substr($char2, 1, $tnum-1)){
		return ("-".$char2."---", 1);
	}elsif(substr($char1, 1, $tnum-1) eq substr($char2, 0, $tnum-1)){
		return ("---".$char2."-", 1);
	}
	$r = 0;
	for($i=0; $i<=$#sen1; $i++){
		if($sen1[$i] eq $sen2[$i]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - $r;
	if($mismatch <= 1){
		return ("--".$char2."--", 1);
	}
	return($char2, 0);
}

sub clusterchar2{
	my $char1 = shift; my $char2 = shift;
	my $tnum = length($char1);
	my @sen1 = split(//, $char1); my @sen2 = split(//, $char2); my $i = 0; my $r = 0; my $mismatch = 0; my $flag = 0;
	if($char1 eq $char2){
		return ("--".$char2."--", 1);
	}elsif(substr($char1, 0, $tnum-1) eq substr($char2, 1, $tnum-1)){
		return ("-".$char2."---", 1);
	}elsif(substr($char1, 1, $tnum-1) eq substr($char2, 0, $tnum-1)){
		return ("---".$char2."-", 1);
	}
	$r = 0;
	for($i=0; $i<=$#sen1; $i++){
		if($sen1[$i] eq $sen2[$i]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - $r;
	if($mismatch <= 1){
		return ("--".$char2."--", 1);
	}
	if(substr($char1, 0, $tnum-2) eq substr($char2, 2, $tnum-2)){
		return ($char2."----", 1);
	}elsif(substr($char1, 2, $tnum-2) eq substr($char2, 0, $tnum-2)){
		return ("----".$char2, 1);
	}
	$r = 0;
	for($i=0; $i<=$#sen1-1; $i++){
		if($sen1[$i] eq $sen2[$i+1]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - 1 - $r;
	if($mismatch <= 1){
		return ("-".$char2."---", 1);
	}
	$r = 0;
	for($i=1; $i<=$#sen1; $i++){
		if($sen1[$i] eq $sen2[$i-1]){
			$r = $r + 1;
		}
	}
	$mismatch = length($char1) - 1 - $r;
	if($mismatch <= 1){
		return ("---".$char2."-", 1);
	}
	return($char2, 0);
}

sub build_motif{
	my $rep_seq = shift; my $cont = shift; my $ALPHA = shift;
	$rep_seq = "--".$rep_seq."--";
	my @cons = @{$cont}; my @sen1 = split(//, $rep_seq); my @sen2 = split(//, $ALPHA);
	my %minf = (); 
	my $i=0; my $j=0; my $tnum = $cons[1];
	for($j=0; $j<=$#sen1; $j++){
		if($sen1[$j] eq $sen2[0]){
			$minf{$j} = [$tnum*0.25,$tnum*0.25,$tnum*0.25,$tnum*0.25];
		}elsif($sen1[$j] eq $sen2[1]){
			$minf{$j} = [$tnum,0,0,0];
		}elsif($sen1[$j] eq $sen2[2]){
			$minf{$j} = [0,$tnum,0,0];
		}elsif($sen1[$j] eq $sen2[3]){
			$minf{$j} = [0,0,$tnum,0];
		}else{
			$minf{$j} = [0,0,0,$tnum];
		}
	}
	for($i=2; $i<=$#cons; $i=$i+2){
		@sen1 = split(//, $cons[$i]);
		$tnum = $cons[$i+1];
		for($j=0; $j<=$#sen1; $j++){
			if($sen1[$j] eq $sen2[0]){
				${$minf{$j}}[0] = ${$minf{$j}}[0] + $tnum*0.25;
				${$minf{$j}}[1] = ${$minf{$j}}[1] + $tnum*0.25;
				${$minf{$j}}[2] = ${$minf{$j}}[2] + $tnum*0.25;
				${$minf{$j}}[3] = ${$minf{$j}}[3] + $tnum*0.25;
			}elsif($sen1[$j] eq $sen2[1]){
				${$minf{$j}}[0] = ${$minf{$j}}[0] + $tnum;
			}elsif($sen1[$j] eq $sen2[2]){
				${$minf{$j}}[1] = ${$minf{$j}}[1] + $tnum;
			}elsif($sen1[$j] eq $sen2[3]){
				${$minf{$j}}[2] = ${$minf{$j}}[2] + $tnum;
			}else{
				${$minf{$j}}[3] = ${$minf{$j}}[3] + $tnum;
			}
		}		
	}
	$tnum = $cons[0];
	return (\%minf, $tnum);
}

sub build_motif2{
	my $rep_seq = shift; my $cont = shift; my $ALPHA1 = shift; my $ALPHA2 = shift;
	#$rep_seq = "--".$rep_seq."--";
	my @Sen = split(/\|/, $rep_seq);
	my @sent1 = split(//, "--".$Sen[0]."--"); my @sent2 = split(//, "--".$Sen[1]."--");
	my @cons = @{$cont}; my @sen1 = split(//, $ALPHA1); my @sen2 = split(//, $ALPHA2);
	my %minf1 = (); my %minf2 = ();
	my $i=0; my $j=0; my $tnum = $cons[1];
	for($j=0; $j<=$#sent1; $j++){
		if($sent1[$j] eq $sen1[0]){
			$minf1{$j} = [$tnum*0.25,$tnum*0.25,$tnum*0.25,$tnum*0.25];
		}elsif($sent1[$j] eq $sen1[1]){
			$minf1{$j} = [$tnum,0,0,0];
		}elsif($sent1[$j] eq $sen1[2]){
			$minf1{$j} = [0,$tnum,0,0];
		}elsif($sent1[$j] eq $sen1[3]){
			$minf1{$j} = [0,0,$tnum,0];
		}else{
			$minf1{$j} = [0,0,0,$tnum];
		}
	}
	for($j=0; $j<=$#sent2; $j++){
		if($sent2[$j] eq $sen2[0]){
			$minf2{$j} = [$tnum*0.5,$tnum*0.5];
		}elsif($sent2[$j] eq $sen2[1]){
			$minf2{$j} = [$tnum,0];
		}else{
			$minf2{$j} = [0,$tnum];
		}
	}
	for($i=2; $i<=$#cons; $i=$i+2){
		@Sen = split(/\|/, $cons[$i]);
		@sent1 = split(//, $Sen[0]); @sent2 = split(//, $Sen[1]); $tnum = $cons[$i+1];
		for($j=0; $j<=$#sent1; $j++){
			if($sent1[$j] eq $sen1[0]){
				${$minf1{$j}}[0] = ${$minf1{$j}}[0] + $tnum*0.25;
				${$minf1{$j}}[1] = ${$minf1{$j}}[1] + $tnum*0.25;
				${$minf1{$j}}[2] = ${$minf1{$j}}[2] + $tnum*0.25;
				${$minf1{$j}}[3] = ${$minf1{$j}}[3] + $tnum*0.25;
			}elsif($sent1[$j] eq $sen1[1]){
				${$minf1{$j}}[0] = ${$minf1{$j}}[0] + $tnum;
			}elsif($sent1[$j] eq $sen1[2]){
				${$minf1{$j}}[1] = ${$minf1{$j}}[1] + $tnum;
			}elsif($sent1[$j] eq $sen1[3]){
				${$minf1{$j}}[2] = ${$minf1{$j}}[2] + $tnum;
			}else{
				${$minf1{$j}}[3] = ${$minf1{$j}}[3] + $tnum;
			}
		}
		for($j=0; $j<=$#sent2; $j++){
			if($sent2[$j] eq $sen2[0]){
				${$minf2{$j}}[0] = ${$minf2{$j}}[0] + $tnum*0.5;
				${$minf2{$j}}[1] = ${$minf2{$j}}[1] + $tnum*0.5;
			}elsif($sent2[$j] eq $sen2[1]){
				${$minf2{$j}}[0] = ${$minf2{$j}}[0] + $tnum;
			}else{
				${$minf2{$j}}[1] = ${$minf2{$j}}[1] + $tnum;
			}
		}
	}
	$tnum = $cons[0];
	return (\%minf1, \%minf2, $tnum);
}

sub kmer_cal3{
	my $Inf = shift; 
	my %inf = %{$Inf}; my %kmer_seq = (); my %kmer_loc = ();
	my @sent = (); my @psign1 = (0)x(16); my @usign1 = (0)x(16); my @sen1 = (); my @sen2 = (); my @sen3 = ();
	my $i = 0; my $j = 0; my $r = 0; 
	my $key = ""; my $seq = ""; my $subseq = ""; my $sta = 0; my $end = 0; my $sna = ""; my $Num = 0;
	foreach $key (keys %inf){
		$seq = ${$inf{$key}}[1];
		$subseq = substr($seq, 5, 6);
		if(exists $kmer_seq{$subseq}){
			$kmer_seq{$subseq} = $kmer_seq{$subseq} + 1;
			push(@{$kmer_loc{$subseq}}, $seq);
		}else{
			$kmer_seq{$subseq} = 1;
			$kmer_loc{$subseq} = [$seq];
		}
	}
	foreach $key ( sort{$kmer_seq{$b} <=> $kmer_seq{$a}} keys %kmer_seq){
		push(@sent, $key);
	}
	#print $#sent+1,"\n";
	$key = $sent[0];
	@sen1 = @{$kmer_loc{$key}};
	#print $#sen1+1,"\n";
	for($i=0; $i<=$#sen1; $i++){
		$Num = $Num + 1;
		@sen2 = split(//, $sen1[$i]);
		for($j=0; $j<=$#sen2; $j++){
			if($sen2[$j] eq "P"){
				$psign1[$j] = $psign1[$j] + 1;
			}else{
				$usign1[$j] = $usign1[$j] + 1;
			}
		}
	}
	for($r=1; $r<=$#sent; $r++){
		$sna = $sent[$r];
		#print $key,"\t",$sna,"\n";
		my ($sna2, $flag) = &clusterchar1($key, $sna);
		if($flag == 1){
			@sen3 = split(//, $sna2);
			$sta = 0;
			for($i=0; $i<=3; $i++){
				if($sen3[$i] eq "-"){
					$sta = $sta + 1;
				}
			}
			@sen1 = @{$kmer_loc{$sna}};
			#print $#sen1+1,"\n";
			$sta = $sta - 2;
			for($i=0; $i<=$#sen1; $i++){
				$Num = $Num + 1;
				@sen2 = split(//, $sen1[$i]);
				for($j=0; $j<=$#sen2; $j++){
					#$j = $j + $sta;
					if((0<=$j + $sta)&&($j + $sta<=$#psign1)){
						if($sen2[$j + $sta] eq "P"){
							$psign1[$j] = $psign1[$j] + 1;
						}else{
							$usign1[$j] = $usign1[$j] + 1;
						}
					}
				}
			}
		}
	}
	return (\@psign1, \@usign1, $Num);
}

sub build_str_motif{
	my $rep_seq = shift; my $cont = shift; my $Kmer_loc = shift; my $Data_Inf = shift; my $protein_name = shift;
	my $flank_len = 5;
	my %kmer_loc = %{$Kmer_loc}; my %data_inf = %{$Data_Inf};
	my %minf1 = (); my %minf2 = (); my %str_inf = ();
	my @cons = @{$cont}; my @sen1 = (); my @sen2 = ();
	my $i = 0; my $j = 0; my $r = 0; my $k = 0; my $num = 0; my $num1 = 0; my $num2 = 0;
	my @psign1 = (0)x(length($rep_seq) + 2*$flank_len); my @usign1 = (0)x(length($rep_seq) + 2*$flank_len);
	my @psign2 = (0)x(length($rep_seq) + 2*$flank_len); my @usign2 = (0)x(length($rep_seq) + 2*$flank_len);
	#open(OUT1, ">", "strtmp/".$rep_seq."str_test.txt");
	for($i=0; $i<=$#cons; $i++){
		#print $cons[$i],"|\t";
		@sen1 = @{$kmer_loc{$cons[$i]}};
		for($j=0; $j<=$#sen1; $j++){
			#print $sen1[$j],"\t";
			$num = $num + 1;
			#print $num,"\t";
			@sen2 = split(/\_/, $sen1[$j]);
			if(!exists $data_inf{$sen2[0]}){
				print $cons[$i],"|",$sen2[0],"\n";
				next;
			}
			my @ics = split(/\|/, ${$data_inf{$sen2[0]}}[3]);
			my @psent = split(/\|/, ${$data_inf{$sen2[0]}}[4]);
			my $max_sign = max(@psent); my $min_sign = min(@psent);
			my $ave_ics = 0; my $n_ics = 0;
			
			my $sta = $sen2[1]; my $end = $sen2[2]; my $key = $sen2[0]; 
			my $str_char = &get_str(${$data_inf{$key}}[1], ${$data_inf{$key}}[3], $protein_name);
			my $str_char1 = substr($str_char, $sta - 5, $end - $sta + 1 + 10);

			my $flag = 0; my $fnum = ""; 
			if(length($str_char1) == 16){
				my @str_sent = split(//, $str_char1);
				for($r=0; $r<=$#str_sent; $r++){
					if($str_sent[$r] eq "."){
						$fnum = $fnum."U";
					}else{
						$fnum = $fnum."P";
					}
				}
				$str_inf{$num} = [$str_char1, $fnum];
			}
		}
		#last;
	}
	#close OUT1;
	#print "strcuture site: ",$num,"\n";
	#$k = keys %str_inf;
	#print $k,"\n";
	my ($Psign1, $Usign1, $NUM) = kmer_cal2(\%str_inf);
	@psign1 = @{$Psign1}; @usign1 = @{$Usign1};
	for($i=0; $i<=$#psign1; $i++){
		#$minf{$i} = [$psign[$i], $usign[$i]];
		if($psign1[$i] + $usign1[$i] != 0){
			#my @prob = ($psign[$i]/($psign[$i] + $usign[$i]), $usign[$i]/($psign[$i] + $usign[$i]));
			#my $Height = &cal_entropy(\@prob);
			$minf1{$i} = [$psign1[$i]/($psign1[$i] + $usign1[$i]), $usign1[$i]/($psign1[$i] + $usign1[$i])];
			#$minf{$i} = [$psign[$i]/($psign[$i] + $usign[$i])*$Height, $usign[$i]/($psign[$i] + $usign[$i])*$Height];
		}else{
			$minf1{$i} = [0.5, 0.5];
			#$minf{$i} = [0, 0];
		}
	}
	for($i=0; $i<=$#psign2; $i++){
		#$minf{$i} = [$psign[$i], $usign[$i]];
		if($psign2[$i] + $usign2[$i] != 0){
			#my @prob = ($psign[$i]/($psign[$i] + $usign[$i]), $usign[$i]/($psign[$i] + $usign[$i]));
			#my $Height = &cal_entropy(\@prob);
			$minf2{$i} = [$psign2[$i]/($psign2[$i] + $usign2[$i]), $usign2[$i]/($psign2[$i] + $usign2[$i])];
			#$minf{$i} = [$psign[$i]/($psign[$i] + $usign[$i])*$Height, $usign[$i]/($psign[$i] + $usign[$i])*$Height];
		}else{
			$minf2{$i} = [0.5, 0.5];
			#$minf{$i} = [0, 0];
		}
	}
	return (\%minf1, \%minf2, $num1, $NUM);
}


sub cal_entropy{
	my $sen = shift;
	my @sent = @{$sen};
	my $i = 0; my $Sum = sum(@sent); my $Entropy = 0;
	for($i=0; $i<=$#sent; $i++){
		$sent[$i] = $sent[$i]/$Sum;
		if($sent[$i] > 0){
			$Entropy = $Entropy - $sent[$i]*log($sent[$i])/log(2);
		}
	}
	return(1 - $Entropy);
}

sub get_str{
	my $seqref = shift; my $ics = shift; my $protien = shift;
	my $i; my $j; my $r;
	my @sen1 = split(//, $seqref); my @sen2 = split(/\|/, $ics);
	my $tmp_seq_file = $protien."_tmp_seq_file.txt"; my $tmp_shape_file = $protien."_tmp_shape_file.txt"; 
	open(SEQ, ">", $tmp_seq_file);
	open(SHAPE, ">", $tmp_shape_file);
	for($i=0; $i<=$#sen1; $i++){
		print SEQ $sen1[$i];
	}
	print SEQ "\n";
	close SEQ;
	$j = 1;
	for($i = 0; $i<=$#sen2; $i++){
		$j = $i + 1;
		if($sen2[$i] < 0){
			print SHAPE $j,"\t-1\n";
		}else{
			print SHAPE $j,"\t",$sen2[$i]*2,"\n";
		}
	}
	close SHAPE;
	#my $str_res = `/Share2/home/zhangqf/usr/ViennaRNA-2.2.3/bin/RNAfold --noPS −−shapeMethod="Dm8b−0.7" --shape=tmp_shape_file.txt < tmp_seq_file.txt`;
	my $str_res = `RNAfold --noPS −−shapeMethod="Dm8b−0.7" --shape=$tmp_shape_file < $tmp_seq_file`;
	my @sent1 = split(/\n/,$str_res);
	my $exa_seq = $sent1[0];
	my @sent2 = split(/\s/,$sent1[1]); 
	my $exa_str = $sent2[0];
	#my @sent3 = split(/\|/, $sna);
	#$inf_str{$sna} = [$sent3[3], $exa_seq, $exa_str];
	return($exa_str);
}

sub get_str2{
	my $seqref = shift; my $ics = shift; my $protien = shift;
	my %inf_seq = %{$seqref}; my %inf_ics = %{$ics}; 
	my %inf_str = ();
	my $sen; my $sna; my $sen1; my $seq;
	my $i; my $j; my $r; my $count = 0;
	my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = ();
	foreach $sna (keys %inf_seq){
		@sen1 = @{$inf_seq{$sna}};
		@sen2 = @{$inf_ics{$sna}};
		my $tmp_seq_file = $protien."_tmp_seq_file.txt"; my $tmp_shape_file = $protien."_tmp_shape_file.txt";
		open(SEQ, ">", $tmp_seq_file);
		open(SHAPE, ">", $tmp_shape_file);
		for($i=0; $i<=$#sen1; $i++){
			print SEQ $sen1[$i];
		}
		print SEQ "\n";
		close SEQ;
		$j = 1;
		for($i = 0; $i<=$#sen2; $i++){
			if($sen2[$i] eq "NULL"){
				print SHAPE $j,"\t-1\n";
			}else{
				print SHAPE $j,"\t",$sen2[$i]*2,"\n";
			}
		}
		close SHAPE;
		#my $str_res = `/Share2/home/zhangqf/usr/ViennaRNA-2.2.3/bin/RNAfold --noPS −−shapeMethod="Dm8b−0.7" --shape=tmp_shape_file.txt < tmp_seq_file.txt`;
		my $str_res = `RNAfold --noPS −−shapeMethod="Dm8b−0.7" --shape=$tmp_shape_file < $tmp_seq_file`;
		my @sent1 = split(/\n/,$str_res);
		my $exa_seq = $sent1[0];
		my @sent2 = split(/\s/,$sent1[1]); 
		my $exa_str = $sent2[0];
		my @sent3 = split(/\|/, $sna);
		$inf_str{$sna} = [$sent3[3], $exa_seq, $exa_str];
	}
	return (\%inf_str);
}

sub mismatch{
	my $char1 = shift; my $char2 = shift;
	my @sen1 = split(//, $char1); my @sen2 = split(//, $char2); my $i = 0; my $r = 0;
	for($i=0; $i<=$#sen1; $i++){
		if($sen1[$i] ne $sen2[$i]){
			$r = $r + 1;
		}
	}
	return($r);
}

exit;