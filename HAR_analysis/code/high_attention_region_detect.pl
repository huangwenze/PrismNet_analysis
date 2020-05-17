#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Parallel::ForkManager;


my $usage = "This script is to get high attention region from attention signal file.
usage: $0
";


my $pred_data_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/pre_data/";
#H9_win.txt
my $pred_saliency_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/datagenome61/";

my $pred_high_region = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/high_signal_region5/";

my $sen = ""; my $sen1 = ""; my $sen2 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $file1 = ""; my $file2 = ""; my $file3 = ""; my $file4 = "";
my %inf = (); my %prot = (); my %ginf = (); my %func = (); my %clipdb = ();
#my %inf1 = (); my %inf2 = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $len = 100; my $flag = 0;
my $sta = 0; my $end = 0; $len = 0; my $posi = 0; my $bscore = 0; my $cell_line = "";


opendir(TEMPDIR, $pred_saliency_path) or die "can't open it:$pred_saliency_path";
#ALKBH5_HEK293_5v_binary_-1_binary_icbind_1_pu_ana61_HepG2.txt
#ALKBH5_HEK293_5v_binary_-1_binary_icbind_1_pu_ana61_K562.txt
#ALKBH5_HEK293_5v_binary_-1_binary_icbind_1_seq_ana61_HepG2.txt
#ALKBH5_HEK293_5v_binary_-1_binary_icbind_1_seq_ana61_K562.txt

my @Dir = readdir TEMPDIR;
my @fa1 = grep /_pu_ana61_/, @Dir;
closedir TEMPDIR;

@fa1 = sort {$a cmp $b} @fa1;
my $proc = 40;
my $mp = Parallel::ForkManager->new($proc);
for($i=0; $i<=$#fa1; $i++){
	$mp->start and next;
	@sen1 = split(/\./, $fa1[$i]); 
	@sen2 = split(/_/, $sen1[0]); 
	my $infile1 = $pred_saliency_path.$fa1[$i];
	my $infile2 = $pred_data_path.$sen2[-1]."_win.txt";
	my $HHdir = $pred_high_region."str2_9_10/";
	`mkdir -p ${HHdir}`;
	my $outfile3 = $HHdir.$sen2[0]."_".$sen2[1]."_".$sen2[-1]."_high_region.txt";
	$HHdir = $pred_high_region."str2_8_9/";
	`mkdir -p ${HHdir}`;
	my $outfile4 = $HHdir.$sen2[0]."_".$sen2[1]."_".$sen2[-1]."_high_region.txt";
	
	$HHdir = $pred_high_region."seq1_9_10/";
	`mkdir -p ${HHdir}`;
	my $outfile1 = $HHdir.$sen2[0]."_".$sen2[1]."_".$sen2[-1]."_high_region.txt";
	$HHdir = $pred_high_region."seq1_8_9/";
	`mkdir -p ${HHdir}`;
	my $outfile2 = $HHdir.$sen2[0]."_".$sen2[1]."_".$sen2[-1]."_high_region.txt";
	
	my ($INF1, $INF2) = &bind_site_seq_str($infile1, $infile2, 20, 0.9, 1.0);
	my %inf1 = %{$INF1}; my %inf2 = %{$INF2};
	open(OUT1, ">", $outfile1);
	foreach $key (sort {$a cmp $b} keys %inf1){
		print OUT1 $key,"\t",${$inf1{$key}}[0],"\n";
	}
	close OUT1;
	open(OUT1, ">", $outfile3);
	foreach $key (sort {$a cmp $b} keys %inf2){
		print OUT1 $key,"\t",${$inf2{$key}}[0],"\n";
	}
	close OUT1;
	($INF1, $INF2) = &bind_site_seq_str($infile1, $infile2, 20, 0.8, 0.9);
	%inf1 = %{$INF1}; %inf2 = %{$INF2};
	open(OUT2, ">", $outfile2);
	foreach $key (sort {$a cmp $b} keys %inf1){
		print OUT2 $key,"\t",${$inf1{$key}}[0],"\n";
	}
	close OUT2;
	open(OUT2, ">", $outfile4);
	foreach $key (sort {$a cmp $b} keys %inf2){
		print OUT2 $key,"\t",${$inf2{$key}}[0],"\n";
	}
	close OUT2;
	#if($i >= 1){
	#	last;
	#}
	#last;
	$mp->finish;
}
$mp->wait_all_children;

sub bind_site_seq_str{
	my $infile1 = shift; my $infile2 = shift; my $len = shift; my $minscore = shift; my $maxscore = shift;
	#my ($score2) = &max_per_str($infile1, $len, $per, $minscore, $maxscore);
	my %inf1 = (); my %inf2 = (); my %site_inf1 = (); my %site_inf2 = ();
	my $key = ""; my $sen = "";
	open(FILE1, $infile1)||die("open $infile1 error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		$key = $sent1[0];
		my @seqsent = (); my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			for(my $i=0; $i<=100; $i++){
				my $num = max($sent2[$i*5], $sent2[$i*5+1], $sent2[$i*5+2], $sent2[$i*5+3]);
				push(@seqsent, $num);
				push(@strsent, $sent2[$i*5+4]);			
			}
			my $i = 0; my $j = 0; my $index = 0; my $max1 = -1; my $max2 = -1; my $sum1 = 0; my $sum2 = 0; 
			my $sta1 = -1; my $end1 = -1; my $sta2 = -1; my $end2 = -1; my $tsum1 = 0; my $tsum2 = 0;
			for($i=0; $i<=$#strsent-$len+1; $i++){
				$sum1 = 0; $sum2 = 0;
				for($j=0; $j<$len; $j++){
					$sum2 = $sum2 + $strsent[$i+$j];
				}
				if($sum2 > $max2){
					$max2 = $sum2;
					$sta2 = $i;
					$end2 = $i + $len - 1;
				}
				for($j=0; $j<$len; $j++){
					$sum1 = $sum1 + $seqsent[$i+$j];
				}
				if($sum1 > $max1){
					$max1 = $sum1;
					$sta1 = $i;
					$end1 = $i + $len - 1;
				}
			}
			#my $sna = $key."_".$sta;
			#$site_inf{$sna} = [$end, $tsum1, $tsum2];
			if($max1 > 0){
				${$site_inf1{$key}}{$sta1} = [$end1, $max1];
			}
			if($max2 > 0){
				${$site_inf2{$key}}{$sta2} = [$end2, $max2];
			}
		}
	}
	close FILE1;
	open(FILE2, $infile2)||die("open $infile2 error!\n");
	#A       ENST00000582829|241|341 CTACACTTAACTTTTCCGGCTTCGGGTTTA
	my $i = 0;
	while($sen = <FILE2>){
		if(exists $site_inf1{$i}){
			my @sent1 = split(/\t/, $sen);
			my @sent2 = split(/\|/, $sent1[1]);
			foreach $key (keys %{$site_inf1{$i}}){
				#$inf{$sent1[1]} = $site_inf{$i};
				my @sent3 = @{${$site_inf1{$i}}{$key}};
				my $sna = $sent2[0]."_".($sent2[1] + $key)."_".($sent2[1] + $sent3[0]);
				$inf1{$sna} = [$sent3[1]];
			}
			foreach $key (keys %{$site_inf2{$i}}){
				#$inf{$sent1[1]} = $site_inf{$i};
				my @sent3 = @{${$site_inf2{$i}}{$key}};
				my $sna = $sent2[0]."_".($sent2[1] + $key)."_".($sent2[1] + $sent3[0]);
				$inf2{$sna} = [$sent3[1]];
			}
		}
		$i = $i + 1;
	}
	close FILE2;
	return(\%inf1, \%inf2);
}


sub bind_site_str{
	my $infile1 = shift; my $infile2 = shift; my $len = shift; my $per = shift; my $minscore = shift; my $maxscore = shift;
	my ($score2) = &max_per_str($infile1, $len, $per, $minscore, $maxscore);
	my %inf = (); my %site_inf = ();
	my $key = ""; my $sen = "";
	open(FILE1, $infile1)||die("open $infile1 error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		$key = $sent1[0];
		my @seqsent = (); my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			for(my $i=0; $i<=100; $i++){
				#my $num = max($sent2[$i*4], $sent2[$i*4+1], $sent2[$i*4+2], $sent2[$i*4+3]);
				#push(@seqsent, $num);
				push(@strsent, $sent2[$i*5+4]);			
			}
		}
		my $i = 0; my $j = 0; my $index = 0; my $maxn = 0; my $sum1 = 0; my $sum2 = 0; my $sta = -1; my $end = -1; my $tsum1 = 0; my $tsum2 = 0;
		for($i=0; $i<=$#strsent-$len+1; $i++){
			$sum1 = 0; $sum2 = 0;
			for($j=0; $j<$len; $j++){
				#$sum1 = $sum1 + $seqsent[$i+$j];
				$sum2 = $sum2 + $strsent[$i+$j];
			}
			if($sum2 > $score2){
				if($end >= $i){
					$end = $i + $len - 1;
				}elsif($end > 0){
					my $sna = $key."_".$sta;
					#$site_inf{$sna} = [$end, $tsum1, $tsum2];
					${$site_inf{$key}}{$sta} = [$end, $tsum2];
					$sta = $i;
					$end = $i + $len - 1;
					$tsum2 = $sum2;
				}else{
					$sta = $i;
					$end = $i + $len - 1;
					$tsum2 = $sum2;
				}
			}
		}
		if($end > 0){
			my $sna = $key."_".$sta;
			#$site_inf{$sna} = [$end, $tsum1, $tsum2];
			${$site_inf{$key}}{$sta} = [$end, $tsum2];
		}
	}
	close FILE1;
	open(FILE2, $infile2)||die("open $infile2 error!\n");
	#A       ENST00000582829|241|341 CTACACTTAACTTTTCCGGCTTCGGGTTTA
	my $i = 0;
	while($sen = <FILE2>){
		if(exists $site_inf{$i}){
			my @sent1 = split(/\t/, $sen);
			my @sent2 = split(/\|/, $sent1[1]);
			foreach $key (keys %{$site_inf{$i}}){
				#$inf{$sent1[1]} = $site_inf{$i};
				my @sent3 = @{${$site_inf{$i}}{$key}};
				my $sna = $sent2[0]."_".($sent2[1] + $key)."_".($sent2[1] + $sent3[0]);
				$inf{$sna} = [$sent3[1]];
			}
		}
		$i = $i + 1;
	}
	close FILE2;
	return(\%inf);
}

sub bind_site_seq{
	my $infile1 = shift; my $infile2 = shift; my $len = shift; my $per = shift; my $minscore = shift; my $maxscore = shift;
	my ($score1) = &max_per_seq($infile1, $len, $per, $minscore, $maxscore);
	my %inf = (); my %site_inf = ();
	my $key = ""; my $sen = "";
	open(FILE1, $infile1)||die("open $infile1 error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		$key = $sent1[0];
		my @seqsent = (); #my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			for(my $i=0; $i<=100; $i++){
				my $num = max($sent2[$i*4], $sent2[$i*4+1], $sent2[$i*4+2], $sent2[$i*4+3]);
				push(@seqsent, $num);
				#push(@strsent, $sent2[$i*5+4]);			
			}
		}
		my $i = 0; my $j = 0; my $index = 0; my $maxn = 0; my $sum1 = 0; my $sum2 = 0; my $sta = -1; my $end = -1; my $tsum1 = 0; my $tsum2 = 0;
		for($i=0; $i<=$#seqsent-$len+1; $i++){
			$sum1 = 0; $sum2 = 0;
			for($j=0; $j<$len; $j++){
				$sum1 = $sum1 + $seqsent[$i+$j];
				#$sum2 = $sum2 + $strsent[$i+$j];
			}
			if($sum1 > $score1){
				if($end >= $i){
					$end = $i + $len - 1;
				}elsif($end > 0){
					my $sna = $key."_".$sta;
					#$site_inf{$sna} = [$end, $tsum1, $tsum2];
					${$site_inf{$key}}{$sta} = [$end, $tsum1];
					$sta = $i;
					$end = $i + $len - 1;
					$tsum1 = $sum1;
				}else{
					$sta = $i;
					$end = $i + $len - 1;
					$tsum1 = $sum1;
				}
			}
		}
		if($end > 0){
			my $sna = $key."_".$sta;
			#$site_inf{$sna} = [$end, $tsum1, $tsum2];
			${$site_inf{$key}}{$sta} = [$end, $tsum1];
		}
	}
	close FILE1;
	open(FILE2, $infile2)||die("open $infile2 error!\n");
	#A       ENST00000582829|241|341 CTACACTTAACTTTTCCGGCTTCGGGTTTA
	my $i = 0;
	while($sen = <FILE2>){
		if(exists $site_inf{$i}){
			my @sent1 = split(/\t/, $sen);
			my @sent2 = split(/\|/, $sent1[1]);
			foreach $key (keys %{$site_inf{$i}}){
				#$inf{$sent1[1]} = $site_inf{$i};
				my @sent3 = @{${$site_inf{$i}}{$key}};
				my $sna = $sent2[0]."_".($sent2[1] + $key)."_".($sent2[1] + $sent3[0]);
				$inf{$sna} = [$sent3[1]];
			}
		}
		$i = $i + 1;
	}
	close FILE2;
	return(\%inf);
}

sub max_per_seq{
	my $infile = shift; my $len = shift; my $per = shift; my $minscore = shift; my $maxscore = shift;
	#my $list = shift; my $len = shift;
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = "";
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my @total1 = (); my @total2 = ();
	open(FILE1, $infile)||die("open $infile error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	$sna = 0;
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		my @seqsent = (); my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			for($i=0; $i<=100; $i++){
				my $num = max($sent2[$i*4], $sent2[$i*4+1], $sent2[$i*4+2], $sent2[$i*4+3]);
				push(@seqsent, $num);
				#push(@strsent, $sent2[$i*5+4]);			
			}
		}
		for($i=0; $i<=$#seqsent-$len+1; $i++){
			push(@total1, sum(@seqsent[$i..($i+$len-1)]));
		}
	}
	close FILE1;
	@total1 = sort {$b <=> $a} @total1;
	my $boun1 = $total1[int(($#total1 + 1)*$per)-1];
	return ($boun1);
}

sub max_per_str{
	my $infile = shift; my $len = shift; my $per = shift; my $minscore = shift; my $maxscore = shift;
	#my $list = shift; my $len = shift;
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = "";
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my @total1 = (); my @total2 = ();
	open(FILE1, $infile)||die("open $infile error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	$sna = 0;
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		my @seqsent = (); my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			for($i=0; $i<=100; $i++){
				#my $num = max($sent2[$i*5], $sent2[$i*5+1], $sent2[$i*5+2], $sent2[$i*5+3]);
				#push(@seqsent, $num);
				push(@strsent, $sent2[$i*5+4]);			
			}
		}
		for($i=0; $i<=$#strsent-$len+1; $i++){
			push(@total2, sum(@strsent[$i..($i+$len-1)]));
		}
	}
	close FILE1;
	#@total1 = sort {$b <=> $a} @total1;
	@total2 = sort {$b <=> $a} @total2;
	#my $boun1 = $total1[int(($#total1 + 1)*$per)-1];
	my $boun2 = $total2[int(($#total2 + 1)*$per)-1];
	return ($boun2);
}

sub max_per_seq_str{
	my $infile = shift; my $len = shift; my $per = shift; my $minscore = shift; my $maxscore = shift;
	#my $list = shift; my $len = shift;
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = "";
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my @total1 = (); my @total2 = ();
	open(FILE1, $infile)||die("open $infile error!\n");
	#184     0.812   0.000,0.000,0.000,0.000,-0.000,0.000,0.000,0.000,0.000,
	$sna = 0;
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		my @seqsent = (); my @strsent = ();
		if(($minscore < $sent1[1])&&($sent1[1] <= $maxscore)){
			my @sent2 = split(/\,/, $sent1[2]);
			
			for($i=0; $i<=100; $i++){
				my $num = max($sent2[$i*5], $sent2[$i*5+1], $sent2[$i*5+2], $sent2[$i*5+3]);
				push(@seqsent, $num);
				push(@strsent, $sent2[$i*5+4]);			
			}
		}
		for($i=0; $i<=$#seqsent-$len+1; $i++){
			push(@total1, sum(@seqsent[$i..($i+$len-1)]));
		}
		for($i=0; $i<=$#strsent-$len+1; $i++){
			push(@total2, sum(@strsent[$i..($i+$len-1)]));
		}
	}
	close FILE1;
	@total1 = sort {$b <=> $a} @total1;
	@total2 = sort {$b <=> $a} @total2;
	my $boun1 = $total1[int(($#total1 + 1)*$per)-1];
	my $boun2 = $total2[int(($#total2 + 1)*$per*2)-1];
	return ($boun1, $boun2);
}

exit;