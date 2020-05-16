#!/usr/bin/perl -w
use strict;
use Math::Complex;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;


my $usage = "This script is to search the riboSNitches according to the SNP file and structurally structurally variable site reault.
usage: $0
";

my $ref_trx = "/150T/zhangqf2/huangwz/Gencode/hg38_transcriptome_std_simple.fa";
my $trx_corr = "/150T/zhangqf2/huangwz/Gencode/hg38_trx_corr.txt";


my $ref_snp_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/SNP_trx_site/";
my $path_dyn_str = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/dynamic_str_region05/combine_region/";


my $str_snp_out_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/ribosnitch_site/";
my $nonstr_snp_out_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/nonstr_snp/";



my $TRX_seq = &get_trx_seq($ref_trx);
my %inf_seq = %{$TRX_seq};


my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; 
my %inf = (); my %minf = (); my %ginf = (); my %func = (); my %clipdb = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $len = 100; #my $flag1 = 0; my $flag2 = 0;
my $sta = 0; my $end = 0; $len = 0; my $cell_line = "";

my @cell = ("HepG2","K562","Hela","HEK293T","HEK293","H9");

#my $trx_len = &get_trx_seq($ref_trx);
#my %inf_len = %{$trx_len};  my %uniq = (); my %uniq2 = (); 

opendir(TEMPDIR, $path_dyn_str) or die "can't open it:$path_dyn_str";
#H9_HEK293T_dystrsite.txt

my @Dir = readdir TEMPDIR;
my @fa1 = grep /dystrsite\.txt$/, @Dir;
closedir TEMPDIR;

open(OUT2, ">", "snp_str_overlap_number.txt");
#my $mp = Parallel::ForkManager->new($proc);
for($i=0; $i<=$#fa1; $i++){
	#$mp->start and next;
	my $file1 = $fa1[$i];
	@sent1 = split(/\_/, $file1);
	$file1 = $path_dyn_str.$file1;
	#$file3 = $path_dyn_str.$sent1[1]."_".$sent1[2]."_dystrsite.txt";
	my $file2 = $ref_snp_path.$sent1[0]."_SNPs_trxsite.bed";
	my $file3 = $ref_snp_path.$sent1[1]."_SNPs_trxsite.bed";
	my $file4 = $str_snp_out_path.$sent1[0]."_".$sent1[1]."_dystrsite_SNP.txt";
	my $file5 = $nonstr_snp_out_path.$sent1[0]."_".$sent1[1]."_nonstrsite_SNP.txt";
	my %dsite1 = (); my %dsite2 = (); my %dstr1 = ();
	open(OUT1, ">", $file4);
	open(OUT3, ">", $file5);
	open(FILE1, $file2)||die("open $file2 error!\n");
	while($sen = <FILE1>){
		#ENST00000491635	1629	1629	G|A|1233.77|0.500|chr10|1000771|1000772|+
		#ENST00000360803	1211	1211	G|A|1233.77|0.500|chr10|1000771|1000772|+
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		@sen2 = split(/\|/, $sen1[3]);
		$sna = $sen1[0]."|".$sen1[1]."|".$sen1[2];
		my $sna1 = $sen2[0]."|".$sen2[1]."|".$sen2[3];
		my $sna2 = $sen2[4]."|".$sen2[5]."|".$sen2[6]."|".$sen2[7];
		$dsite1{$sna} = [$sna2, $sna1, $sent1[0]];
	}
	close FILE1;
	open(FILE1, $file3)||die("open $file3 error!\n");
	while($sen = <FILE1>){
		#ENST00000468709	2149	2149	G|A|1233.77|0.500|chr10|1000771|1000772|+
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		@sen2 = split(/\|/, $sen1[3]);
		$sna = $sen1[0]."|".$sen1[1]."|".$sen1[2];
		my $sna1 = $sen2[0]."|".$sen2[1]."|".$sen2[3];
		my $sna2 = $sen2[4]."|".$sen2[5]."|".$sen2[6]."|".$sen2[7];
		if(exists $dsite1{$sna}){
			if($sna1 eq $dsite1{$sna}[1]){
				delete($dsite1{$sna});
			}else{
				$dsite1{$sna}[1] = $dsite1{$sna}[1]."|".$sna1;
				$dsite1{$sna}[2] = $dsite1{$sna}[2]."|".$sent1[1];
			}
			#$dsite1{$sna}[0] = $dsite1{$sna}[0]."|".$sen1[3];
			#$dsite1{$sna}[1] = $dsite1{$sna}[1]."|".$sent1[1];
		}else{
			$dsite1{$sna} = [$sna2, $sna1, $sent1[1]];
			#$dsite1{$sna} = [$sen1[3], $sent1[1]];
		}
	}
	foreach $key (keys %dsite1){
		@sen1 = split(/\|/, $key);
		if(exists $dsite2{$sen1[0]}){
			push(@{$dsite2{$sen1[0]}}, $sen1[1]);
		}else{
			$dsite2{$sen1[0]} = [$sen1[1]];
		}
	}
	close FILE1;
	
	my $n1 = 0; my $n2 = 0; my $n3 = 0; my $n4 = 0; my $n5 = keys %dsite1; my %dover = (); my $n6 = 0;
	open(FILE1, $file1)||die("open $file1 error!\n");
	#ENST00000341423	963	990	1.47073689013
	while($sen = <FILE1>){
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		my $flag = 0; my $flag2 = 0; my $out_str = "";
		$sta = int(rand(length($inf_seq{$sen1[0]}) - ($sen1[2] - $sen1[1])) + 1); $end = $sta + ($sen1[2] - $sen1[1]);
		if(exists $dsite2{$sen1[0]}){
			@sen2 = @{$dsite2{$sen1[0]}};
			for($j=0; $j<=$#sen2; $j++){
				if(($sen1[1] <= $sen2[$j]) && ($sen2[$j] <= $sen1[2])){
					$flag = 1;
					$sna = $sen1[0]."|".$sen2[$j]."|".$sen2[$j];
					print OUT1 $sen,"\t",$sen1[0],"|",$sen2[$j],"|",${$dsite1{$sna}}[1],"|",${$dsite1{$sna}}[2],"|",${$dsite1{$sna}}[0],"\n";
					$dover{$sna} = 1;
					#ENST00000341423	963	990	1.47073689013	ENST00000341423|980|T|C|48.77|K562	G|C##ENST00000635625_448_495_2.12679430129|ENST00000635625_401_561_H9|0.773
				}
				if(($sta <= $sen2[$j]) && ($sen2[$j] <= $end)){
					$flag2 = 1;
					#ENST00000341423	963	990	1.47073689013	ENST00000341423|980|T|C|48.77|K562	G|C##ENST00000635625_448_495_2.12679430129|ENST00000635625_401_561_H9|0.773
				}
			}
		}
		if($flag == 1){
			$n1 = $n1 + 1;
		}
		$n2 = $n2 + 1;
		if($flag2 == 1){
			$n3 = $n3 + 1;
		}
		$n4 = $n4 + 1;
	}
	$n6 = keys %dover;
	foreach $key (keys %dsite1){
		if(exists $dover{$key}){
		}else{
			@sen1 = split(/\|/, $key);
			print OUT3 $sen1[0],"\t",$sen1[1],"\t",$sen1[2],"\t",$sen1[0],"|",$sen1[1],"|",${$dsite1{$key}}[1],"|",${$dsite1{$key}}[2],"|",${$dsite1{$key}}[0],"\n";
		}
	}
	print OUT2 $sent1[0]."\t".$sent1[1],"\t",$n1,"\t",$n2,"\t",$n6,"\t",$n5,"\t",$n3,"\t",$n4,"\n";
	close FILE1;
	close OUT1;
	close OUT3;
	#$mp->finish;
	#last;
}
#$mp->wait_all_children;
close OUT2;

sub random_trx_regoin{
	my $trx_inf = shift; my $len = shift; my $number = shift; 
	my %tlen = %{$trx_inf};
	my %random_reg = ();
	my $sta = 0; my $end = 0; my $r = 0; 
	foreach my $key (keys %tlen){
		$r = $r + 1;
		$sta = int(rand($tlen{$key} - $len - 10) + 5);
		$end = $sta + $len - 1;
		$sna = $key."|".$sta."|".$end;
		$random_reg{$sna} = 1;
		if($r >= $number){
			last;
		}
	}
	return (\%random_reg);
}

sub get_exon_seq{
	my $trx_id = shift; my $sta = shift; my $end = shift; my $T_inf = shift; my $T_seq = shift;
	my %trx_inf = %{$T_inf}; my %trx_seq = %{$T_seq}; my $res = "";
	#ENST00000451998	315_439	ENSG00000075234.16|SE|0.2829463723|0.0014985015|chr22|46292203|46292328|+
	my @sen2 = @{$trx_inf{$trx_id}};
	my @sen3 = split(/\,/,$sen2[4]);
	$res = length($trx_seq{$trx_id})."\t".$sen2[3];
	for($i=1; $i<$#sen3; $i++){
		#@sent1 = split(/\_/,$sen1[1]);
		my @sent2 = split(/\-/,$sen3[$i]);
		#if(($sent2[0] < ($sta + $end)/2) && ($sent2[1] > ($sta + $end)/2)){
		if((abs($sent2[0] - $sta) < 3) && (abs($sent2[1] - $end) < 3)){
			$res = $res."\t".$sen3[$i-1].",".$sen3[$i].",".$sen3[$i+1];
			my @sent1 = split(/\-/,$sen3[$i-1]);
			$res = $res."\t".substr($trx_seq{$trx_id}, $sent1[0] - 1, $sent1[1] - $sent1[0] + 1);
			@sent1 = split(/\-/,$sen3[$i]);
			$res = $res."\t".substr($trx_seq{$trx_id}, $sent1[0] - 1, $sent1[1] - $sent1[0] + 1);
			@sent1 = split(/\-/,$sen3[$i+1]);
			$res = $res."\t".substr($trx_seq{$trx_id}, $sent1[0] - 1, $sent1[1] - $sent1[0] + 1);
		}
	}
	return ($res);
}

sub read_splice_data{
	my $splice_file = shift; my $T_inf = shift;
	my $sen1 = ""; my %trx_inf = %{$T_inf};
	my @sen2 = (); my @sent1 = (); my @sent2 = (); my $flag = 0;
	open(FILE2, $splice_file)||die("open $splice_file error!\n");
	#chr1    100345025       100345183       +       ENSG00000079335.18|AF|0.36626841619999995|0.03996004    100345026       100345183       ENST00000635056.1       2       159
	#chr1    100345025       100345183       +       ENSG00000079335.18|AF|0.36626841619999995|0.03996004    100345026       100345183       ENST00000446055.1       251     408
	#chr1    100351734       100351808       +       ENSG00000079335.18|AF|0.24619165420000003|0.0455794206  100351735       100351808       ENST00000455467.5       2       75
	#chr1    10400376        10400393        +       ENSG00000142657.20|A3|0.0810068099|0.0157342657 10400377        10400393        ENST00000477958.5       131     147
	#chr1    10400376        10400393        +       ENSG00000142657.20|A3|0.0810068099|0.0157342657 10400393        10400393        ENST00000465632.5       208     208
	my %tmp = (); my %skip = (); my %inclu = ();
	my $numb1 = 0; my $numb2 = 0;
	while($sen1 = <FILE2>){
		chomp($sen1);
		@sen2 = split(/\t/,$sen1);
		@sent1 = split(/\|/,$sen2[4]);
		@sent2 = split(/\./,$sen2[7]);
		if((($sen2[1]+1) ne $sen2[5]) || ($sen2[2] ne $sen2[6]) || ($sent1[1] ne "SE")){
			next;
		}

		$sna = $sent2[0]."_".$sen2[8]."_".$sen2[9];
		$tmp{$sna} = [$sent1[2], $sent1[3], $sen2[4]."|".$sen2[0]."|".$sen2[1]."|".$sen2[2]."|".$sen2[3]];
		if(($sent1[2] < 0)&&($sent1[3] < 0.05)){
			$numb1 = $numb1 + 1;
			if(exists $skip{$sent2[0]}){
				push(@{$skip{$sent2[0]}}, $sen2[8], $sen2[9]);
			}else{
				$skip{$sent2[0]} = [$sen2[8], $sen2[9]];
			}
		}elsif(($sent1[2] > 0)&&($sent1[3] < 0.05)){
			$numb2 = $numb2 + 1;
			if(exists $inclu{$sent2[0]}){
				push(@{$inclu{$sent2[0]}}, $sen2[8], $sen2[9]);
			}else{
				$inclu{$sent2[0]} = [$sen2[8], $sen2[9]];
			}
		}
	}
	print $numb1,"\t",$numb2,"\n";
	close FILE2;
	return(\%tmp, \%skip, \%inclu);
}

sub get_snp_inf{
	my $trx_corr = shift;
	my $sen = "";
	my %inf = ();
	open(FILE1, $trx_corr)||die("open $trx_corr error!\n");
	#1	1014042	475283	G	A	.	.	AF_ESP=0.00546;AF_EXAC=0.00165;AF_TGP=0.00619;
	while($sen = <FILE1>){
		chomp($sen);
		my @sen1 = split(/\t/, $sen);
		my $sna = "chr".$sen1[0]."|".$sen1[1]."|".$sen1[1];
		$inf{$sna} = $sen1[3]."|".$sen1[4]."|".$sen1[7];
	}
	close FILE1;
	return (\%inf);
}

sub get_trx_inf{
	my $trx_corr = shift;
	my $sen = "";
	my %inf = ();
	open(FILE1, $trx_corr)||die("open $trx_corr error!\n");
	#ENST00000456328.2	DDX11L1=ENSG00000223972.5	processed_transcript	1657	1-359,360-468,469-1657	1	1657
	while($sen = <FILE1>){
		chomp($sen);
		my @sen1 = split(/\t/, $sen);
		my @sent1 = split(/\./, $sen1[0]);
		my @sent2 = split(/\=/, $sen1[1]);
		my @sent3 = split(/\./, $sent2[1]);
		my $sna = $sent1[0];
		$inf{$sna} = [$sent2[0], $sent3[0], $sen1[2], $sen1[3], $sen1[4], $sen1[5], $sen1[6]];
		#gene name, gene id, gene type, length, splicing, 3UTR, 5UTR
	}
	close FILE1;
	return (\%inf);
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