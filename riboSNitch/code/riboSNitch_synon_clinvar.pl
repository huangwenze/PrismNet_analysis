#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;


my $outfile = $ARGV[0];

my $usage = "This script is to search the riboSNitch associated with clinvar mutation
usage: $0 <outfile>
";
die $usage if $#ARGV<0;

my $ref_ribosnitch_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/riboSNitch/snp_str/";
my $ref_snp_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/SNP_trx_site/";
#HEK293T_SNPs_trx.bed
#Hela_SNPs_trx.bed

my $ref_trx = "/150T/zhangqf2/huangwz/Gencode/hg38_transcriptome_std_simple.fa";
my $trx_seq = &get_trx_seq($ref_trx);
my %seql = %{$trx_seq};

#chr10	1000772	1000772	+	G|A|303.77	1000773	1000772	ENST00000360803.8	1212	1211
#chr10	100150741	100150741	-	C|T|1356.77	100150742	100150741	ENST00000421367.6	5145	5144

#my $ref_trx = "/Share2/home/zhangqf/huangwz/Gencode/gencode.v26.transcripts.std.fa";
#my $ref_trx_corr = "/Share2/home/zhangqf/huangwz/Gencode/hg38_trx_corr.txt";
my $ref_trx_corr = "/150T/zhangqf2/huangwz/Gencode/hg38.transCoor.utr.bed";

my $high_signal_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/high_signal_region5/";

my $broad_bind_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/broad_binding_region/";
#my $Trx_seq = &get_trx_seq($ref_trx);
#my %sequ = %{$Trx_seq};
#my $Trx_cor = &get_trx_inf($ref_trx_corr);
#my %corr = %{$Trx_cor};

my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $file1 = ""; my $file2 = ""; my $file3 = ""; my $file4 = "";
my %inf = (); my %prot = (); my %ginf = (); my %func = (); my %clipdb = ();
my %inf1 = (); my %inf2 = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $len = 100; my $flag = 0;
my $sta = 0; my $end = 0; $len = 0; my $posi = 0; my $bscore = 0; my $cell_line = "";

my @score_region = ("seq1_5_6","seq1_6_7","seq1_7_8","seq1_8_9","seq1_9_10","str2_5_6","str2_6_7","str2_7_8","str2_8_9","str2_9_10");


#my $ribosNitch_synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/disease_SNP/riboSNitch/Synon_SNP/";
#my $nonstr_synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/disease_SNP/background/Synon_SNP/";

my $ribosNitch_synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/ribosnitch/Synon_SNP/";
my $nonstr_synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/background/Synon_SNP/";

my $ExAC_single_file = "/150T/zhangqf2/huangwz/ExAC/ExAC.singleton.bed";
my $ExAC_common_file = "/150T/zhangqf2/huangwz/ExAC/ExAC.common.bed";

my $clinvar_snp_info_file = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/Disease/Disease/Clinvar/clinvar_20190504.vcf";

my $SNP_inf = &get_snp_inf($clinvar_snp_info_file);
my %snp_inf = %{$SNP_inf};


my %minf1 = (); my %tinf1 = (); my %minf2 = (); my %tinf2 = ();

opendir(TEMPDIR, $ribosNitch_synon_path) or die "can't open it:$ribosNitch_synon_path";
#HEK293T_HEK293_dystrsite_synon.txt

#HEK293T_HEK293_nonstrsite_synon.txt

my @Dir = readdir TEMPDIR;
my @fa1 = grep /dystrsite_synon\.txt$/, @Dir;
closedir TEMPDIR;

open(OUT, ">", $outfile);
for($i=0; $i<=$#fa1; $i++){
	@sen1 = split(/_/, $fa1[$i]); 
	my $infile1 = $ribosNitch_synon_path.$fa1[$i];
	my $infile2 = $nonstr_synon_path.$sen1[0]."_".$sen1[1]."_nonstrsite_synon.txt";
	open(FILE1, $infile1)||die("open $infile1 error!\n");
	my $n1 = 0; my $n2 = 0; my $n3 = 0; my $n4 = 0;
	%minf1 = (); %minf2 = ();
	while($sen = <FILE1>){
		#ENST00000341421	1763	1784	1.43194552969	ENST00000341421|1768|C|T|1719.77|C|T|445.77|HEK293T|HEK293	ENST00000341421|1768|C|T|chr3|42209771|42209771|+	Synon
		#ENST00000341421	1763	1784	1.43194552969	ENST00000341421|1768|C|T|1719.77|C|T|445.77|HEK293T|HEK293	ENST00000341421|1768|C|T|chr3|42209771|42209771|+	Synon
		chomp($sen);
		@sent1 = split(/\t/, $sen);
		@sent2 = split(/\|/, $sent1[4]);
		@sent3 = split(/\|/, $sent1[5]);
		#$sid = $sent3[4]."|".$sent3[5]."|".$sent2[2]."|".$sent2[3];
		#$minf{$sid} = 1;
		$sid = $sent2[0]."|".$sent2[1]."|".$sent2[2]."|".$sent2[3];
		$minf1{$sid} = $sent3[4]."|".$sent3[6]."|".$sent2[2]."|".$sent2[3];
		$tinf1{$sid} = $sent3[4]."|".$sent3[6]."|".$sent2[2]."|".$sent2[3];
	}
	close FILE1;
	open(FILE1, $infile2)||die("open $infile2 error!\n");
	while($sen = <FILE1>){
		#ENST00000354574	2695	2695	ENST00000354574|2695|T|C|266.77|H9	ENST00000354574|2696|A|G|chr12|109932493|109932493|-	Synon
		#ENST00000490537	792	792	ENST00000490537|792|A|C|422.77|A|C|691.77|H9|HEK293T	ENST00000490537|792|A|C|chr13|37003588|37003588|+	Synon
		chomp($sen);
		@sent1 = split(/\t/, $sen);
		@sent2 = split(/\|/, $sent1[3]);
		@sent3 = split(/\|/, $sent1[4]);
		#$sid = $sent3[4]."|".$sent3[5]."|".$sent2[2]."|".$sent2[3];
		#$minf{$sid} = 1;
		$sid = $sent2[0]."|".$sent2[1]."|".$sent2[2]."|".$sent2[3];
		$minf2{$sid} = $sent3[4]."|".$sent3[6]."|".$sent2[2]."|".$sent2[3];
		$tinf2{$sid} = $sent3[4]."|".$sent3[6]."|".$sent2[2]."|".$sent2[3];
	}
	close FILE1;
	my $RECORD = &get_SNV_clinvar_bind(\%minf1, \%minf2, $SNP_inf);
	my %record = %{$RECORD};
	foreach $key (sort {$a cmp $b} keys %record){
		my %tmp1 = %{$record{$key}};
		foreach my $ke1 (sort {$a cmp $b} keys %tmp1){
			print OUT $sen1[0],"\t",$sen1[1],"\t",$key,"\t",$ke1,"\t",${$tmp1{$ke1}}[0],"\t",${$tmp1{$ke1}}[1],"\t",${$tmp1{$ke1}}[2],"\t",${$tmp1{$ke1}}[3],"\n";
		}
	}
	#last;
}

my %tinf3 = ();
foreach $key (keys %tinf2){
	if(!exists $tinf1{$key}){
		$tinf3{$key} = $tinf2{$key};
	}
}

my $RECORD = &get_SNV_clinvar_bind(\%tinf1, \%tinf3, $SNP_inf);
my %record = %{$RECORD};
foreach $key (sort {$a cmp $b} keys %record){
	my %tmp1 = %{$record{$key}};
	foreach my $ke1 (sort {$a cmp $b} keys %tmp1){
		print OUT "Total1\tTotal2\t",$key,"\t",$ke1,"\t",${$tmp1{$ke1}}[0],"\t",${$tmp1{$ke1}}[1],"\t",${$tmp1{$ke1}}[2],"\t",${$tmp1{$ke1}}[3],"\n";
	}
}

close OUT;


sub get_SNV_clinvar_bind{
	my $str_SNV = shift; my $unstr_SNV = shift; my $SNP_INF = shift;
	my %tinf1 = %{$str_SNV}; my %tinf2 = %{$unstr_SNV}; my %snp_inf = %{$SNP_INF};
	#my @region1 = ("seq1_5_6","seq1_6_7","seq1_7_8","seq1_8_9","seq1_9_10","str2_5_6","str2_6_7","str2_7_8","str2_8_9","str2_9_10");
	my @region1 = ("only_seq1_5_6","only_seq1_6_7","only_seq1_7_8","only_seq1_8_9","only_seq1_9_10","seq1_str2_5_6","seq1_str2_6_7","seq1_str2_7_8","seq1_str2_8_9","seq1_str2_9_10","only_str2_5_6","only_str2_6_7","only_str2_7_8","only_str2_8_9","only_str2_9_10");
	my @region01 = ("only_seq1_5_6","only_seq1_6_7","only_seq1_7_8","only_seq1_8_9","only_seq1_9_10");
	my @region02 = ("seq1_str2_5_6","seq1_str2_6_7","seq1_str2_7_8","seq1_str2_8_9","seq1_str2_9_10");
	my @region03 = ("only_str2_5_6","only_str2_6_7","only_str2_7_8","only_str2_8_9","only_str2_9_10");	
	my %ninf1 = (); my %ninf2 = (); my %ninf3 = (); my %ninf4 = ();
	my $n1 = 0; my $n2 = 0; my $n3 = 0; my $n4 = 0;
	my %record = ();
	foreach my $sid (keys %tinf1){
		my @sent1 = split(/\|/, $sid);
		if(exists $snp_inf{$tinf1{$sid}}){
			$n1 = $n1 + 1;
			if(exists $ninf1{$sent1[0]}){
				push(@{$ninf1{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf1{$sent1[0]} = [$sent1[1]];
			}
		}else{
			$n2 = $n2 + 1;
			if(exists $ninf2{$sent1[0]}){
				push(@{$ninf2{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf2{$sent1[0]} = [$sent1[1]];
			}
		}
	}
	foreach my $sid (keys %tinf2){
		my @sent1 = split(/\|/, $sid);
		if(exists $snp_inf{$tinf2{$sid}}){
			$n3 = $n3 + 1;
			if(exists $ninf3{$sent1[0]}){
				push(@{$ninf3{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf3{$sent1[0]} = [$sent1[1]];
			}
		}else{
			$n4 = $n4 + 1;
			if(exists $ninf4{$sent1[0]}){
				push(@{$ninf4{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf4{$sent1[0]} = [$sent1[1]];
			}
		}
	}
	#$record{"Total"} = [$n1, $n2, $n3, $n4];
	#print OUT "Total1\tTotal2\t",$n1,"\t",$n2,"\t",$n3,"\t",$n4,"\n";
	
	my $high_signal_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/high_signal_region5/";
	my $high_signal_div_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis/high_signal/high_signal_region5/HRR_total/";
	#for(my $r=0; $r<=$#region1; $r=$r+5){
	#my $high_signal_div_path = $high_signal_path."/".$region1[$r];
	my @sent2 = split(/\_/, "HRR_total"); #only_seq1_5_6
	my @sent3 = ("HRR_total");
	opendir(TEMPDIR, $high_signal_div_path) or die "can't open it:$high_signal_div_path";
	#ALKBH5_HEK293_HepG2_high_region.txt
	#ALKBH5_HEK293_K562_high_region.txt
	my @Dir01 = readdir TEMPDIR;
	my @fa01 = grep /_K562_high_region\.txt$/, @Dir01;
	closedir TEMPDIR;
	for(my $k=0; $k<=$#fa01; $k++){
		my @filen = split(/_/, $fa01[$k]); 
		my ($r1, $r2) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf1);
		my ($r3, $r4) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf2);
		my ($r5, $r6) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf3);
		my ($r7, $r8) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf4);
		#print OUT "Total1\t",$region1[$r],"\t",$r1,"\t",$r3,"\t",$r5,"\t",$r7,"\n";
		my $sna1 = $sent2[0]."_".$sent2[1];
		my $sna2 = $filen[0]."_".$filen[1];
		${$record{$sna1}}{$sna2} = [$r1, $r3, $r5, $r7];
	}
	#last;
	#print OUT "Total2\t",$score_region[$r],"\t",$r2,"\t",$r4,"\t",$r6,"\t",$r8,"\n";
	#}
	#my ($r1, $r2) = &overlap_high_signal_region3($high_signal_path, \%ninf1);
	#my ($r3, $r4) = &overlap_high_signal_region3($high_signal_path, \%ninf2);
	#my ($r5, $r6) = &overlap_high_signal_region3($high_signal_path, \%ninf3);
	#my ($r7, $r8) = &overlap_high_signal_region3($high_signal_path, \%ninf4);
	#$record{"unbind"} = [$n1 - $r1, $n2 - $r3, $n3 - $r5, $n4 - $r7];
	return(\%record);
}

sub get_SNV_clinvar_bind01{
	my $str_SNV = shift; my $unstr_SNV = shift; my $SNP_INF = shift;
	my %tinf1 = %{$str_SNV}; my %tinf2 = %{$unstr_SNV}; my %snp_inf = %{$SNP_INF};
	#my @region1 = ("seq1_5_6","seq1_6_7","seq1_7_8","seq1_8_9","seq1_9_10","str2_5_6","str2_6_7","str2_7_8","str2_8_9","str2_9_10");
	my @region1 = ("only_seq1_5_6","only_seq1_6_7","only_seq1_7_8","only_seq1_8_9","only_seq1_9_10","seq1_str2_5_6","seq1_str2_6_7","seq1_str2_7_8","seq1_str2_8_9","seq1_str2_9_10","only_str2_5_6","only_str2_6_7","only_str2_7_8","only_str2_8_9","only_str2_9_10");
	my @region01 = ("only_seq1_5_6","only_seq1_6_7","only_seq1_7_8","only_seq1_8_9","only_seq1_9_10");
	my @region02 = ("seq1_str2_5_6","seq1_str2_6_7","seq1_str2_7_8","seq1_str2_8_9","seq1_str2_9_10");
	my @region03 = ("only_str2_5_6","only_str2_6_7","only_str2_7_8","only_str2_8_9","only_str2_9_10");	
	my %ninf1 = (); my %ninf2 = (); my %ninf3 = (); my %ninf4 = ();
	my $n1 = 0; my $n2 = 0; my $n3 = 0; my $n4 = 0;
	my %record = ();
	foreach my $sid (keys %tinf1){
		my @sent1 = split(/\|/, $sid);
		if(exists $snp_inf{$tinf1{$sid}}){
			$n1 = $n1 + 1;
			if(exists $ninf1{$sent1[0]}){
				push(@{$ninf1{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf1{$sent1[0]} = [$sent1[1]];
			}
		}else{
			$n2 = $n2 + 1;
			if(exists $ninf2{$sent1[0]}){
				push(@{$ninf2{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf2{$sent1[0]} = [$sent1[1]];
			}
		}
	}
	foreach my $sid (keys %tinf2){
		my @sent1 = split(/\|/, $sid);
		if(exists $snp_inf{$tinf2{$sid}}){
			$n3 = $n3 + 1;
			if(exists $ninf3{$sent1[0]}){
				push(@{$ninf3{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf3{$sent1[0]} = [$sent1[1]];
			}
		}else{
			$n4 = $n4 + 1;
			if(exists $ninf4{$sent1[0]}){
				push(@{$ninf4{$sent1[0]}}, $sent1[1]);
			}else{
				$ninf4{$sent1[0]} = [$sent1[1]];
			}
		}
	}
	#$record{"Total"} = [$n1, $n2, $n3, $n4];
	#print OUT "Total1\tTotal2\t",$n1,"\t",$n2,"\t",$n3,"\t",$n4,"\n";
	
	for(my $r=0; $r<=$#region1; $r=$r+5){
		my $high_signal_div_path = $high_signal_path."/".$region1[$r];
		my @sent2 = split(/\_/, $region1[$r]); #only_seq1_5_6
		my @sent3 = @region1[$r..($r+4)];
		opendir(TEMPDIR, $high_signal_div_path) or die "can't open it:$high_signal_div_path";
		#ALKBH5_HEK293_HepG2_high_region.txt
		#ALKBH5_HEK293_K562_high_region.txt
		my @Dir01 = readdir TEMPDIR;
		my @fa01 = grep /_K562_high_region\.txt$/, @Dir01;
		closedir TEMPDIR;
		for(my $k=0; $k<=$#fa01; $k++){
			my @filen = split(/_/, $fa01[$k]); 
			my ($r1, $r2) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf1);
			my ($r3, $r4) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf2);
			my ($r5, $r6) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf3);
			my ($r7, $r8) = &overlap_high_signal_region($high_signal_path, \@sent3, $fa01[$k], \%ninf4);
			#print OUT "Total1\t",$region1[$r],"\t",$r1,"\t",$r3,"\t",$r5,"\t",$r7,"\n";
			my $sna1 = $sent2[0]."_".$sent2[1];
			my $sna2 = $filen[0]."_".$filen[1];
			${$record{$sna1}}{$sna2} = [$r1, $r3, $r5, $r7];
		}
		#last;
		#print OUT "Total2\t",$score_region[$r],"\t",$r2,"\t",$r4,"\t",$r6,"\t",$r8,"\n";
	}
	#my ($r1, $r2) = &overlap_high_signal_region3($high_signal_path, \%ninf1);
	#my ($r3, $r4) = &overlap_high_signal_region3($high_signal_path, \%ninf2);
	#my ($r5, $r6) = &overlap_high_signal_region3($high_signal_path, \%ninf3);
	#my ($r7, $r8) = &overlap_high_signal_region3($high_signal_path, \%ninf4);
	#$record{"unbind"} = [$n1 - $r1, $n2 - $r3, $n3 - $r5, $n4 - $r7];
	return(\%record);
}


sub site_in_region{
	my $Site = shift; my $INF = shift;
	my %site = %{$Site}; my %inf1 = %{$INF};
	my $flag = 0; my $i = 0; my $j = 0; my $r = 0;
	my $n1 = 0; my $n2 = 0;
	foreach my $ke1 (keys %site){
		my @sent1 = @{$site{$ke1}};
		if(exists $inf1{$ke1}){
			for($r=0; $r<=$#sent1; $r++){
				my @sent2 = @{$inf1{$ke1}};
				$flag = 0;
				for($j=0; $j<=$#sent2; $j++){
					my @sent3 = split(/\_/, $sent2[$j]);
					if(($sent3[0] <= $sent1[$r])&&($sent1[$r] <= $sent3[1])){
						$flag = 1;
						last;
					}
				}
				if($flag == 1){
					$n1 = $n1 + 1;
				}
			}
		}
		$n2 = $n2 + $#sent1 + 1;
	}
	return($n1, $n2);
}


sub overlap_high_signal_region{
	my $Path1 = shift; my $Reglist = shift; my $filen = shift; my $mut_site = shift;
	my @reglist = @{$Reglist};
	my @sent2 = split(/_/, $filen);
	my %sinf = ();	
	for(my $i=0; $i<=$#reglist; $i++){
		my $file1 = $Path1."/".$reglist[$i]."/".$filen;
		#my $file2 = $Path1."/".$reglist[$i]."/".$sent2[0]."_".$sent2[1]."_K562_high_region.txt";
	
		#print $Path1,"\t",$#Dir2,"\t",$Dir2[0],"\t",$Dir2[1],"\t",$Dir2[2],"\t",$Dir2[3],"\t",$#fa2,"\n";
		if(-e $file1){
			open(FILE1, $file1)||die("open $file1 error!\n");
			#ENST00000529473_2564_2573
			#ENST00000617316_107_116
			while($sen = <FILE1>){
				chomp($sen);
				my @sent1 = split(/\t/, $sen);
				if(exists $sinf{$sent1[0]}){
					push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
				}else{
					$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
				}
			}
			close FILE1;
		}
		#if(-e $file2){
		#	open(FILE1, $file2)||die("open $file2 error!\n");
		#	#ENST00000529473_2564_2573
		#	#ENST00000617316_107_116
		#	while($sen = <FILE1>){
		#		chomp($sen);
		#		my @sent1 = split(/\_/, $sen);
		#		if(exists $sinf{$sent1[0]}){
		#			push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
		#		}else{
		#			$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
		#		}
		#	}
		#	close FILE1;
		#}
	}
	my ($n1, $n2) = &site_in_region($mut_site, \%sinf);
	return($n1, $n2);
}


sub overlap_high_signal_region1{
	my $Path1 = shift; my $filen = shift; my $mut_site = shift;
	my @sent2 = split(/_/, $filen); 
	my $file1 = $Path1."/".$filen;
	my $file2 = $Path1."/".$sent2[0]."_".$sent2[1]."_K562_high_region.txt";
	my %sinf = ();
	#print $Path1,"\t",$#Dir2,"\t",$Dir2[0],"\t",$Dir2[1],"\t",$Dir2[2],"\t",$Dir2[3],"\t",$#fa2,"\n";
	open(FILE1, $file1)||die("open $file1 error!\n");
	#ENST00000529473_2564_2573
	#ENST00000617316_107_116
	while($sen = <FILE1>){
		chomp($sen);
		my @sent1 = split(/\_/, $sen);
		if(exists $sinf{$sent1[0]}){
			push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
		}else{
			$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
		}
	}
	close FILE1;
	if(-e $file2){
		open(FILE1, $file2)||die("open $file2 error!\n");
		#ENST00000529473_2564_2573
		#ENST00000617316_107_116
		while($sen = <FILE1>){
			chomp($sen);
			my @sent1 = split(/\_/, $sen);
			if(exists $sinf{$sent1[0]}){
				push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
			}else{
				$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
			}
		}
		close FILE1;
	}
	my ($n1, $n2) = &site_in_region($mut_site, \%sinf);
	return($n1, $n2);
}


sub overlap_high_signal_region2{
	my $Path1 = shift; my $mut_site = shift;
	my %sinf = ();
	#my @region1 = ("seq1_5_6","seq1_6_7","seq1_7_8","seq1_8_9","seq1_9_10","str2_5_6","str2_6_7","str2_7_8","str2_8_9","str2_9_10");
	opendir(TEMPDIR, $Path1) or die "can't open it:$Path1";
	#ALKBH5_HEK293_HepG2_bindsite.txt
	#ALKBH5_HEK293_K562_bindsite.txt
	my @Dir = readdir TEMPDIR;
	my @fa1 = grep /bindsite\.txt$/, @Dir;
	closedir TEMPDIR;
	
	foreach my $key (@fa1){
		my $file1 = $Path1."/".$key;
		#print $Path1,"\t",$#Dir2,"\t",$Dir2[0],"\t",$Dir2[1],"\t",$Dir2[2],"\t",$Dir2[3],"\t",$#fa2,"\n";
		open(FILE1, $file1)||die("open $file1 error!\n");
		#ENST00000579197	4801	4901	0.561
		#ENST00000438436	1061	1161	0.796
		while($sen = <FILE1>){
			chomp($sen);
			my @sent1 = split(/\t/, $sen);
			if(exists $sinf{$sent1[0]}){
				push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
			}else{
				$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
			}
		}
		close FILE1;
	}
	my ($n1, $n2) = &site_in_region($mut_site, \%sinf);
	return($n1, $n2);
}


sub overlap_high_signal_region3{
	my $Path1 = shift; my $mut_site = shift;
	my %sinf = ();
	my @region1 = ("seq1_5_6","seq1_6_7","seq1_7_8","seq1_8_9","seq1_9_10","str2_5_6","str2_6_7","str2_7_8","str2_8_9","str2_9_10");
	foreach my $Dir (@region1){
		my $file1 = $Path1."/".$Dir."_high_signal_region.txt";
		#print $Path1,"\t",$#Dir2,"\t",$Dir2[0],"\t",$Dir2[1],"\t",$Dir2[2],"\t",$Dir2[3],"\t",$#fa2,"\n";
		open(FILE1, $file1)||die("open $file1 error!\n");
		#ENST00000000233	27	120
		#ENST00000000233	169	280
		while($sen = <FILE1>){
			chomp($sen);
			my @sent1 = split(/\t/, $sen);
			if(exists $sinf{$sent1[0]}){
				push(@{$sinf{$sent1[0]}}, $sent1[1]."_".$sent1[2]);
			}else{
				$sinf{$sent1[0]} = [$sent1[1]."_".$sent1[2]];
			}
		}
		close FILE1;
	}
	my ($n1, $n2) = &site_in_region($mut_site, \%sinf);
	return($n1, $n2);
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
		my $sna = "chr".$sen1[0]."|".$sen1[1]."|".$sen1[3]."|".$sen1[4];
		$inf{$sna} = $sen1[3]."|".$sen1[4]."|".$sen1[7];
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
			$inf{$sen2[0]} = length($sen);
		}
	}
	close FILE3;
	return \%inf;
}

exit;

