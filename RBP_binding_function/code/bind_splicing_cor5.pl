#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Statistics::Basic qw(:all);


my $cell_line = $ARGV[0]; 
my $RBP = $ARGV[1];
my $outfile1 = $ARGV[2];

my $usage = "This script is to correlation between RBP differential binding and alternative splicing of transcripts.
usage: $0 <cell_line> <RBP> <outfile1>
";
die $usage if $#ARGV<2;


my $trx_shape_file = "/150T/zhangqf2/huangwz/total_smart_icshape/new_smartSHAPE_0412/".$cell_line."_smartSHAPE.out";
my %trx_exp = (); 
open(FILE2, $trx_shape_file)||die("open $trx_shape_file error!\n");
#ENST00000318041.13	1728	*	NULL	NULL
while(my $sen = <FILE2>){
	chomp($sen);
	my @sent1 = split(/\t/, $sen);
	my @sent2 = split(/\./, $sent1[0]);
	$trx_exp{$sent2[0]} = 1;
}
close FILE2;

my $splicing_file1 = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/SRSF1_splicing/TRA2_data_raw/TRA2_events_".$cell_line."_trx.bed";

my $CLIP_path = "/150T/zhangqf2/huangwz/RBP_CLIP/human_trx_clip/";

my $path1 = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/data/";

my $pred_res_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/datagenome_pu/";

#my $path1_out_pre = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/new_figure3/model_pred_region/";
my $path1_out_pre = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/datagenome_pu_region/";

my $ref_trx_corr = "/150T/zhangqf2/huangwz/Gencode/hg38.transCoor.utr.bed";


my @cell = ("H9","HEK293","HEK293T","Hela","HepG2","K562");


opendir(TEMPDIR, $path1_out_pre) or die "can't open it:$path1_out_pre";
#ALKBH5_HEK293_H9_pu_predsite.txt
#ALKBH5_HEK293_HEK293T_pu_predsite.txt
my @Dir = readdir TEMPDIR;
#my $prex = "_HepG2_pu_predsite.txt";
#my @fa1 = grep /${prex}$/, @Dir;
my @fa1 = grep /_${cell_line}_pu_predsite\.txt$/, @Dir;
#my @fa2 = grep /_HepG2_K562_ics$/, @Dir;
#my @fa = (@fa1, @fa2);
closedir TEMPDIR;

my %hpred_file = ();
foreach my $ke (@fa1){
	my @sent1 = split(/\_/, $ke);
	$hpred_file{$sent1[0]} = $ke; 
}

opendir(TEMPDIR, $CLIP_path) or die "can't open it:$CLIP_path";
#CLIPDB20001_ALKBH5_HEK293_trx.bed

@Dir = readdir TEMPDIR;
#my @fa1 = grep /${RBP}/, @Dir;
@fa1 = grep /${cell_line}\_trx\.bed/, @Dir;
closedir TEMPDIR;

my %hclip_file = ();
foreach my $ke (@fa1){
	my @sent1 = split(/\_/, $ke);
	if(exists $hclip_file{$sent1[1]}){
		push(@{$hclip_file{$sent1[1]}}, $ke);
	}else{
		$hclip_file{$sent1[1]} = [$ke]; 
	}
}

#my $TE_inf = &read_TE_file($TE_file1, \%trx_exp);
#my %te_inf = %{$TE_inf};

my $SPLIC_inf = &read_splicing_file($splicing_file1, \%trx_exp);
my %splic_inf = %{$SPLIC_inf};


my $pred_bind_file = $path1_out_pre.$hpred_file{$RBP};

my $Pred_inf = &get_icbind_peak4($pred_bind_file); my %icbind_inf = %{$Pred_inf};

my ($cor1, $cor2, $Com_over1) = &overlap_reg01($Pred_inf, $SPLIC_inf);

my %com_over1 = %{$Com_over1};

open(OUT1, ">", $outfile1);
foreach my $ke (keys %com_over1){
	print OUT1 $ke,"\t",join("\t", @{$com_over1{$ke}}),"\n";
}
close OUT1;

print "Result: ",$RBP,"\t",$cell_line,"\t",$cor1,"\t",$cor2,"\n";




sub get_icbind_peak2{
	my $pred_sam = shift; my $icbind_peak = shift;
	my %inf = ();
	my $sen = ""; my $sen1 = ""; my $sen2 = "";
	open(FILE1, $pred_sam)||die("open $pred_sam error!\n");
	#A	ENST00000582829|241|341	CTACACTTAACTTTT
	#A	ENST00000582829|261|361	TTCGGGTTTAAGACG
	open(FILE2, $icbind_peak)||die("open $icbind_peak error!\n");
	#0.021
	#0.009
	my %minf = (); my %ginf = ();
	my $sid = ""; my $sta = 0; my $end = 0; my $bscore = 0;
	my $i = 0;
	while($sen2 = <FILE2>){
		$sen1 = <FILE1>;
		chomp($sen1);
		chomp($sen2);
		my @sent1 = split(/\t/, $sen1);
		my @sent2 = split(/\|/, $sent1[1]);
		my $sid = $sent2[0]; my $sta = $sent2[1]; my $end = $sent2[2]; my $bscore = $sen2;
		if($bscore > 0.5){
			if(exists $minf{$sid}){
				#$minf{$sid} = $minf{$sid} + $bscore;
				if($bscore > $minf{$sid}){
					$minf{$sid} = $bscore;
				}
			}else{
				$minf{$sid} = $bscore;
			}
		}
	}
	close FILE1;
	close FILE2;
	return (\%minf);
}


sub get_icbind_peak3{
	my $pred_sam = shift; #my $icbind_peak = shift;
	my %inf = ();
	my $sen = ""; my $sen1 = ""; my $sen2 = "";
	open(FILE1, $pred_sam)||die("open $pred_sam error!\n");
	#ENST00000582829	261	361	0.829
	#ENST00000582829	501	601	0.567
	#ENST00000582829	701	801	0.974
	my %minf = (); my %ginf = ();
	my $sid = ""; my $sta = 0; my $end = 0; my $bscore = 0;
	my $i = 0;
	while($sen1 = <FILE1>){
		chomp($sen1);
		my @sent2 = split(/\t/, $sen1);
		my $sid = $sent2[0]; my $sta = $sent2[1]; my $end = $sent2[2]; my $bscore = $sent2[3];
		if($bscore > 0.999){
			$bscore = 0.999;
		}elsif($bscore < 0.001){
			$bscore = 0.001;
		}
		if($bscore > 0.5){
			#$bscore = -log(1/$bscore - 1);
			if(exists $minf{$sid}){
				push(@{$minf{$sid}}, $bscore);
			}else{
				$minf{$sid} = [$bscore];
			}
		}
	}
	foreach my $ke (keys %minf){
		my @sent1 = @{$minf{$ke}};
		#$ginf{$ke} = sum(@sent1)/($#sent1 + 1);
		$ginf{$ke} = sum(@sent1);
		#$ginf{$ke} = max(@sent1);
	}
	close FILE1;
	return (\%ginf);
}


sub get_icbind_peak4{
	my $pred_sam = shift; #my $icbind_peak = shift;
	my %inf = ();
	my $sen = ""; my $sen1 = ""; my $sen2 = "";
	open(FILE1, $pred_sam)||die("open $pred_sam error!\n");
	#ENST00000582829	261	361	0.829
	#ENST00000582829	501	601	0.567
	#ENST00000582829	701	801	0.974
	my %minf = (); my %ginf = ();
	my $sid = ""; my $sta = 0; my $end = 0; my $bscore = 0;
	my $i = 0;
	while($sen1 = <FILE1>){
		chomp($sen1);
		my @sent2 = split(/\t/, $sen1);
		my $sid = $sent2[0]; my $sta = $sent2[1]; my $end = $sent2[2]; my $bscore = $sent2[3];
		if($bscore > 0.999){
			$bscore = 0.999;
		}elsif($bscore < 0.001){
			$bscore = 0.001;
		}
		if($bscore > 0.5){
			#$bscore = -log(1/$bscore - 1);
			if(exists $minf{$sid}){
				push(@{$minf{$sid}}, $sta, $end, $bscore);
			}else{
				$minf{$sid} = [$sta, $end, $bscore];
			}
		}
	}
	#foreach my $ke (keys %minf){
	#	my @sent1 = @{$minf{$ke}};
	#	#$ginf{$ke} = sum(@sent1)/($#sent1 + 1);
	#	$ginf{$ke} = sum(@sent1);
	#	#$ginf{$ke} = max(@sent1);
	#}
	close FILE1;
	return (\%minf);
}


sub overlap_reg01{
	my $Reg1 = shift; my $Reg2 = shift; 
	my %inf1 = %{$Reg1}; my %inf2 = %{$Reg2};
	my %com = (); my %spec = ();
	my @func1 = (); my @bind_inf1 = (); my @func2 = (); my @bind_inf2 = ();
	foreach my $ke (keys %inf1){
		if(exists $inf2{$ke}){
			my @sent1 = @{$inf1{$ke}}; my @sent2 = @{$inf2{$ke}};
			for(my $i=0; $i<=$#sent1; $i=$i+3){
				my $flag = 0;
				for(my $j=0; $j<=$#sent2; $j=$j+5){
					if( abs(($sent1[$i] + $sent1[$i+1])/2 - $sent2[$j]) < ($sent1[$i+1] - $sent1[$i])/2 + 100 ){
						#$flag = 1;
						my $sna = $ke."_".$sent1[$i]."_".$sent1[$i+1];
						$com{$sna} = [$sent1[$i+2], $ke."_".$sent2[$j]."_".$sent2[$j+1]."_".$sent2[$j+2], $sent2[$j+2], $sent2[$j+3], $sent2[$j+4], "5UTR"];
						push(@func1, ($sent2[$j+3]+$sent2[$j+4])/2);
						push(@bind_inf1, $sent1[$i+2]);
						#last;
					}elsif(abs(($sent1[$i] + $sent1[$i+1])/2 - $sent2[$j+1]) < ($sent1[$i+1] - $sent1[$i])/2 + 100){
						#$flag = 1;
						my $sna = $ke."_".$sent1[$i]."_".$sent1[$i+1];
						$com{$sna} = [$sent1[$i+2], $ke."_".$sent2[$j]."_".$sent2[$j+1]."_".$sent2[$j+2], $sent2[$j+2], $sent2[$j+3], $sent2[$j+4], "3UTR"];
						push(@func2, ($sent2[$j+3]+$sent2[$j+4])/2);
						push(@bind_inf2, $sent1[$i+2]);
						#last;
					}
				}
			}
		}
	}
	my $cor1 = correlation([@func1], [@bind_inf1]);
	my $cor2 = correlation([@func2], [@bind_inf2]);
	return($cor1, $cor2, \%com);
}


sub read_clip_file{
	my $clip_file = shift;
	my %inf = (); my %ginf = ();
	open(FILE1, $clip_file)||die("open $clip_file error!\n");
	while(my $sen = <FILE1>){
		#chr1	634517	634568	+	ENSG00000198744.5_0_8#AUH#eCLIP,CLIPper#ENCFF816BPV#K562#0.000134856683388#569918|38	634518	634568	ENST00000416718.2	143	193
		#chr1	827087	827117	-	ENSG00000225880.4_0_3#AUH#eCLIP,CLIPper#ENCFF816BPV#K562#0.0288697458512#762481|15	827088	827117	ENST00000473798.1	406	435
		chomp($sen);
		my @sent1 = split(/\t/, $sen);
		my @sent2 = split(/\|/, $sent1[4]);
		my @sent3 = split(/\./, $sent1[7]);
		$sent2[1] = log($sent2[1] + 1);
		if(exists $inf{$sent3[0]}){
			push(@{$inf{$sent3[0]}}, $sent2[1]);
		}else{
			$inf{$sent3[0]} = [$sent2[1]];
		}
	}
	close FILE1;
	foreach my $ke (keys %inf){
		my @sent1 = @{$inf{$ke}};
		#$ginf{$ke} = sum(@sent1)/($#sent1 + 1);
		$ginf{$ke} = sum(@sent1);
		#$ginf{$ke} = max(@sent1);
	}
	#my @Value = values %inf;
	#my $min = min(@Value); 
	return(\%ginf);
}


sub read_splicing_file{
	my $splice_file = shift; my $EXP = shift;
	my %exp_inf = %{$EXP};
	open(FILE2, $splice_file)||die("open $splice_file error!\n");
	#chrX	100636191	100636608	-	ENSG00000000003.14|A5|0.06015437353828982|0.05557552747633143	100636192	100636608	ENST00000496771.5	82	498
	#chrX	100636191	100636608	-	ENSG00000000003.14|A5|0.06015437353828982|0.05557552747633143	100636608	100636608	ENST00000373020.8	199	199
	my %tmp = (); my %exon = (); my %skip = (); my %inclu = ();
	my $numb1 = 0; my $numb2 = 0;
	while(my $sen1 = <FILE2>){
		chomp($sen1);
		my @sen2 = split(/\t/,$sen1);
		my @sent1 = split(/\|/,$sen2[4]);
		my @sent2 = split(/\./,$sen2[7]);
		if((($sen2[1]+1) ne $sen2[5]) || ($sen2[2] ne $sen2[6])){
			next;
		}
		if((exists $exp_inf{$sent2[0]})&&($sent1[2]=~/[0-9\-\.]+/)&&($sent1[3]=~/[0-9\-\.]+/)){
			my $sna = $sent2[0]."_".$sen2[8]."_".$sen2[9];
			$tmp{$sna} = [$sent1[1], $sent1[2], $sent1[3]];
			if(exists $exon{$sent2[0]}){
				push(@{$exon{$sent2[0]}}, $sen2[8], $sen2[9], $sent1[1], $sent1[2], $sent1[3]);
			}else{
				$exon{$sent2[0]} = [$sen2[8], $sen2[9], $sent1[1], $sent1[2], $sent1[3]];
			}
		}
	}
	close FILE2;
	return(\%exon);
}


sub read_TE_file{
	my $TE_file = shift; my $EXP = shift;
	my %inf = (); my %exp_inf = %{$EXP};
	open(FILE1, $TE_file)||die("open $TE_file error!\n");
	#ENST00000592274.1	0.197461112611
	#ENST00000216962.8	0.877877024834
	while(my $sen = <FILE1>){
		chomp($sen);
		my @sen1 = split(/\t/, $sen);
		my @sent1 = split(/\./, $sen1[0]);
		my $sna = $sent1[0];
		if(exists $exp_inf{$sna}){
			$inf{$sna} = $sen1[1];
		}
		#trx id, TE
	}
	close FILE1;
	return(\%inf);
}


sub read_HL_file{
	my $HL_file = shift;
	my %inf = ();
	open(FILE1, $HL_file)||die("open $HL_file error!\n");
	#"HEK293_RNA half-life (min) rep1" "HEK293_RNA half-life (min) rep2"
	#"A1BG" 91.4258577065151 67.5981915356002
	#"A2LD1" 45.2599993330771 85.2998924659824
	my $sen = <FILE1>;
	while($sen = <FILE1>){
		chomp($sen);
		$sen=~s/\"//g;
		my @sen1 = split(/\s/, $sen);
		#@sent1 = split(/\./, $sen1[0]);
		my $sna = $sen1[0];
		$inf{$sna} = ($sen1[1] + $sen1[2])/2;
		#gene name, half_life
	}
	close FILE1;
	return(\%inf);
}


sub get_trx_inf{
	my $trx_corr = shift;
	my $sen = "";
	my %inf = ();
	open(FILE1, $trx_corr)||die("open $trx_corr error!\n");
	#ENST00000332831.3	OR4F16=ENSG00000273547.1	protein_coding	939	1-939	937-939	0,939,3
	#ENST00000420190.6	SAMD11=ENSG00000187634.11	protein_coding	1578	1-1021,1022-1113,1114-1295,1296-1346,1347-1471,1472-1561,1562-1578	1-504	504,1578,0
	#ENST00000437963.5	SAMD11=ENSG00000187634.11	protein_coding	387	1-40,41-132,133-314,315-365,366-387	1-40,41-60	60,387,0
	#ENST00000618181.4	SAMD11=ENSG00000187634.11	protein_coding	2179	1-60,61-152,153-334,335-459,460-549,550-690,691-769,770-1269,1270-1394,1395-1505,1506-2179	1-60,61-80,1749-2179	80,2179,431
	while($sen = <FILE1>){
		chomp($sen);
		my @sen1 = split(/\t/, $sen);
		my @sent1 = split(/\./, $sen1[0]);
		my @sent2 = split(/\=/, $sen1[1]);
		my @sent3 = split(/\./, $sent2[1]);
		my @sent4 = split(/\,/, $sen1[-1]);
		my $sna = $sent1[0];
		$inf{$sna} = [$sent2[0], $sent3[0], $sen1[2], $sen1[3], $sen1[4], $sent4[0], $sent4[1], $sent4[2]];
		#trx id => gene name, gene id, gene type, length, splicing, 3UTR, CDS, 5UTR
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
