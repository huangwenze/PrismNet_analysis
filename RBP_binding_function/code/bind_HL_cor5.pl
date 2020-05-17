#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;
use Statistics::Basic qw(:all);


my $HL_file1 = $ARGV[0];
my $RBP = $ARGV[1];
my $clip_file = $ARGV[2];
my $outfile1 = $ARGV[3];

my $usage = "This script is to correlation between RBP binding and half life of transcripts.
usage: $0 <clip_file1> <cell> <clip_file> <outfile1>
";
die $usage if $#ARGV<3;


my $trx_shape_file = "/150T/zhangqf2/huangwz/total_smart_icshape/new_smartSHAPE_0412/HEK293_smartSHAPE.out";
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


my $CLIP_record_file = "/150T/zhangqf2/huangwz/RBP_CLIP/human.RBP.CLIP.Piranha.HEK293.combined.tbl";

$HL_file1 = "/150T/zhangqf2/huangwz/structure_motif/functional/Half_life/".$HL_file1;

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
my $prex = "_HEK293_pu_predsite.txt";
my @fa1 = grep /_HEK293_pu_predsite\.txt$/, @Dir;
closedir TEMPDIR;

my %hpred_file = ();
foreach my $ke (@fa1){
	my @sent1 = split(/\_/, $ke);
	$hpred_file{$sent1[0]} = $ke; 
}


my %hclip_file = ();
open(FILE1, $CLIP_record_file)||die("open $CLIP_record_file error!\n");
#CLIPDB20002	ALKBH5	PAR-CLIP,Piranha	GSE38201,GSM936506	HEK293
#CLIPDB20012	C17ORF85	PAR-CLIP,Piranha	GSE38201,GSM936507	HEK293

while(my $sen = <FILE1>){
	chomp($sen);
	my @sent1 = split(/\t/, $sen);
	my $file1 = $sent1[0]."_".$sent1[1]."_".$sent1[4]."_trx.bed";
	if(exists $hclip_file{$sent1[1]}){
		push(@{$hclip_file{$sent1[1]}}, $file1);
	}else{
		$hclip_file{$sent1[1]} = [$file1]; 
	}	
}
close FILE1;

my $TRX_INF = &get_trx_inf($ref_trx_corr);

my $TE_inf = &read_HL_file($HL_file1, $TRX_INF, \%trx_exp);
my %te_inf = %{$TE_inf};

open(OUT2, ">", $outfile1);

my $pred_bind_file = $path1_out_pre.$hpred_file{$RBP};
my $Pred_inf = &get_icbind_peak3($pred_bind_file); my %icbind_inf = %{$Pred_inf};
my $clip_file1 = $CLIP_path.$clip_file; 
my $CLIP_inf = &read_clip_file($clip_file1); my %clip_inf = %{$CLIP_inf}; 

my @func_score1 = (); my @func_score2 = (); my @func_score3 = (); my @func_score4 = ();

foreach my $ke1 (keys %te_inf){
	if(exists $clip_inf{$ke1}){
		push(@func_score1, log($te_inf{$ke1}));
		print OUT2 $RBP,"\tCLIP\t",$clip_file,"\t",log($te_inf{$ke1}),"\t",$clip_inf{$ke1},"\n";
	}else{
		push(@func_score2, log($te_inf{$ke1}));
		print OUT2 $RBP,"\tCLIP\t",$clip_file,"\t",log($te_inf{$ke1}),"\t-1\n";
	}
}

foreach my $ke2 (keys %te_inf){
	if(exists $icbind_inf{$ke2}){
		push(@func_score3, log($te_inf{$ke2}));
		print OUT2 $RBP,"\ticbind\t",$clip_file,"\t",log($te_inf{$ke2}),"\t",$icbind_inf{$ke2},"\n";
	}else{
		push(@func_score4, log($te_inf{$ke2}));
		print OUT2 $RBP,"\ticbind\t",$clip_file,"\t",log($te_inf{$ke2}),"\t-1\n";
	}
}

print $RBP,"\taverage_compare\t",$#func_score1+1,"\t",sum(@func_score1)/($#func_score1+1),"\t",$#func_score2+1,"\t",sum(@func_score2)/($#func_score2+1),"\t",$#func_score3+1,"\t",sum(@func_score3)/($#func_score3+1),"\t",$#func_score4+1,"\t",sum(@func_score4)/($#func_score4+1),"\n";

close OUT2;

sub get_icbind_peak{
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
		if((exists $minf{$sid})&&($minf{$sid} < $bscore)){
			$minf{$sid} = $bscore;
		}elsif(! exists $minf{$sid}){
			$minf{$sid} = $bscore;
		}
	}
	close FILE1;
	close FILE2;
	return (\%minf);
}


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
			$bscore = -log(1/$bscore - 1);
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
		#$ginf{$ke} = sum(@sent1);
		$ginf{$ke} = max(@sent1);
	}
	close FILE1;
	return (\%ginf);
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
		#$ginf{$ke} = sum(@sent1);
		$ginf{$ke} = max(@sent1);
	}
	#my @Value = values %inf;
	#my $min = min(@Value); 
	return(\%ginf);
}

sub read_TE_file{
	my $TE_file = shift;
	my %inf = ();
	open(FILE1, $TE_file)||die("open $TE_file error!\n");
	#ENST00000592274.1	0.197461112611
	#ENST00000216962.8	0.877877024834
	while(my $sen = <FILE1>){
		chomp($sen);
		my @sen1 = split(/\t/, $sen);
		my @sent1 = split(/\./, $sen1[0]);
		my $sna = $sent1[0];
		$inf{$sna} = $sen1[1];
		#trx id, TE
	}
	close FILE1;
	return(\%inf);
}


sub read_HL_file{
	my $HL_file = shift; my $TRX_coor = shift; my $EXP = shift;
	my %inf = (); my %trx_coor = %{$TRX_coor}; my %ginf = (); my %exp_inf = %{$EXP};
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
	foreach my $ke (keys %trx_coor){
		my $sna = ${$trx_coor{$ke}}[0];
		if((exists $inf{$sna}) && (exists $exp_inf{$ke})){
			$ginf{$ke} = log($inf{$sna}+1);
		}
	}
	return(\%ginf);
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
