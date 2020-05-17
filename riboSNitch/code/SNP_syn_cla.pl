#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

my $outfile = $ARGV[0];

my $usage = "This script is to filter synonymous mutation SNP
usage: $0 <outfile>
";
die $usage if $#ARGV<0;

my $ref_ribosnitch_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/ribosnitch_site/";
my $ref_snp_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/SNP_trx_site/";


my $ref_trx = "/Share2/home/zhangqf/huangwz/Gencode/gencode.v26.transcripts.std.fa";
#my $ref_trx_corr = "/Share2/home/zhangqf/huangwz/Gencode/hg38_trx_corr.txt";
my $ref_trx_corr = "/150T/zhangqf2/huangwz/Gencode/hg38.transCoor.utr.bed";

my $Trx_seq = &get_trx_seq($ref_trx);
my %sequ = %{$Trx_seq};
my $Trx_cor = &get_trx_inf($ref_trx_corr);
my %corr = %{$Trx_cor};

my $sen = ""; my $sen1 = ""; my $seq = ""; my $sna = ""; my $file = ""; my $file1 = ""; my $file2 = ""; my $file3 = ""; my $file4 = "";
my %inf = (); my %minf = (); my %prot = (); my %ginf = (); my %func = (); my %clipdb = ();
my %inf1 = (); my %inf2 = ();
my @sen = (); my @sen1 = (); my @sen2 = (); my @sen3 = (); my @sent1 = (); my @sent2 = (); my @sent3 = (); my @sent4 = (); my @sequ = ();
my $sid = ""; my $key; my $ics; my $ssid = "";
my $i = 0; my $j = 0; my $r = 0; my $num = 0; my $tnum = 0; my $k1; my $len = 100; my $flag = 0;
my $sta = 0; my $end = 0; $len = 0; my $posi = 0; my $bscore = 0; my $cell_line = "";

my $synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/ribosnitch/Synon_SNP/";
my $non_synon_path = "/150T/zhangqf2/huangwz/structure_motif/human_dataset6/ex_analysis2/ribosnitch2/ribosnitch/non_Synon_SNP/";

opendir(TEMPDIR, $ref_ribosnitch_path) or die "can't open it:$ref_ribosnitch_path";
#H9_HEK293T_dystrsite.txt

my @Dir = readdir TEMPDIR;
my @fa1 = grep /dystrsite_SNP\.txt$/, @Dir;
closedir TEMPDIR;

open(OUT, ">", $outfile);

for($i=0; $i<=$#fa1; $i++){
	%inf1 = (); %inf2 = ();
	@sen1 = split(/_/, $fa1[$i]); 
	print $fa1[$i],"\n";
	my $infile = $ref_ribosnitch_path.$fa1[$i];
	my $outfile1 = $synon_path.$sen1[0]."_".$sen1[1]."_dystrsite_synon.txt";
	my $outfile2 = $non_synon_path.$sen1[0]."_".$sen1[1]."_dystrsite_nonsynon.txt";
	$file1 = $ref_snp_path.$sen1[0]."_SNPs_trx.bed";
	$file2 = $ref_snp_path.$sen1[1]."_SNPs_trx.bed";
	
	open(FILE1, $file1)||die("open $file1 error!\n");
	while($sen = <FILE1>){
		#chr10	1000772	1000772	+	G|A|303.77	1000773	1000772	ENST00000360803.8	1212	1211
		#chr10	100150741	100150741	-	C|T|1356.77	100150742	100150741	ENST00000421367.6	5145	5144
		chomp($sen);
		@sent1 = split(/\t/, $sen);
		@sent2 = split(/\./, $sent1[7]);
		$inf1{$sent2[0]."|".$sent1[9]} = $sent1[0]."|".$sent1[1]."|".$sent1[2]."|".$sent1[3];
	}
	close FILE1;
	
	open(FILE2, $file2)||die("open $file2 error!\n");
	while($sen = <FILE2>){
		#chr10	1000772	1000772	+	G|A|303.77	1000773	1000772	ENST00000360803.8	1212	1211
		#chr10	100150741	100150741	-	C|T|1356.77	100150742	100150741	ENST00000421367.6	5145	5144
		chomp($sen);
		@sent1 = split(/\t/, $sen);
		@sent2 = split(/\./, $sent1[7]);
		$inf1{$sent2[0]."|".$sent1[9]} = $sent1[0]."|".$sent1[1]."|".$sent1[2]."|".$sent1[3];
	}
	close FILE2;
	
	open(FILE1, $infile)||die("open $infile error!\n");
	#H9_HEK293T_dystrsite_SNP.txt
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	#ENST00000341421	1763	1784	1.43194552969	ENST00000341421|1768|C|T|1719.77|C|T|445.77|HEK293T|HEK293
	open(OUT1, ">", $outfile1);
	open(OUT2, ">", $outfile2);
	while($sen = <FILE1>){
		chomp($sen);
		@sent1 = split(/\t/, $sen);
		@sent2 = split(/\|/, $sent1[4]);
		if(exists $sequ{$sent1[0]}){
			$seq = $sequ{$sent1[0]};
			my $ref_char = ""; my $mut_char = "";
			if(exists $inf1{$sent2[0]."|".$sent2[1]}){
				@sent3 = split(/\|/, $inf1{$sent2[0]."|".$sent2[1]});
				if($sent3[3] eq "+"){
					$posi = $sent2[1];
					$ref_char = $sent2[2];
					$mut_char = $sent2[3];
					#if($sent2[2] eq substr($seq, $sent2[1]-1, 1)){
					#	print OUT1 $sen,"\t",substr($seq, $sent2[1]-3, 5),"\t",$inf1{$sent2[0]."|".$sent2[1]},"\n";
					#}else{
					#	print OUT2 $sen,"\t",substr($seq, $sent2[1]-3, 5),"\t",$inf1{$sent2[0]."|".$sent2[1]},"\n";
					#}				
				}elsif($sent3[3] eq "-"){
					$posi = $sent2[1];
					$ref_char = &sup_char($sent2[2]);
					$mut_char = &sup_char($sent2[3]);				
					#if($sent2[2] eq &sup_char(substr($seq, $sent2[1], 1))){
					#	print OUT1 $sen,"\t",substr($seq, $sent2[1]-2, 5),"\t",$inf1{$sent2[0]."|".$sent2[1]},"\n";
					#}else{
					#	print OUT2 $sen,"\t",substr($seq, $sent2[1]-2, 5),"\t",$inf1{$sent2[0]."|".$sent2[1]},"\n";
					#}				
				}
				my @trx_corr = @{$corr{$sent2[0]}};
				if(($trx_corr[2] ne "protein_coding")||($posi <= $trx_corr[5])||($trx_corr[3] - $trx_corr[7] < $posi)){
					print OUT1 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tSynon\n";
					$inf2{$sent1[0]."|".$posi."|".$ref_char."|".$mut_char} = 1;
				}else{
					my $codep = ($posi - $trx_corr[5])%3;
					my $Char1 = ""; my $Char2 = ""; $flag = 0;
					if($codep == 1){
						$Char1 = substr($seq, $posi-1, 3);
						my @c1 = split(//, $Char1);
						if($c1[0] eq $ref_char){
							$c1[0] = $mut_char;
							$Char2 = join("", @c1);
						}else{
							print OUT2 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tERROR|",$c1[0],"|",$ref_char,"\n";
							$flag = 1;
						}
					}elsif($codep == 2){
						$Char1 = substr($seq, $posi-2, 3);
						my @c1 = split(//, $Char1);
						if($c1[1] eq $ref_char){
							$c1[1] = $mut_char;
							$Char2 = join("", @c1);
						}else{
							print OUT2 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tERROR|",$c1[1],"|",$ref_char,"\n";
							$flag = 1;
						}				
					}else{
						$Char1 = substr($seq, $posi-3, 3);
						my @c1 = split(//, $Char1);
						if($c1[2] eq $ref_char){
							$c1[2] = $mut_char;
							$Char2 = join("", @c1);
						}else{
							print OUT2 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tERROR|",$c1[2],"|",$ref_char,"\n";
							$flag = 1;
						}
					}
					if($flag == 0){
						if(&syn_muta($Char1, $Char2) == 1){
							print OUT1 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tSynon\n";
							$inf2{$sent1[0]."|".$posi."|".$ref_char."|".$mut_char} = 1;
						}elsif(&syn_muta($Char1, $Char2) == 2){
							print $sent1[0]."|".$posi."|".$ref_char."|".$mut_char,"\t",$Char1,"\t",$Char2,"\n";
						}else{
							print OUT2 $sen,"\t",$sent1[0],"|",$posi,"|",$ref_char,"|",$mut_char,"|",$inf1{$sent2[0]."|".$sent2[1]},"\tnon-Synon\n";
						}
					}
				}
			}
		}
	}
	close FILE1;
	close OUT1;
	close OUT2;
	$r = keys %inf2;
	print OUT $sen1[0],"\t",$sen1[1],"\t",$r,"\n";
	#last;
}
close OUT;

sub sup_char{
	my $Char1 = shift;
	if($Char1 eq "A"){
		return "T";
	}elsif($Char1 eq "T"){
		return "A";
	}elsif($Char1 eq "C"){
		return "G";
	}else{
		return "C";
	}
}

sub syn_muta{
	my $Char1 = shift; my $Char2 = shift;
	$Char1 =~s/T/U/g; $Char2 =~s/T/U/g;
	my %Amino = (
	"UUU" => "F", "UUC" => "F", "UUA" => "L", "UUG" => "L", "UCU" => "S", "UCC" => "S", "UCA" => "S", "UCG" => "S",
	"UAU" => "Y", "UAC" => "Y", "UAA" => "", "UAG" => "", "UGU" => "C", "UGC" => "C", "UGA" => "", "UGG" => "W",
	"CUU" => "L", "CUC" => "L", "CUA" => "L", "CUG" => "L", "CCU" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
	"CAU" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q", "CGU" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
	"AUU" => "I", "AUC" => "I", "AUA" => "I", "AUG" => "M", "ACU" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
	"AAU" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K", "AGU" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
	"GUU" => "V", "GUC" => "V", "GUA" => "V", "GUG" => "V", "GCU" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
	"GAU" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E", "GGU" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
	);
	if((!exists $Amino{$Char1})||(!exists $Amino{$Char2})){
		#print $Char1,"\t",$Char2,"\n";
		return(2);
	}
	if($Amino{$Char1} eq $Amino{$Char2}){
		return(1);
	}else{
		return(0);
	}
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
		#gene name, gene id, gene type, length, splicing, 3UTR, CDS, 5UTR
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