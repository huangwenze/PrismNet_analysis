#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

my $prot_cell = $ARGV[0];


my $usage = "This script is to compare the motif similarity of each RBP model and combine the similar motifs.
usage: $0 <prot_cell>

example: perl motif_summary.pl SRSF1_HepG2
";
die $usage if $#ARGV<0;

my $summary_file = $prot_cell."_summary.txt";
my $meme_file = $prot_cell."_motif_10_seq.meme";
my $seq_file = $prot_cell."_motif_10_seq.txt";
my $str_file = $prot_cell."_motif_10_str.txt";

my $tmeme_file = $prot_cell."_top10_motif_10_seq.meme";
my $tmeme_file2 = $prot_cell."_top10_motif_10_seq2.meme";

my $sum_out = $prot_cell."_combine_summary.txt";
my $seqstr_out = $prot_cell."_topmotif.txt";
#RBFOX2_mes_summary.txt
#RBFOX2_mes_motif_10_seq.txt
#RBFOX2_mes_motif_10_str.txt
#RBFOX2_mes_motif_10_seq.meme

`head -n 240 $meme_file > $tmeme_file`;
`cp $tmeme_file $tmeme_file2`;
`tomtom -o $prot_cell $tmeme_file $tmeme_file2`;

my $motif_similar = $prot_cell."/tomtom.txt";
my $Sinf = &read_summary($summary_file);
my $Seq_inf = &read_seq_count($seq_file);
my $Str_inf = &read_str_count($str_file);

my %sinf = %{$Sinf}; my %seq_inf = %{$Seq_inf}; my %str_inf = %{$Str_inf};
my %finf = (); my %fsinf = ();

my $Motif_com = &read_tomtom($motif_similar, $prot_cell, \%sinf);
my %motif_com = %{$Motif_com}; 
my $num = 0; my $pnum = 0; my $unum = 0;
foreach my $key (keys %sinf){
	my $count = (${$sinf{$key}}[1] =~ s/U/U/g);
	if($count >= 4){
		$unum = $unum + ${$sinf{$key}}[2];
	}else{
		$pnum = $pnum + ${$sinf{$key}}[2];
	}
	$num = $num + ${$sinf{$key}}[2];
}

foreach my $key ( sort{$a<=>$b} keys %motif_com){
	my @sen = @{$motif_com{$key}};
	#print $key,"\t",join("|", @sen),"\n";
}

foreach my $key ( sort{$a<=>$b} keys %motif_com){
	my @sen = @{$motif_com{$key}};
	#print ">",$key,"\t",join("|",@sen),"\n";
	if($#sen == -1){
		$finf{$key} = $sinf{$key};
		$fsinf{$key} = $seq_inf{$key};
		${$fsinf{$key}}{4} = ${$str_inf{$key}}{0};
		${$fsinf{$key}}{5} = ${$str_inf{$key}}{1};
	}else{
		$finf{$key} = $sinf{$key};
		$fsinf{$key} = $seq_inf{$key};
		${$fsinf{$key}}{4} = ${$str_inf{$key}}{0};
		${$fsinf{$key}}{5} = ${$str_inf{$key}}{1};
		for(my $i=0; $i<=$#sen; $i++){
			my @sent1 = split(/\|/, $sen[$i]);
			my $shf = -$sent1[1];
			${$finf{$key}}[2] = ${$finf{$key}}[2] + ${$sinf{$sent1[0]}}[2];
			for(my $i=0; $i<=3; $i++){
				for(my $j=max(0-$shf, 0); $j<=min(9-$shf, 9); $j++){
					${${$fsinf{$key}}{$i}}[$j] = ${${$fsinf{$key}}{$i}}[$j] + ${${$seq_inf{$sent1[0]}}{$i}}[$j+$shf];
				}
			}
			for(my $i=0; $i<=1; $i++){
				for(my $j=max(0-$shf, 0); $j<=min(9-$shf, 9); $j++){
					${${$fsinf{$key}}{$i+4}}[$j] = ${${$fsinf{$key}}{$i+4}}[$j] + ${${$str_inf{$sent1[0]}}{$i}}[$j+$shf];
				}
			}
		}
	}
}

open(OUT1, ">", $sum_out);
open(OUT2, ">", $seqstr_out);
print OUT1 $num,"\t",$unum,"\t",$unum/$num,"\t",$pnum,"\t",$pnum/$num,"\n";
foreach my $key (sort{${$finf{$b}}[2] <=> ${$finf{$a}}[2]} keys %finf){
	print OUT1 ${$finf{$key}}[0],"|",${$finf{$key}}[1],"\t",${$finf{$key}}[2],"\t",${$finf{$key}}[2]/$num,"\n";
	for(my $i=0; $i<=5; $i++){
		for(my $j=0; $j<=9; $j++){
			print OUT2 ${${$fsinf{$key}}{$i}}[$j],"\t";
		}
		print OUT2 "\n";
	}
}

close OUT1;
close OUT2;

sub read_summary{
	my $fdata_file = shift; 
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; 
	my @sen1 = (); my @sen2 = ();
	my %inf = ();
	open(FILE1, $fdata_file)||die("open $fdata_file error!\n");
	#2
	#UGCAUG|UUUUUU	157
	#UGCAUG|PPPPPP	119
	$sen = <FILE1>;
	while($sen = <FILE1>){
		$i = $i + 1;
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		@sen2 = split(/\|/, $sen1[0]);
		$inf{$i} = [$sen2[0], $sen2[1], $sen1[1]];
		if($i >= 10){
			last;
		}
	}
	close FILE1;
	return(\%inf);
}

sub mismatch{
	my $char1 = shift; my $char2 = shift;
	my @sen1 = split(//, $char1); my @sen2 = split(//, $char2); 
	my $i = 0; my $r = 0;
	for($i=0; $i<=$#sen1; $i++){
		if($sen1[$i] ne $sen2[$i]){
			$r = $r + 1;
		}
	}
	return($r);
}

sub read_seq_count{
	my $fdata_file = shift; 
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; 
	my @sen1 = (); my @sen2 = ();
	my %inf = ();
	open(FILE1, $fdata_file)||die("open $fdata_file error!\n");
	#39.2500	31.0000	35.0000	0.0000	0.0000	157.0000	0.0000	8.2500	33.2500	39.2500	
	#39.2500	45.0000	15.0000	0.0000	157.0000	0.0000	0.0000	8.2500	42.2500	39.2500
	while($sen = <FILE1>){
		$i = $i + 1;
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		${$inf{$i}}{0} = [@sen1[0..9]];
		for($j=1; $j<=3; $j++){
			$sen = <FILE1>;
			chomp($sen);
			@sen1 = split(/\t/, $sen);
			${$inf{$i}}{$j} = [@sen1[0..9]];
		}
		if($i >= 10){
			last;
		}
	}
	close FILE1;
	return(\%inf);
}

sub read_str_count{
	my $fdata_file = shift; 
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; 
	my @sen1 = (); my @sen2 = ();
	my %inf = ();
	open(FILE1, $fdata_file)||die("open $fdata_file error!\n");
	#39.2500	31.0000	35.0000	0.0000	0.0000	157.0000	0.0000	8.2500	33.2500	39.2500	
	#39.2500	45.0000	15.0000	0.0000	157.0000	0.0000	0.0000	8.2500	42.2500	39.2500
	while($sen = <FILE1>){
		$i = $i + 1;
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		${$inf{$i}}{0} = [@sen1[0..9]];
		for($j=1; $j<=1; $j++){
			$sen = <FILE1>;
			chomp($sen);
			@sen1 = split(/\t/, $sen);
			${$inf{$i}}{$j} = [@sen1[0..9]];
		}
		if($i >= 10){
			last;
		}
	}
	close FILE1;
	return(\%inf);
}

sub read_tomtom{
	my $fdata_file = shift; my $prot_name = shift; my $Sinf = shift;
	my $i = 0; my $j = 0; my $r = 0;
	my $sen = ""; my $Pattern = $prot_name."_";
	my @sen1 = (); my @sen2 = (); my @sent1 = (); my @sent2 = ();
	my %inf = (); my %uniq = (); my %sinf = %{$Sinf};
	$inf{1} = []; $uniq{1} = 1;
	open(FILE1, $fdata_file)||die("open $fdata_file error!\n");
	##Query ID	Target ID	Optimal offset	p-value	E-value	q-value	Overlap	Query consensus	Target consensus	Orientation
	#RBFOX2_mes1	RBFOX2_mes1	0	2.32559e-10	2.32559e-09	4.65117e-09	10	ACTGCATGTA	ACTGCATGTA	+
	#RBFOX2_mes3     RBFOX2_mes9     -1      0.00661106      0.0661106       0.0440737       9       AACATGTTCA      AAATGTGCCA      +
	#RBFOX2_mes1     RBFOX2_mes4     1       0.0902425       0.902425        0.150404        9       ACTGCATGTA      AAATGCTTGA      -
	#RBFOX2_mes3     RBFOX2_mes2     2       0.230382        2.30382 0.394583        8       AACATGTTCA      TCTGCATGCT      +
	$sen = <FILE1>;
	while($sen = <FILE1>){
		#$i = $i + 1;
		chomp($sen);
		@sen1 = split(/\t/, $sen);
		if(($sen1[5] < 0.05) && ($sen1[9] eq "+") &&($sen1[1] ne $sen1[0])){
			$sen1[0]=~s/$prot_name/$Pattern/g;
			$sen1[1]=~s/$prot_name/$Pattern/g;
			@sent1 = split(/\_/, $sen1[0]);
			@sent2 = split(/\_/, $sen1[1]);
			if(&mismatch(${$sinf{$sent1[2]}}[1], ${$sinf{$sent2[2]}}[1]) <= 1){
				if((!exists $uniq{$sent1[2]})&&(!exists $uniq{$sent2[2]})){
					if(exists $inf{$sent1[2]}){
						push(@{$inf{$sent1[2]}}, $sent2[2]."|".$sen1[2]."|".$sen1[7]."|".$sen1[8]);
						$uniq{$sent1[2]} = 1;
						$uniq{$sent2[2]} = 1;
					}else{
						$inf{$sent1[2]} = [$sent2[2]."|".$sen1[2]."|".$sen1[7]."|".$sen1[8]];
						$uniq{$sent1[2]} = 1;
						$uniq{$sent2[2]} = 1;
					}
				}
			}
		}
	}
	for($i=1; $i<=10; $i++){
		if((!exists $uniq{$i}) && (!exists $inf{$i})){
			$inf{$i} = [];
		}
	}
	
	close FILE1;
	return(\%inf);
}

exit;
