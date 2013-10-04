use strict;
use File::Basename;

#depends on the template (zero based).
my @MET_POSITION=(0,35,54,63,71,75,81,84,94,102,143,148);
#insert size
my $ISIZE=151;
# quality score cutoff
my $QSCORE_CUTOFF=20;
my %BIN=(
	cutadapt=>"$ENV{HOME}/.local/bin/cutadapt",
	flash=>"$ENV{HOME}/src/FLASH-1.2.7/flash"
);
my %CADAPT_PARAM=(
	left=>{
		a=>['ATCTCGTATGCCGTCTTCTGCTTG','AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
		q=> 20,
		O=>6
	},
	right=>{
		a=>['AGATCTCGGTGGTCGCCGTATCATT','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
		q=> 20,
		O=>6
	},
	foxp3_left=>{
		g=>'TATAGATTATGTTTTTATAT',
		m=>171
	},
	foxp3_right=>{
		a=>'GAATTGGTTGTTTTGTTTTGTAGTAG'
	}
);



##prepare directories

my $BASE_DIR=shift | '/stats/oliver/BISSEQ/test/';
my $FQ1=shift || '/dunwich/scratch/olly/DRAINBOW/bis_test/fastq/B790001_P3_D2.1.fq.gz';
my $FQ2=shift || '/dunwich/scratch/olly/DRAINBOW/bis_test/fastq/B790001_P3_D2.2.fq.gz';;

if(!$BASE_DIR || !-e $FQ1 || !-e $FQ2){
	print "Incorrect parameters\n";
	exit;
}

if(! -e $BASE_DIR){
	mkdir($BASE_DIR);
	for my $d(qw/remove_il_adaptor manual_trim stitch final_trim results/){
		mkdir($BASE_DIR.'/'.$d);
	}
}else{
	print "$BASE_DIR already exists aborting\n";
}





foreach my $lr(keys %CADAPT_PARAM){
	my $cmd = $BIN{cutadapt}." ";
	foreach my $p(keys %{$CADAPT_PARAM{$lr}}){
		my $pv = $CADAPT_PARAM{$lr}->{$p};
		$cmd.=ref($pv) eq 'ARRAY'?join(" ",map{"-$p $_ "}@{$pv}):"-$p $pv ";
	}
	$CADAPT_PARAM{$lr}{cmd}=$cmd;
}

##################################
##TRIM ILLUMINA ADAPTOR SEQUENCE##
##################################

##trim LH side ill
my $lhfile="$BASE_DIR/remove_il_adaptor/".basename($FQ1,'.fq.gz');
my $cmd =  $CADAPT_PARAM{left}{cmd}." $FQ1 2> $lhfile.cut.stats.out | gzip - > $lhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd`;
##trim RH side ill
my $rhfile="$BASE_DIR/remove_il_adaptor/".basename($FQ2,'.fq.gz');
$cmd =  $CADAPT_PARAM{right}{cmd}." $FQ2 2> $lhfile.cut.stats.out | gzip - > $rhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd`;

#########################
##STITCH READS TOGETHER##
#########################

my $sfile = basename($lhfile,'.1');
$cmd="$BIN{flash} $lhfile.trimmed.fq.gz $rhfile.trimmed.fq.gz -z -o $sfile -d $BASE_DIR/stitch/ >$BASE_DIR/stitch/$sfile.stats.txt"; 
print "$cmd\n";
`$cmd`;

#############################
##TRIM UP TO FIRST METH SITE#
#############################

my $ltfile = "$BASE_DIR/final_trim/$sfile.LTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_left}{cmd}." $BASE_DIR/stitch/$sfile.extendedFrags.fastq.gz 2>$ltfile.cut.stats.out | gzip - > $ltfile.fq.gz";
print "$cmd\n";
`$cmd`;

#############################
##TRIM REMAINING ADAPTOR 5' #
#############################

my $rtfile = "$BASE_DIR/final_trim/$sfile.LRTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_right}{cmd}." $ltfile.fq.gz 2>$rtfile.cut.stats.out | gzip - > $rtfile.fq.gz";
print "$cmd\n";
`$cmd`;

###########################
#COUNT METHYLATION STATUS##
###########################

open(IN,"gunzip -c $rtfile.fq.gz |") || die "Cannot open $rtfile.fq.gz\n";
open(OUT,"> $BASE_DIR/results/$sfile.out") || die "Cannot open $BASE_DIR/results/$sfile.out\n";

my %res;
my @len;
##read in fastq formatted file
while(<IN>){
	chomp;
	if($. % 4==1){
		@len=();
	}elsif($. % 4 ==2 ){
		push @len,$_;
	}elsif($. % 4 ==0){
		push @len,$_;
		my $call = call_met($len[0],$len[1],$ISIZE,\@MET_POSITION,$QSCORE_CUTOFF) if @len;
		$res{$call}++ if $call;
	}
}

#sort and output results seems quicker than UNIX 
my $total;
foreach my $r(sort{$b->[1] <=> $a->[1]}map{[$_,$res{$_}]}keys %res){
	print OUT join("\t",$sfile,@$r)."\n";
	$total+=$r->[1];
}
print "Total $total reads processed for $sfile\n";


##call methylation status at known bp's
sub call_met{
	my ($s,$q,$isize,$met_index,$qscore)=@_;
	return "MISMATCH" if(length($s) ne length($q));
	my @sequence = split(//,$s);
	my @quality = map{ord($_)-33}split(//,$q);
	my $slength=length($s);
	return "INSERT LENGTH" if $slength != $isize;
	my $qflag=1;
	for my $i(@$met_index){
		$qflag=0,last if $quality[$i]<$qscore;
		next if $sequence[$i] eq 'C';
		$qflag=0,last if $sequence[$i] ne 'T';
		$qflag=0,last if $sequence[$i+1] ne 'G' or $quality[$i+1]<$qscore;
	}
	return 'FAILED QC' unless $qflag;
	return  join("",@sequence[@$met_index]);
}
