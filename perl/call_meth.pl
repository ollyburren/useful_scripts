use strict;
use File::Basename;
use Config::IniFiles;
use Getopt::Long;
use Data::Dumper;

## SET THIS TO DEBUG THE CALLING ALGO ON A TEST DATASET
## SHOULD BE ZERO UNDER NORMAL CIRCUMSTANCES
my $TEST_MET_CALL_ONLY=0;
## USED IN CALL METH WHEN POLY T's are PRESENT TO STOP INFINITE RECURSION
## IN PRACTICE SHOULDN'T EVER BE REACHED.
my $MAX_RECURSION_LEVEL=10; 

my $USAGE=<<EOL;
$0: Count methylation sites in a set of paired end reads.
This software given a set set of paired end FASTQ formatted files will attempt
to QC and call methylation status at a series of bases based on parameters defined
in the given ini file.

$0 --base_dir analysis_path --f1 forward.fq.gz -f2 reverse.fq.gz --ini inifile.ini --gene gene

	MANDATORY PARAMETERS:
		base_dir|b:	path to a location to store intermiediate and final analysis files
		f1:	path to forward set of reads in gzipped fastq format
		f2: path to reverse set of reads in gzipped fastq format
		ini: path to ini file containing configuration options
		gene|g: name of gene to analysed, there must be the [GENE], [GENE FOXP3_(left|right) 
			sections defined in the ini file..
			
	OPTIONAL:
		help|h: print this message.
		
	Original method Chris Penkett, Software Olly Burren.
		
Last updated 15/10/2013
EOL

####################
#OPTIONS PROCESSING#
####################


my ($base_dir,$fq1,$fq2,$gene,$ini,$help);

GetOptions(	
	"base_dir|b=s" => \$base_dir,
	"f1=s" => \$fq1,
	"f2=s" => \$fq2,
	"gene|g=s"=>\$gene,
	"ini=s"=>\$ini,
	"help|h"=>\$help);

my $ABORT_FLAG=0;
if(!$base_dir){
	print "[ERROR] Require --base_dir/-b option: Location for analysis files to be written\n";
	$ABORT_FLAG++;
}elsif(!$fq1){
	print "[ERROR] Require --f1 option: Path to fastaq file with forward PE reads to be analysed\n";
	$ABORT_FLAG++;
}elsif(!$fq2){
	print "[ERROR] Require --f2 option: Path to fastaq file with reverse PE reads to be analysed\n";
	$ABORT_FLAG++;
}elsif(!$ini){
	print "[ERROR] Require --ini option: Path to ini\n";
	$ABORT_FLAG++;
}elsif(!$gene){
	print "[ERROR] Require gene\n";
	$ABORT_FLAG++;
}

exit(1) if $ABORT_FLAG;

if(! -e $ini){
	print "[ERROR] Cannot find $ini ini file\n";
	$ABORT_FLAG++;
}elsif(! -e $fq1){
	print "[ERROR] Cannot find $fq1 forward FASTAQ file\n";
	$ABORT_FLAG++;
}elsif(! -e $fq2){
	print "[ERROR] Cannot find $fq2 reverse FASTAQ file\n";
	$ABORT_FLAG++;
}



my $cfg = new Config::IniFiles( -file => $ini);

#check that sections for gene given exist

if(!$cfg->SectionExists($gene) || 
		!$cfg->SectionExists("CUTADAPT ${gene}_left") || 
		!$cfg->SectionExists("CUTADAPT ${gene}_right")){
	print "[ERROR] Cannot find all conf sections for $gene\n";
	print "[ERROR] Expecting [$gene],[CUTADAPT ${gene}_left] and [CUTADAPT ${gene}_right]\n";
	$ABORT_FLAG++;
}

exit(1) if $ABORT_FLAG;

#################################
##PARSE IN SETTING FROM CFG FILE#
#################################

##BINARIES

my %BIN=(
	cutadapt=>$cfg->val('GENERAL','cutadapt_bin'),
	flash=>$cfg->val('GENERAL','flash_bin'),
	sickle=>$cfg->val('GENERAL','sickle_bin')
);

# quality score cutoff
my $QSCORE_CUTOFF=$cfg->val('GENERAL','qscore_cutoff');

my @MET_POSITION;
# methylation positions
foreach my $met($cfg->val($gene,'met_pos')){
	push @MET_POSITION,[split(",",$met)];
}

#die Dumper(@MET_POSITION)."\n";
my $trim_flag = $cfg->val($gene,'trim') || 'sickle';
#insert size
my $ISIZE=$cfg->val($gene,'insert_size');

#cutadapt parameters
my %CADAPT_PARAM=(
	left=>{
		a=>$cfg->val('CUTADAPT ilmn_left','a'),
		#q=> $cfg->val('CUTADAPT ilmn_left','q'),
		O=>$cfg->val('CUTADAPT ilmn_left','O'),
		m=>$cfg->val('CUTADAPT ilmn_left','m')
	},
	right=>{
		a=>$cfg->val('CUTADAPT ilmn_right','a'),
		#q=> $cfg->val('CUTADAPT ilmn_right','q'),
		O=>$cfg->val('CUTADAPT ilmn_right','O'),
		m=>$cfg->val('CUTADAPT ilmn_left','m')
	},
	foxp3_left=>{
		g=>$cfg->val("CUTADAPT ${gene}_left",'g'),
		m=>$cfg->val("CUTADAPT ${gene}_left",'m')
	},
	foxp3_right=>{
		a=>$cfg->val("CUTADAPT ${gene}_right",'a'),
	}
);


## using cut adapt q parameter and sickle can cause
## problems as sickle does not cope with reads with 0 bases in
unless($trim_flag eq 'sickle'){
	$CADAPT_PARAM{left}{q}=$cfg->val('CUTADAPT ilmn_left','q');
	$CADAPT_PARAM{right}{q}=$cfg->val('CUTADAPT ilmn_right','q');
}






###############################
##PREPARE ANALYSIS DIRECTORIES#
###############################


if(! -e $base_dir){
	mkdir($base_dir);
}
$base_dir.='/'.$gene;
if(! -e $base_dir){
	mkdir($base_dir);
	for my $d(qw/remove_il_adaptor manual_trim stitch final_trim results sickle pre_trim log/){
		mkdir($base_dir.'/'.$d);
	}
}else{
	print "$base_dir already exists overwriting\n";
}

foreach my $lr(keys %CADAPT_PARAM){
	my $cmd = $BIN{cutadapt}." ";
	foreach my $p(keys %{$CADAPT_PARAM{$lr}}){
		my $pv = $CADAPT_PARAM{$lr}->{$p};
		$cmd.=ref($pv) eq 'ARRAY'?join(" ",map{"-$p $_ "}@{$pv}):"-$p $pv ";
	}
	$CADAPT_PARAM{$lr}{cmd}=$cmd;
}



my %RC;

push @{$RC{RAW}},count_fastqreads($fq1);
push @{$RC{RAW}},count_fastqreads($fq2);

######################
##PRE TRIM SEQUENCES##
######################

my $cmd;
## some sequences need manually trimming this is done by setting trim option to 
## an integer 
if($trim_flag=~/[0-9]+/){
	$cmd = "zcat $fq1 | cut -c1-$trim_flag | gzip - >  $base_dir/pre_trim/".basename($fq1);
	print "$cmd\n";
	`$cmd` unless $TEST_MET_CALL_ONLY;
	$fq1 = "$base_dir/pre_trim/".basename($fq1);
	$cmd = "zcat $fq2 | cut -c1-$trim_flag | gzip - >  $base_dir/pre_trim/".basename($fq2);
	print "$cmd\n";
	`$cmd` unless $TEST_MET_CALL_ONLY;
	$fq2 = "$base_dir/pre_trim/".basename($fq2);
	$trim_flag='none'; ## make sure we don't sickle as well
	push @{$RC{PRE_TRIM}},count_fastqreads($fq1);
	push @{$RC{PRE_TRIM}},count_fastqreads($fq2);
}

##################################
##TRIM ILLUMINA ADAPTOR SEQUENCE##
##################################

##trim LH side ill
my $lhfile="$base_dir/remove_il_adaptor/".basename($fq1,'.fq.gz');
$cmd =  $CADAPT_PARAM{left}{cmd}." $fq1 2> $lhfile.cut.stats.out | gzip - > $lhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd` unless $TEST_MET_CALL_ONLY;
##trim RH side ill
my $rhfile="$base_dir/remove_il_adaptor/".basename($fq2,'.fq.gz');
$cmd =  $CADAPT_PARAM{right}{cmd}." $fq2 2> $rhfile.cut.stats.out | gzip - > $rhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd` unless $TEST_MET_CALL_ONLY;

push @{$RC{TRIM_IL}},count_fastqreads("$lhfile.trimmed.fq.gz");
push @{$RC{TRIM_IL}},count_fastqreads("$rhfile.trimmed.fq.gz");

##############################################
##TRIM POOR QUALITY SEQUENCE AT END OF READS##
##############################################

my $lsfile = "$base_dir/sickle/".basename($lhfile).".trimmed.fq";
my $rsfile = "$base_dir/sickle/".basename($rhfile).".trimmed.fq";
if($trim_flag eq 'sickle'){
	my $sickle_fail = "$base_dir/sickle/".basename($lhfile,".1");
	my $sickle_out = "$base_dir/sickle/".basename($lhfile,".1.fq.gz").".sickle.out";
	$cmd = "$BIN{sickle} pe -f $lhfile.trimmed.fq.gz -r $rhfile.trimmed.fq.gz -t sanger -o $lsfile - -p $rsfile -q $QSCORE_CUTOFF -s $sickle_fail -x > $sickle_out";
	print "$cmd\n";
	`$cmd` unless $TEST_MET_CALL_ONLY;
	`gzip -f $rsfile $lsfile`;
	$lsfile.='.gz';
	$rsfile.='.gz';
	push @{$RC{SICKLE}},count_fastqreads("$lsfile");
	push @{$RC{SICKLE}},count_fastqreads("$rsfile");
}elsif($trim_flag eq 'none'){ ## NOOP
	unlink($lsfile) if -e $lsfile;
	$cmd = "ln -s $lhfile.trimmed.fq.gz $lsfile";
	print "$cmd\n";
	`$cmd` unless $TEST_MET_CALL_ONLY;
	unlink($rsfile) if -e $rsfile;
	$cmd = "ln -s $rhfile.trimmed.fq.gz  $rsfile";
	print "$cmd\n";
	`$cmd` unless $TEST_MET_CALL_ONLY;
}

#########################
##STITCH READS TOGETHER##
#########################

my $sfile = basename($lhfile,'.1');
#$cmd="$BIN{flash} $lhfile.trimmed.fq.gz $rhfile.trimmed.fq.gz -z -o $sfile -d $base_dir/stitch/ >$base_dir/stitch/$sfile.stats.txt";
$cmd="$BIN{flash} $lsfile $rsfile -z -o $sfile -d $base_dir/stitch/ >$base_dir/stitch/$sfile.stats.txt";
print "$cmd\n";
`$cmd` unless $TEST_MET_CALL_ONLY;
push @{$RC{STITCH}},count_fastqreads("$base_dir/stitch/$sfile.extendedFrags.fastq.gz");

#############################
##TRIM UP TO FIRST METH SITE#
#############################

my $ltfile = "$base_dir/final_trim/$sfile.LTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_left}{cmd}." $base_dir/stitch/$sfile.extendedFrags.fastq.gz 2>$ltfile.cut.stats.out | gzip - > $ltfile.fq.gz";
print "$cmd\n";
`$cmd` unless $TEST_MET_CALL_ONLY;
push @{$RC{LEFT_TRIM}},count_fastqreads("$ltfile.fq.gz");

#############################
##TRIM REMAINING ADAPTOR 5' #
#############################

my $rtfile = "$base_dir/final_trim/$sfile.LRTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_right}{cmd}." $ltfile.fq.gz 2>$rtfile.cut.stats.out | gzip - > $rtfile.fq.gz";
print "$cmd\n";
`$cmd` unless $TEST_MET_CALL_ONLY;
push @{$RC{RIGHT_TRIM}},count_fastqreads("$rtfile.fq.gz");

###########################
#COUNT METHYLATION STATUS##
###########################

open(IN,"gunzip -c $rtfile.fq.gz |") || die "Cannot open $rtfile.fq.gz\n";
open(OUT,"> $base_dir/results/$sfile.out") || die "Cannot open $base_dir/results/$sfile.out\n";
open(LOG,"> $base_dir/log/$sfile.log") || die "Cannot open $base_dir/log/$sfile.log\n";
open(COUNT,"> $base_dir/log/$sfile.count") || die "Cannot open $base_dir/log/$sfile.count\n";

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
		my @s = split(//,$len[0]);
		my @q =  map{ord($_)-33}split(//,$len[1]);
		my $call = call_met(\@s,\@q,$ISIZE,\@MET_POSITION,$QSCORE_CUTOFF,0) if @len;
		$res{$call}++ if $call;
	}
}

#my $total;
my $pass_total=0;
my $error_total=0;
foreach my $r(sort{$b->[1] <=> $a->[1]}map{[$_,$res{$_}]}keys %res){
	if($r->[0] =~s/^\[LOG\]//){
		print LOG join("\t",$sfile,$gene,@$r)."\n";
		$error_total+=$r->[1];
	}else{
		print OUT join("\t",$sfile,$gene,@$r)."\n";
		$pass_total+=$r->[1];
	}
	#$total+=$r->[1];
}
print LOG join("\t",$sfile,$gene,'PASSED',$ISIZE,$pass_total)."\n";
push @{$RC{CALLED_READS}},$pass_total;
push @{$RC{OTHER_READS}},$error_total;

my $prev_count = $RC{RAW}->[0];
foreach my $a(qw/SICKLE STITCH LEFT_TRIM CALLED_READS OTHER_READS/){
	next unless $RC{$a};
	if($a=~/READS$/){
		print COUNT join("\t","$sfile",$gene,$a,@{$RC{$a}})."\n";
		next;
	}
	my @adj;
	foreach my $c(@{$RC{$a}}){
		push @adj,$prev_count-$c;
	}
	#print COUNT join("\t","$sfile",$gene,$a,'RAWC',@{$RC{$a}})."\n";
	print COUNT join("\t","$sfile",$gene,"${a}_REMOVED",@adj)."\n";
	$prev_count = $RC{$a}->[0];
}
close(OUT);
close(LOG);
close(COUNT);

##call methylation status at known bp's
##use recursion to cope with poly T's 
sub call_met{
	my ($s,$q,$isize,$met_aref,$qscore,$recur_count)=@_;
	## dead mans handle to make sure that recursion doesn't go insane !
	return "[LOG]MAX_RECURSION_LEVEL_REACHED" if $recur_count > $MAX_RECURSION_LEVEL;
	$isize-=$recur_count;
	my $met_index=$met_aref->[$recur_count] || return "[LOG]INSERT_LENGTH\t".scalar(@$s);
	return "[LOG]MISMATCH" if @$s != @$q;
	
	##RECURSION TO COPE WITH POLY T's 
	$recur_count++;
	return call_met($s,$q,$isize,$met_aref,$qscore,$recur_count) if @$s != $isize;
	my $qflag=1;
	for my $i(@$met_index){
		$qflag=0,last if $q->[$i]<$qscore;
		next if $s->[$i] eq 'C';
		$qflag=0,last if $s->[$i] ne 'T';
		$qflag=0,last if $s->[$i+1] ne 'G' or $q->[$i+1]<$qscore;
	}
	return "[LOG]FAILED_QC\tNA" unless $qflag;
	my @sequence = @$s;
	return  join("",@sequence[@$met_index]);
}

sub count_fastqreads{
	my $file = shift;
	my $line = `zcat $file | wc -l`;
	return $line/4;
}
