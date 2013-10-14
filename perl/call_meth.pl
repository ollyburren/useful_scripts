use strict;
use File::Basename;
use Config::IniFiles;
use Getopt::Long;

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
		
Last updated 14/10/2013
EOL

####################
#OPTIONS PROCESSING#
####################

##MONDAY CONVERT TO USING GETOPTS THEN ADD CD3 PARAMETERS TO INI FILE.

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
	flash=>$cfg->val('GENERAL','flash_bin')
);

# quality score cutoff
my $QSCORE_CUTOFF=$cfg->val('GENERAL','qscore_cutoff');

# methylation positions
my @MET_POSITION=split(",",$cfg->val($gene,'met_pos'));
	
#insert size
my $ISIZE=$cfg->val($gene,'insert_size');

#cutadapt parameters
my %CADAPT_PARAM=(
	left=>{
		a=>$cfg->val('CUTADAPT ilmn_left','a'),
		q=> $cfg->val('CUTADAPT ilmn_left','q'),
		O=>$cfg->val('CUTADAPT ilmn_left','O')
	},
	right=>{
		a=>$cfg->val('CUTADAPT ilmn_right','a'),
		q=> $cfg->val('CUTADAPT ilmn_right','q'),
		O=>$cfg->val('CUTADAPT ilmn_right','O')
	},
	foxp3_left=>{
		g=>$cfg->val("CUTADAPT ${gene}_left",'g'),
		m=>$cfg->val("CUTADAPT ${gene}_left",'m')
	},
	foxp3_right=>{
		a=>$cfg->val("CUTADAPT ${gene}_right",'a'),
	}
);


###############################
##PREPARE ANALYSIS DIRECTORIES#
###############################


if(! -e $base_dir){
	mkdir($base_dir);
}
$base_dir.='/'.$gene;
if(! -e $base_dir){
	mkdir($base_dir);
	for my $d(qw/remove_il_adaptor manual_trim stitch final_trim results/){
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

##################################
##TRIM ILLUMINA ADAPTOR SEQUENCE##
##################################

##trim LH side ill
my $lhfile="$base_dir/remove_il_adaptor/".basename($fq1,'.fq.gz');
my $cmd =  $CADAPT_PARAM{left}{cmd}." $fq1 2> $lhfile.cut.stats.out | gzip - > $lhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd`;
##trim RH side ill
my $rhfile="$base_dir/remove_il_adaptor/".basename($fq2,'.fq.gz');
$cmd =  $CADAPT_PARAM{right}{cmd}." $fq2 2> $lhfile.cut.stats.out | gzip - > $rhfile.trimmed.fq.gz"; 
print $cmd."\n";
`$cmd`;

#########################
##STITCH READS TOGETHER##
#########################

my $sfile = basename($lhfile,'.1');
$cmd="$BIN{flash} $lhfile.trimmed.fq.gz $rhfile.trimmed.fq.gz -z -o $sfile -d $base_dir/stitch/ >$base_dir/stitch/$sfile.stats.txt"; 
print "$cmd\n";
`$cmd`;

#############################
##TRIM UP TO FIRST METH SITE#
#############################

my $ltfile = "$base_dir/final_trim/$sfile.LTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_left}{cmd}." $base_dir/stitch/$sfile.extendedFrags.fastq.gz 2>$ltfile.cut.stats.out | gzip - > $ltfile.fq.gz";
print "$cmd\n";
`$cmd`;

#############################
##TRIM REMAINING ADAPTOR 5' #
#############################

my $rtfile = "$base_dir/final_trim/$sfile.LRTRIM";
my $cmd =  $CADAPT_PARAM{foxp3_right}{cmd}." $ltfile.fq.gz 2>$rtfile.cut.stats.out | gzip - > $rtfile.fq.gz";
print "$cmd\n";
`$cmd`;

###########################
#COUNT METHYLATION STATUS##
###########################

open(IN,"gunzip -c $rtfile.fq.gz |") || die "Cannot open $rtfile.fq.gz\n";
open(OUT,"> $base_dir/results/$sfile.out") || die "Cannot open $base_dir/results/$sfile.out\n";

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
