use strict;
use warnings;

## uncomment if BioPerl is not available in PERL5LIB path
#use lib 
use Bio::SeqIO;
use Data::Dumper;

my $match_dir = 'location of psl files delimiter is ;';

my $seqout =  Bio::SeqIO->new( -format => 'Fasta', -fh => \*STDOUT);

## number of iterations to run if testing
my $count=2;

while(<DATA>){
	next if /^#/;
  chomp;
  my $id = $_;
	my $nurl = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?DB_NAME=EMBL&id=$id&style=raw";
  open(EMBL,"curl -s \"$nurl\" |");
  my $stream = Bio::SeqIO->new(-fh => \*EMBL,-format=>'embl');
  my $seq = $stream->next_seq();
  my %genes;
  my $seqname;
  foreach my $f($seq->get_SeqFeatures()){
  		$seqname =  $f->seq_id;
     print $f->gff_string()."\n" if $f->primary_tag =~/source|gene/;
  }
  parsePSL("$match_dir/$seqname-A.psl",$seqname,'probeA');
  parsePSL("$match_dir/$seqname-B.psl",$seqname,'probeB');
  ## uncomment to test
  #last unless $count--;
}

sub parsePSL{
	my ($file,$acc,$pt) = @_;
	open(PSL,"cut -d\\; -f1,9-11,21 $file |");
	while(<PSL>){
		chomp;
		my @vals=split(";",$_);
		next unless $vals[0] == $vals[3];
		$vals[4]=~s/,//;
		next if $vals[4]=~/,/; ## this mean there is a mismatch on haplotype !
		my ($start,$end);
		if($vals[1] eq '+'){
			$start=$vals[4];
			$end=$vals[4]+$vals[0]-1;
		}else{
			$start=$vals[4]-$vals[0];
			$end=$vals[4];
		}
			
		print join("\t",$acc,'blat',$pt,$start,$end,$vals[1],'.','.',"id=$vals[2]")."\n"
	}
}
	
__DATA__
#AC006293
#AC011501
#AL133414
AY320039
GU182338
GU182339
GU182340
GU182341
GU182342
GU182343
GU182344
GU182345
GU182346
GU182347
GU182348
GU182349
GU182350
GU182351
GU182352
GU182353
GU182354
GU182355
GU182356
GU182357
GU182358
GU182359
GU182360
GU182361
#GU182362
