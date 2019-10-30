use Cwd 'abs_path';

$tissue = "";
$dir = './';
for($i=0;$i<= $#ARGV; $i++){
    if($ARGV[$i] eq "-d" || $ARGV[$i] eq "-dir"){
        $dir = $ARGV[++$i];
        next;
    }
    if($ARGV[$i] eq "-v" || $ARGV[$i] eq "-vcf"){
        $vcf = $ARGV[++$i];
        next;
    }
    if($ARGV[$i] eq "-t" || $ARGV[$i] eq "-tissue"){
        $tissue = $ARGV[++$i];
        next;
    }

}

if(! -e $vcf){
    print STDERR "Error: VCF header is not specified ... \n";
    exit(1);
}


$dir = abs_path($dir);;
@files = <$dir/*>;
$count=0;
foreach $f (@files){
    process_fm($f);
}

open FILE, "zcat $vcf |";
while(<FILE>){
    next if /^\s*\#/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    if(defined($snp{$data[2]})){
        print "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$snp{$data[2]}->{info}\n";
    }
}


sub process_fm{
    my ($f) = @_;
    $f =~/$dir\/(\S+?)\./;

    my $gene = $1;
    my %cluster;
    open FILE, "grep \\\{ $f | ";
    while(<FILE>){
        s/\{//;
        s/\}//;
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $cluster{$data[0]} = "\[$data[2]:$data[1]\]";
    }

    open FILE, "grep \\\(\\\( $f | ";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        next if $data[4] == -1;
        next if $data[2] < 1e-4;
        my $info = "$gene:$data[4]\@$tissue\=$data[2]".$cluster{$data[4]};
        if(!defined($snp{$data[1]})){
            $data[1] =~ /chr(\S+)\_(\d+)\_(\S+)\_(\S+)\_b38/;
            $map{$1}->{$2} = $data[1];
            $snp{$data[1]}->{header} = "chr$1\t$2\t$data[1]\t$3\t$4";
            $snp{$data[1]}->{info} = "$info";
        }else{
            $snp{$data[1]}->{info} .="|".$info;
        }


    }

}


