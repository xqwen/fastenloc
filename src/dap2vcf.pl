#!/usr/bin/perl
use Cwd 'abs_path';

$tissue = "";
$dir = './';
$thresh = 0;

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

    if($ARGV[$i] eq "-p"){
        $thresh = $ARGV[++$i];
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
        next if $data[2] < $thresh;
        $cluster{$data[0]} = "\[$data[2]:$data[1]\]";
    }

    open FILE, "grep \\\(\\\( $f | ";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        next if $data[4] == -1;
        next if !defined($cluster{$data[4]});

        next if $data[2] < 1e-4;
        my $info = "$gene:$data[4]\@$tissue\=$data[2]".$cluster{$data[4]};
        $id = "$gene\_$data[1]";
        if(!defined($snp{$id})){
            $snp{$id}->{info} = "$info";
        }else{
            $snp{$id}->{info} .="|".$info;
        }


    }

}


