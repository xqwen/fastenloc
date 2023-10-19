open FILE, "enloc.sig.out";
while(<FILE>){

    chomp;
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;

    $rcd{$data[0]} = sprintf "$data[-2] $data[-1]";
}

open FILE, "compare/enloc.sig.out";
while(<FILE>){
    chomp;
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    next if !defined ($rcd{$data[0]});
    print "$data[0]\t$rcd{$data[0]}\t$data[-2] $data[-1]\n";
}
