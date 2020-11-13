open FILE, "sort -grk6 $ARGV[0] | ";
while(<FILE>){
    
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    $data[0]=~/(\S+)\:\d+/;
    $gene = $1;
    next if defined $rcd{$gene};

    print "$gene\t$data[-1]\n";
    $rcd{$gene} = 1;
}


