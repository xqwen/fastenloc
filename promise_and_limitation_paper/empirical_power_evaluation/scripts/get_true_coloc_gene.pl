open FILE, "sim_data/sim.truth.rst";
while(<FILE>){

    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;

    next if $data[-1] == 0 || $data[-2] == 0;
    $data[0]=~/g(\d+)/;
    $rcd{$1} = 1;
}

foreach (sort {$a <=> $b} keys %rcd){
    print "$_\n";
}
