open FILE, "sim_data/sim.truth.rst";
while(<FILE>){
    
    chomp;
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    next if $data[-1] == 0 || $data[-2] == 0;
    $rcd{$data[0]} =  sprintf "%9.5f %9.5f", $data[-2], $data[-1];

}


open FILE, "fastenloc_out/sim.true_prior.enloc.snp.out";
while(<FILE>){
    
    chomp;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    next if !defined($rcd{$data[1]});
    $cluster{$data[1]} = $data[0];
    $rcp{$data[0]} = 0;
    $scp{$data[1]} = $data[-1];
}

open FILE, "fastenloc_out/sim.true_prior.enloc.sig.out";
while(<FILE>){

    chomp;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;

    next if !defined($rcp{$data[0]});
    $rcp{$data[0]} = $data[-1];
    $eqtl{$data[0]} = $data[2];
    $gwas{$data[0]} = $data[4];
}

printf "coloc_snp\tsig_cluster\tSCP\tRCP\teQTL_eff\tGWAS_eff\teQTL_spip\tGWAS_spip\n";
foreach $s (keys %rcd){

    my $c = $cluster{$s};
    if (!defined($cluster{$s})){
        $c = "NA";
    }
    printf "%s\t%s\t\t%7.3e\t%7.3e\t%s\t%7.3e\t%7.3e\n", $s, $c, $scp{$s}, $rcp{$c}, $rcd{$s}, $eqtl{$c}, $gwas{$c};
}
