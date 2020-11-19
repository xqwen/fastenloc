sub get_thresh {

    my $alpha = 0.05;
    my ($d) = @_;
    my @data = sort {$b <=> $a} @{$d};
    my $sum = 0;
    my $count = 0;
    my $thresh = 1;
    foreach $v (@data){
        $sum += (1-$v);
        $count++;
        last if($sum/$count > $alpha);
        $thresh = $v;
    }

    return $thresh;
}

 
open FILE, "awk \'\{print \$6\}\' fastenloc_out/sim.enloc.sig.out |";
@rcp = <FILE>;
chomp @rcp;
$rcp_thresh = get_thresh(\@rcp);


open FILE,  "awk \'\{print \$2\}\' summary/sim.gwas_analysis.summary |";
@gwas = <FILE>;
chomp @gwas;
$gwas_thresh = get_thresh(\@gwas);


open FILE,  "awk \'\{print \$2\}\' summary/sim.eqtl_analysis.summary |";
@eqtl = <FILE>;
chomp @eqtl;
$eqtl_thresh = get_thresh(\@eqtl);

#print "$rcp_thresh\n";
#print "$gwas_thresh\n";
#print "$eqtl_thresh\n";

printf "SNP                  cluster type    eqtl_ratio  gwas_ratio  rcp\n";
open FILE, "results/fastenloc.est_prior.summary";
while(<FILE>){

    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    

    my $type = 0;
    if($data[3]>= $rcp_thresh){
        $type = 3;
    }else{

        if($data[6]>= $eqtl_thresh && $data[7] >= $gwas_thresh){
            $type = 2;
        }elsif ($data[6] < $eqtl_thresh && $data[7] <  $gwas_thresh){
            $type = 0;
        }else{
            $type = 1;
        }
    }
    
    my ($er, $gr) = parse_data($data[0], $type, $data[1]);
    printf "%16s  %10s  %d   %7.2f  %7.2f   %7.3f\n", $data[0], $data[1], $type, $er, $gr,$data[3];

}


sub parse_data{

    my($snpid, $type, $sid) = @_;
    $snpid=~/(g\d+)\_(rs\d+)/;
    $gene = $1;
    $snp =  $2;

    my $eqtl_ratio = 0;
    my $gwas_ratio = 0;
    my $ecluster;
    my $gcluster;

    my $epip;
    my $lead_epip;

    $out = `grep \\\(\\\( dap_eqtl/$gene.fm.rst | grep $snp`;
    chomp $out;
    my @data = split /\s+/, $out;
    shift @data until $data[0]=~/^\S/;
    $ecluster = $data[4];
    $epip = $data[2];

    if($ecluster != -1){
        open DATA, "grep \\\(\\\( dap_eqtl/$gene.fm.rst |";
        while(<DATA>){
            my @d2 = split /\s+/, $_;
            shift @d2 until $d2[0]=~/^\S/;
            if($d2[4] == $ecluster){
                $eqtl_ratio = $epip/$d2[2];
                $lead_epip = $d2[2];
                last;
            }
        }

    }


    $out = `grep \\\(\\\( dap_gwas/$gene.fm.rst | grep $snp`;
    chomp $out;
    @data = split /\s+/, $out;
    shift @data until $data[0]=~/^\S/;
    $gcluster = $data[4];
    $gpip = $data[2];

    if($gcluster != -1){
        open DATA, "grep \\\(\\\( dap_gwas/$gene.fm.rst |";
        while(<DATA>){
            my @d2 = split /\s+/, $_;
            shift @d2 until $d2[0]=~/^\S/;
            if($d2[4] == $gcluster){
                $gwas_ratio = $gpip/$d2[2];
                $lead_gpip = $d2[2];
                last;
            }
        }

    }

    #printf "%16s  %10s  %d   %7.2f  %7.2f\n", $snpid, $sid, $type, $eqtl_ratio, $gwas_ratio;
    return ($eqtl_ratio, $gwas_ratio);   
}



