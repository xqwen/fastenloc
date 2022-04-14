#include "controller.h"
#include <stdlib.h>
#include <string>
#include <string.h>

#define LENGTH 1024
#define UNDEF -99999999

int show_banner(){

    fprintf(stderr, "\t\t==================================================================\n\n");
    fprintf(stderr, "\t\t                     fastENLOC (v1.0)                             \n\n");
    fprintf(stderr, "\t\t                Release date:  July, 2020                         \n\n"); 
    fprintf(stderr, "\t\t==================================================================\n\n\n");
    return 1;
}

int print_usage(){

    fprintf(stderr, "\nUsage: fastenloc -eqtl eqtl_file -gwas gwas_file [-tissue tissue_name] [-total_variants total_number_of_gwas_variants] [-thread number_of_processing_thread] [-s shrinkage_param] [-prefix output_prefix] \n\n");
    return 1;
}

int main(int argc, char **argv){


    char eqtl_file[LENGTH];
    char gwas_file[LENGTH];
    char tissue[LENGTH];
    char prefix[LENGTH];

    memset(eqtl_file, 0, LENGTH);
    memset(gwas_file, 0, LENGTH);
    memset(tissue, 0, LENGTH);
    memset(prefix, 0, LENGTH);

    int set_enrich_p = 0;
    int set_enrich_a = 0;

    double a0 = UNDEF;
    double a1 = UNDEF;

    double p1 = UNDEF;
    double p2 = UNDEF;
    double p12 = UNDEF;

    double output_thresh = 1e-4;

    int ImpN = 25;
    int nthread = 1;
    int enrich_est_only = 0;
    double shrinkage = -1;

    int total_snp = 0;

    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-e")==0 || strcmp(argv[i], "-eqtl")==0){
            strcpy(eqtl_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-gwas")==0){
            strcpy(gwas_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-tissue")==0){
            strcpy(tissue,argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-total_variant")==0 || strcmp(argv[i], "-total_variants")==0 || strcmp(argv[i], "-tv")==0){
            total_snp = atoi(argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-imp")==0 || strcmp(argv[i], "-impute")==0){
            ImpN = atoi(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-shrinkage")==0 || strcmp(argv[i], "-s")==0){
            shrinkage = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-thread")==0){
            nthread = atoi(argv[++i]);
            continue;
        }
        if(strcmp(argv[i], "-a0")==0){
            a0 = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-a1")==0){
            a1 = atof(argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-p1")==0){
            p1 = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-p2")==0){
            p2 = atof(argv[++i]);
            continue;
        }
        if(strcmp(argv[i], "-p12")==0){
            p12 = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-prefix")==0){
            strcpy(prefix,argv[++i]);
            continue;
        }
        
        if(strcmp(argv[i], "--all")==0 || strcmp(argv[i], "--output_all")==0){
            output_thresh = 0;
            continue;
        }

        if(strcmp(argv[i], "--est")==0 || strcmp(argv[i], "--est_only")==0 || strcmp(argv[i], "--enrich")==0 || strcmp(argv[i], "--enrich_only")==0){
            enrich_est_only = 1;
            continue;
        }

        fprintf(stderr,"Error: unknown command option \'%s\'\n", argv[i]); 
        continue;


    }
    
    show_banner();

    // check required options

    if(strlen(eqtl_file)==0){
        fprintf(stderr, "Error: molecular QTL annotation file is unspecified ... exit\n");
        print_usage();
        exit(1);
    }

    if(strlen(gwas_file)==0){
        fprintf(stderr, "Error: GWAS fine-mapping  file is unspecified ... exit\n");
        print_usage();
        exit(2);
    }


    // print command line options
    
    fprintf(stderr, "\nParameters and options:\n\n");
    fprintf(stderr, "Input files:\n");
    fprintf(stderr, "    * Molecular qtl annotation file: %s\n", eqtl_file);
    fprintf(stderr, "    * GWAS fine-mapping file: %s\n", gwas_file);
    if(strlen(tissue) !=0){
        fprintf(stderr, "    * Tissue specified: %s\n", tissue);
    }


    fprintf(stderr, "\nEnrichment parameters:\n");
    if(p1 != UNDEF || a1 != UNDEF){
        fprintf(stderr, "    * Specified, skip estimation procedure\n");
    }else{
        fprintf(stderr, "    * Rounds of multiple imputation: %d\n", ImpN);
        double shp = 1.0;
        if(shrinkage != -1){
            shp = shrinkage;
        }
        fprintf(stderr, "    * Shrinkage parameter: %.1f\n", shp);
    }


    fprintf(stderr, "\nMiscsellaneous options:\n");
    fprintf(stderr, "    * Total GWAS variants: ");
    if(total_snp == 0){
        fprintf(stderr, "unspecified, use GWAS file input\n");
    }else{
        fprintf(stderr, "%d\n", total_snp);
    }

    fprintf(stderr, "    * Simultaneous running threads: %d\n", nthread);
    


    fprintf(stderr, "\nOutput options:\n");
    if(strlen(prefix) !=0)
        fprintf(stderr, "    * Output file prefix: %s\n", prefix);
    fprintf(stderr, "    * RCP and SCP output threshold: %.1e\n", output_thresh);


    fprintf(stderr, "\n\n");





    controller con;

    con.set_imp_num(ImpN);
    con.set_snp_size(total_snp);
    con.set_thread(nthread);
    con.set_prefix(prefix);
    con.set_output_thresh(output_thresh);

    // default shrinkage
    double pv = 1;
    // no shrinkage
    if(shrinkage == 0){
        pv = -1;
    }
    // user-spcified valid shrinkage: 1/prior_variance
    if(shrinkage >0){
        pv = 1/shrinkage;
    }

    con.set_prior_variance(pv);

    








    con.load_eqtl(eqtl_file, tissue);
    con.load_gwas_torus(gwas_file);


    if(a0 != UNDEF && a1 != UNDEF){
        con.set_enrich_params(a0,a1);
    }else if(p1 != UNDEF && p2 != UNDEF && p12!= UNDEF){
        con.set_enrich_params(p1,p2,p12);
    }else{
        con.enrich_est();
    }
    if(!enrich_est_only){
        con.compute_coloc_prob();
    }
}   
