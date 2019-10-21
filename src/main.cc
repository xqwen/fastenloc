#include "controller.h"
#include <stdlib.h>
#include <string>
#include <string.h>

#define LENGTH 1024
#define UNDEF -99999999


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

    int ImpN = 25;
    int nthread = 1;
    double pc = 1;

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


        if(strcmp(argv[i], "-total_variant")==0){
            total_snp = atoi(argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-imp")==0 || strcmp(argv[i], "-impute")==0){
            ImpN = atoi(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-pc")==0 || strcmp(argv[i], "-pseudo_count")==0){
            pc = atof(argv[++i]);
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


        fprintf(stderr,"Error: unknown command option \'%s\'\n", argv[i]); 
        continue;


    }


    controller con;

    con.set_imp_num(ImpN);
    con.set_snp_size(total_snp);
    con.set_thread(nthread);
    con.set_prefix(prefix);
    con.set_pseudo_count(pc);
    


    con.load_eqtl(eqtl_file, tissue);
    con.load_gwas_torus(gwas_file);


    if(a0 != UNDEF && a1 != UNDEF){
        con.set_enrich_params(a0,a1);
    }else if(p1 != UNDEF && p2 != UNDEF && p12!= UNDEF){
        con.set_enrich_params(p1,p2,p12);
    }else{
        con.enrich_est();
    }

    con.compute_coloc_prob();
}   
