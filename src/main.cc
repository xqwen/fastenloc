#include "controller.h"
#include <stdlib.h>
#include <string>
#include <string.h>

#define LENGTH 1024



int main(int argc, char **argv){


    char eqtl_file[LENGTH];
    char gwas_file[LENGTH];
    char tissue[LENGTH];
    char prefix[LENGTH];

    memset(eqtl_file, 0, LENGTH);
    memset(gwas_file, 0, LENGTH);
    memset(tissue, 0, LENGTH);
    memset(prefix, 0, LENGTH);

    int ImpN = 25;
    int nthread = 1;

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

        if(strcmp(argv[i], "-imp")==0 || strcmp(argv[i], "-impute")==0){
            ImpN = atoi(argv[++i]);
            continue;
        }
                                    
        if(strcmp(argv[i], "-thread")==0){
            nthread = atoi(argv[++i]);
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
    con.set_thread(nthread);
    con.set_prefix(prefix);

    con.load_eqtl(eqtl_file, tissue);
    con.load_gwas_torus(gwas_file);
    
    con.enrich_est();
    con.compute_coloc_prob();
}   
