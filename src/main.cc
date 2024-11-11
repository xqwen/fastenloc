#include "controller.h"
#include <stdlib.h>
#include <string>
#include <string.h>
#include <math.h>

#define LENGTH 1024
#define UNDEF -99999999

int show_banner(){

    fprintf(stderr, "\t\t==================================================================\n\n");
    fprintf(stderr, "\t\t                    FastENLOC (v3.1)                        \n\n");
    fprintf(stderr, "\t\t              release date: November, 2024            \n\n"); 
    fprintf(stderr, "\t\t==================================================================\n\n\n");
    fprintf(stderr, "\n\n");
    return 1;
}

int print_usage(){

    fprintf(stderr, "\nUsage: fastenloc -eqtl eqtl_file -gwas gwas_file [-tissue tissue_name] [-total_variants total_number_of_gwas_variants] [-thread number_of_processing_thread] [-s shrinkage_param] [-prefix output_prefix] \n\n");
    return 1;
}

int main(int argc, char **argv){


    char eqtl_file[LENGTH];
    char gwas_file[LENGTH];

    char combined_summary_file[LENGTH];

    char tissue[LENGTH];
    char prefix[LENGTH];

    memset(eqtl_file, 0, LENGTH);
    memset(gwas_file, 0, LENGTH);
    memset(combined_summary_file, 0, LENGTH);

    memset(tissue, 0, LENGTH);
    memset(prefix, 0, LENGTH);

    int set_enrich_p = 0;
    int set_enrich_a = 0;

    int outlier_control = 1;

    double a0 = UNDEF;
    double a1 = UNDEF;

    double p1 = UNDEF;
    double p2 = UNDEF;
    double p12 = UNDEF;


    double sdY_eqtl = UNDEF;
    double sdY_gwas = UNDEF;
    
    double binary_eqtl = 0;
    double binary_gwas = 0; 

    int sum_conversion = 0; 

    double cap_a1 = -UNDEF;

    double output_thresh = 1e-4;

    int ImpN = 25;
    int nthread = 1;
    int enrich_est_only = 0;
    double shrinkage = -1;

    int total_snp = 0;
    int set_warning = 0;
    int set_error = 0;

    int gwas_format = 1;
    int coloc_prob_option = 2;  // default using exact algorithm

    unsigned long random_seed = 0;

    show_banner();
  
    /// Taking user command-line input
   
    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-e")==0 || strcmp(argv[i], "-eqtl")==0){
            strcpy(eqtl_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-gwas")==0){
            strcpy(gwas_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-go")==0){
            strcpy(gwas_file,argv[++i]);
            gwas_format = 2;
            continue;
        }

         if(strcmp(argv[i], "-gs")==0){
            strcpy(gwas_file,argv[++i]);
            gwas_format = 3;
            continue;
        }

        if(strcmp(argv[i], "-sum")==0 || strcmp(argv[i], "-summary")==0){
            strcpy(combined_summary_file,argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "--conv") ==0){
            sum_conversion = 1;
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


        if(strcmp(argv[i], "--coloc_default_prior") == 0){
            p1 = 1e-4;
            p2 = 1e-4;
            p12 = 1e-5;
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

        if(strcmp(argv[i], "-cap_a1")==0){
            cap_a1 = atof(argv[++i]);
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
        
        if(strcmp(argv[i], "-seed")==0){
            random_seed = (unsigned long)atoi(argv[++i]);
            continue;
        }
        
        if(strcmp(argv[i], "--all")==0 || strcmp(argv[i], "--output_all")==0){
            output_thresh = 0;
            continue;
        }

        if(strcmp(argv[i], "--approx")==0 || strcmp(argv[i], "--apprx") ==0 || strcmp(argv[i], "--legacy") == 0) { 
            coloc_prob_option = 1;
            continue;
        }

        if(strcmp(argv[i], "--exact")==0 ){
            coloc_prob_option = 2;
            continue;
        }


        if(strcmp(argv[i], "--est")==0 || strcmp(argv[i], "--est_only")==0 || strcmp(argv[i], "--enrich")==0 || strcmp(argv[i], "--enrich_only")==0){
            enrich_est_only = 1;
            continue;
        }



        if(strcmp(argv[i], "-sdy_e")==0 || strcmp(argv[i], "-sdy_eqtl")==0){
            sdY_eqtl = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-sdy_g")==0 || strcmp(argv[i], "-sdy_gwas")==0){
            sdY_gwas = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "--binary_e")==0 || strcmp(argv[i], "--binary_eqtl")==0){
            binary_eqtl = 1;
            continue;
        }


        if(strcmp(argv[i], "--binary_g")==0 || strcmp(argv[i], "--binary_gwas")==0){
            binary_gwas = 1;
            continue;
        }





        fprintf(stderr,"\nError: unknown command option \'%s\'\n", argv[i]); 
        set_error = 1;
        continue;


    }
    

    if(set_error==1){
        print_usage();
        exit(1);
    }



    // pre-processing/convertion of user input

    // default shrinkage
    double pv = 1;
    // no shrinkage
    if (shrinkage == 0)
    {
        pv = -1;
    }
    // user-spcified valid shrinkage: 1/prior_variance
    if (shrinkage > 0)
    {
        pv = 1 / shrinkage;
    }


    // start processing

    controller con;
    
    //  Enrichment Analysis Option

    con.set_a1_cap(cap_a1);
    int run_enrich = 1;
    if (p1 != UNDEF && p2 != UNDEF && p12 != UNDEF)
    {
        fprintf(stderr,"Applying user-specified colocalization priors, skipping enrichment analysis\n\n");
        con.set_enrich_params(p1, p2, p12);
        run_enrich = 0;
    }else{

        // total_snp is required input now
        
        if(total_snp == 0){
            fprintf(stderr, "Error: number of total GWAS variants is unspecified\n\n");
            exit(2);
        }

        con.set_snp_size(total_snp);

        if (a0 != UNDEF && a1 != UNDEF)
        {
            fprintf(stderr,"Applying user-specified colocalization priors, skipping enrichment analysis\n\n");
            con.set_enrich_params(a0, a1);
            run_enrich = 0; 
        
        }else{
            // setup enrichment analysis required parameters
            con.set_imp_num(ImpN);
            con.set_random_seed(random_seed);
            con.set_outlier_control(outlier_control);
            con.set_prior_variance(pv);
        }

    }


    // input data file format options 

    if (strlen(combined_summary_file) > 0)
    {
        // BF computation options
        vector<double> prior_vec1;
        if (sdY_eqtl != UNDEF)
            prior_vec1.push_back(pow(0.15 * sdY_eqtl, 2));
        else if (binary_eqtl)
            prior_vec1.push_back(0.04);

        vector<double> prior_vec2;
        if (sdY_gwas != UNDEF)
            prior_vec2.push_back(pow(0.15 * sdY_gwas, 2));
        else if (binary_gwas)
            prior_vec2.push_back(0.04);

        con.set_abf_piror_vec(prior_vec1, prior_vec2);
        con.load_combined_summary(combined_summary_file);
        con.init_pip();

    } else if (strlen(eqtl_file) == 0){
        // no eqtl fine-mapping info specified
        fprintf(stderr, "Error: molecular QTL annotation file is unspecified \n\n");
        exit(2);
    }else{
         con.load_eqtl(eqtl_file, tissue);

        if (strlen(gwas_file) == 0)
        {
            fprintf(stderr, "Error: GWAS fine-mapping file is unspecified\n\n");
            exit(2);
        }

        if (gwas_format == 1)
            con.load_gwas(gwas_file);

        if (gwas_format == 2)
            con.load_gwas_torus(gwas_file);

        if (gwas_format == 3)
        {
            vector<double> prior_vec1;
            vector<double> prior_vec2;
            con.set_abf_piror_vec(prior_vec1, prior_vec2);
            con.load_gwas_summary(gwas_file);
        }
    }

    // set coloc computational option
    con.set_coloc_prob_option(coloc_prob_option);

    // set output options
    con.set_prefix(prefix);
    con.set_output_thresh(output_thresh);

    // start enrichment analysis
    if (run_enrich)
    {
        con.enrich_est();
    }

    if (!enrich_est_only)
    {
        con.compute_coloc_prob();
    }
    fprintf(stderr, "\n\n\nFastENLOC analysis is completed\n\n");
}
