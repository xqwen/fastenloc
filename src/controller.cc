#include "controller.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <omp.h>
#include <zlib.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Input Processing ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

///// eQTL with fine-mapping information

void controller::load_eqtl(char *eqtl_file, char *tissue)
{

    fprintf(stderr, "Processing eQTL annotations ... \n");

    istringstream ins;

    string target_tissue = string(tissue);
    string chr;
    string pos;
    string snp_id; // important
    string allele1;
    string allele2;
    string content;

    string delim = "|";

    int format_status = 0;

    vector<string> lines = readLines(string(eqtl_file));

    for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it) 
    {

        string line = *it;
        ins.clear();
        ins.str(line);

        if (ins >> chr >> pos >> snp_id >> allele1 >> allele2 >> content)
        {

            format_status = 1;

            // passing content
            size_t pos = 0;
            string token;

            do
            {
                pos = content.find(delim);
                token = content.substr(0, pos);
                size_t pos1 = token.find("@");
                size_t pos2 = token.find("=");
                size_t pos3 = token.find("[");
                size_t pos4 = token.find("]");
                string sig_id;
                string tissue_type;
                if (pos1 != string::npos && pos2 - pos1 > 1)
                {
                    tissue_type = token.substr(pos1 + 1, pos2 - pos1 - 1);
                }

                sig_id = token.substr(0, pos1);

                // printf("sig id = %s %d %d\n",sig_id.c_str(),pos2,pos1);
                string pip = token.substr(pos2 + 1, pos3 - pos2 - 1);
                string sig_pip = token.substr(pos3 + 1, pos4 - pos3 - 1);
                // processing

                if (strlen(tissue) == 0 || tissue_type.compare(target_tissue) == 0)
                {

                    if (snp_index.find(snp_id) == snp_index.end())
                    {
                        snp_index[snp_id] = snp_vec.size();
                        snp_vec.push_back(snp_id);
                    }

                    if (eqtl_sig_index.find(sig_id) == eqtl_sig_index.end())
                    {
                        eqtl_sig_index[sig_id] = eqtl_vec.size();
                        sigCluster cluster;
                        cluster.id = sig_id;
                        cluster.cpip = 0;
                        eqtl_vec.push_back(cluster);
                    }

                    int index = eqtl_sig_index[sig_id];
                    double val = atof(pip.c_str());
                    eqtl_vec[index].cpip += val;
                    eqtl_vec[index].snp_vec.push_back(snp_id);
                    eqtl_vec[index].pip_vec.push_back(val);

                    // std::cout << snp_id<< "    "<< sig_id<<"  "<< tissue_type<<"   "<< pip << " "<< sig_pip<< std::endl;
                }

                content.erase(0, pos + delim.length());
            } while (pos != string::npos);
        }
    }

    if (format_status == 0)
    {
        fprintf(stderr, "\nError: unexpected format in eQTL annotation file \"%s\".\n\n", eqtl_file);
        exit(1);
    }

    double sum = 0;
    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        sum += eqtl_vec[i].cpip;
    }

    if (sum == 0)
    {
        fprintf(stderr, "\nError: no eQTL annotated in \"%s\". \n\n", eqtl_file);
        exit(1);
    }

    fprintf(stderr, "read in %d SNPs, %d eQTL signal clusters, %.1f expected eQTLs\n\n", int(snp_vec.size()), int(eqtl_vec.size()), sum);

    P_eqtl = sum;
}




///// GWAS with fine-mapping information


void controller::load_gwas(char *gwas_file, char *tissue)
{

    fprintf(stderr, "Processing complex trait data ... \n");

    istringstream ins;

    string chr;
    string pos;
    string snp_id; // important
    string allele1;
    string allele2;
    string content;

    string delim = "|";

    double gwas_sum = 0;
    int gwas_count = 0;
    int format_status = 0;

    vector<string> lines = readLines(string(gwas_file));
    for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it) 
    {

        string line = *it;
        ins.clear();
        ins.str(line);

        if (ins >> chr >> pos >> snp_id >> allele1 >> allele2 >> content)
        {
            format_status = 1;
            // passing content
            size_t pos = 0;
            string token;

            do
            {
                pos = content.find(delim);
                token = content.substr(0, pos);
                size_t pos1 = token.find("@");
                size_t pos2 = token.find("=");
                size_t pos3 = token.find("[");
                size_t pos4 = token.find("]");
                string sig_id;
                string tissue_type;
                if (pos1 != string::npos && pos2 - pos1 > 1)
                {
                    tissue_type = token.substr(pos1 + 1, pos2 - pos1 - 1);
                }

                sig_id = token.substr(0, pos1);

                // printf("sig id = %s %d %d\n",sig_id.c_str(),pos2,pos1);
                string pip = token.substr(pos2 + 1, pos3 - pos2 - 1);
                string sig_pip = token.substr(pos3 + 1, pos4 - pos3 - 1);
                // processing

                //

                if (snp_index.find(snp_id) == snp_index.end())
                {
                    snp_index[snp_id] = snp_vec.size();
                    snp_vec.push_back(snp_id);
                }

                // new GWAS locus
                if (snp2gwas_locus.find(snp_id) == snp2gwas_locus.end())
                {
                    snp2gwas_locus[snp_id] = sig_id;
                }
                else
                {
                    snp2gwas_locus[snp_id] += "_" + sig_id;
                }

                if (gwas_sig_index.find(sig_id) == gwas_sig_index.end())
                {
                    gwas_sig_index[sig_id] = gwas_vec.size();
                    sigCluster cluster;
                    cluster.id = sig_id;
                    cluster.cpip = 0;
                    gwas_vec.push_back(cluster);
                }

                int index = gwas_sig_index[sig_id];
                double val = atof(pip.c_str());
                gwas_vec[index].cpip += val;
                gwas_vec[index].snp_vec.push_back(snp_id);
                gwas_vec[index].pip_vec.push_back(val);
                gwas_sum += val;
                gwas_count++;

                content.erase(0, pos + delim.length());
            } while (pos != string::npos);
        }
    }
    if (format_status == 0)
    {
        fprintf(stderr, "\nError: unexpected format in complex trait file \"%s\".\n\n", gwas_file);
        exit(1);
    }

    if (gwas_sum == 0)
    {
        fprintf(stderr, "\nError: no complex trait associations detected in \"%s\". \n\n", gwas_file);
        exit(1);
    }

    gwas_pip_vec = vector<double>(snp_vec.size(), 0.0);

    for (int i = 0; i < gwas_vec.size(); i++)
    {
        for (int j = 0; j < gwas_vec[i].snp_vec.size(); j++)
        {
            string snp = gwas_vec[i].snp_vec[j];
            gwas_pip_vec[snp_index[snp]] = gwas_vec[i].pip_vec[j];
        }
    }

    if (total_snp < snp_vec.size())
        total_snp = snp_vec.size();

    fprintf(stderr, "read in %d SNPs (eQTL+gwas), %d GWAS loci, %.1f expected hits\n\n", int(snp_vec.size()), int(gwas_vec.size()), gwas_sum);

    // estimated marginal priors
    P_gwas = gwas_sum / total_snp;
    P_eqtl = P_eqtl / total_snp;
}


/// combined GWAS eQTL input with only summary statistics (bhat, se_bhat)


void controller::load_combined_summary(char* summary_input)
{

    fprintf(stderr, "Processing eQTL + GWAS combined summary statistics input ... \n");

    istringstream ins;

    string snp_id; // important
    string sig_id;
    double bhat_e;
    double se_e;
    double bhat_g;
    double se_g;

    vector<string> lines = readLines(string(summary_input));
    for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it)
    {
        string line = *it;
        ins.clear();
        ins.str(line);

        if (ins >> sig_id >> snp_id >> bhat_e >> se_e >> bhat_g >> se_g)
        {


            if (snp_index.find(snp_id) == snp_index.end())
            {
                snp_index[snp_id] = snp_vec.size();
                snp_vec.push_back(snp_id);
            }

            if (eqtl_sig_index.find(sig_id) == eqtl_sig_index.end())
            {
                eqtl_sig_index[sig_id] = eqtl_vec.size();
                sigCluster eqtl_cluster;
                eqtl_cluster.set_prior_var_vec(eqtl_abf_prior_vec);
                eqtl_cluster.id = sig_id;
                eqtl_cluster.cpip = 0;
                eqtl_vec.push_back(eqtl_cluster);
            }

            int index = eqtl_sig_index[sig_id];
            eqtl_vec[index].snp_vec.push_back(snp_id);
            eqtl_vec[index].compute_BF(bhat_e, se_e);


            if (snp2gwas_locus.find(snp_id) == snp2gwas_locus.end())
            {
                snp2gwas_locus[snp_id] = sig_id;
            }
            else
            {
                snp2gwas_locus[snp_id] += "_" + sig_id;
            }

            if (gwas_sig_index.find(sig_id) == gwas_sig_index.end())
            {
                gwas_sig_index[sig_id] = gwas_vec.size();
                sigCluster gwas_cluster;
                gwas_cluster.set_prior_var_vec(gwas_abf_prior_vec);
                gwas_cluster.id = sig_id;
                gwas_cluster.cpip = 0;
                gwas_vec.push_back(gwas_cluster);
            }

            index = gwas_sig_index[sig_id];
            gwas_vec[index].snp_vec.push_back(snp_id);
            gwas_vec[index].compute_BF(bhat_g, se_g);


        }
    }

    // legacy GWAS format
    fprintf(stderr, "read in %d SNPs, %d loci \n\n", int(snp_vec.size()), int(gwas_vec.size()) );

    total_snp = snp_vec.size();


}





///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Data Pre-processing //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////





// for summary statistics input, after setting P_eqtl and p_gwas
void controller::init_pip(){

   if(P_eqtl <0){ 
        P_eqtl = torus_estimate(eqtl_vec);
        printf("Estimated eQTL frequency = %9.3e\n", P_eqtl);
   }

    if(P_gwas < 0) {
        P_gwas = torus_estimate(gwas_vec);
        printf("Estimated GWAS frequency = %9.3e\n", P_gwas);
    }



    for(int i=0;i<eqtl_vec.size();i++){
        eqtl_vec[i].compute_pip(P_eqtl);
    }


    for(int i=0;i<gwas_vec.size();i++){
        gwas_vec[i].compute_pip(P_gwas);
    }



    gwas_pip_vec = vector<double>(snp_vec.size(), 0.0);

    for (int i = 0; i < gwas_vec.size(); i++)
    {
        for (int j = 0; j < gwas_vec[i].snp_vec.size(); j++)
        {
            string snp = gwas_vec[i].snp_vec[j];
            gwas_pip_vec[snp_index[snp]] = gwas_vec[i].pip_vec[j];
        }
    }



}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Enrichment Specification and Estimation ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Set enrichment parameters without estimation procedure

void controller::set_enrich_params(double a0, double a1)
{
    if (a1 > cap_a1)
    {
        fprintf(stderr, "\n\n*** estimated a1 (%.3f) exceeds the user-specified cap value, set a1 = %.3f ***\n", a1, cap_a1);
        a1 = cap_a1;
    }
    pi1_e = exp(a0 + a1) / (1 + exp(a0 + a1));
    pi1_ne = exp(a0) / (1 + exp(a0));
}

// Set enrichment parameters without estimation procedure: p_1 = p(gwas only), p_2 = p(eqtl only), p12 = p( colocalization )

void controller::set_enrich_params(double p1, double p2, double p12)
{

    double a0 = log(p1 / (1 - p1 - p2 - p12));
    double a1 = log(p12 * (1 - p1 - p2 - p12) / (p1 * p2));
    fprintf(stderr, "converting enrichment parameters: \n");
    fprintf(stderr, "%10s   %7.3f\n", "Intercept", a0);
    fprintf(stderr, "%10s   %7.3f\n\n", "Enrichment", a1);

    P_gwas = p1 + p12;
    P_eqtl = p2 + p12;
    set_enrich_params(a0, a1);
}


// Estimate enrichment parameters by multiple imputation


void controller::enrich_est()
{

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    if (seed != 0)
    {
        gsl_rng_set(r, seed);
        fprintf(stderr, "Use user-specified random seed: %lu\n\n", seed);
    }

    vector<double> a0_vec = vector<double>(ImpN, 0.0);
    vector<double> v0_vec = vector<double>(ImpN, 0.0);
    vector<double> a1_vec = vector<double>(ImpN, 0.0);
    vector<double> v1_vec = vector<double>(ImpN, 0.0);

    vector<double> eqtl_sample_vec = vector<double>(ImpN, 0.0);


#pragma omp parallel for num_threads(nthread)
    for (int k = 0; k < ImpN; k++)
    {


        vector<int> eqtl_sample = vector<int>(snp_vec.size(), 0);

        for (int i = 0; i < eqtl_vec.size(); i++)
        {
            // impute/sample eQTN
            int index = eqtl_vec[i].impute_qtn(r);

            if (index >= 0)
            {
                string snp = eqtl_vec[i].snp_vec[index];
                eqtl_sample[snp_index[snp]] = 1;
                eqtl_sample_vec[k] += 1;
            }
        }
        vector<double> rst = run_EM(eqtl_sample);
#pragma omp critical
        {
            a0_vec[k] = rst[0];
            a1_vec[k] = rst[1];
            v0_vec[k] = rst[2];
            v1_vec[k] = rst[3];
        }
        fprintf(stderr, "Imputation round %2d is completed\n", k + 1);
    }

    // MI post-processing

    // Detect and remove outliers from MI results

    if (outlier_control)
    {

        vector<double> a1_shrink_temp_vec;
        double mean = 0;
        double sd = 0;
        for (int k = 0; k < ImpN; k++)
        {
            double post_a1 = (a1_vec[k] * prior_variance) / (prior_variance + v1_vec[k]);
            a1_shrink_temp_vec.push_back(post_a1);
            mean += post_a1;
            sd += pow(post_a1, 2);
        }

        mean = mean / ImpN;
        sd = sqrt(sd / ImpN - pow(mean, 2));

        int remove_count = 0;
        vector<double> a0_vec_temp;
        vector<double> a1_vec_temp;
        vector<double> v0_vec_temp;
        vector<double> v1_vec_temp;

        for (int k = 0; k < ImpN; k++)
        {
            double dist = fabs(a1_shrink_temp_vec[k] - mean) / sd;
            if (dist > 3)
            {
                remove_count++;
                continue;
            }
            else
            {

                a1_vec_temp.push_back(a1_vec[k]);
                v1_vec_temp.push_back(v1_vec[k]);
                a0_vec_temp.push_back(a0_vec[k]);
                v0_vec_temp.push_back(v0_vec[k]);
            }
        }

        if (remove_count > 0)
        {

            fprintf(stderr, "\n*** %d outlying MI results removed ***\n", remove_count);
            ImpN -= remove_count;

            a1_vec = a1_vec_temp;
            a0_vec = a0_vec_temp;
            v1_vec = v1_vec_temp;
            v0_vec = v0_vec_temp;
        }
    }

    fprintf(stderr, "\nEffective MI numbers: %d", ImpN);
    // main processing

    string mi_file = prefix + string("enloc.mi.out");
    FILE *fd_mi = fopen(mi_file.c_str(), "w");
    fprintf(fd_mi, "%7s\t%7s\t%7s\t%7s\n", "a0", "a1", "p_eqtl", "p_gwas");

    a0_est = 0;
    a1_est = 0;
    double a1_shrink_est = 0;
    double var1_shrink = 0;

    double var0 = 0;
    double var1 = 0;
    vector<double> a1_shrink_vec;
    vector<double> v1_shrink_vec;
    for (int k = 0; k < ImpN; k++)
    {
        a0_est += a0_vec[k];
        a1_est += a1_vec[k];
        var0 += v0_vec[k];
        var1 += v1_vec[k];

        double post_var = 1.0 / (1.0 / prior_variance + 1 / (v1_vec[k]));
        double post_a1 = (a1_vec[k] * prior_variance) / (prior_variance + v1_vec[k]);

        a1_shrink_vec.push_back(post_a1);
        v1_shrink_vec.push_back(post_var);
        a1_shrink_est += post_a1;
        var1_shrink += post_var;

        double p_eqtl_mi = eqtl_sample_vec[k] / total_snp;
        double p_gwas_mi = p_eqtl_mi * exp(a0_vec[k] + a1_vec[k]) / (1 + exp(a0_vec[k] + a1_vec[k])) + (1 - p_eqtl_mi) * exp(a0_vec[k]) / (1 + exp(a0_vec[k]));

        fprintf(fd_mi, "%7.3f\t%7.3f\t\t%7.3e\t%7.3e\n", a0_vec[k], a1_shrink_vec[k], p_eqtl_mi, p_gwas_mi);
    }

    fclose(fd_mi);
    a0_est = a0_est / ImpN;
    a1_est = a1_est / ImpN;
    a1_shrink_est = a1_shrink_est / ImpN;

    double bv0 = 0;
    double bv1 = 0;
    double bv1_shrink = 0;

    for (int k = 0; k < ImpN; k++)
    {

        bv0 += pow(a0_vec[k] - a0_est, 2);
        bv1 += pow(a1_vec[k] - a1_est, 2);
        bv1_shrink += pow(a1_shrink_vec[k] - a1_shrink_est, 2);
    }

    bv0 = bv0 / (ImpN - 1);
    bv1 = bv1 / (ImpN - 1);
    bv1_shrink = bv1_shrink / (ImpN - 1);
    var0 = var0 / ImpN;
    var1 = var1 / ImpN;
    var1_shrink = var1_shrink / ImpN;

    double sd0 = sqrt(var0 + bv0 * (ImpN + 1) / ImpN);
    double sd1 = sqrt(var1 + bv1 * (ImpN + 1) / ImpN);
    double sd1_shrink = sqrt(var1_shrink + bv1_shrink * (ImpN + 1) / ImpN);

    a0_est = log(P_gwas / (1 + P_eqtl * exp(a1_est) - P_eqtl - P_gwas));

    // apply double shrinkage to a1 estimate
    a1_shrink_est = prior_variance * a1_shrink_est / (prior_variance + pow(sd1_shrink, 2));
    sd1_shrink = sqrt(1.0 / (1.0 / prior_variance + 1 / pow(sd1_shrink, 2)));

    string enrich_file = prefix + string("enloc.enrich.out");
    FILE *fd = fopen(enrich_file.c_str(), "w");
    fprintf(fd, "%25s   %7.3f     %7s\n", "Intercept", a0_est, "-");
    fprintf(fd, "%25s   %7.3f     %7.3f\n", "Enrichment (no shrinkage)", a1_est, sd1);
    fprintf(fd, "%25s   %7.3f     %7.3f\n", "Enrichment (w/ shrinkage)", a1_shrink_est, sd1_shrink);

    // use shrinkage estimate
    a1_est = a1_shrink_est;

    double p1 = (1 - P_eqtl) * exp(a0_est) / (1 + exp(a0_est));
    double p2 = P_eqtl / (1 + exp(a0_est + a1_est));
    double p12 = P_eqtl * exp(a0_est + a1_est) / (1 + exp(a0_est + a1_est));

    fprintf(fd, "\n\n## Alternative (coloc) parameterization: p1 = %7.3e, p2 = %7.3e, p12 = %7.3e\n\n", p1, p2, p12);

    fclose(fd);

    set_enrich_params(a0_est, a1_est);

    // pi1_e = exp(a0_est+a1_est)/(1+exp(a0_est + a1_est));
    // pi1_ne =  exp(a0_est)/(1+exp(a0_est));
    // printf("%7.3e  %7.3e     %7.3e\n", pi1_e, pi1_ne, pi1);

    fprintf(stderr, "\nEnrichment analysis is completed\n");
    fprintf(stderr, "pi1 = %7.3e  pi1_e = %7.3e pi1_ne = %7.3e\np_eqtl = %7.3e\n", P_gwas, pi1_e, pi1_ne, P_eqtl);

    gsl_rng_free(r);
    return;
}






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Computing and reporting colocalization probabilities ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////// Main entry control function

void controller::compute_coloc_prob()
{

    if (coloc_prob_option == 1)
        compute_coloc_prob_legacy();
    else
        compute_coloc_prob_exact();
}



////// exact computation algorithm implemented in version 3

void controller::compute_coloc_prob_exact()
{

    fprintf(stderr, "\nComputing colocalization probabilities ... \n\n");

    string snp_file = prefix + string("enloc.snp.out");
    string sig_file = prefix + string("enloc.sig.out");
    string gen_file = prefix + string("enloc.gene.out");

    FILE *fd1 = fopen(snp_file.c_str(), "w");
    FILE *fd2 = fopen(sig_file.c_str(), "w");
    FILE *fd3 = fopen(gen_file.c_str(), "w");

    fprintf(fd1, "Signal\tSNP\tPIP_qtl\tPIP_gwas_marginal\tPIP_gwas_qtl_prior\tSCP\n");
    fprintf(fd2, "Signal\tNum_SNP\tCPIP_qtl\tCPIP_gwas_marginal\tCPIP_gwas_qtl_prior\tRCP\tLCP\n");
    fprintf(fd3, "Gene\t\tGRCP\tGLCP\n");


    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        eqtl_vec[i].coloc_prob = 0;
        eqtl_vec[i].locus_coloc_prob = 0;

        vector<double> gprob_vec_m; // marginal gwas pips
        double gprob_cpip_m = 0;    // marginal gwas cpip

        vector<double> glog10bf_vec; // gwas BF

        int p = eqtl_vec[i].snp_vec.size(); // number of SNPs

        //printf("%s: %d\n",eqtl_vec[i].id.c_str(), p);


        for (int k = 0; k < p; k++)
        {
            string snp = eqtl_vec[i].snp_vec[k];
            double d = gwas_pip_vec[snp_index[snp]];
            gprob_vec_m.push_back(d);
            gprob_cpip_m += d;
        }



        // renormalizing GWAS marginal PIP if necessary (advoid prob overflow due to LD difference)
        if (gprob_cpip_m > 1 - 1e-5)
        {
            // renormalize
            for (int k = 0; k < p; k++)
            {
                gprob_vec_m[k] = (1 - 1e-5) * (gprob_vec_m[k] / gprob_cpip_m);
            }

            gprob_cpip_m = 1 - 1e-5;
        }

        double locus_epip = 0;

        // recovering GWAS BFs from fine-mapping results
        // Important: the contrast is the scenario that all SNPs in the cluster has zero effects
        // The procedure is based on DAP-1 approximation

        double log10_MNC = log10(1 - P_gwas ) - log10(1 - gprob_cpip_m);



        for (int k = 0; k < p; k++)
        {
            double gpost = gprob_vec_m[k];

            if (gpost == 0)
                gpost = 1e-10;
            
            
            double log10bf = log10(gpost) - log10(P_gwas) + log10_MNC;
            glog10bf_vec.push_back(log10bf);
            locus_epip += eqtl_vec[i].pip_vec[k];
        }
        if (locus_epip > 1 - 1e-8)
            locus_epip = 1 - 1e-8;

        // computing p(d,r | E, G) by considering all configureations within a single signal cluster (DAP-1 approximation for both traits)

        double NC = 0;

        // case 1: null-null scenario
        NC += pow(10, p * log10(1 - pi1_ne) + log10(1 - locus_epip));

        // case 2: gwas only-eqtl null
        vector<double> g_only_prob_vec;
        for (int k = 0; k < p; k++)
        {
            double prob = pow(10, log10(pi1_ne) + (p - 1) * log10(1 - pi1_ne) + glog10bf_vec[k] + log10(1 - locus_epip));
            NC += prob;
            g_only_prob_vec.push_back(prob);
        }

        // case 3: eqtl only-gwas null
        for (int k = 0; k < p; k++)
        {
            NC += pow(10, log10(1 - pi1_e) + (p - 1) * log10(1 - pi1_ne) + log10(eqtl_vec[i].pip_vec[k]));
        }

        // case 4: gwas-eqtl combination
        vector<vector<double> > coloc_config; // coloc_config[m][n] indicates m-th SNP is the GWAS hit and n-th SNP is the causal eQTL
                                             // initialization
        for (int k = 0; k < p; k++)
        {
            coloc_config.push_back(vector<double>(p, 0));
        }

        for (int m = 0; m < p; m++)
        {

            for (int n = 0; n < p; n++)
            {
                double log10_prior;
                if (m == n)
                {
                    log10_prior = log10(pi1_e) + (p - 1) * log10(1 - pi1_ne);
                }
                else
                {
                    log10_prior = log10(pi1_ne) + log10(1 - pi1_e) + (p - 2) * log10(1 - pi1_ne);
                }
                double prob = pow(10, log10_prior + glog10bf_vec[m] + log10(eqtl_vec[i].pip_vec[n]));
                NC += prob;
                coloc_config[m][n] = prob;
            }
        }

        // collecting results
        vector<double> scp_vec;    // snp-level colocalization probabilities; sum of it is RCP
        vector<double> gpip_e_vec; // gwas pips informed by eQTL annotation
        double RCP = 0;            // region-wise SNP-level colocalization probability
        double LCP = 0;            // locus-level colocalization probability

        for (int m = 0; m < p; m++)
        {
            for (int n = 0; n < p; n++)
            {
                double prob = coloc_config[m][n] / NC;
                LCP += prob;

                if (m == n)
                {
                    RCP += prob;
                    scp_vec.push_back(prob);
                }
            }
        }


        eqtl_vec[i].coloc_prob = RCP;
        eqtl_vec[i].locus_coloc_prob = LCP;

        fprintf(stderr, "Locus %s:  RCP = %7.3e  LCP = %7.3e\n", eqtl_vec[i].id.c_str(), RCP, LCP);

        double gprob_cpip_e = 0; // GWAS cluster/locus cpip with eQTL annotation

        for (int m = 0; m < p; m++)
        {
            // marginalize over (p+1) eQTL configurations for each GWAS SNP
            double gprob_e = g_only_prob_vec[m] / NC;
            for (int n = 0; n < p; n++)
            {
                gprob_e += coloc_config[m][n] / NC;
            }
            gpip_e_vec.push_back(gprob_e);

            gprob_cpip_e += gprob_e;
        }

        /////////////////////////////////////////////// Reporting ////////////////////////////////////

        double max_scp = 0;
        string max_snp;

        for (int k = 0; k < eqtl_vec[i].snp_vec.size(); k++)
        {
            string snp = eqtl_vec[i].snp_vec[k];
            string locus_id = eqtl_vec[i].id + "(@)" + snp2gwas_locus[snp];
            if (scp_vec[k] >= output_thresh)
                fprintf(fd1, "%15s   %15s   %7.3e %7.3e    %7.3e      %7.3e\n", locus_id.c_str(), snp.c_str(), eqtl_vec[i].pip_vec[k], gprob_vec_m[k], gpip_e_vec[k], scp_vec[k]);
            if (scp_vec[k] >= max_scp)
            {
                max_scp = scp_vec[k];
                max_snp = snp;
            }
        }

        string locus_id = eqtl_vec[i].id + "(@)" + snp2gwas_locus[max_snp];

        if (RCP >= output_thresh || LCP >= output_thresh)
        {
            string locus_id = eqtl_vec[i].id + "(@)" + snp2gwas_locus[max_snp];
            fprintf(fd2, "%15s   %4d  %7.3e %7.3e    %7.3e      %7.3e\t%7.3e\n", locus_id.c_str(), int(eqtl_vec[i].snp_vec.size()), eqtl_vec[i].cpip, gprob_cpip_m, gprob_cpip_e, RCP, LCP);
        }
    }

    // gene-level quantification
    map<string, double> GLCP;
    map<string, double> GRCP;

    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        string loc_id = eqtl_vec[i].id;
        int pos = loc_id.find(":");
        string gene_id = loc_id.substr(0, pos);
        if (GLCP.find(gene_id) == GLCP.end())
        {
            GLCP[gene_id] = 1;
            GRCP[gene_id] = 1;
        }

        GLCP[gene_id] *= (1 - eqtl_vec[i].locus_coloc_prob);
        GRCP[gene_id] *= (1 - eqtl_vec[i].coloc_prob);
    }

    map<string, double>::iterator it = GRCP.begin();
    while (it != GRCP.end())
    {
        string gene = it->first;
        double v1 = 1 - it->second;
        double v2 = 1 - GLCP[gene];
        fprintf(fd3, "%s\t\t%7.3e\t%7.3e\n", gene.c_str(), v1, v2);
        it++;
    }

    fclose(fd1);
    fclose(fd2);
    fclose(fd3);
}








/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Utility functions //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// File parsing utility

vector<string> controller::readLines(const string &filename)
{

    vector<string> lines;

    // Check if file is gzip compressed
    ifstream file(filename, ios::binary);
    if (!file)
    {
        cerr << "Error: cannot open file: " << filename << endl;
        exit(1);
    }

    unsigned char magic[2];
    file.read(reinterpret_cast<char *>(magic), 2);
    file.close();

    if (magic[0] == 0x1F && magic[1] == 0x8B)
    { // gzip magic numbers
        gzFile gz = gzopen(filename.c_str(), "rb");
        if (!gz)
        {
            cerr << "Error: cannot open gzip file: " << filename << endl;
            exit(2);
        }

        char buffer[4096];
        while (true)
        {
            if (gzgets(gz, buffer, sizeof(buffer)) == NULL)
                break;
            lines.push_back(string(buffer));
        }
        gzclose(gz);
    }
    else
    {
        ifstream txtFile(filename);
        if (!txtFile)
        {
            cerr << "Error: cannot open text file: " << filename << endl;
            return lines;
        }

        string line;
        while (getline(txtFile, line))
        {
            lines.push_back(line);
        }
        txtFile.close();
    }

    return lines;
}





///////////// Find MLE of the enrichment parameters given an MI eQTL sample ////////////////////

vector<double> controller::run_EM(vector<int> &eqtl_sample)
{

    // starting point
    double a0 = log(P_gwas / (1 - P_gwas));
    double a1 = 0;

    double var0;
    double var1;

    double r1 = exp(a0 + a1);    // prior odds for eQTL
    double r0 = exp(a0);         // prior odds for non-eQTL
    double rm = P_gwas / (1 - P_gwas); // marginal prior odds

    while (1)
    {
        double pseudo_count = 1.0;
        double e0g0 = pseudo_count * (1 - P_gwas) * (1 - P_eqtl);
        double e0g1 = pseudo_count * (1 - P_eqtl) * P_gwas;
        double e1g0 = pseudo_count * (1 - P_gwas) * P_eqtl;
        double e1g1 = pseudo_count * P_gwas * P_eqtl;

        for (int i = 0; i < snp_vec.size(); i++)
        {

            double val = gwas_pip_vec[i];
            if (val == 1)
                val = 1 - 1e-8;
            // posterior ratio
            val = val / (1 - val); // convert to BF
                                   // val/r_null is marginal likelihood/bayes factor
            if (eqtl_sample[i] == 0)
            {
                val = r0 * (val / rm); // posterior odds
                val = val / (1 + val); // updated posterior with current prior given eqtl = 0

                // bookkeeping
                e0g1 += val;
                e0g0 += 1 - val;
            }

            if (eqtl_sample[i] == 1)
            {
                val = r1 * (val / rm); // posterior odds
                val = val / (1 + val); // updated posterior with current prior given eqtl = 1

                // bookkeeping
                e1g1 += val;
                e1g0 += 1 - val;
            }
        }

        e0g0 += total_snp - (e0g0 + e0g1 + e1g0 + e1g1);

        double a1_new = log(e1g1 * e0g0 / (e1g0 * e0g1));

        //printf("EM:  (%f\t%f\t%f\t%f)\t\t%f\t%f\n", e0g0, e1g1, e1g0, e0g1,a1_new, a0);
        if (fabs(a1_new - a1) < 0.01)
        {
            a1 = a1_new;
            var1 = (1.0 / e0g0 + 1.0 / e1g0 + 1.0 / e1g1 + 1.0 / e0g1);
            var0 = (1.0 / e0g1 + 1.0 / e0g0);
            break;
        }

        a1 = a1_new;
        a0 = log(e0g1 / e0g0);
        // a0 = log((e0g1+1)/(e0g0+1));

        r0 = exp(a0);
        r1 = exp(a0 + a1);
    }
    vector<double> av;
    av.push_back(a0);
    av.push_back(a1);
    av.push_back(var0);
    av.push_back(var1);
    return (av);
}





////  Setting various controller options



void controller::set_abf_piror_vec(vector<double> & eqtl_prior_vec, vector<double> & gwas_prior_vec){

    vector<double> default_prior_vec;
    
    default_prior_vec.push_back(0.1);
    default_prior_vec.push_back(0.2);
    default_prior_vec.push_back(0.4);
    default_prior_vec.push_back(0.8);
    default_prior_vec.push_back(1.6);

    if(eqtl_prior_vec.size()>0){
        eqtl_abf_prior_vec = eqtl_prior_vec;
    }else
        eqtl_abf_prior_vec = default_prior_vec;

    if(gwas_prior_vec.size()>0){
        gwas_abf_prior_vec = gwas_prior_vec;
    }else   
        gwas_abf_prior_vec = default_prior_vec;

}



void controller::set_prefix(char *str)
{

    if (strlen(str) == 0)
    {
        prefix = string("");
    }
    else
    {
        prefix = string(str) + string(".");
    }
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Legacy Code (Non-default anymore but still working with explicit argument /////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



//////// 1. input function


// load GWAS data in TORUS format 


void controller::load_gwas_torus(char *gwas_file)
{

    fprintf(stderr, "Processing complex trait data ... \n");

    istringstream ins;

    string snp_id; // important
    string sig_id;
    double prior;
    double posterior;

    string delim = "|";

    double gwas_sum = 0;
    int gwas_count = 0;
    int format_status = 0;

    vector<string> lines = readLines(string(gwas_file));
    for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it)
    {
        string line = *it;
        ins.clear();
        ins.str(line);

        if (ins >> snp_id >> sig_id >> posterior)
        {
            format_status = 1;

            if (snp_index.find(snp_id) == snp_index.end())
            {
                snp_index[snp_id] = snp_vec.size();
                snp_vec.push_back(snp_id);
            }

            if (snp2gwas_locus.find(snp_id) == snp2gwas_locus.end())
            {
                snp2gwas_locus[snp_id] = sig_id;
            }
            else
            {
                snp2gwas_locus[snp_id] += "_" + sig_id;
            }

            if (gwas_sig_index.find(sig_id) == gwas_sig_index.end())
            {
                gwas_sig_index[sig_id] = gwas_vec.size();
                sigCluster cluster;
                cluster.id = sig_id;
                cluster.cpip = 0;
                gwas_vec.push_back(cluster);
            }

            int index = gwas_sig_index[sig_id];
            gwas_vec[index].cpip += posterior;
            gwas_vec[index].snp_vec.push_back(snp_id);
            gwas_vec[index].pip_vec.push_back(posterior);
            gwas_sum += posterior;
            gwas_count++;
        }
    }
    if (format_status == 0)
    {
        fprintf(stderr, "\nError: unexpected format in complex trait file \"%s\".\n\n", gwas_file);
        exit(1);
    }

    if (gwas_sum == 0)
    {
        fprintf(stderr, "\nError: no complex trait associations detected in \"%s\". \n\n", gwas_file);
        exit(1);
    }

    gwas_pip_vec = vector<double>(snp_vec.size(), 0.0);

    for (int i = 0; i < gwas_vec.size(); i++)
    {
        for (int j = 0; j < gwas_vec[i].snp_vec.size(); j++)
        {
            string snp = gwas_vec[i].snp_vec[j];
            gwas_pip_vec[snp_index[snp]] = gwas_vec[i].pip_vec[j];
        }
    }

    if (total_snp < snp_vec.size())
        total_snp = snp_vec.size();


    fprintf(stderr, "read in %d SNPs (eQTL+gwas), %d GWAS loci, %.1f expected hits\n\n", int(snp_vec.size()), int(gwas_vec.size()), gwas_sum);

    // estimated marginal priors
    P_gwas = gwas_sum / total_snp;
    P_eqtl = P_eqtl / total_snp;
}






/////// 2. compute colocalization probabilities

// Approximate computation based on the algorithms presented in PLOS Genetics 2016; first implemented in version 1 & 2

void controller::compute_coloc_prob_legacy()
{

    fprintf(stderr, "\nComputing colocalization probabilities ... \n\n");

    string snp_file = prefix + string("enloc.snp.out");
    string sig_file = prefix + string("enloc.sig.out");
    string gen_file = prefix + string("enloc.gene.out");

    FILE *fd1 = fopen(snp_file.c_str(), "w");
    FILE *fd2 = fopen(sig_file.c_str(), "w");
    FILE *fd3 = fopen(gen_file.c_str(), "w");

    fprintf(fd1, "Signal\tSNP\tPIP_qtl\tPIP_gwas_marginal\tPIP_gwas_qtl_prior\tSCP\n");
    fprintf(fd2, "Signal\tNum_SNP\tCPIP_qtl\tCPIP_gwas_marginal\tCPIP_gwas_qtl_prior\tRCP\tLCP\n");
    fprintf(fd3, "Gene\t\tGRCP\tGLCP\n");


    double r1 = pi1_e / (1 - pi1_e);
    double r0 = pi1_ne / (1 - pi1_ne);

    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        eqtl_vec[i].coloc_prob = 0;
        eqtl_vec[i].locus_coloc_prob = 0;
        vector<double> gprob_vec_m;
        double gprob_m = 0;

        for (int k = 0; k < eqtl_vec[i].snp_vec.size(); k++)
        {
            string snp = eqtl_vec[i].snp_vec[k];
            double d = gwas_pip_vec[snp_index[snp]];
            gprob_vec_m.push_back(d);
            gprob_m += d;
        }

        // renormalize if cluster gwas cpip >= 1
        if (gprob_m > 1 - 1e-5)
        {
            for (int k = 0; k < eqtl_vec[i].snp_vec.size(); k++)
            {
                gprob_vec_m[k] = (1 - 1e-5) * (gprob_vec_m[k] / gprob_m);
            }

            gprob_m = 1 - 1e-5;
        }

        double nc = ((1 - P_gwas) / P_gwas) / (1 - gprob_m);
        double sum_bf = nc - (1 - P_gwas) / P_gwas;
        vector<double> bf_vec;

        // for consolidated locus-level colocalization
        double locus_gpip = gprob_m;
        double locus_epip = 0;
        // lcp def done

        for (int k = 0; k < eqtl_vec[i].snp_vec.size(); k++)
        {
            // lcp comp
            locus_epip += eqtl_vec[i].pip_vec[k];
            // lcp comp done
            double bf = gprob_vec_m[k] * nc;
            bf_vec.push_back(bf);
        }

        // check bf_vec values with exact method

        double gwas_cpip = 0;
        double max_scp = 0;
        string max_snp;

        // Adjust the GWAS probabilities with eQTL priors

        for (int k = 0; k < eqtl_vec[i].snp_vec.size(); k++)
        {

            string snp = eqtl_vec[i].snp_vec[k];

            // Scenario 1:  SNP-level colocalization

            double prob = bf_vec[k] / ((1 - pi1_e) / pi1_e + ((1 - pi1_e) / pi1_e) * (pi1_ne / (1 - pi1_ne)) * (sum_bf - bf_vec[k]) + bf_vec[k]);
            double p_coloc = prob * eqtl_vec[i].pip_vec[k];
            double gprob_e = p_coloc;
            eqtl_vec[i].coloc_prob += p_coloc;
            eqtl_vec[i].locus_coloc_prob += p_coloc;

            // Scenario 2: No eQTLs, single GWAS hit in the cluster

            prob = bf_vec[k] / ((1 - pi1_ne) / pi1_ne + sum_bf);
            gprob_e += prob * (1 - eqtl_vec[i].cpip);

            // Scenario 3: Locus-level colocalization -- j-th SNP is the eQTL, and k-th SNP is the GWAS
            for (int j = 0; j < eqtl_vec[i].snp_vec.size(); j++)
            {
                if (j == k)
                    continue;
                prob = bf_vec[k] / (((1 - pi1_ne) / pi1_ne) + (sum_bf - bf_vec[j]) + (pi1_e / (1 - pi1_e)) * ((1 - pi1_ne) / pi1_ne) * bf_vec[j]);
                gprob_e += prob * eqtl_vec[i].pip_vec[j];
                eqtl_vec[i].locus_coloc_prob += prob * eqtl_vec[i].pip_vec[j];
            }

            string locus_id = eqtl_vec[i].id + "(@)" + snp2gwas_locus[snp];
            if (p_coloc >= output_thresh)
                fprintf(fd1, "%15s   %15s   %7.3e %7.3e    %7.3e      %7.3e\n", locus_id.c_str(), snp.c_str(), eqtl_vec[i].pip_vec[k], gprob_vec_m[k], gprob_e, p_coloc);
            gwas_cpip += gprob_e;
            if (p_coloc >= max_scp)
            {
                max_scp = p_coloc;
                max_snp = snp;
            }
        }

        if (eqtl_vec[i].coloc_prob >= output_thresh || eqtl_vec[i].locus_coloc_prob >= output_thresh)
        {
            string locus_id = eqtl_vec[i].id + "(@)" + snp2gwas_locus[max_snp];
            fprintf(fd2, "%15s   %4d  %7.3e %7.3e    %7.3e      %7.3e\t%7.3e\n", locus_id.c_str(), int(eqtl_vec[i].snp_vec.size()), eqtl_vec[i].cpip, gprob_m, gwas_cpip, eqtl_vec[i].coloc_prob, eqtl_vec[i].locus_coloc_prob);
        }
    }

    // gene-level quantification
    map<string, double> GLCP;
    map<string, double> GRCP;

    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        string loc_id = eqtl_vec[i].id;
        int pos = loc_id.find(":");
        string gene_id = loc_id.substr(0, pos);
        if (GLCP.find(gene_id) == GLCP.end())
        {
            GLCP[gene_id] = 1;
            GRCP[gene_id] = 1;
        }

        GLCP[gene_id] *= (1 - eqtl_vec[i].locus_coloc_prob);
        GRCP[gene_id] *= (1 - eqtl_vec[i].coloc_prob);
    }

    map<string, double>::iterator it = GRCP.begin();
    while (it != GRCP.end())
    {
        string gene = it->first;
        double v1 = 1 - it->second;
        double v2 = 1 - GLCP[gene];
        fprintf(fd3, "%s\t\t%7.3e\t%7.3e\n", gene.c_str(), v1, v2);
        it++;
    }

    fclose(fd1);
    fclose(fd2);
    fclose(fd3);
}









//////////////////////////////////////////////////////////////////////////////
///////////// Code under active development  /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



void controller::load_eqtl_summary(char *eqtl_file, char *tissue)
{

    fprintf(stderr, "Processing eQTL annotations (summary statistics)... \n");

    istringstream ins;

    string target_tissue = string(tissue);

    string chr;
    string pos;
    string snp_id; // important
    string allele1;
    string allele2;
    string content;

    string delim = "|";

    int format_status = 0;

    vector<string> lines = readLines(string(eqtl_file));
    for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it)
    {

        string line = *it;

        ins.clear();
        ins.str(line);

        if (ins >> chr >> pos >> snp_id >> allele1 >> allele2 >> content)
        {

            format_status = 1;

            // passing content
            size_t pos = 0;
            string token;

            do
            {
                pos = content.find(delim);
                token = content.substr(0, pos);
                size_t pos1 = token.find("@");
                size_t pos2 = token.find("=");
                size_t pos3 = token.find("(");
                size_t pos4 = token.find(")");
                string sig_id;
                string tissue_type;
                if (pos1 != string::npos && pos2 - pos1 > 1)
                {
                    tissue_type = token.substr(pos1 + 1, pos2 - pos1 - 1);
                }

                sig_id = token.substr(0, pos1);

                // printf("sig id = %s %d %d\n",sig_id.c_str(),pos2,pos1);
                string bhat = token.substr(pos2 + 1, pos3 - pos2 - 1);
                string se = token.substr(pos3 + 1, pos4 - pos3 - 1);
                // processing

                if (strlen(tissue) == 0 || tissue_type.compare(target_tissue) == 0)
                {

                    if (snp_index.find(snp_id) == snp_index.end())
                    {
                        snp_index[snp_id] = snp_vec.size();
                        snp_vec.push_back(snp_id);
                    }

                    if (eqtl_sig_index.find(sig_id) == eqtl_sig_index.end())
                    {
                        eqtl_sig_index[sig_id] = eqtl_vec.size();
                        sigCluster cluster;
                        cluster.id = sig_id;
                        cluster.cpip = 0;
                        eqtl_vec.push_back(cluster);
                    }

                    int index = eqtl_sig_index[sig_id];

                    // std::cout << snp_id<< "    "<< sig_id<<"  "<< tissue_type<<"   "<< pip << " "<< sig_pip<< std::endl;
                }

                content.erase(0, pos + delim.length());
            } while (pos != string::npos);
        }
    }

    if (format_status == 0)
    {
        fprintf(stderr, "\nError: unexpected format in eQTL annotation file \"%s\".\n\n", eqtl_file);
        exit(1);
    }

    double sum = 0;
    for (int i = 0; i < eqtl_vec.size(); i++)
    {
        sum += eqtl_vec[i].cpip;
    }

    if (sum == 0)
    {
        fprintf(stderr, "\nError: no eQTL annotated in \"%s\". \n\n", eqtl_file);
        exit(1);
    }

    fprintf(stderr, "read in %d SNPs, %d eQTL signal clusters, %.1f expected eQTLs\n\n", int(snp_vec.size()), int(eqtl_vec.size()), sum);

    P_eqtl = sum;
}

void controller::load_gwas_summary(char *gwas_file, char *tissue)
{
}



// EM algorithm to estimate marginal prior, for summary statistics input only

double controller::torus_estimate(vector<sigCluster> & sig_vec){

    // initalize pi value
    double pi = 1.0/total_snp;
    while(1){

        // E step
        double sum = 0;
        for(int i=0;i<sig_vec.size();i++){
            sum += sig_vec[i].get_cpip(pi);
        }
        

        // M step
        double new_pi = sum/total_snp;
        if(fabs(new_pi - pi)/pi < 1e-2)
            return new_pi;
        
        pi = new_pi;
    }
        





}
    
    





