using namespace std;

#include "sigCluster.h"
#include <vector>
#include <map>
#include <string>

class controller
{

private:
    vector<sigCluster> eqtl_vec;
    vector<sigCluster> gwas_vec;
    map<string, int> eqtl_sig_index;
    map<string, int> gwas_sig_index;

    map<string, string> snp2gwas_locus;

    map<string, int> snp_index;
    vector<string> snp_vec;
    
    vector<double> gwas_pip_vec;
    
    // vector<double> gwas_abf_vec;

    string prefix;

    int ImpN;
    unsigned long seed; // random seed for MI
    int nthread;

    int outlier_control;

    int total_snp;

    double pi1_e; // gwas prior for eqtl
    double pi1_ne; // gwas prior for non-eqtl

    double a0_est;
    double a1_est;

    double cap_a1;

    // for enrichment prior
    double prior_variance;
    double P_eqtl;
    double P_gwas;

    // alternative parameterization
    double p1;
    double p2;
    double p12;

    // threshold value to output signal/snp coloc probs
    double output_thresh;
    int coloc_prob_option;

    int use_sum_stat;
    int dump_sum_conversion;


private:

    vector<double> gwas_abf_prior_vec;
    vector<double> eqtl_abf_prior_vec;
    



public:



    controller(){
        // default initialization
        P_gwas = -1;
        P_eqtl = -1;
        ImpN = 25;

        use_sum_stat = 0;

    }
    void set_coloc_prob_option(int option)
    {
        coloc_prob_option = option;
    }

    void set_a1_cap(double cap)
    {
        cap_a1 = cap;
    }

    void set_imp_num(int imp)
    {
        ImpN = imp;
    }

    void set_snp_size(int size)
    {
        total_snp = size;
    }

    void set_outlier_control(int control)
    {
        outlier_control = control;
    }

    void set_prior_variance(double pv)
    {
        prior_variance = pv;
    }

    void set_thread(int thread)
    {
        nthread = thread;
    }

    void set_random_seed(unsigned long rseed){
        seed = rseed;
    }

    void set_sum_conversion(int conv){
        dump_sum_conversion = conv;
    }


    void set_prefix(char *str);

    void set_enrich_params(double p1, double p2, double p12);
    void set_enrich_params(double a0, double a1);

    void load_eqtl(char *eqtl_file, char *tissue);
    void load_eqtl_summary(char *eqtl_file, char *tissue);


    void load_gwas(char *gwas_file, char *tissue = 0);
    void load_gwas_summary(char *gwas_file, char *tissue = 0);

    void load_combined_summary(char *summary_file);
    void set_abf_piror_vec(vector<double> & eqtl_prior_vec, vector<double> & gwas_prior_vec);

    void load_gwas_torus(char *gwas_file);

    void set_output_thresh(double value)
    {
        output_thresh = value;
    }

    void enrich_est();
    void compute_coloc_prob();


    void init_pip(); // pre-processing for summary statistics input
    double torus_estimate(vector<sigCluster> & sig_vec); 

private:

    void compute_coloc_prob_exact();
    void compute_coloc_prob_legacy();

    vector<string> readLines(const string& filename); 

    vector<double> run_EM(vector<int> &eqtl_sample);
};
