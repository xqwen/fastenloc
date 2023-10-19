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

    string prefix;

    int ImpN;
    unsigned long seed; // random seed for MI
    int nthread;

    int outlier_control;

    int total_snp;

    double pi1;
    double pi1_e;
    double pi1_ne;

    double a0_est;
    double a1_est;

    double cap_a1;

    // for enrichment prior
    double prior_variance;
    double P_eqtl;
    double P_gwas;

    // threshold value to output signal/snp coloc probs
    double output_thresh;

    int coloc_prob_option;

public:
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


    void set_prefix(char *str);

    void set_enrich_params(double p1, double p2, double p12);
    void set_enrich_params(double a0, double a1);

    void load_eqtl(char *eqtl_file, char *tissue);
    void load_gwas(char *gwas_file, char *tissue = 0);
    void load_eqtl_summary(char *eqtl_file, char *tissue);

    void load_gwas_torus(char *gwas_file);

    void set_output_thresh(double value)
    {
        output_thresh = value;
    }

    void enrich_est();
    void compute_coloc_prob();

private:
    void compute_coloc_prob_exact();
    void compute_coloc_prob_default();

    vector<string> readLines(const string& filename); 

private:
    vector<double> run_EM(vector<int> &eqtl_sample);
};
