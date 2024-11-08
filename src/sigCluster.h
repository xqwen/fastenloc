#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace std;
class sigCluster {

    public:

        string gene;
        string id;  //signal id of the gene
        double cpip;
        vector<double> pip_vec;
        vector<string> snp_vec;
        vector<double> coloc_vec;    

        vector<double> log10bf_vec; // for summary statistics input
        vector<double> grid_vec; // prior effect size vector for BF computation

        double coloc_prob; 
        double locus_coloc_prob; 

    private:
        double *pip_prob;
        vector<double> prior_var_vec;

    public:


        sigCluster(){
            pip_prob = 0;
            coloc_prob = locus_coloc_prob = cpip = 0;
        }



        int impute_qtn(const gsl_rng *r);
    
        void compute_BF(double bhat, double se);

        void compute_pip(double pi1);

        double get_cpip(double pi1);

        void set_prior_var_vec(vector<double> & w2_vec){
            prior_var_vec = w2_vec;
        }

    private:
    
        double log10_weighted_sum(vector<double> &vec, vector<double> &wts);



};
