using namespace std;

#include "sigCluster.h"
#include <vector>
#include <map>
#include <string>

class controller {

    private:
        
        vector<sigCluster> eqtl_vec;
        vector<sigCluster> gwas_vec;
        map<string, int> eqtl_sig_index;
        map<string, int> gwas_sig_index;
        
        map<string, int> snp_index;
        vector<string> snp_vec;
        vector<double> gwas_pip_vec;

        string prefix;

        int ImpN;
        int nthread;
        
        double pi1;
        double pi1_e;
        double pi1_ne;
        
        double a0_est;
        double a1_est;

    public:

        void set_imp_num(int imp){
            ImpN = imp;
        }

        void set_thread(int thread){
            nthread = thread;
        }

        void set_prefix(char *str);
        
        void load_eqtl(char *eqtl_file, char *tissue);
        void load_gwas_torus(char *gwas_file);

        void enrich_est();
        void compute_coloc_prob();

    private:

        vector<double> run_EM(vector<int> & eqtl_sample);
};
