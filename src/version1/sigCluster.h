using namespace std;
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

class sigCluster {

    public:

        string gene;
        string id;  //signal id of the gene
        double cpip;
        vector<double> pip_vec;
        vector<string> snp_vec;
        vector<double> coloc_vec;    

        double coloc_prob;

    private:
        double *pip_prob;

    public:

        sigCluster(){
            pip_prob = 0;
        }

        int impute_qtn(const gsl_rng *r);

};
