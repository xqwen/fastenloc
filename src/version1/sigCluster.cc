#include "sigCluster.h"
#include <gsl/gsl_randist.h>

int sigCluster::impute_qtn(const gsl_rng *r){

    if(pip_prob==0){
        if(pip_vec.size() == 0)
            return -1;
        else{
            pip_prob = new double[pip_vec.size()+1];
            for(int i=0;i<pip_vec.size();i++){
                pip_prob[i] = pip_vec[i];
            }
            pip_prob[pip_vec.size()] = 1 - cpip;
        }   
    }
    size_t K = pip_vec.size()+1;
    unsigned int *sample = new unsigned int[pip_vec.size()+1];
    gsl_ran_multinomial(r,  K, 1, pip_prob, sample);
    for(int i=0;i<pip_vec.size();i++){
        if(sample[i] == 1)
            return i;
    }

    return -1;

}
