#include "sigCluster.h"
#include <math.h>
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




void sigCluster::compute_BF(double beta, double se_beta){

    if(se_beta == 0){
        log10bf_vec.push_back(0);
        return;
    }
    
    int size = prior_var_vec.size();

    double z2 = pow((beta/se_beta), 2.0);
    double v2 = pow(se_beta, 2.0);
    vector<double> rstv;
    vector<double> wtv(size,1.0/double(size));
    for(int i=0;i<size;i++){
        double w2 = prior_var_vec[i];
        double val = 0.5*log(v2/(v2+w2)) + 0.5*z2*(w2/(v2+w2));
        rstv.push_back(val/log(10));
    }

    double rst = log10_weighted_sum(rstv,wtv); 
    log10bf_vec.push_back(rst);

    return;

}


void sigCluster::compute_pip(double pi1){

    if(log10bf_vec.size()==0)
        return;
    
    vector<double> bf_vec = log10bf_vec;
    bf_vec.push_back(0);

    vector<double> pv(log10bf_vec.size(), pi1);
    pv.push_back(1-pi1);

    double log10sum = log10_weighted_sum(bf_vec, pv);

    //printf("log10sum = %f\n", log10sum);

    cpip = 0;

    //double max = 0;
    for(int i=0;i<log10bf_vec.size();i++){
        double pip = pow(10, log10(pi1)+log10bf_vec[i] - log10sum);
        pip_vec.push_back(pip);
        /*
        if(pip>max)
            max = pip;
        */
        cpip += pip;        
    }

    //printf("cpip = %f  %f\n\n", cpip, max);

    return;



}



double sigCluster::get_cpip(double pi1){

    if(log10bf_vec.size()==0)
        return 0;
    
    vector<double> bf_vec = log10bf_vec;
    bf_vec.push_back(0);

    vector<double> pv(log10bf_vec.size(), pi1);
    pv.push_back(1-pi1);

    double log10sum = log10_weighted_sum(bf_vec, pv);

    //printf("log10sum = %f\n", log10sum);

    double rst = 0;

    //double max = 0;
    for(int i=0;i<log10bf_vec.size();i++){
        double pip = pow(10, log10(pi1)+log10bf_vec[i] - log10sum);
        rst += pip;        
    }

    //printf("cpip = %f  \n", cpip);

    return rst;
}









double sigCluster::log10_weighted_sum(vector<double> &vec, vector<double> &wts){


    double max = vec[0];
    for(size_t i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(size_t i=0;i<vec.size();i++){
        sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}
