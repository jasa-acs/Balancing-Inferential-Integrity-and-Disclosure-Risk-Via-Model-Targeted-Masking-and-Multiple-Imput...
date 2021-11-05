#include <iostream>
#include <ctime>
#include "scythestat/rng/mersenne.h"
#include "scythestat/distributions.h"
#include "scythestat/ide.h"
#include "scythestat/la.h"
#include "scythestat/matrix.h"
#include "scythestat/rng.h"
#include "scythestat/smath.h"
#include "scythestat/stat.h"
#include "scythestat/optimize.h"
#include "scythestat/defs.h"

using namespace scythe;
using namespace std;

mersenne myrng;

void update_x_cont(Matrix<double> &temp_x, const int K, const int n, Matrix<double> x_pseudo_sum, //n X 1 matrix: summed across K multiple pseudo-copies
                     Matrix<double> &x_mean, double &temp_var_w, double &temp_var_x,
                     Matrix<double> analy_w, Matrix<double> analy_w_mean, Matrix<double> analy_cov_x, double analy_beta_x, Matrix<int> ind_syn){
    // update x
  
    for (int i=0; i < n; ++i){
        if(ind_syn(i)==1){
            double var_tilde = 1.0/(K/temp_var_w + 1.0/temp_var_x);
            
            double mu_tilde = (x_pseudo_sum(i)/temp_var_w + x_mean(i)/temp_var_x + (analy_w(i) - analy_w_mean(i) + analy_cov_x(i)*analy_beta_x)*analy_beta_x)*var_tilde;
            
            temp_x(i) = myrng.rnorm(mu_tilde, sqrt(var_tilde));
            
        }
    }
    
}
//update M_X and M_masking for W
void update_M_XandMasking_cont(Matrix<double> temp_x, const int n, const int K,//never update M_data
                         Matrix<double> x_pseudo, //n X K K multiple pseudo-copies
                         Matrix<double> &x_mean, double &temp_var_pseudo, double &temp_var_x, Matrix<double> M_cov, int p, Matrix<int> ind_syn){//M_cov covariates for X_cont in M_X
    
    // update var_w
    double rate = 0.0;
    for (int i=0; i < n; ++i){
        if(ind_syn(i)==1){
            Matrix<double> x_i(1,K,true,temp_x(i));
            
            rate = rate + ((x_pseudo(i,_) - x_i)*t(x_pseudo(i,_) - x_i))(0);
        }
    }
    temp_var_pseudo = 1.0/myrng.rgamma(0.01 + 0.5*sum(ind_syn)*K, 0.01 + 0.5*rate); //everyone has pseudo copies
    
    // update prior model: regression x ~ cov, where temp_var_x is the variance
    
    Matrix<double> M_temp1 = eye(p)%(1.0/1000.0);
    
    Matrix<double> M_temp2(p,1,true,0.0);
    
    for (int i=0; i < n; ++i){
        
        M_temp1 += (1.0/temp_var_x)%kronecker(t(M_cov(i,_)),M_cov(i,_));
        
        M_temp2 += (1.0/temp_var_x)%t(M_cov(i,_))%temp_x(i);
        
    }//end of loop over i
    
    Matrix<double> Sigma = invpd(M_temp1);
    
    Matrix<double> Mu = Sigma*M_temp2;
    
    Matrix<double> temp_beta = myrng.rmvnorm(Mu, Sigma);
    
    x_mean = M_cov*temp_beta;
    
    Matrix<double> M_resid = temp_x - x_mean;
    
    temp_var_x = 1.0/myrng.rgamma(0.01 + 0.5*n, 0.01 + 0.5*(t(M_resid)*M_resid)(0));

}

void update_x_binary(Matrix<double> &temp_x, int n, Matrix<double> wx_mean, Matrix<double> analy_w_resid, double analy_beta,
                     Matrix<double> x_pseudo, double &alpha, Matrix<int> ind_syn){
    
    for (int i=0; i < n; ++i){
        if(ind_syn(i)==1){
        
            double prob_1 = dnorm(wx_mean(i),0.0,1.0)*exp(-0.5*(analy_w_resid(i) - 1.0*analy_beta)*(analy_w_resid(i)-1.0*analy_beta))*exp(-0.5*(x_pseudo(i) - alpha)*(x_pseudo(i)-alpha));
            
            double prob_0 = dnorm(-wx_mean(i),0.0,1.0)*exp(-0.5*analy_w_resid(i)*analy_w_resid(i))*exp(-0.5*(x_pseudo(i) - 0.0)*(x_pseudo(i)-0.0));
            
            prob_1 = prob_1/(prob_1+prob_0);
            
            temp_x(i) = myrng.rbern(prob_1);//C++ automatically converts int to double
        };
    };
 
}

//update M_X and M_masking for W

void update_XandMasking_binary(Matrix<double> temp_x, const int n, Matrix<double> &wx_mean, Matrix<double> M_cov, int p, Matrix<double> &temp_wx,
                           Matrix<double> x_pseudo, double &alpha, Matrix<int> ind_syn){// M_cov in M_X
  
    //int p = M_cov.ncols();
    double pseudo_sum=0.0;
    double x_sum = 0.0;
    for (int i=0; i < n; ++i){

        if (temp_x(i)==1.0){
            temp_wx(i) = myrng.rtbnorm_combo(wx_mean(i),1,0,50);
            if(ind_syn(i)==1){
                pseudo_sum = pseudo_sum + x_pseudo(i);
                x_sum = x_sum + temp_x(i);
            }
        }else{
            temp_wx(i) = myrng.rtanorm_combo(wx_mean(i),1,0,50);
        };
        
    }//end of loop over i
   
    alpha= myrng.rnorm(pseudo_sum*(1.0/x_sum),sqrt(1.0/x_sum));
    
    Matrix<double> M_temp1 = eye(p)%(1.0/1000.0);
    
    Matrix<double> M_temp2(p,1,true,0.0);
    
    for (int i=0; i < n; ++i){
        
        M_temp1 += kronecker(t(M_cov(i,_)),M_cov(i,_));
        
        M_temp2 += t(M_cov(i,_))%temp_wx(i);
        
    }//end of loop over i
    
    Matrix<double> Sigma = invpd(M_temp1);
    
    Matrix<double> Mu = Sigma*M_temp2;
    
    Matrix<double> temp_beta = myrng.rmvnorm(Mu, Sigma);
    
    wx_mean = M_cov*temp_beta;
}

// update analysis model Z_resp ~ 1 + Z_cov + sensitive X

void update_M_analy(Matrix<double> M_o, const int n, const int p, Matrix<double>&W_mean, Matrix<double> &temp_W,
                   Matrix<double> M_cov, Matrix<double> &temp_beta){
    for (int i=0; i < n; ++i){
        if (M_o(i)==1.0){
            temp_W(i) = myrng.rtbnorm_combo(W_mean(i),1,0,50);
        }else{
            temp_W(i) = myrng.rtanorm_combo(W_mean(i),1,0,50);
        }
    }//end of loop over i
    
    Matrix<double> M_temp1 = eye(p)%(1.0/1000.0);
    
    Matrix<double> M_temp2(p,1,true,0.0);
    
    for (int i=0; i < n; ++i){
        
        M_temp1 += kronecker(t(M_cov(i,_)),M_cov(i,_));
        
        M_temp2 += t(M_cov(i,_))%temp_W(i);
        
    }//end of loop over i
    
    Matrix<double> Sigma = invpd(M_temp1);
    
    Matrix<double> Mu = Sigma*M_temp2;
    
    temp_beta = myrng.rmvnorm(Mu, Sigma);
    
    W_mean = M_cov*temp_beta;
    
}

extern "C"{

void MCMCjoint(const int* seed, const int* mcmc, const int* burnin, const int* draw, const int* n, const int* K, const int* p, const int* p_syn, const double* sum_syn, const int* ind_syn, const double* pseudo, const double* pseudo_cont, const double* data, double* temp_x, double* temp_beta_analy,
    double* temp_wx_mean, double* temp_x_mean, double* temp_analy_w_mean,
    double* save_beta_analy, double* save_x){
    
    myrng.initialize(seed[0]);
    
    int thin = draw[0]/mcmc[0];
    
    Matrix<double> M_pseudo(n[0],p_syn[0],pseudo); // note: x_pseudo is n X p_syn: pseudo-copies; for continuous, x_pseudo[,0] is sum over K multiple copies
    
     Matrix<double> M_ind_syn(n[0],1,ind_syn);
    // initialize
    
    // for imputer's model: defined sequentially
    
    Matrix<double> M_data(n[0], p[0], data);

    //to save perturbed variables: continuous age + binary 1,2,3
    
    Matrix<double> M_x(n[0],p_syn[0],temp_x); //if ind_syn[i] = 0; M_x is never updated; initialize using original un-perturbed values
   
    // analysis model:
    
    Matrix<double> M_analy_beta(n[0],p[0]-1,temp_beta_analy); //probit model Z_resp ~ 1 + Z_cov + sensitive X
    
    Matrix<double> M_analy_w(n[0],1,true,0.0);
    
    Matrix<double> M_analy_w_mean(n[0],1,temp_analy_w_mean);
    
    // specifically for binary: wx_mean ~ N(x_mean, 1) + binary x=1 if wx_mean > 0
    
    Matrix<double> wx_mean(n[0],p_syn[0]-1,temp_wx_mean);
    
    Matrix<double> M_wx(n[0],p_syn[0]-1,true,0.0);

    // specifically for continuous:
    
    Matrix<double> M_pseudo_cont(n[0],K[0],pseudo_cont);
    
    Matrix<double> x_mean(n[0],1,temp_x_mean);
    
    //initialize
    
    double var_pseudo = 1.0;
    
    double var_x = 1.0;
    
    int tot_iter = burnin[0] + draw[0];
        
    Matrix<double> x_cont(n[0],1,true,0.0);
    
    Matrix<double> x_binary(n[0],1,true,0.0);
    
    Matrix<double> wx_mean_binary(n[0],1,true,0.0);
    
    Matrix<double> M_alpha(p_syn[0]-1,1,true,2.0);
    
    Matrix<double> wx_binary(n[0],1,true,0.0);
    double alpha=2.0;
    
    for (int iter = 0; iter < tot_iter; ++iter){

        update_M_analy(M_data(_,0), n[0], p[0]-1, M_analy_w_mean, M_analy_w, M_data(0,1,n[0]-1,(p[0]-1)), M_analy_beta);

        update_M_XandMasking_cont(M_data(_,p[0]-p_syn[0]), n[0], K[0], M_pseudo_cont, x_mean, var_pseudo, var_x, M_data(0,1,n[0]-1,(p[0]-p_syn[0]-1)),
                                  (p[0]-p_syn[0]-1), M_ind_syn);
        
        for (int j = 0; j < (p_syn[0]-1); ++j){
            
            wx_mean_binary=wx_mean(_,j);
 
            update_XandMasking_binary(M_data(_,p[0]-p_syn[0]+1+j),n[0],wx_mean_binary, M_data(0,1,n[0]-1,(p[0]-p_syn[0]+j)),(p[0]-p_syn[0]+j),wx_binary,M_pseudo(_,j+1), alpha, M_ind_syn);
            
            M_alpha(j) = alpha;
            
            wx_mean(_,j) = wx_mean_binary;

        };

        if ((iter >= burnin[0]) && (iter%thin==0)) {
            
            int g = (iter-burnin[0])/thin;
            
            //obtain X_pub
            
            update_x_cont(x_cont, K[0], n[0], M_pseudo(_,0), x_mean, var_pseudo, var_x,
                          M_analy_w, M_analy_w_mean, M_data(_,p[0]-p_syn[0]), M_analy_beta(p[0]-p_syn[0]-1), M_ind_syn);
            
            for (int i=0; i < n[0]; ++i){
                if(M_ind_syn(i)==1){
                    M_x(i,0) = x_cont(i);
                };
            };
            
            for (int j = 0; j < (p_syn[0]-1); ++j){
            
                Matrix<double> analy_w_resid = M_analy_w - M_analy_w_mean + M_data(_,p[0]-p_syn[0]+1+j)%M_analy_beta(p[0]-p_syn[0]+j);
                
                alpha=M_alpha(j);

                update_x_binary(x_binary, n[0], wx_mean(_,j), analy_w_resid, M_analy_beta(p[0]-p_syn[0]+j), M_pseudo(_,j+1), alpha, M_ind_syn);
                
                for (int i=0; i < n[0]; ++i){
                    if(M_ind_syn(i)==1){
                        M_x(i,j+1) = x_binary(i);
                    };
                };
                
            }
            
            for (int j = 0; j < (p[0]-1); ++j){
                
                int idx = j*mcmc[0] + g;
                
                save_beta_analy[idx] = M_analy_beta(j);
                
            }
            
            for (int j = 0; j < n[0]*p_syn[0]; ++j){
                
                int idx = j*mcmc[0] + g;
                
                save_x[idx] = M_x(j);
                
            }
            

        }
    
    } // end of loop over MCMC iterations

} //end of MCMCjoint
    
}//end of extern "C"
