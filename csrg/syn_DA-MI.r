## this file contain the main functions to generate synthetic data sets using DA-MI
## it is to be called by run_DA_MI.r

## current setup ##
setup=paste("csrg_",names(sen_true)[1],"=",idx_4,"_",names(sen_true)[2],"=",idx_1,"_",
names(sen_true)[3],"=",idx_2,"_",names(sen_true)[4],"=",idx_3,sep="")
print(setup)

M_data=data.frame(O,cov,sen_true)
print(str(M_data))

syn_ind=rep(1,n)
M_analy=glm(O~.-1,data=M_data,family=binomial(probit))
temp_beta_analy=coef(M_analy)
summary(M_analy)

temp_analy_w_mean=predict(M_analy,type="link")
#plot(temp_analy_w_mean,O)

tmp=colnames(cov)[1]

for (j in 2:dim(cov)[2]){
	
	tmp=paste(tmp,"+",colnames(cov)[j],sep="")
}

M_X = lm(as.formula(paste(names(sen_true)[1],"~",tmp,"-1",sep="")),data=M_data)
temp_x_mean=predict(M_X,type="response")

temp_wx_mean=matrix(NA,nrow=n,ncol=p_syn-1)

for(v in 2:p_syn){

	tmp=paste(tmp,"+",colnames(sen_true)[v-1],sep="")
	
	M_X = glm(as.formula(paste(names(sen_true)[v],"~",tmp,"-1",sep="")),data=M_data,family=binomial(probit))
	
	temp_wx_mean[,v-1]=predict(M_X,type="link")

}

## MCMC setup ##

M.max=3
draw_R=3 #thin=draw_R/mcmc_R
mcmc_R=M.max
burnin_R=100 #total_iter=burnin+draw

#####################################
start_time=Sys.time()
output <- .C("MCMCjoint",
seed=as.integer(2021),
mcmc = as.integer(mcmc_R),
burnin = as.integer(burnin_R),
draw = as.integer(draw_R),
n = as.integer(n),
K = as.integer(KK),
p = as.integer(dim(M_data)[2]), 
p_syn = as.integer(p_syn),
sum_syn = as.double(n),
ind_syn = as.integer(rep(1,n)),# indicator if the copy is synthesized
pseudo = as.double(W),
pseudo_cont = as.double(W_cont),
data = as.double(as.matrix(M_data)), # 1 + all data 
temp_x = as.double(as.matrix(sen_true)), 
temp_beta_analy = as.double(temp_beta_analy), 
temp_wx_mean = as.double(temp_wx_mean),
temp_x_mean = as.double(temp_x_mean),
temp_analy_w_mean = as.double(temp_analy_w_mean),
save_beta_analy=as.double(rep(0,mcmc_R*(dim(M_data)[2]-1))), 
save_x=as.double(rep(0,mcmc_R*n*p_syn))
)
time <- Sys.time()-start_time
print(time)
beta <- array(output$save_beta_analy,dim=c(mcmc_R,dim(M_data)[2]-1))
print(round(apply(beta,2,mean),3))

x_syn <- array(output$save_x,dim=c(mcmc_R,n,p_syn))

csrg_syn=list()

for(l in 1:mcmc_R){
	csrg_syn[[l]]=csrg_dat	
	csrg_syn[[l]][,1:p_syn]=as.matrix(x_syn[l,,])
}






##########################################################

