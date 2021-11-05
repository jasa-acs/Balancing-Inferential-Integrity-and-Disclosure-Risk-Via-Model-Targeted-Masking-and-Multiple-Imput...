rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

fdir <- "working_directory"
rdir<-paste(fdir,"R_code/",sep="")
cdir<-paste(fdir,"C++_code",sep="")

result_dir<-paste(fdir,"Result_simu/",sep="")
source(paste(rdir,"gen_data.r",sep=""))
source(paste(rdir,"for_summary.r",sep=""))
source(paste(rdir,"fun_for_summary_simu.r",sep=""))

setwd(cdir)
system ("R CMD SHLIB joint_multivariate_sensitive_as_covariate.cpp") 
dyn.load(file.path(cdir,paste("joint_multivariate_sensitive_as_covariate",.Platform$dynlib.ext,sep="")))

### setup ###

run_no=as.numeric(Sys.getenv('SGE_TASK_ID'));
c1=c(1,3,5)
b1=c(1,3,5)
prop=c(0.5,1)#proportion of perturbation
icount=0
setup_matrix=matrix(NA,nrow=2000,ncol=4)

for (idx1 in c1){
	for (idx2 in b1){
		for (idx3 in prop){
			for (data_no in c(1:100)){  
				icount=icount+1
				setup_matrix[icount,1]=idx1
				setup_matrix[icount,2]=idx2
				setup_matrix[icount,3]=idx3	
                setup_matrix[icount,4]=data_no
			}
		}			
	}
}
setup_matrix[1:icount,]
icount
var_w=setup_matrix[run_no,1]
sep=setup_matrix[run_no,2]
p=setup_matrix[run_no,3]
data_no=setup_matrix[run_no,4]
KK=15

## current setup ##

setup=paste("cont=",var_w,"_binary=",sep,"_prop=",p,"_data_",data_no,sep="")
print(setup)

simu=gen_data(n=500,prop=p,seed.no=data_no)
simu_w=gen_pseudo(x=simu$x,KK=KK,pseudo_noise=simu$pseudo_noise,sep=sep,var_w=var_w,seed.no=data_no)
table(simu$simu_data$O)
M_data=simu$simu_data
head(M_data)
n=dim(M_data)[1]
p_syn=dim(simu$x)[2]
ind_syn=simu$ind_syn
M_analy=glm(O~.-1,data=M_data,family=binomial(probit))
temp_beta_analy=coef(M_analy)
summary(M_analy)

temp_analy_w_mean=predict(M_analy,type="link")
cov=M_data[,2:3]

tmp=colnames(cov)[1]

for (j in 2:dim(cov)[2]){
	
	tmp=paste(tmp,"+",colnames(cov)[j],sep="")
}

M_X = lm(as.formula(paste("x.1 ~",tmp,"-1",sep="")),data=M_data)
temp_x_mean=predict(M_X,type="response")

temp_wx_mean=matrix(NA,nrow=n,ncol=p_syn-1)

for(i in 2:p_syn){

	tmp=paste(tmp,"+ x.",i-1,sep="")
	
	M_X = glm(as.formula(paste("x.",i,"~",tmp,"-1",sep="")),data=M_data,family=binomial(probit),
	control=glm.control(epsilon=1e-200,maxit=2000500, trace=F))
	
	temp_wx_mean[,i-1]=predict(M_X,type="link")

}
summary(M_X)

## MCMC setup ##
M.max=50
draw_R=1500
mcmc_R=M.max
burnin_R=5000

#####################################	
	
output <- .C("MCMCjoint",
seed=as.integer(2019),
mcmc = as.integer(mcmc_R),
burnin = as.integer(burnin_R),
draw = as.integer(draw_R),
n = as.integer(n),
K = as.integer(KK),
p = as.integer(dim(M_data)[2]), 
p_syn = as.integer(p_syn),
sum_syn = as.double(sum(ind_syn)),
ind_syn = as.integer(ind_syn),# indicator if the copy is synthesized
pseudo = as.double(simu_w$W),
pseudo_cont = as.double(simu_w$W_cont),
data = as.double(as.matrix(M_data)), # 1 + all data 
temp_x = as.double(as.matrix(simu$x)), 
temp_beta_analy = as.double(temp_beta_analy), 
temp_wx_mean = as.double(temp_wx_mean),
temp_x_mean = as.double(temp_x_mean),
temp_analy_w_mean = as.double(temp_analy_w_mean),
save_beta_analy=as.double(rep(0,mcmc_R*(dim(M_data)[2]-1))), 
save_x=as.double(rep(0,mcmc_R*n*p_syn))
)


x_syn <- array(output$save_x,dim=c(mcmc_R,n,p_syn))
syn.list=list(syn=x_syn)
save(syn.list,file=paste(result_dir,"syn_",setup,"_covariate.Rdata",sep=""))
