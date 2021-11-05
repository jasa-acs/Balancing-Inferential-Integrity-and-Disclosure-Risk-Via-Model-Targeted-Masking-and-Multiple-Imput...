

library(synthpop)
library(Hmisc)

fdir <- "working folder"
run_no=as.numeric(Sys.getenv('SGE_TASK_ID'));

rdir<-paste(fdir,"R_code/",sep="")
result_dir<-paste(fdir,"Result_simu/",sep="")
manu_dir <-paste(fdir,"latex/jasa/table/",sep="")

source(paste(rdir,"gen_data.r",sep=""))
source(paste(rdir,"for_summary.r",sep=""))
source(paste(rdir,"fun_for_summary_simu.r",sep=""))

run_no=as.numeric(Sys.getenv('SGE_TASK_ID'));

prop=c(0.5,1)#proportion of perturbation
icount=0
setup_matrix=matrix(NA,nrow=2000,ncol=4)

for (idx3 in prop){
	for (data_no in c(1:100)){  
		icount=icount+1
		setup_matrix[icount,3]=idx3	
        setup_matrix[icount,4]=data_no
	}
}			

icount
p=setup_matrix[run_no,3]
data_no=setup_matrix[run_no,4]


## current setup ##
simu=gen_data(n=500,prop=p,seed.no=data_no)
table(simu$simu_data$O)

M_data=simu$simu_data
head(M_data)
n=dim(M_data)[1]
p_syn=dim(simu$x)[2]
ind_syn=simu$ind_syn

## analysis using unperturbed data
p_analy=dim(M_data)[2]-1
M_analy=glm(O~.-1,data=M_data,family=binomial(probit))
temp_beta_analy=coef(M_analy)
summary(M_analy)
ci0=round(confint.default(M_analy),3)[2:p_analy,]
est_analy=round(coef(M_analy),3)[2:p_analy]
ci_analy=paste("(",ci0[,1],", ",ci0[,2],")",sep="")
est_analy;ci_analy

for(threshold in c(0.1,0.5)){
	
	original_dat=M_data[,-2][,c(3,4,1,2)] # remove 1's and in the order of z_resp/o, z_cov, x_binary, x_cont
	head(original_dat)
    name_syn=colnames(original_dat[,1:2])
    original_dat[,2]=as.factor(original_dat[,2])
    
    # CART 
    
	res11=gen_syn_synthpop_simu(indx=1,original_dat,M=3,name_syn,threshold,p,ind_syn) 
	res12=gen_syn_synthpop_simu(indx=1,original_dat,M=5,name_syn,threshold,p,ind_syn)
	res13=gen_syn_synthpop_simu(indx=1,original_dat,M=20,name_syn,threshold,p,ind_syn)
	
	# Norm + Logit

	res21=gen_syn_synthpop_simu(indx=2,original_dat,M=3,name_syn,threshold,p,ind_syn) 
	res22=gen_syn_synthpop_simu(indx=2,original_dat,M=5,name_syn,threshold,p,ind_syn)
	res23=gen_syn_synthpop_simu(indx=2,original_dat,M=20,name_syn,threshold,p,ind_syn)
		
	res.list=list(res11=res11,res12=res12,res13=res13,res21=res21,res22=res22,res23=res23)
	
	save(res.list,file=paste(result_dir,"synthpop_result_","prop=",p,"_threshold=",threshold,"_data_",data_no,".Rdata",sep=""))
	
}
warnings()


