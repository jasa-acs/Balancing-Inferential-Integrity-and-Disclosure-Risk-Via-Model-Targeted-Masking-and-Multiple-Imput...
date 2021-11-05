#########################
## using synthpop  ###
#########################
rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

library(synthpop)
library(Hmisc)

fdir <- "/simu/"
run_no=as.numeric(Sys.getenv('SGE_TASK_ID'));

rdir<-paste(fdir,"R_code/",sep="")
data_dir<-paste(fdir,"Data/",sep="")
sum_dir<-paste(fdir,"summary/",sep="")
result_dir<-paste(fdir,"Result_simu/",sep="")
n=500
T_rep=100

### proposed methods  ###

N_prop=c(0.5,1)
M=c(3,5,20)
threshold=c(0.1,0.5)
name_prop=c("P=50%","P=100%")
N_method=c("CART","norm+logit")
risk_syn<-array(NA, dim=c(length(N_prop),length(threshold),length(M),T_rep,length(N_method), 2))
est_beta=mse_beta=var_beta=sdMI_beta=sd_beta=sdMI_beta=cov_beta=lenCI_beta=array(NA,dim=c(length(N_prop),length(threshold),length(M),T_rep,length(N_method),3))
true_beta_analy=0.5*c(1,1,1)

for (k in 1:length(N_prop)){
	for (l in 1:length(threshold)){	
		for (data_no in 1:T_rep){	
			
			p=N_prop[k]
			
			load(paste(result_dir,"synthpop_result_","prop=",p,"_threshold=",threshold[l],"_data_",data_no,".Rdata",sep=""))
			
			for (method in 1:length(N_method)){
				
          	for (m in 1:length(M)){
            			
            	if(m==1 & method==1){
            		res=res.list$res11
            	}else if(m==2 & method==1){
            		res=res.list$res12
            	}else if(m==3 & method==1){
            		res=res.list$res13
            	}
            	
            	if(m==1 & method==2){
            		res=res.list$res21
            	}else if(m==2 & method==2){
            		res=res.list$res22
            	}else if(m==3 & method==2){
            		res=res.list$res23
            	}
            	
            	mse_beta[k,l,m,data_no,method,] = (as.numeric(res$res_inf$EST) - true_beta_analy)^2
            	est_beta[k,l,m,data_no,method,] = res$res_inf$EST
        
			    sd_beta[k,l,m,data_no,method,] = res$res_inf$sd_beta
			            
			    sdMI_beta[k,l,m,data_no,method,] = res$res_inf$MIsd_beta
			    
			    var_beta[k,l,m,data_no,method,] = res$res_inf$var_beta
			            
			    sdMI_beta[k,l,m,data_no,method,] = res$res_inf$MIsd_beta
			             	
			    cov_beta[k,l,m,data_no,method,] = res$res_inf$cov_beta
			                
			    lenCI_beta[k,l,m,data_no,method,] = res$res_inf$lenCI_beta*0.5
			
			    risk_syn[k,l,m,data_no,method,] = c(res$res_risk$true_rate,res$res_risk_avg$true_rate)
		            	
                
            }
 
		}
	}
}
}				
                

p_syn=2

for(method in 1:length(N_method)){
	
	for(idx in 1:(p_syn+1)){
		
		rmse=apply(mse_beta[,1,,1:T_rep,method,idx],c(1,2),function(x){sqrt(mean(x))})
		sdMC=apply(est_beta[,1,,1:T_rep,method,idx],c(1,2),function(x){sqrt(var(x))})
		
		sdMI=round(apply(sdMI_beta[,1,,1:T_rep,method,idx],c(1,2),mean),3)
		sd=round(sqrt(apply(var_beta[,1,,1:T_rep,method,idx],c(1,2),mean)),3)
		sd1=round(apply(sd_beta[,1,,1:T_rep,method,idx],c(1,2),mean),3)
		lenCI=round(apply(lenCI_beta[,1,,1:T_rep,method,idx],c(1,2),mean),6)
		cov=round(apply(cov_beta[,1,,1:T_rep,method,idx],c(1,2),mean),3)
		print(cov)
		dimnames(rmse)=dimnames(sdMC)=dimnames(sdMI)=dimnames(sd)=dimnames(sdMI)=dimnames(sd)=dimnames(lenCI)=dimnames(cov)=list(c("p=0.5","p=1"),c("M=3","M=5","M=20"))
		
		inf.list=list(rmse=rmse,cov=cov,sdMC=sdMC,sd=sd,sdMI=sdMI,lenCI=lenCI)
		save(inf.list,file=paste(sum_dir,"sum_x.",idx,"_synthpop_",method,".Rdata",sep=""))
		print(c(method,idx))
		print(rmse)
	
	}
}

risk=apply(risk_syn,c(1,2,3,5,6),mean)

dimnames(risk)=list(c("P=0.5","P=1"),c(paste("threshold=",threshold[1],sep=""),paste("threshold=",threshold[2],sep="")),c("M=3","M=5","M=20"),N_method,c("MI","mean/mode(MI)"))
risk.list=list(risk=risk)
save(risk.list,file=paste(sum_dir,"risk_synthpop.Rdata",sep=""))

#### NEW ####
for(method in 1:length(N_method)){
	sd=round(apply(sd_beta[,1,,1:T_rep,method,],c(1,2),mean),3)
	lenCI=round(apply(lenCI_beta[,1,,1:T_rep,method,],c(1,2),mean),6)
	rmse=round(apply(mse_beta[,1,,1:T_rep,method,],c(1,2),function(x){sqrt(mean(x))}),6)
	dimnames(sd)=dimnames(lenCI)=dimnames(rmse)=list(c("P=0.5","P=1"),c("M=3","M=5","M=20"))
	print(sd)
	print(lenCI)
	inf.list=list(sd=sd,lenCI=lenCI,rmse=rmse)
	save(inf.list,file=paste(sum_dir,"sum_synthpop_",method,"_new.Rdata",sep=""))
}

