rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

fdir <- "/simu/"
rdir<-paste(fdir,"R_code/",sep="")
data_dir<-paste(fdir,"Data/",sep="")
sum_dir<-paste(fdir,"summary/",sep="")
result_dir<-paste(fdir,"Result_simu/",sep="")
n=500
T_rep=100

### proposed methods  ###

N_varw=c(1,3,5)
N_sep=c(1,3,5)
N_T=c(10)
N_prop=c(0.5,1)
M=c(3,5,20)
threshold=c(0.1,0.5)
name_varw=c(expression(sigma[e]^2==N_varw[1]),expression(sigma[e]^2==N_varw[2]),expression(sigma[e]^2==N_varw[3]))
name_sep=c(expression(alpha==N_sep[1]),expression(alpha==N_sep[2]),expression(alpha==N_sep[3]))

name_prop=c("P=50%","P=100%")

risk_syn<-array(NA, dim=c(length(N_varw),length(N_sep),length(N_prop),length(threshold),length(M),T_rep, 2))
est_beta=mse_beta=var_beta=sdMI_beta=sd_beta=sdMI_beta=cov_beta=lenCI_beta=array(NA, dim=c(length(N_varw),length(N_sep),length(N_prop),length(threshold),length(M),T_rep, 3))
true_beta_analy=0.5*c(1,1,1)

for (i in 1:length(N_varw)){
	for (j in 1:length(N_sep)){
		for (k in 1:length(N_prop)){
			for (l in 1:length(threshold)){	
				for (data_no in 1:T_rep){	
					print(c(i,j,k,l))
					var_w=N_varw[i]	
					sep=N_sep[j]	
					p=N_prop[k]
					
					setup=paste("cont=",var_w,"_binary=",sep,"_prop=",p,"_data_",data_no,sep="")
					load(paste(result_dir,"DA_",setup,"_threshold=",threshold[l],".Rdata",sep=""))

            		for (m in 1:length(M)){
            			if(m==1){
            				res=res.list$res1
            			}else if(m==2){
            				res=res.list$res2
            			}else if(m==3){
            				res=res.list$res3
            			}
            			mse_beta[i,j,k,l,m,data_no,] = (as.numeric(res$res_inf_hybrid$EST)-true_beta_analy)^2
		            	
		            	est_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$EST
		            
		            	sdMI_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$MIsd_beta
		            	
		            	var_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$var_beta
		             	
		            	cov_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$cov_beta
		            			                		        
		            	sd_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$sd_beta
		            	lenCI_beta[i,j,k,l,m,data_no,] = res$res_inf_hybrid$lenCI_beta*0.5	
		            	risk_syn[i,j,k,l,m,data_no,] = c(res$res_risk_hybrid$true_rate,res$res_risk_avg_hybrid$true_rate)
		            }
            	}
			}		
		}
	}
}

####risk ###
risk=apply(risk_syn,c(1,2,3,4,5,7),mean)

dimnames(risk)=list(name_varw,name_sep,c("P=0.5","P=1"),c(paste("threshold=",threshold[1],sep=""),paste("threshold=",threshold[2],sep="")),c("M=3","M=5","M=20"),c("MI","mean/mode(MI)"))
risk.list=list(risk=risk)
save(risk.list,file=paste(sum_dir,"risk_DA_new.Rdata",sep=""))

#### inference ####

sd=round(apply(sd_beta[,,,1,,1:T_rep,],c(1,2,3,4),mean),3) #standardized differences
lenCI=round(apply(lenCI_beta[,,,1,,1:T_rep,],c(1,2,3,4),mean),6) #percentage 95% CI Overlap
rmse=apply(mse_beta[,,,1,,1:T_rep,],c(1,2,3,4),function(x){sqrt(mean(x))})#rmse

dimnames(sd)=dimnames(lenCI)=dimnames(rmse)=list(name_varw,name_sep,c("P=0.5","P=1"),c("M=3","M=5","M=20"))
print(sd)
print(lenCI)
print(rmse)
inf.list=list(sd=sd,lenCI=lenCI,rmse=rmse)
save(inf.list,file=paste(sum_dir,"sum_DA_new.Rdata",sep=""))

