

rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

cdir<-manu_dir<-data_dir<-rdir<-fdir<-"/csrg/"
result_dir<-paste(fdir,"Result_CSRG/",sep="")

############## load the data #####

csrg=read.csv(paste(data_dir,"/csrg.csv",sep=""))
csrg=csrg[,-1]
head(csrg)

source(paste(rdir,"/csrg_MI.r",sep="")) # code to impute missing data
# note that m=5 is specified in csrg_MI.r
print(res_MI,3) # the benchmark MI results without synthesis for workdisability 
print(res_ILD_MI,3) # the benchmark MI results without synthesis for ILD

### load shared C++ library ###
setwd(cdir)
system ("R CMD SHLIB joint_multivariate_sensitive_as_covariate.cpp") 
dyn.load(file.path(cdir,paste("joint_multivariate_sensitive_as_covariate",.Platform$dynlib.ext,sep="")))


#### setup to generate W in MI-DA approach ###

sep=c(6,6,6)
var_w=1
KK=20 
prop=1 # proportion of perturbation
idx_1=sep[1]
idx_2=sep[2]
idx_3=sep[3]
idx_4=var_w

## for each of the m=5 MI data set obtained above 
## generate r=3 synthetic data sets

res=list()

for(s in 1:m){
	
	csrg_dat=csrg_MI[[s]][,1:11]
	
	source(paste(rdir,"/generate_w.r",sep="")) # 
	
	source(paste(rdir,"/syn_DA-MI.r",sep="")) # 

    res[[s]]=csrg_syn;

}

## for the final m*r complete synthetic data sets produced from workdisability analysis,
## impute missing in one-year ILD analysis, leading to m*r*r2 complete synthetic data sets for ILD analysis

r2=5 #no. of imputation for one-year ILD
source(paste(rdir,"/combine.r",sep="")) # 

## Note: only calculate the disclosure risk for those without missing sensitive variables
## original_dat1 are those without missing sensitive variables

source(paste(rdir,"/disclose_risk.r",sep="")) # 

r=3 #no. of synthesis in DA-MI

# threshold =5

sum.list=prop_risk_inf_missing(m,r,r2,res,original_dat1,ind_sen,res_MI,
threshold=5/sd(csrg$Age_raw),result_dir,setup,id_ILD_ind,csrg$OneYearILD,res_ILD_MI);
save(sum.list,file=paste(result_dir,"result_",setup,"_K=",KK,"_r=",r,"_5.Rdata",sep=""))

print(time)

# threshold =3
sum.list=prop_risk_inf_missing(m,r,r2,res,original_dat1,ind_sen,res_MI,
threshold=3/sd(csrg$Age_raw),result_dir,setup,id_ILD_ind,csrg$OneYearILD,res_ILD_MI);
save(sum.list,file=paste(result_dir,"result_",setup,"_K=",KK,"_r=",r,"_3.Rdata",sep=""))














