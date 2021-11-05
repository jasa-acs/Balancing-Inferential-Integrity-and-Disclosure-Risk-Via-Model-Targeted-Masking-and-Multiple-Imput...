rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

library(synthpop)

manu_dir<-data_dir<-rdir<-fdir<-"/csrg/"
result_dir<-paste(fdir,"Result_CSRG/",sep="")


############## load the data #####

csrg=read.csv(paste(data_dir,"/csrg.csv",sep=""))
csrg=csrg[,-1]
head(csrg)

source(paste(rdir,"/csrg_MI.r",sep="")) # code to impute missing data/seed is set

print(res_MI,3) # the benchmark MI results without synthesis for workdisability 
print(res_ILD_MI,3) # the benchmark MI results without synthesis for ILD

source(paste(rdir,"/syn_synthpop.r",sep="")) # load synthpop_syn_missing definition

# (1) threshold d=3 when defining a match on age

# note m=5 and r2=5:defined in csrg_MI.r

# reg: norm+logit

res_reg=synthpop_syn_missing(csrg,csrg_MI,m=m,ind_sen,method=1,threshold=3,r=3,r2=r2,res_MI,res_ILD_MI);
	
## CART

res_cart=synthpop_syn_missing(csrg,csrg_MI,m=m,ind_sen,method=2,threshold=3,r=3,r2=r2,res_MI,res_ILD_MI);

sum_synthpop.list=list(res_reg=res_reg,res_cart=res_cart)

save(sum_synthpop.list,file=paste(result_dir,"synthpop_d=3.Rdata",sep=""))


# (2) threshold d=5 when defining a match on age

# reg: norm+logit

res_reg=synthpop_syn_missing(csrg,csrg_MI,m=m,ind_sen,method=1,threshold=5,r=3,r2=r2,res_MI,res_ILD_MI);
	
## CART
res_cart=synthpop_syn_missing(csrg,csrg_MI,m=m,ind_sen,method=2,threshold=5,r=3,r2=r2,res_MI,res_ILD_MI);

sum_synthpop.list=list(res_reg=res_reg,res_cart=res_cart)

save(sum_synthpop.list,file=paste(result_dir,"synthpop_d=5.Rdata",sep=""))











