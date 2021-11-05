### this file is to reproduce the main tables: Tables 3-5 in the paper

rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()


library(Hmisc)#for latex

manu_dir<-data_dir<-rdir<-fdir<-"/csrg/"
result_dir<-paste(fdir,"Result_CSRG/",sep="")


###############
### Table 3 ###
###############

tab_risk=matrix(0,nrow=3,ncol=4)
colnames(tab_risk)=rep(c("\\bf MI","\\bf mean/mode(MI)"),2)
rownames(tab_risk)=c("\\bf DA-MI","\\bf CART","\\bf Norm+Logit")

# threshold d=3
load(paste(result_dir,"synthpop_d=3.Rdata",sep=""))
tab_risk[2,1:2]=round(100*c(sum_synthpop.list$res_cart[[2]]$true_rate,sum_synthpop.list$res_cart[[3]]$true_rate),3)
tab_risk[3,1:2]=round(100*c(sum_synthpop.list$res_reg[[2]]$true_rate,sum_synthpop.list$res_reg[[3]]$true_rate),3)

# threshold d=5
load(paste(result_dir,"synthpop_d=5.Rdata",sep=""))
tab_risk[2,3:4]=round(100*c(sum_synthpop.list$res_cart[[2]]$true_rate,sum_synthpop.list$res_cart[[3]]$true_rate),3)
tab_risk[3,3:4]=round(100*c(sum_synthpop.list$res_reg[[2]]$true_rate,sum_synthpop.list$res_reg[[3]]$true_rate),3)

# threshold d=3
load(paste(result_dir,"result_csrg_Age_std=1_Gender=6_NewWhite=6_MoreThanHighSchool=6_K=20_r=3_3.Rdata",sep=""))
tab_risk[1,1:2]=round(100*c(sum.list[[3]]$true_rate,sum.list[[4]]$true_rate),3)

# threshold d=5
load(paste(result_dir,"result_csrg_Age_std=1_Gender=6_NewWhite=6_MoreThanHighSchool=6_K=20_r=3_5.Rdata",sep=""))
tab_risk[1,3:4]=round(100*c(sum.list[[3]]$true_rate,sum.list[[4]]$true_rate),3)

latex(tab_risk,title="", na.blank=TRUE, table.env=FALSE, col.just=rep("c",dim(tab_risk)[2]),center="none", 
cgroup = c("\\bf threshold $d=3$","\\bf threshold $d=5$"), n.cgroup=c(2,2),
file=paste(manu_dir,"csrg_simu_risk.tex",sep=""))


### note: Tables 4-5 report the results when the threshold d=5

## load the results obtained using DA-MI approach
load(paste(result_dir,"result_csrg_Age_std=1_Gender=6_NewWhite=6_MoreThanHighSchool=6_K=20_r=3_5.Rdata",sep=""))

## load the results obtained using reg:norm+logit and cart approaches implemented in Synthpop

load(paste(result_dir,"synthpop_d=5.Rdata",sep=""))

## load the benchmark results 
load(paste(result_dir,"csrg_benchmark.Rdata",sep="")) 


mypaste=function(x){#needed for creating Tables 4 and 5
	x=round(x,3)
	paste("(",x[1],", ",x[2],")",sep="")
}


###############
### Table 4 ###
###############

N_var=c("Age","Gender","Race","Education","Diffuse","Centro Positive","Disease Duration","FVC","DLCO","TLC")

N_method=c("\\bf Unperturbed MI",
"\\bf DA-MI","\\bf CART","\\bf Norm+Logit",
"\\bf Unperturbed MI",
"\\bf DA-MI","\\bf CART","\\bf Norm+Logit")

tab=matrix(NA,nrow=2*length(N_var),ncol=8)
colnames(tab)=N_method
rownames(tab)=rep(c("Estimate","95\\% CI"),length(N_var))

iest=seq(1,19,by=2)
ici=seq(2,20,by=2)

## Work Disability 
tab[iest,1]=round(csrg.list$res_MI[,1],3)
tab[ici,1]=apply(csrg.list$res_MI[,3:4],1,mypaste)

tab[iest,2]=round(sum.list[[1]]$res_reiter[,1],3) 
tab[ici,2]=apply(sum.list[[1]]$res_reiter[,2:3],1,mypaste)

tab[iest,3]=round(sum_synthpop.list$res_cart[[1]]$res_reiter[,1],3) 
tab[ici,3]=apply(sum_synthpop.list$res_cart[[1]]$res_reiter[,2:3],1,mypaste)

tab[iest,4]=round(sum_synthpop.list$res_reg[[1]]$res_reiter[,1],3) 
tab[ici,4]=apply(sum_synthpop.list$res_reg[[1]]$res_reiter[,2:3],1,mypaste)

## one-year ILD

tab[iest,5]=round(csrg.list$res_ILD_MI$res[,1],3)
tab[ici,5]=apply(csrg.list$res_ILD_MI$res[,3:4],1,mypaste)

tab[iest,6]=round(sum.list[[2]]$res[,1],3) 
tab[ici,6]=apply(sum.list[[2]]$res[,3:4],1,mypaste)

tab[iest,7]=round(sum_synthpop.list$res_cart[[4]]$res[,1],3) 
tab[ici,7]=apply(sum_synthpop.list$res_cart[[4]]$res[,3:4],1,mypaste)

tab[iest,8]=round(sum_synthpop.list$res_reg[[4]]$res[,1],3) 
tab[ici,8]=apply(sum_synthpop.list$res_reg[[4]]$res[,3:4],1,mypaste)

latex(tab,title="", na.blank=TRUE, table.env=FALSE, col.just=rep("c",dim(tab)[2]),center="none", 
cgroup = c("\\bf Predict Work Disability","\\bf Predict Onset of ILD"), n.cgroup=c(4,4),
rgroup = N_var, n.rgroup=rep(2,length(N_var)),
file=paste(manu_dir,"csrg_simu_est.tex",sep=""))

###############
### Table 5 ###
###############

##### Standard difference + CI overlap

N_method=c("\\bf DA-MI","\\bf CART","\\bf Norm+Logit",
"\\bf DA-MI","\\bf CART","\\bf Norm+Logit")

tab=matrix(NA,nrow=2*length(N_var)+2,ncol=6)
colnames(tab)=N_method
rownames(tab)=rep(c("Std. Diff.","95\\% CI Overlap"),length(N_var)+1)

isd=seq(1,19,by=2)
ici=seq(2,20,by=2)

tmp1=round(sum.list[[1]]$sd_diff,3)
tmp2=round(sum.list[[1]]$overlap_reiter,3)
tab[isd,1]=tmp1
tab[ici,1]=tmp2
tab[21,1]=round(mean(tmp1),3);tab[22,1]=round(mean(tmp2),3)

tmp1=round(sum_synthpop.list$res_cart[[1]]$sd_diff,3)
tmp2=round(sum_synthpop.list$res_cart[[1]]$overlap_reiter,3)
tab[isd,2]=tmp1
tab[ici,2]=tmp2
tab[21,2]=round(mean(tmp1),3);tab[22,2]=round(mean(tmp2),3)

tmp1=round(sum_synthpop.list$res_reg[[1]]$sd_diff,3)
tmp2=round(sum_synthpop.list$res_reg[[1]]$overlap_reiter,3)
tab[isd,3]=tmp1
tab[ici,3]=tmp2
tab[21,3]=round(mean(tmp1),3);tab[22,3]=round(mean(tmp2),3)

### ILD

tmp1=round(sum.list[[2]]$res[,5],3)
tmp2=round(sum.list[[2]]$res[,6],3)
tab[isd,4]=tmp1
tab[ici,4]=tmp2
tab[21,4]=round(mean(tmp1),3);tab[22,4]=round(mean(tmp2),3)

tmp1=round(sum_synthpop.list$res_cart[[4]]$res[,5],3)
tmp2=round(sum_synthpop.list$res_cart[[4]]$res[,6],3)
tab[isd,5]=tmp1
tab[ici,5]=tmp2
tab[21,5]=round(mean(tmp1),3);tab[22,5]=round(mean(tmp2),3)

tmp1=round(sum_synthpop.list$res_reg[[4]]$res[,5],3)
tmp2=round(sum_synthpop.list$res_reg[[4]]$res[,6],3)
tab[isd,6]=tmp1
tab[ici,6]=tmp2
tab[21,6]=round(mean(tmp1),3);tab[22,6]=round(mean(tmp2),3)

latex(tab,title="", na.blank=TRUE, table.env=FALSE, col.just=rep("c",dim(tab)[2]),center="none", 
cgroup = c("\\bf Predict Work Disability","\\bf Predict Onset of ILD"), n.cgroup=c(3,3),
rgroup = c(N_var,"Average"), n.rgroup=rep(2,length(N_var)+1),
file=paste(manu_dir,"csrg_simu_sd_overlap.tex",sep=""))

