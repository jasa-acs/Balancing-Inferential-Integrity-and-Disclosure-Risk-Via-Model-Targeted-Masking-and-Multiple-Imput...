rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))
gc()

### setup
N_varw=c(1,3,5)
N_sep=c(1,3,5)
N_prop=c(0.5,1)
M=c(3,5,20)
threshold=c(0.1,0.5)
name_varw=c(expression(1),expression(3),expression(5))
name_sep=c(expression(alpha==1),expression(alpha==3),expression(alpha==5))
name_prop=c("(1) P=50%","(2) P=100%")
name_prop1=c("P=50%","P=100%")


source("/simu/function_for_plot.r")


########################
## utility ###
########################


load("/simu/summary/sum_DA.Rdata") ##MI-DA
sd_DA=inf.list$sd
lenCI_DA=inf.list$lenCI*100
rmse_DA=inf.list$rmse

load("/simu/summary/sum_synthpop_1_new.Rdata") #CART
sd_syn1=inf.list$sd
lenCI_syn1=inf.list$lenCI*100
rmse_syn1=inf.list$rmse


load("/simu/summary/sum_synthpop_2_new.Rdata") #Norm+Logit
sd_syn2=inf.list$sd
lenCI_syn2=inf.list$lenCI*100
rmse_syn2=inf.list$rmse
range(lenCI_syn1)
range(lenCI_syn2)


#####################
# ratio of rmse: MI-DA to CART 
######################

ratio_1=rmse_DA
for(p in 1:2){
	for(m in 1:3){	
		print(paste("proportion of synthesis:", name_prop1[p], "&", "M=", M[m] ))	
		print(round(rmse_DA[,,p,m]/rmse_syn1[p,m],2))
		ratio_1[,,p,m] = round(rmse_DA[,,p,m]*100/rmse_syn1[p,m],2)
	}
}
###########################
# ratio of rmse: MI-DA to Norm+Logit
###########################
ratio_2=rmse_DA
for(p in 1:2){
	for(m in 1:3){
		print(paste("proportion of synthesis:", name_prop1[p], "&", "M=", M[m] ))		
		print(round(rmse_DA[,,p,m]/rmse_syn2[p,m],2))
		ratio_2[,,p,m] = round(rmse_DA[,,p,m]*100/rmse_syn2[p,m],2)
	}
}



plot_simu0(var_DA1=ratio_1,var_DA2=ratio_2,name="overall",name_var="Ratio of RMSE x 100",shortname_var="ratio",pos1="bottomright",pos2="bottom",legend_loc=1);

### CI overlaps##

plot_simu(var_DA=lenCI_DA,var_syn1=lenCI_syn1,var_syn2=lenCI_syn2,name="overall",name_var="percentage 95% CI Overlap (%)",shortname_var="IO",pos1="bottomright",pos2="bottom",legend_loc=1);

### standardized difference ##

plot_simu(var_DA=sd_DA,var_syn1=sd_syn1,var_syn2=sd_syn2,name="overall",name_var="Standardized Difference",shortname_var="Coef_Diff",pos1="topright",pos2="top",legend_loc=1);



########################
## risk ###
########################

load("/simu/summary/risk_synthpop.Rdata")
risk_synthpop=round(risk.list$risk*100*100)

load("/simu/summary/risk_DA.Rdata")
risk_DA=round(risk.list$risk*100*100)

#threshold=0.1

var_syn1=risk_synthpop[,1,,1,1] 
var_syn2=risk_synthpop[,1,,2,1] 
var_DA=risk_DA[,,,1,,1] 

plot_simu(var_DA,var_syn1,var_syn2,name="threshold=1_MI",name_var="true match rate (%) x 100",shortname_var="rate",pos1="bottomright",pos2="bottom",legend_loc=1);

var_syn1=risk_synthpop[,1,,1,2] 
var_syn2=risk_synthpop[,1,,2,2]
var_DA=risk_DA[,,,1,,2] 

plot_simu(var_DA,var_syn1,var_syn2,name="threshold=1_meanMI",name_var="true match rate (%) x 100",shortname_var="rate",pos1="bottomright",pos2="bottom",legend_loc=2);

#threshold=0.5

var_syn1=risk_synthpop[,2,,1,1]
var_syn2=risk_synthpop[,2,,2,1] 
var_DA=risk_DA[,,,1,,1] 

plot_simu(var_DA,var_syn1,var_syn2,name="threshold=5_MI",name_var="true match rate (%) x 100 ",shortname_var="rate",pos1="topright",pos2="top",legend_loc=1);

var_syn1=risk_synthpop[,2,,1,2] 
var_syn2=risk_synthpop[,2,,2,2] 
var_DA=risk_DA[,,,1,,2]

plot_simu(var_DA,var_syn1,var_syn2,name="threshold=5_meanMI",name_var="true match rate (%) x 100",shortname_var="rate",pos1="topright",pos2="top",legend_loc=2);

