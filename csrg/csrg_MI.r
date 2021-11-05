
## this file is to impute missing by MICE for the loaded csrg data set


##### scaling so that the coefficients are of similar magnitudes
Age_std=scale(csrg$Age_raw)
csrg$Age_std=as.vector(Age_std)


## complete case for work disability analysis 

my_ild_mod = glm(WorkDisabled~Age_std+Gender+NewWhite+MoreThanHighSchool+Diffuse+CentroPositive+DiseaseDuration+
Diffuse+BaseFVC+BaseDLCO+BaseTLC,data=csrg,family=binomial(probit))
summary(my_ild_mod)
res_complete=cbind(coef(my_ild_mod),summary(my_ild_mod)$coef[,2],confint.default(my_ild_mod))
res_complete=res_complete[-1,]

###########################################
# the original data that have missing data
###########################################

dim(csrg)
original_dat=subset(csrg,select=c(Age_std,Gender,NewWhite,MoreThanHighSchool,Diffuse,
CentroPositive,DiseaseDuration,BaseFVC,BaseDLCO, BaseTLC,WorkDisabled))

##############
#### ILD #####
##############

dim(csrg)
id_ILD=which(csrg$BaselineILD==0)
length(id_ILD)
table(is.na(csrg$BaselineILD==0))

id_ILD_ind=1*I(csrg$BaselineILD==0 & !is.na(csrg$BaselineILD))
length(id_ILD_ind)
table(id_ILD_ind)

## focus on those without baseline ILD for one-year ILD analysis
 
csrg_ILD=subset(csrg,BaselineILD==0,select=c(Age_std,Gender,NewWhite,MoreThanHighSchool,Diffuse,
CentroPositive,DiseaseDuration,BaseFVC,BaseDLCO, BaseTLC,OneYearILD))
dim(csrg_ILD)

for(i in 1:dim(csrg_ILD)[2]){
	
	print(colnames(csrg_ILD)[i])
	print(table(is.na(csrg_ILD[,i])))

}

## complete case for ILD analysis 

fit_ILD = glm(OneYearILD ~ Age_std+Gender+NewWhite+MoreThanHighSchool+Diffuse+CentroPositive+DiseaseDuration+
Diffuse+BaseFVC+BaseDLCO+BaseTLC,data=csrg_ILD,family=binomial(probit))
summary(fit_ILD)

res_ILD_complete=cbind(coef(fit_ILD),summary(fit_ILD)$coef[,2],confint.default(fit_ILD))
res_ILD_complete=res_ILD_complete[-1,]


#########
### MI ##
#########

m=5
library(mi)
mdf=missing_data.frame(original_dat)
#show(mdf)
imputations <- mi(mdf, n.iter = 50, n.chains = 1, max.minutes = 20,seed=2020,parallel=FALSE,verbose=FALSE)# 
csrg_MI=complete(imputations,m=m) 
str(csrg_MI[[m]])
head(csrg_MI[[m]])

source(paste(rdir,"/combine.r",sep="")) # code to impute missing
res_MI=combine_MI(M=m,MI_data=csrg_MI,p=10)
res_MI
res_complete


## ind_sen: id in original_dat (csrg:n=698 but more columns) that have missing sensitive variables

ind_sen=NULL

for(i in 1:4){
	ind_sen=c(ind_sen,which(is.na(original_dat[,i])))
}
original_dat1=original_dat[-ind_sen,]
print(head(original_dat1))


# the disclosure risk for those who have non-missing sensitive variables


source(paste(rdir,"/disclose_risk.r",sep="")) 
# using d=3 as cutoff to define a match for age

risk3=risk(original_dat=original_dat1,syn_dat=list(original_dat1),M=1,3/sd(csrg$Age_raw),name_syn=colnames(original_dat1)[1:4],p_syn=4)

# using d=5 as cutoff to define a match for age

risk5=risk(original_dat=original_dat1,syn_dat=list(original_dat1),M=1,5/sd(csrg$Age_raw),name_syn=colnames(original_dat1)[1:4],p_syn=4)
sum(risk5$K)


#################################################
# focus on those without baseline ILD
# among those remaining 475 subjects (csrg_ILD), 
# 79 of them had missing one-year incidence of ILD
# impute missing ILD for these 475 in csrg_MI 
#################################################

str(csrg_MI[[1]])
show(mdf)

tt=list() ## for given MI data set, contain mean var of the estimand combined across m MI data sets in ILD

r2=5
for(i in 1:m){
	
	temp=cbind(csrg_MI[[i]][id_ILD,1:10],csrg_ILD$OneYearILD)
	mdf_ILD=missing_data.frame(temp)
    #show(mdf_ILD)
	imput_ILD <- mi(mdf_ILD, n.iter = 50, n.chains = 1, max.minutes = 20,seed=2020,parallel=FALSE,verbose=FALSE)# 
    csrg_ILD_MI=complete(imput_ILD,m=r2) 
    head(csrg_ILD_MI[[i]])
    tt[[i]]=combine_qbu(M=r2,MI_data=csrg_ILD_MI[c(1:r2)],p=10) 

}
res_ILD_MI=combine_MI_MI(tt,m,r2,res_true=res_ILD_complete);
res_ILD_MI 

csrg.list=list(res_MI=res_MI,res_complete=res_complete,risk3=risk3,risk5=risk5,
res_ILD_complete=res_ILD_complete,res_ILD_MI=res_ILD_MI)
save(csrg.list,file=paste(result_dir,"csrg_benchmark.Rdata",sep=""))

