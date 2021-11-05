
# this file defines functions to produce the results in both work disability and ILD analyses
# for the synthetic data sets obtained using the CART and Norm+logit approaches (implemented in synthpop package)

library(synthpop)

synthpop_syn_missing=function(csrg,csrg_MI,m,ind_sen,method,threshold,r,r2,res_MI,res_ILD_MI){

	source(paste(rdir,"/combine.r",sep="")) # includes functions to combine results
	
	tt=list() ## for given MI data set, contain mean var of the estimand combined across r synthetic data sets
	tt_syn=list() ## for given MI data set, contain combined estimate for r synthetic data sets using Reiter's rule
	
	syn_dat=list()
	t=0 
	
	for(i in 1:m){
		if(method==1){
			start_time <- Sys.time()
			res_syn <- syn(data=csrg_MI[[i]][,1:11],visit.sequence=c(11:5,1:4),
			method=c("norm",rep("logreg",3),rep("",7)),m=r,proper=TRUE,print.flag=FALSE,seed=2000) 	
			time <- Sys.time()-start_time
	        print(time)
				
		}else{
			start_time <- Sys.time()
			res_syn <- syn(data=csrg_MI[[i]][,1:11],visit.sequence=c(11:5,1:4),method=c(rep("cart",4),rep("",7)),
			m=r,proper=TRUE,print.flag=FALSE,seed=2000)
			time <- Sys.time()-start_time
	        print(time)
		}
	
		for(j in 1:r){
			t=t+1
			syn_dat[[t]]=as.data.frame(res_syn$syn[[j]])
		}
		tt[[i]]=combine_qbu(M=r,MI_data=res_syn$syn[c(1:r)],p=10) 
		tt_syn[[i]]=combine_syn(M=r,MI_data=res_syn$syn[c(1:r)],p=10)
	};
	# res_reg_complete=combine_MI_syn(tt,tt_syn,m,r,res_true=res_complete);
	res_reg_MI=combine_MI_syn(tt,tt_syn,m,r,res_true=res_MI);
	res_reg_MI
	
	source(paste(rdir,"/disclose_risk.r",sep="")) # 
	
	###### ILD ############
	
	### MI (with r2) to each of syn_dat[1:mr] ###
	tt_ILD=list() 
	
	## for given MI+Syn data set
	## contain mean var of the estimand combined across r2 MI data sets
	
	for(i in 1:(m*r)){
		
		temp=cbind(syn_dat[[i]][id_ILD,1:10],csrg_ILD$OneYearILD)
		mdf_ILD=missing_data.frame(temp)
	    show(mdf_ILD)
		imput_ILD <- mi(mdf_ILD, n.chains = 3,seed=2021,parallel=TRUE,verbose=F)# 
	    csrg_ILD_MI=complete(imput_ILD,m=r2) 
	    head(csrg_ILD_MI[[1]])
	    tt_ILD[[i]]=combine_qbu(M=r2,MI_data=csrg_ILD_MI[c(1:r2)],p=10) 
	
	}
	
	res_synthpop_ILD_MI=combine_MI_MI(tt_ILD,m*r,r2,res_true=res_ILD_MI$res);
	res_synthpop_ILD_MI
	 
	#### end of ILD #######
	
	## only calculate the disclosure risk for those without missing sensitive variables
	## convert binary factor to 0,1
	syn_dat=lapply(syn_dat,convert_to01)
	
	# original_dat1: those without missing on sensitive
	syn_dat1=lapply(syn_dat,function(x){x[-ind_sen,]})
	
	risk.list=risk(original_dat=original_dat1,
	syn_dat=syn_dat1,
	M=m*r,
	threshold/sd(csrg$Age_raw),
	name_syn=colnames(original_dat1)[1:4],
	p_syn=4)
	
	##############
	
	n=dim(syn_dat1[[1]])[1]
	p=dim(syn_dat1[[1]])[2]
	
	syn_dat1_mean=as.data.frame(apply(array(unlist(syn_dat1),dim=c(n,p,m*r)),c(1,2),mean_or_mode))
	colnames(syn_dat1_mean)=colnames(original_dat1)
	
	risk_mean.list=risk(original_dat=original_dat1,
	syn_dat=list(syn_dat1_mean),
	M=1,
	threshold/sd(csrg$Age_raw),
	name_syn=colnames(original_dat)[1:4],
	p_syn=4)
	
	sum.list=list(res_reg_MI,risk.list,risk_mean.list,res_synthpop_ILD_MI)
	return(sum.list)

}#end of synthpop_syn_missing












