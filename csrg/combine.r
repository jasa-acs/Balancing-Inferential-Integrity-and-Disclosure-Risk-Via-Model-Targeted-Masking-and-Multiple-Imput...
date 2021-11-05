combine_syn<-function(M,MI_data,p){ 
		
	MI_est=lapply(MI_data,function(x){
		dat=x[,1:(p+1)]
		mod<-glm(dat[,p+1] ~ ., data=dat[,1:p], family=binomial(link="probit")) 
		return(summary(mod)$coefficients[2:(p+1),1:2]) #p x 2: (mean,se)
		});
	coef=sapply(MI_est,function(x){x[,1]})
	sd=sapply(MI_est,function(x){x[,2]})
	
	## Reiter's combining rule for synthesis
	est_syn=apply(coef,1,mean)
	B=apply(coef,1,var)
	W=apply(sd^2,1,mean)
	se_syn <- sqrt(W + B/M)
	
	df_syn <-(M-1)*((1+M*W/B)^2)

    cr=qt(p=0.975, df=df_syn)
    lci_syn=est_syn-cr*se_syn
    uci_syn=est_syn+cr*se_syn
    return(cbind(est_syn,se_syn,lci_syn,uci_syn))
}

convert_to01=function(data){
	
	n=dim(data)[1]
	p=dim(data)[2]
	for(j in 1:p){
		
		if(is.factor(data[,j])){
			
			data[,j]=as.numeric(data[,j])-1
		}
	}
	return(data)
	
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

mean_or_mode=function(vec){
	
	if(length(unique(vec))==2){
		return(getmode(vec));
	}else{
		return(mean(vec));
	};
	
}

combine_qbu<-function(M=r,MI_data,p=10){ 
	# by default: the last column in dat contain the response
	MI_est=lapply(MI_data,function(x){
		dat=x[,1:(p+1)]
		mod<-glm(dat[,p+1] ~ ., data=dat[,1:p], family=binomial(link="probit")) 
		return(summary(mod)$coefficients[2:(p+1),1:2]) #p x 2: (mean,se)
		});
	coef=sapply(MI_est,function(x){x[,1]})
	sd=sapply(MI_est,function(x){x[,2]})
	
	## Rubin's combining rule
	q_l=apply(coef,1,mean) # l=1,...,m m: no. MI (Qm)
	b_l=apply(coef,1,var) #(B)
	u_l=apply(sd^2,1,mean) #(W)
	
	se_MI <- sqrt(u_l + b_l * (1 + 1/M))
    df_MI <- (M-1)*(1+u_l/(b_l*(1/M+1)))^2
    cr=qt(p=0.975, df=df_MI)
    lci_MI=q_l-cr*se_MI
    uci_MI=q_l+cr*se_MI
	len_CI=uci_MI-lci_MI
	
    return(cbind(q_l,b_l,u_l,len_CI))
}


combine_MI_MI<-function(tt,m,r,res_true){
	
	# by default, res_true=cbind(est,se,lb,ub)
	
	# based on Rubin (2003) Statistica Neerlandica
	
	Q_bar=apply(sapply(tt,function(x){x[,1]}),1,mean)
	MS_w=apply(sapply(tt,function(x){x[,2]}),1,mean)
	U_bar=apply(sapply(tt,function(x){x[,3]}),1,mean)
	MS_b=r*apply(sapply(tt,function(x){x[,1]}),1,var)
	tmp=(1/r)*(1+1/m)*MS_b
	T_M=U_bar+tmp+(1-1/r)*MS_w
	v_M=(1/(m-1))*(tmp/T_M)^2+(1/(m*r-m))*((1-1/r)*MS_w/T_M)^2
	lb=Q_bar-qt(0.975,1/v_M)*sqrt(T_M)
	ub=Q_bar+qt(0.975,1/v_M)*sqrt(T_M)
		
	## standardized difference: 
	sd_diff=abs(Q_bar-res_true[,1])/res_true[,2] #by default,est,se,lb,ub
	
	## interval overlap measure
	## using Rubin (2003)
	diff=0.5*(apply(cbind(res_true[,4],ub),1,min)- apply(cbind(res_true[,3],lb),1,max))
		
	overlap=diff/(ub-lb)+diff/(res_true[,4]-res_true[,3]) 
	
    return(list(res=cbind(Q_bar,sqrt(T_M),lb,ub,sd_diff,overlap)))
}#end of function


combine_MI_syn<-function(tt,tt_syn,m,r,res_true){ 
	
	# based on Reiter (2004)
	
	q_M=apply(sapply(tt,function(x){x[,1]}),1,mean)
	b_M=apply(sapply(tt,function(x){x[,2]}),1,mean)
	u_M=apply(sapply(tt,function(x){x[,3]}),1,mean)
	B_M=apply(sapply(tt,function(x){x[,1]}),1,var)
	
	T_M=(1+1/m)*B_M-b_M/r+u_M
	v_M=1/(((1+1/m)*B_M)^2/((m-1)*T_M^2)+(b_M/r)^2/(m*(r-1)*(T_M^2)))
	
	if(min(T_M) >0){
		lb=q_M-qt(0.975,v_M)*sqrt(T_M)
		ub=q_M+qt(0.975,v_M)*sqrt(T_M)

	}else{# if T_M is negative
		T_M_adj=(1+1/m)*B_M+u_M 
		v_M_adj=(m-1)*(1+m*u_M/((m+1)*B_M))^2
		lb=q_M-qt(0.975,v_M_adj)*sqrt(T_M_adj)
	    ub=q_M+qt(0.975,v_M_adj)*sqrt(T_M_adj)
	}
	
	res_reiter=cbind(q_M,lb,ub)
	
	est=round(apply(sapply(tt_syn,function(x){x[,1]}),1,mean),3) # consistent with q_M
	
	## standardized difference: 
	sd_diff=abs(est-res_true[,1])/res_true[,2] #by default,est,se,lb,ub
	
	## interval overlap measure
	
	diff=0.5*(apply(cbind(res_true[,4],ub),1,min)- apply(cbind(res_true[,3],lb),1,max))
		
	overlap_reiter=diff/(ub-lb)+diff/(res_true[,4]-res_true[,3]) 
    
    return(list(res_reiter=res_reiter, overlap_reiter=overlap_reiter,sd_diff=sd_diff))

}#end of function


combine_MI<-function(M=5,MI_data=aa,p=10){ 
		
	MI_est=lapply(MI_data,function(x){
		dat=x[,1:(p+1)]
		mod<-glm(dat[,p+1] ~ ., data=dat[,1:p], family=binomial(link="probit")) 
		return(summary(mod)$coefficients[2:(p+1),1:2]) #p x 2: (mean,se)
		});
	coef=sapply(MI_est,function(x){x[,1]})
	sd=sapply(MI_est,function(x){x[,2]})
	
	## Rubin's combining rule
	est_MI=apply(coef,1,mean)
	B=apply(coef,1,var)
	W=apply(sd^2,1,mean)
	se_MI <- sqrt(W + B * (1 + 1/M))
	df_MI <-(M-1)*(1+W/(B*(1/M+1)))^2
    cr=qt(p=0.975, df=df_MI)
    lci_MI=est_MI-cr*se_MI
    uci_MI=est_MI+cr*se_MI
    return(cbind(est_MI,se_MI,lci_MI,uci_MI))
}

NN_masking=function(original_dat1,syn_dat1,M=m*r,threshold){ #find nearest neighbours/return id/replace

	risk.list=risk(original_dat=original_dat1,syn_dat=syn_dat1,M=M,threshold,name_syn=colnames(original_dat1)[1:4],p_syn=4);

    temp=risk.list$true_risk;
    
    if(temp > 0){
    	id_high=c(1:dim(original_dat1)[1])[risk.list$K==1]
    	id_nb=rep(NA,length(id_high))
    	t=0
    	for(i in id_high){
			t=t+1
			a=which(original_dat1[,2]==original_dat1[i,2] & original_dat1[,3]==original_dat1[i,3] 
			& original_dat1[,4]==original_dat1[i,4]);
			a_no_i=a[!a==i];
			id=which.min(abs(original_dat1[a_no_i,1]-original_dat1[i,1]));
			id_nb[t]=a_no_i[id];
			for(j in 1:(m*r)){	
				syn_dat1[[j]][i,1:4]=syn_dat1[[j]][id_nb[t],1:4]
			};
		};
    	
    }; # end of if

	
	return(list(syn_dat1=syn_dat1));

};


prop_risk_inf_missing=function(m,r,r2,res,original_dat1,ind_sen,res_MI,threshold,
                               result_dir,setup,id_ILD_ind,OneYearILD,res_ILD_MI){
	
	# ind_sen: indicator for missing sensitive
	# original_dat1: rm obs with missing sensitive
	# id_ILD: indicator for those without baseline ILD/included in ILD 
	
	syn_dat=list()
	t=0
	for(s in 1:m){
		for(j in 1:r){
			t=t+1
			syn_dat[[t]]=res[[s]][[j]]
		}
	};
	
	syn_dat1=lapply(syn_dat,function(x){x[-ind_sen,]})
	id_ILD1=id_ILD_ind[-ind_sen]
	OneYearILD1=OneYearILD[-ind_sen]
	
    syn_dat_rm=lapply(syn_dat,function(x){x[ind_sen,]}); # those are the ones not considered in risk
    id_ILD_rm=id_ILD_ind[ind_sen]
    OneYearILD_rm=OneYearILD[ind_sen]
     
    ### masking step ###	
    start_time <- Sys.time()
	masking.list=NN_masking(original_dat1,syn_dat1,M=m*r,threshold);	
	time <- Sys.time()-start_time
	print(time)
	syn_dat1=masking.list$syn_dat1	
		
	for(s in 1:m){
		for(j in 1:r){
			t=(s-1)*r+j
			res[[s]][[j]]=rbind(syn_dat1[[t]],syn_dat_rm[[t]])
		};
	};   
	
	# ## summarize the result for given r

    tt=list() ## for given MI data set, contain mean var of the estimand combined across r synthetic data sets
    tt_syn=list() ## for given MI data set, contain combined estimate for r synthetic data sets using Reiter's rule

    for(s in 1:m){
		tt[[s]]=combine_qbu(M=r,MI_data=res[[s]][1:r],p=10) 
		tt_syn[[s]]=combine_syn(M=r,MI_data=res[[s]][1:r],p=10)
    }
    
    n=dim(syn_dat1[[1]])[1]
    p=dim(syn_dat1[[1]])[2]

    syn_dat1_mean=as.data.frame(apply(array(unlist(syn_dat1),dim=c(n,p,m*r)),c(1,2),mean_or_mode))
    colnames(syn_dat1_mean)=colnames(original_dat1)
    
    ## ILD analysis ##
    syn_dat=list()
	t=0
	for(s in 1:m){
		for(j in 1:r){
			t=t+1
			syn_dat[[t]]=res[[s]][[j]]
		}
	};
    tt_ILD=list() 
    
    ## for given mr MI+Syn data set
    
    id_ILD_new=c(id_ILD1,id_ILD_rm)
    OneYearILD_new=c(OneYearILD1,OneYearILD_rm)
    
	for(i in 1:(m*r)){
	
		temp=cbind(syn_dat[[i]][id_ILD_new==1,1:10],OneYearILD_new[id_ILD_new==1])#
		mdf_ILD=missing_data.frame(temp)
	    show(mdf_ILD)
		imput_ILD <- mi(mdf_ILD, n.iter = 50, n.chains = 1, seed=2020, parallel=FALSE,verbose=FALSE)# 
	    csrg_ILD_MI=complete(imput_ILD,m=r2)
	    tt_ILD[[i]]=combine_qbu(M=r2, MI_data=csrg_ILD_MI[c(1:r2)], p=10) 
    }

    res_prop_ILD_MI=combine_MI_MI(tt_ILD,m*r,r2,res_true=res_ILD_MI$res);

    sum.list=list(combine_MI_syn(tt,tt_syn,m,r,res_true=res_MI), # sum.list[[1]]=result for work disability
               res_prop_ILD_MI=res_prop_ILD_MI, # sum.list[[2]]=result for ILD
	           masking.list$risk.list, # sum.list[[3]]=true rate using strategy I
	           risk(original_dat=original_dat1,# sum.list[[4]]=true rate using strategy II
               syn_dat=list(syn_dat1_mean),
               M=1,
               threshold,
               name_syn=colnames(original_dat1)[1:4],
               p_syn=4))
	          
	print(sum.list)
	return(sum.list)
		
}
