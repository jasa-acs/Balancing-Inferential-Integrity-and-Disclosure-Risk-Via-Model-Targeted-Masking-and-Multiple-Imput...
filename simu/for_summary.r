
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

utility_inference<-function(name_syn, p_syn,p,syn,M, original_dat){ 
	
	# response is not perturbed
    # original_dat is a data frame and only includes response O and all covariates
    # syn: n X (p_syn x M) matrices of M-synthetic copies stacked in the order of name_syn
    
    M_analy=glm(O~.,data=original_dat,family=binomial(probit))
    beta_analy=coef(M_analy)[2:p]
    ci_analy=confint.default(M_analy)[2:p,]
    se_analy=summary(M_analy)$coefficients[2:p,2]

    O=original_dat$O
    beta <- array(NA,dim=c(p-1,M,2))
    
    for(i in 1:M){
    	
    	for(j in 1:p_syn){#p_syn
    		
    		original_dat[,name_syn[j]]= syn[,((j-1)*M+i)]
    		
    	}
        mod<-glm(O ~ ., data=original_dat, family=binomial(link="probit")) 
        
        beta[,i,] <- summary(mod)$coefficients[2:p,1:2]
        
    }
    rd=3
    
    CI=EST=var_beta=sd_beta=cov_beta=lenCI_beta=MIvar_beta=MIsd_beta=rep(NA,p-1)
    names(CI)=names(EST)=names(var_beta)=names(sd_beta)=names(cov_beta)=names(lenCI_beta)=names(MIvar_beta)=names(MIsd_beta)=rownames(summary(mod)$coefficients)[2:p]
    
    for (pp in 1:(p-1)){   	
	     res <- reit.t(beta[pp,,1], beta[pp,,2],true_b=0.5) # true_b=beta_analy[pp]
	     lci_beta = res[6]
		 uci_beta = res[7]
		 CI[pp] = paste("(",round(lci_beta,rd),", ",round(uci_beta,rd),")",sep="")		 
		 EST[pp] = paste(round(res[1],rd))
		 var_beta[pp]=res[2]
		 
		 ## standardized beta difference as in Snoke et al. (2018)
		 sd_beta[pp]=abs(res[1]-beta_analy[pp])/se_analy[pp]
        
		 cov_beta[pp]=res[3]
	
		 ## lencCI reports interval overlap measure
		 diff=min(uci_beta,ci_analy[pp,2])- max(lci_beta,ci_analy[pp,1])
		 lenCI_beta[pp]=diff/(ci_analy[pp,2]-ci_analy[pp,1])+diff/res[4]

		 MIvar_beta[pp]=res[5]
		 MIsd_beta[pp]=sqrt(MIvar_beta[pp])
	 
	}
	return(list(CI=CI,EST=EST,var_beta=var_beta,sd_beta=sd_beta,cov_beta=cov_beta,
	lenCI_beta=lenCI_beta,MIvar_beta=MIvar_beta,MIsd_beta=MIsd_beta,beta_analy=beta_analy))

}

is.defined <- function(sym) {
  sym <- deparse(substitute(sym))
  env <- parent.frame()
  exists(sym, env)
}

utility_inference_synthpop<-function(res_syn,name_syn, p_syn,p,M, original_dat){ 
	
	# res_syn: output from running syn in synthpop package
	# response is not perturbed
    # original_dat is a data frame and only includes response O and all covariates
    # syn: n X (p_syn x M) matrices of M-synthetic copies stacked in the order of name_syn
    
    M_analy=glm(O~.,data=original_dat,family=binomial(probit))
    beta_analy=coef(M_analy)[2:p]
    ci_analy=confint.default(M_analy)[2:p,]
    se_analy=summary(M_analy)$coefficients[2:p,2]
    
    beta <- array(NA,dim=c(p-1,M,2))
    
    for(i in 1:M){
    	if(is.defined(res_syn$syn[[i]]$Age_raw)){
    		res_syn$syn[[i]]$Age_raw=res_syn$syn[[i]]$Age_raw/10
    	}
        mod<-glm(O ~ ., data=res_syn$syn[[i]], family=binomial(link="probit"))
        
        beta[,i,] <- summary(mod)$coefficients[2:p,1:2]

    }
    rd=3
    CI=EST=var_beta=sd_beta=cov_beta=lenCI_beta=MIvar_beta=MIsd_beta=rep(NA,p-1)
    names(CI)=names(EST)=names(var_beta)=names(sd_beta)=names(cov_beta)=names(lenCI_beta)=names(MIvar_beta)=names(MIsd_beta)=rownames(summary(mod)$coefficients)[2:p]
    for (pp in 1:(p-1)){
     	
	     res <- reit.t(beta[pp,,1], beta[pp,,2],true_b=0.5)#beta_analy[pp]
	     lci_beta = res[6]
		 uci_beta = res[7]
		 CI[pp] = paste("(",round(lci_beta,rd),", ",round(uci_beta,rd),")",sep="")
		 EST[pp] = paste(round(res[1],rd))
		 var_beta[pp]=res[2]
		
		 ## standardized beta difference as in Snoke et al. (2018)
		 sd_beta[pp]=abs(res[1]-beta_analy[pp])/se_analy[pp]
		 
		 cov_beta[pp]=res[3]
		 
		 #lenCI_beta reports interval overlap measure
		 diff=min(uci_beta,ci_analy[pp,2])- max(lci_beta,ci_analy[pp,1])
		 lenCI_beta[pp]=diff/(ci_analy[pp,2]-ci_analy[pp,1])+diff/res[4]
		 MIvar_beta[pp]=res[5]
		 MIsd_beta[pp]=sqrt(MIvar_beta[pp])
	 
	}
	return(list(CI=CI,EST=EST,var_beta=var_beta,sd_beta=sd_beta,cov_beta=cov_beta,
	lenCI_beta=lenCI_beta,MIvar_beta=MIvar_beta,MIsd_beta=MIsd_beta,beta_analy=beta_analy))
}

# Multiple Imputation combining rules - partially synthetic data
# http://www2.stat.duke.edu/~jerry/Papers/jasa07.pdf

reit.t<-function(results, stderrors,true_b=1){
    m<-length(results)
    betas<-results
    varbetas<-(stderrors)^2
    qbar<-mean(betas)
    ubar<-mean(varbetas)
    b<-var(betas)
    tvar<-ubar + b/m 
    
    MIvar <- (1 + 1/m)*b + ubar
    
    df <- (m-1)*((1+m*ubar/b)^2)
    
    crit<-qt(p=0.975, df=df)
    lci<-qbar - crit * sqrt(tvar)
    uci<-qbar + crit * sqrt(tvar)
 
    cov= 1*(true_b >= lci)*(true_b <= uci)
    
    lenCI = uci - lci
    
    res<-c(qbar, tvar, cov, lenCI, MIvar, lci, uci) #returns estimate, variance, lowerCI, upperCI
    res
}

risk_realdata_synthpop=function(original_dat,res_syn,M,threshold,name_syn,p_syn){
	
n=dim(original_dat)[1]
c.j<-rep(NA,n)
I.j<-rep(NA,n)
p_syn=length(name_syn)

for(i in 1:n){

	# print(i)

	match.prob<-array(NA,dim=c(n,M))
	
	for (j in 1:M){

        ### identify all matching records
        
        match=abs(original_dat[i,name_syn[1]]- res_syn$syn[[j]][,name_syn[1]])<= threshold
        
        for(k in 2:p_syn){
        	
        	match=match*(original_dat[i,name_syn[k]]==res_syn$syn[[j]][,name_syn[k]])
        	
        }
        
        ### if more than one match, intruder would pick one at random
		match.prob[,j] <- ifelse(match==1, 1/sum(match), 0)
	}

    #### calculate P(J=j)
	pr.J<-apply(match.prob,1,mean)

    ### calculate c.j etc

	c.j[i]<-length(pr.J[pr.J==max(pr.J)])
	I.j[i]<-(pr.J[i]==max(pr.J))
}

K<-1*(c.j*I.j==1)
F<-(c.j*(1-I.j)==1)
s<-length(c.j[c.j==1&is.na(c.j)==FALSE])

## true_risk

sum(na.omit(K))

####expected match risk

sum(1/c.j*I.j)

###true match rate

sum(na.omit(K))/n

###false match rate

sum(na.omit(F))/s

return(list(true_risk=sum(na.omit(K)),exp_risk=sum(1/c.j*I.j),true_rate=sum(na.omit(K))/n,false_rate=sum(na.omit(F))/s,K=K,c.j=c.j))

}# end of def


risk_realdata=function(original_dat,syn,M,threshold,name_syn,p_syn){
	
n=dim(original_dat)[1]
c.j<-rep(NA,n)
I.j<-rep(NA,n)
p_syn=length(name_syn)

for(i in 1:n){

	# print(i)

	match.prob<-array(NA,dim=c(n,M))
	
	for (j in 1:M){

        match=abs(original_dat[i,name_syn[1]]- syn[,j])<= threshold
        
        for(k in 2:p_syn){
        	
        	match=match*(original_dat[i,name_syn[k]]==syn[,((k-1)*M+j)])
        	
        }

        ### if more than one match, intruder would pick one at random
		match.prob[,j] <- ifelse(match==1, 1/sum(match), 0)
	}

    #### calculate P(J=j)
	pr.J<-apply(match.prob,1,mean)

    ### calculate c.j etc

	c.j[i]<-length(pr.J[pr.J==max(pr.J)])
	I.j[i]<-(pr.J[i]==max(pr.J))
}

K<-1*(c.j*I.j==1)
F<-(c.j*(1-I.j)==1)
s<-length(c.j[c.j==1&is.na(c.j)==FALSE])

## true_risk

sum(na.omit(K))

####expected match risk

sum(1/c.j*I.j)

###true match rate

sum(na.omit(K))/n

###false match rate

sum(na.omit(F))/s

return(list(true_risk=sum(na.omit(K)),exp_risk=sum(1/c.j*I.j),true_rate=sum(na.omit(K))/n,false_rate=sum(na.omit(F))/s,K=K,c.j=c.j))

}# end of def

