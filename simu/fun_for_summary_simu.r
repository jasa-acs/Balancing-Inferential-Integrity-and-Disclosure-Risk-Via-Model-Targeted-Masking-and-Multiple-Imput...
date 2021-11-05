
gen_syn_synthpop_simu=function(indx,original_dat,M,name_syn,threshold,prop,ind_syn){
	
	#syn_copies: n X (p_syn x M) in the order of name_syn
	p=dim(original_dat)[2]
	p_syn=length(name_syn)
	
    # by default,visit.sequence=1:ncol(original_dat)
    p_rest=dim(original_dat)[2]-length(name_syn)
    if(indx==2){#Norm+Logit
    	
    	res_syn <- syn(data=original_dat,visit.sequence=c(p:(p_syn+1),1:p_syn),method=c("norm",rep("logreg",p_syn-1),rep("",p_rest)),m=50,proper=TRUE,
    	print.flag=FALSE,seed=3612) 

    }else if(indx==1){ ## CART
    	
    	res_syn <- syn(data=original_dat,visit.sequence=c(p:(p_syn+1),1:p_syn),method=c(rep("cart",p_syn),rep("",p_rest)),m=50,proper=TRUE,
    	print.flag=FALSE,seed=3611)
    	
    }
    
    if(prop<1){
    	
    	for(m in 1:M){
    		res_syn$syn[[m]][ind_syn==0,]=original_dat[ind_syn==0,]	
    	}
    }
	
	res_inf=utility_inference_synthpop(res_syn,name_syn, p_syn,p=dim(original_dat)[2],M, original_dat)# return a list of CI 
	print(res_inf)
	
	res_risk=risk_realdata_synthpop(original_dat,res_syn,M,threshold,name_syn,p_syn)

	syn_avg=matrix(NA,nrow=nrow(original_dat),ncol=p_syn)
	colnames(syn_avg)=name_syn
	syn_copies=NULL
	
	for(j in 1:p_syn){
		temp=as.matrix(sapply(res_syn$syn,function(a){as.numeric(a[[name_syn[j]]])})[,1:M])#by default "0"->1 and "1"->2
		
		if(j==1){
    		syn_avg[,j]=as.matrix(apply(temp,1,mean))
    	}else{
    		syn_avg[,j]=as.matrix(apply(temp-1,1,getmode))
    	} 
    	syn_copies=cbind(syn_copies,temp-1)
    	
    }
        
    res_risk_avg=risk_realdata(original_dat,syn=syn_avg,M=1,threshold,name_syn,p_syn)
    res_risk
    res_risk_avg
    return(list(res_inf=res_inf,res_risk=res_risk,res_risk_avg=res_risk_avg))
	
}#end of function

proposed_risk_utility=function(result_dir,M_data,syn.list,M,threshold){# return a list of utility measures and risks
	
	#name_syn: names of sensitive variables in the order of continuous and then categorical
    #original_dat: original un-perturbed data frame
	#K: a vector of numbers of repeated pseudo-copies for each of the to-be-perturbed sensitive variables 
	#var_w: the variances of the known added noises to the pseudo-copies
	#threshold: to define a match for continuous variables
	
	original_dat=M_data[,-2][,c(3,4,1,2)] # remove 1's and in the order of z_resp/o, z_cov, x_binary, x_cont
	head(original_dat)
    name_syn=colnames(original_dat[,1:2])
    original_dat[,2]=as.factor(original_dat[,2])

	n=dim(original_dat)[1]
	
	p_syn <- length(name_syn)
    
    n = dim(original_dat)[1]
    
    p = dim(original_dat)[2] # total number of coefficients including intercept 
	
	# stack synthetic copies
	
	syn=matrix(NA,nrow=n,ncol=p_syn*M)
	syn0=matrix(NA,nrow=n,ncol=M)
	syn_avg=matrix(NA,nrow=n,ncol=p_syn)
    
	M.max=dim(syn.list$syn)[1]
	dim(syn.list$syn)
	
	for(j in 1:p_syn){ # syn_avg in the order of name_syn
	    
	    syn[,((j-1)*M+1):(j*M)]=t(syn.list$syn[1:M,,j])
    	
    	if(j==1){
    		syn0=syn[,((j-1)*M+1):(j*M)]
    		syn_avg[,j]=apply(syn[,((j-1)*M+1):(j*M)],1,mean)
    	}else{

    		syn_avg[,j]=apply(syn[,((j-1)*M+1):(j*M)],1,getmode)
    	} 
	    					    	
    }
    a=NULL
    for(i in 1:length(name_syn)){a=c(a,paste(name_syn[i],c(1:M),sep=""))}
    
    colnames(syn)=a
    
    print(table(original_dat[,2],syn_avg[,2]))
        
    syn_org=syn
    res_inf=utility_inference(name_syn, p_syn,p,syn=cbind(syn0,syn[,(M+1):(p_syn*M)]), M, original_dat)# return a list of CI and EST
		
    res_risk=risk_realdata(original_dat,syn,M,threshold,name_syn,p_syn)
    
    res_risk_avg=risk_realdata(original_dat,syn_avg,M=1,threshold,name_syn,p_syn)
        
    res.list=list(res_inf=res_inf,res_risk=res_risk,res_risk_avg=res_risk_avg)
    
	risk=1*I(res.list$res_risk$K==1 | res.list$res_risk_avg$K==1)	
	
	for(i in c(0,1)){
	 
		    id=c(1:n)[risk==1&syn_avg[,2]==i]
            print(original_dat[id,name_syn])
            
			temp_id=c(1:n)[syn_avg[,2]==i]
            temp=syn_avg[syn_avg[,2]==i,]		

            for(ind in id){
            	if(is.null(dim(temp))){
	            	x=abs(syn_avg[ind,1]-temp[1])
	            	ii=temp_id
	            }else{
	            	x=abs(syn_avg[ind,1]-temp[,1])
	         		ii=temp_id[x==min(x[x!=min(x)])]	
	         	}
                syn[ind,]=syn[ii,]
                syn_avg[ind,]=syn_avg[ii,]
                
	        }
	}

	res_inf_hybrid=utility_inference(name_syn, p_syn,p,syn=cbind(syn0,syn[,(M+1):(p_syn*M)]), M, original_dat)# return a list of CI and EST
    
    res_risk_hybrid=risk_realdata(original_dat,syn,M,threshold,name_syn,p_syn)
    
    res_risk_avg_hybrid=risk_realdata(original_dat,syn=syn_avg,M=1,threshold,name_syn,p_syn)
        
    res_final.list=list(res_inf=res_inf,res_risk=res_risk,res_risk_avg=res_risk_avg,
    res_inf_hybrid=res_inf_hybrid,
    res_risk_hybrid=res_risk_hybrid,res_risk_avg_hybrid=res_risk_avg_hybrid,syn_org=syn_org,syn=syn,prop=sum(risk)/n)

    return(res_final.list)	

}#end of function



