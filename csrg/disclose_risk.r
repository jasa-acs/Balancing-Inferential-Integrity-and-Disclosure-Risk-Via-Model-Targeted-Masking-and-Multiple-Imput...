

risk=function(original_dat,syn_dat,M,threshold,name_syn,p_syn){ 
	#note: syn_dat is a list
	
n=dim(original_dat)[1]
c.j<-rep(NA,n)
I.j<-rep(NA,n)
p_syn=length(name_syn)

for(i in 1:n){

	# print(i)

	match.prob<-array(NA,dim=c(n,M))
	
	for (j in 1:M){

        ### identify all matching records
        
        match=abs(original_dat[i,name_syn[1]] - syn_dat[[j]][,name_syn[1]])<= threshold
        
        for(k in 2:p_syn){
        	
        		match=match*(original_dat[i,name_syn[k]]==syn_dat[[j]][,name_syn[k]])
        	
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

return(list(true_risk=sum(na.omit(K)),exp_risk=sum(1/c.j*I.j),true_rate=sum(na.omit(K))/n,false_rate=sum(na.omit(F))/s,K=K))

}# end of def