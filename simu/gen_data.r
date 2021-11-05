
runif_cor=function(n,p){
         
	X1  <- runif(n)-1/2
	Z  <- runif(n)-1/2
	U1  <- rbinom(n,1,p)  
	X2  <- X1*U1 + Z*(1-U1)
    return(2*cbind(X1,X2))
}

gen_data=function(n,prop,seed.no){
	set.seed(seed.no)
	T.max=20
	beta=0.5*c(1,1,1) # fixed for now
	z_cov=runif(n,-2,2)
	p_syn=2
	u=runif_cor(n=500,p=0.5)
	#u[,1]=runif(n,-2,2)
	#u[,2]=runif(n,-2,2)
	x=matrix(NA,nrow=n,ncol=p_syn)
	x[,1]=rnorm(n,z_cov^2+u[,1],rep(1,n))
	tmp=exp(z_cov^2+x[,1]^2+x[,1]*z_cov+u[,2])
	x[,2]=rbinom(n,1,tmp/(1+tmp))
	
	## logit 
	#tmp=exp(cbind(z_cov,x)%*%beta)
	#z_resp=rbinom(n,1,tmp/(1+tmp))
	
	## probit
	w.mean=cbind(z_cov,x)%*%beta
	w=rnorm(n,w.mean,rep(1,n))
	z_resp=1*I(w>0)
	
	table(z_resp)
	ind_syn=rbinom(n,1,prop)
	
	T.max=20
	pseudo_noise=matrix(NA,nrow=n,ncol=T.max+(p_syn-1))
	for(j in 1:p_syn){
		pseudo_noise[,1:T.max]=matrix(rnorm(n*T.max,0,1),nrow=n,ncol=T.max);
		pseudo_noise[,T.max+j-1]=rnorm(n,0,1);
	}
			
	simu_data=data.frame(O=z_resp,z_cov=cbind(1,z_cov),x=x)
	simu.list=list(simu_data=simu_data,x=x,ind_syn=ind_syn,pseudo_noise=pseudo_noise)
	return(simu.list)

}#end of function definition


gen_pseudo=function(x,KK,pseudo_noise,sep,var_w,seed.no){
	set.seed(seed.no)	
	n=dim(x)[1]
	p_syn=dim(x)[2]
	T.max=20
	W=matrix(NA,nrow=n,ncol=p_syn)
	W_cont=matrix(NA,nrow=n,ncol=T.max)
	
	## generate pseudo-copies according to idx1 to idx4
	
	for (i in 1:n){
		for (j in 1:p_syn){
			# independent errors to different sensitive variables			
			if(j>1){
				W[i,j]=x[i,j]*sep[j-1] + pseudo_noise[i,(T.max+j-1)]
			}else{
				e=pseudo_noise[i,1:T.max]
				W_cont[i,]=rep(x[i,j],T.max) + sqrt(var_w)*e
			}
		
		}
		
	}
	W[,1]=apply(W_cont[,1:KK],1,sum)
	return(list(W=W,W_cont=W_cont))
}





