
n=dim(csrg_dat)[1]
sen_true=data.frame(Age_std=csrg_dat$Age_std,
Gender=1*I(csrg_dat$Gender==1),
NewWhite=1*I(csrg_dat$NewWhite==1),
MoreThanHighSchool=1*I(csrg_dat$MoreThanHighSchool==1))
head(sen_true)
cov=data.frame(1,
Diffuse=1*I(csrg_dat$Diffuse==1),
CentroPositive=1*I(csrg_dat$CentroPositive==1),
DiseaseDuration=csrg_dat$DiseaseDuration,
BaseFVC=csrg_dat$BaseFVC,
BaseDLCO=csrg_dat$BaseDLCO,
BaseTLC=csrg_dat$BaseTLC)
table(csrg_dat$Diffuse)


O=1*I(csrg_dat$WorkDisabled==1)
p_syn=dim(sen_true)[2]
W=matrix(NA,nrow=n,ncol=p_syn)
T.max=50
W_cont=matrix(NA,nrow=n,ncol=T.max)

## generate pseudo-copies according to idx1 to idx4
set.seed(2020);

for (i in 1:n){
	for (j in 1:p_syn){
		# independent errors to different sensitive variables					
		if(j>1){
			t=rnorm(1);
			W[i,j]=as.numeric(sen_true[i,j])*sep[j-1] + t
		}else{
			e=rnorm(T.max,0,1)
			W_cont[i,]=rep(sen_true[i,j],T.max) + sqrt(var_w)*e
		}
	}
}
W[,1]=apply(W_cont[,1:KK],1,sum)

### no missing data pattern

ind_obs=matrix(1,nrow=n,ncol=p_syn)

dat=list(csrg_dat=csrg_dat,n=n,sen_true=sen_true,prop=1,O=O,cov=cov)





