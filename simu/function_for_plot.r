plot_simu0=function(var_DA1,var_DA2,name,name_var,shortname_var,pos1,pos2,legend_loc){
fdir <- "working directory"
manu_dir<-paste(fdir,"latex/",sep="")

mysize=1.7
color=c("turquoise4","orange","royalblue")

mylty=c("solid","dashed","dotted","dotdash","longdash","twodash")
y.min=min(c(var_DA1,var_DA2))
y.max=max(c(var_DA1,var_DA2))

if(shortname_var=="IO"){
	
	y.max=100;
}

#twodash

for(m in c(1,2,3)){
	
	pdf(file=paste(manu_dir,"/plot/simu_",shortname_var,"_m=",M[m],"_",name,".pdf",sep=""),height=6,width=15)

	par(mfrow=c(1,length(N_prop)),mar=c(4.5,4.5,1.5,1.2)) ## let j==1 first 

	for(p in c(1:length(N_prop))){
		  plot(N_varw,c(var_DA1[1,1,p,m],var_DA1[2,1,p,m],var_DA1[3,1,p,m]),type="l",lwd=3,col=color[1],xlab=expression(sigma[e]^2),ylab=name_var,ylim=c(y.min,y.max),main=name_prop[p],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,xaxt = "n")

		points(N_varw,c(var_DA1[1,1,p,m],var_DA1[2,1,p,m],var_DA1[3,1,p,m]),pch=15,cex=2,col=color[1])
		axis(1, at=N_varw, labels=name_varw,cex.axis=mysize)
	    	    lines(N_varw,c(var_DA1[1,2,p,m],var_DA1[2,2,p,m],var_DA1[3,2,p,m]),lwd=3,col=color[2],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        points(N_varw,c(var_DA1[1,2,p,m],var_DA1[2,2,p,m],var_DA1[3,2,p,m]),pch=16,cex=2,col=color[2])
           lines(N_varw,c(var_DA1[1,2,p,m],var_DA1[2,2,p,m],var_DA1[3,2,p,m]),lwd=3,col=color[2],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        points(N_varw,c(var_DA1[1,3,p,m],var_DA1[2,3,p,m],var_DA1[3,3,p,m]),pch=17,cex=2,col=color[3])
        lines(N_varw,c(var_DA1[1,3,p,m],var_DA1[2,3,p,m],var_DA1[3,3,p,m]),lwd=3,col=color[3],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        #abline(h=var_syn1[p,m],lty="dotted",lwd=3)
	    #abline(h=var_syn2[p,m],lty="dashed",lwd=3)


       	points(N_varw,c(var_DA2[1,1,p,m],var_DA2[2,1,p,m],var_DA2[3,1,p,m]),pch=15,cex=2,col=color[1])
		axis(1, at=N_varw, labels=name_varw,cex.axis=mysize)
	    	    lines(N_varw,c(var_DA2[1,1,p,m],var_DA2[2,1,p,m],var_DA2[3,1,p,m]),lwd=3,col=color[1],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="dotted")
        
        points(N_varw,c(var_DA2[1,2,p,m],var_DA2[2,2,p,m],var_DA2[3,2,p,m]),pch=16,cex=2,col=color[2])
           lines(N_varw,c(var_DA2[1,2,p,m],var_DA2[2,2,p,m],var_DA2[3,2,p,m]),lwd=3,col=color[2],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="dotted")
        
        points(N_varw,c(var_DA2[1,3,p,m],var_DA2[2,3,p,m],var_DA2[3,3,p,m]),pch=17,cex=2,col=color[3])
        lines(N_varw,c(var_DA2[1,3,p,m],var_DA2[2,3,p,m],var_DA2[3,3,p,m]),lwd=3,col=color[3],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="dotted")

        if (p==legend_loc ){
		legend(pos1,legend=name_sep,col=color[1:3],pch=c(15:17),cex=mysize,pt.cex=2,bty="n")
		legend(pos2,legend=c("DA-MI to CART","DA-MI to Norm+Logit"),seg.len=2.5,lty=c("solid","dotted"),bty="n",lwd=3,cex=mysize,pt.cex=mysize)
	    }	
	}
	dev.off()

}# end of m
}#end of function

plot_simu=function(var_DA,var_syn1,var_syn2,name,name_var,shortname_var,pos1,pos2,legend_loc){
fdir <- "/working directory/"
manu_dir<-paste(fdir,"latex/",sep="")

mysize=1.7
color=c("turquoise4","orange","royalblue")

mylty=c("solid","dashed","dotted","dotdash","longdash","twodash")
y.min=min(c(var_DA,var_syn1,var_syn2))
y.max=max(c(var_DA,var_syn1,var_syn2))

if(shortname_var=="IO"){
	
	y.max=100;
}

#twodash

for(m in c(1,2,3)){
	
	pdf(file=paste(manu_dir,"/plot/simu_",shortname_var,"_m=",M[m],"_",name,".pdf",sep=""),height=6,width=15)

	par(mfrow=c(1,length(N_prop)),mar=c(4.5,4.5,1.5,1.2)) ## let j==1 first 

	for(p in c(1:length(N_prop))){
		  plot(N_varw,c(var_DA[1,1,p,m],var_DA[2,1,p,m],var_DA[3,1,p,m]),type="l",lwd=3,col=color[1],xlab=expression(sigma[e]^2),ylab=name_var,ylim=c(y.min,y.max),main=name_prop[p],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,xaxt = "n")
		points(N_varw,c(var_DA[1,1,p,m],var_DA[2,1,p,m],var_DA[3,1,p,m]),pch=15,cex=2,col=color[1])
		axis(1, at=N_varw, labels=name_varw,cex.axis=mysize)
	    	    lines(N_varw,c(var_DA[1,2,p,m],var_DA[2,2,p,m],var_DA[3,2,p,m]),lwd=3,col=color[2],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        points(N_varw,c(var_DA[1,2,p,m],var_DA[2,2,p,m],var_DA[3,2,p,m]),pch=16,cex=2,col=color[2])
           lines(N_varw,c(var_DA[1,2,p,m],var_DA[2,2,p,m],var_DA[3,2,p,m]),lwd=3,col=color[2],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        points(N_varw,c(var_DA[1,3,p,m],var_DA[2,3,p,m],var_DA[3,3,p,m]),pch=17,cex=2,col=color[3])
        lines(N_varw,c(var_DA[1,3,p,m],var_DA[2,3,p,m],var_DA[3,3,p,m]),lwd=3,col=color[3],cex.main=mysize,cex.lab=mysize,cex.sub=mysize,cex.axis=mysize,lty="solid")
        
        abline(h=var_syn1[p,m],lty="dotted",lwd=3)
	    abline(h=var_syn2[p,m],lty="dashed",lwd=3)
       
        if (p==legend_loc ){
		
		legend(pos1,legend=name_sep,col=color[1:3],pch=c(15:17),cex=mysize,pt.cex=2,bty="n")
		legend(pos2,legend=c("DA-MI","CART","Norm+Logit"),seg.len=2.5,
		lty=c("solid","dotted","dashed"),bty="n",lwd=3,cex=mysize,pt.cex=mysize)

			
		}	
	}
	dev.off()

}# end of m
}#end of function