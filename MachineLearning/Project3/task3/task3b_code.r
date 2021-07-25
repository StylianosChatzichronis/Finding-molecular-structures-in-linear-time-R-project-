#the code for calculating EM must first be run for initializing the means!


#use of stats library
# mixtools library
library("mixtools")
library("stats")
library("mvtnorm")

#EM algorithm

pdf('em.pdf')

#--------------(1)-------Initialize the means μk , covariances Σk and mixing coefficients πk-------- 

	
	#in the following comment there are the means in case the k-means algorithm is not run.
		
	#mean_Ar = matrix(ncol=3,nrows=2);
	#mean_Ar[[1,1]] = 0.5137380;mean_Ar[[1,2]] = 0.005703245;
	#mean_Ar[[2,1]] = -0.5077973;mean_Ar[[2,2]] = 0.520423435;
	#mean_Ar[[3,1]] = -0.5120889;mean_Ar[[3,2]] = -0.532858622;

	k=3;
	#construct the covariance matrices
	#add each covariance matrix to a node in the list S
	
	a = matrix(nrow=3,ncol=2)
	S = rep(list(a),3)
	for(t in 1:k){
		a= which(rnk_matrix==t)
		tr1=matrix(nrow=length(a),ncol=2)
		counter1=1;
		for (i in a){
			tr1[counter1,1]=dataset[i,1]
			tr1[counter1,2]=dataset[i,2]
			counter1= counter1+1
		}
		S[[t]]=cov(tr1)
	}
		
	N=length(dataset[,1])
	#array pk -> mixing coefficients 
	
	pk = matrix(nrow=3,ncol=1)
	for(t in 1:k){
		pk[t,] = length(which(rnk_matrix==t))/N
	}

	



par(mfrow=c(2,2))
 
plot(dataset[,1],dataset[,2],col="green",xlim=c(-1.6,1.6),ylim=c(-1.6,1.6),panel.first = grid(),xlab="x",ylab="y",pch=20,main="Initialized data from kmeans")
color_=rgb(0,0,0,alpha=0.6) 
points(mean_Ar[,1],mean_Ar[,2],col=color_,pch=4)
ellipse(mean_Ar[1,], S[[1]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)
ellipse(mean_Ar[2,], S[[2]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)
ellipse(mean_Ar[3,], S[[3]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)



#-----evaluate the initial value of the log likelihood-------


likelihood=0;



Nk =matrix(nrow=k,ncol=1)
for(t in 1:k){
	Nk[t,1] =  length(which(rnk_matrix==t));	
}

#evaluate the likelihood using the formula
likelihood=0;
for(i in 1:N){
	sum1=0;
	for(t in 1:k){
	
	b=S[t];
 	b = matrix(unlist(b),ncol = 2, byrow = TRUE);	
sum1=sum1 +pk[t,]*(  ((2*pi)^(-2/2)) )*(  (det(b)) ^(-0.5) )*exp(-0.5*t(as.vector(dataset[i,]-mean_Ar[t,]))%*%(solve(b))%*%(as.vector(dataset[i,]-mean_Ar[t,])));
	}
	likelihood = likelihood +log(sum1);
}

print(likelihood)


#----------------(2)----------------------------E step. Evaluate the responsibilities using the current parameter values-------


flag=0;
#rg is the responsibility gamma 
rg=matrix(nrow=N,ncol=k)


lh<-c(likelihood)
if(likelihood>0 ){
new_likelihood =likelihood-10000000;
}else{
new_likelihood =likelihood+10000000000
}






while(abs(new_likelihood-likelihood)>1){
if(flag==0){
flag=1;
new_likelihood =likelihood;
}



likelihood=new_likelihood;


	
	
for(i in 1:N){	
for(t in 1:k){	
		sum1=0;				

		for(j in 1:k){ 
		
			b=S[[j]];
 			b = matrix(unlist(b),ncol=2, byrow=TRUE);
			a=pk[j,1]*dmvnorm(dataset[i,], mean_Ar[j,],b);
			sum1 = sum1 + a;
		}

		b=S[[t]];
 		b = matrix(unlist(b),ncol = 2, byrow = TRUE);
	 a=pk[t,]*dmvnorm(dataset[i,], mean_Ar[t,],b);
	rg[i,t]=a/sum1;
	}
}


#----------------(3)----------M step. Re-estimate the parameters using the current responsibilities------------------------------------


#find new mean
for(t in 1:k){	
	sum1=c(0,0);
	for(i in 1:N){
		sum1= sum1+rg[i,t]*dataset[i,]
	}
	mean_Ar[t,]= sum1/Nk[t,1] 		
}


#find new Nk
for(j in 1:k){
Nk[j,1]=0;
for(i in 1:N){
	Nk[j,1] =  Nk[j,1] + rg[i,j]
}
}

#find new covariance matrix
for(t in 1:k){	
	sum1=matrix(c(0,0,0,0),nrow=2,ncol=2);
	for(i in 1:N){
		sum1= sum1+rg[i,t]*(dataset[i,]- mean_Ar[t,])%*%t(dataset[i,]- mean_Ar[t,])
	}
	S[[t]]= sum1/Nk[t,1] 		
}









#find new mixture coefficient
for(t in 1:k){	
pk[t,]=Nk[t,1]/N;
}


#(4)-------------------------Evaluate the log likelihood--------------


#evaluate the likelihood using the formula
new_likelihood=0;
for(i in 1:N){
sum1=0;
	for(t in 1:k){
	
	b=S[t];
 	b = matrix(unlist(b),ncol = 2, byrow = TRUE);	
sum1=sum1 +pk[t,]*(  (  (2*pi)^(-2/2)))* (  (det(b)) ^(-0.5) )*exp(-0.5*t(as.vector(dataset[i,]-mean_Ar[t,]))%*%(solve(b))%*%(as.vector(dataset[i,]-mean_Ar[t,])));
	}
	new_likelihood = new_likelihood +log(sum1);
}

print(new_likelihood)
lh= c(lh,new_likelihood)



#-----------------plots--------------------------------------------

#---------------------plot classes-----------------------------

plot(0,0,col="white",xlim=c(-1.6,1.6),ylim=c(-1.7,1.7),pch=20,panel.first = grid(),xlab="x",ylab="y",main=toString(likelihood))
for (i in 1:N){


color_=rgb(0.9,0,0,alpha=rg[i,1]) 
points(dataset[i,1],dataset[i,2],col=color_,pch=20)

color_=rgb(0,0.9,0,alpha=rg[i,2]) 
points(dataset[i,1],dataset[i,2],col=color_,pch=20)

color_=rgb(0,0,0.9,alpha=rg[i,3]) 
points(dataset[i,1],dataset[i,2],col=color_,pch=20)


}
color_=rgb(0,0,0,alpha=0.4) 
points(mean_Ar[,1],mean_Ar[,2],col=color_,pch=4)

ellipse(mean_Ar[1,], S[[1]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)
ellipse(mean_Ar[2,], S[[2]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)
ellipse(mean_Ar[3,], S[[3]], alpha = 0.32, npoints = 250,col=color_, newplot = FALSE, draw = TRUE)
}










#plot two-standard deviation
#---------------end of while------if there is convergence in likelihood stop









dev.off()

print(mean_Ar)
print(S)
print(pk)

