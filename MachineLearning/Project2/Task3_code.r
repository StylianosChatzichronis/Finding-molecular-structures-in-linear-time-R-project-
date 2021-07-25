pdf('last.pdf')
#Task a

#read csv file in the working directory
gauss_R = read.csv(file="/home/stelios/training_set_1.csv", header=TRUE, sep=";") 
gauss_G = read.csv(file="/home/stelios/training_set_2.csv", header=TRUE, sep=";")  
gauss_B = read.csv(file="/home/stelios/training_set_3.csv", header=TRUE, sep=";")  

typeof(gauss_B) #type of gauss_R,gauss_B,gauss_G is list

label_class = c("r","g","b");

mean_R = c(0,0)
mean_G = c(0,0)
mean_B = c(0,0)
label_gauss = c("-","-","-");


Nk = c(0,0,0);
#finding mean of each class


#first Mean

label_gauss[1] = levels(gauss_R[[3]][1])
length_of_gauss = lengths(gauss_R);
Nk[1] = length_of_gauss[[1]][1];

for (i in 1:Nk[1]){
 mean_R[1] =  mean_R[1] + gauss_R[[1]][i];
 mean_R[2]  =  mean_R[2] + gauss_R[[2]][i];
}

mean_R[1]= mean_R[1]/length_of_gauss[[1]][1];
mean_R[2]= mean_R[2]/length_of_gauss[[2]][1];


#second Mean

label_gauss[2] = levels(gauss_G[[3]][1])
length_of_gauss = lengths(gauss_G)
Nk[2] = length_of_gauss[[1]][1];

for (i in 1:Nk[2]){
 mean_G[1] =  mean_G[1] + gauss_G[[1]][i]
 mean_G[2]  =  mean_G[2] + gauss_G[[2]][i]
}

mean_G[1]= mean_G[1]/length_of_gauss[[1]][1]
mean_G[2]= mean_G[2]/length_of_gauss[[2]][1]


#third Mean

label_gauss[3] = levels(gauss_B[[3]][1])
length_of_gauss = lengths(gauss_B);
Nk[3] = length_of_gauss[[1]][1];

for (i in 1:Nk[3]){
 mean_B[1] =  mean_B[1] + gauss_B[[1]][i];
 mean_B[2]  =  mean_B[2] + gauss_B[[2]][i];
}

mean_B[1]= mean_B[1]/length_of_gauss[[1]][1];
mean_B[2]= mean_B[2]/length_of_gauss[[2]][1];



#all means in a list
means_of_gaussians = list(mean_R,mean_G,mean_B) #All means in a list 

cov_r=0;
cov_g=0;
cov_b=0;


# calculate the covariance matrix of the r class

length_of_gauss = lengths(gauss_R);
 
temp_a = matrix(c(0,0),nrow=1,ncol=2);
for (i in 1:Nk[1]){			
	temp_a =  matrix(c(gauss_R[[1]][i]-mean_R[1], gauss_R[[2]][i]-mean_R[2]),nrow=2,ncol=1);	
	temp_trans = t(temp_a)
	cov_r = cov_r + temp_a%*%temp_trans	
}
cov_r = cov_r/length_of_gauss[[1]][1];



# calculate the covariance matrix of the g class
length_of_gauss = lengths(gauss_G);
 
vector_a = c(0,0);
for (i in 1:Nk[2]){			
	vector_a =  c(gauss_G[[1]][i]-mean_G[1], gauss_G[[2]][i]-mean_G[2])	
	vector_trans = t(vector_a)
	cov_g = cov_g + vector_a%*%vector_trans	
}
cov_g = cov_g/length_of_gauss[[1]][1];


# calculate the covariance matrix of the b class
length_of_gauss = lengths(gauss_B);
 
vector_a = c(0,0);
for (i in 1:Nk[3]){			
	vector_a =  c(gauss_B[[1]][i]-mean_B[1], gauss_B[[2]][i]-mean_B[2])
	vector_trans = t(vector_a)
	cov_b = cov_b + vector_a%*%vector_trans	
}
cov_b = cov_b/length_of_gauss[[1]][1];
  
N = Nk[1] + Nk[2] + Nk[3]
S = (Nk[1]/N)*cov_r + (Nk[2]/N)*cov_g + (Nk[3]/N)*cov_b




#----------------------------------------------------------------
#task b





#calculating the prior
Pc_r = Nk[1] / N 
Pc_g = Nk[2] / N
Pc_b = Nk[3] / N


#calculating posterior and classification training set 1
vector_a = c(0,0);
max_o = 0;
temp=0;
classification_r= matrix(nrow=3000,ncol=1) #classification of the first dataset
post_R = matrix(nrow=3000,ncol=3) #posterior for R dataset

for (i in 1:Nk[1]){			

	vector_a =  c(gauss_R[[1]][i]-mean_R[1], gauss_R[[2]][i]-mean_R[2])
	vector_trans = t(vector_a)
	vector_trans= t(vector_a)
	dens_r = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	

	
	
	
	vector_a =  c(gauss_R[[1]][i]-mean_G[1], gauss_R[[2]][i]-mean_G[2])
	vector_trans = t(vector_a)
	dens_g = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	
			
	vector_a =  c(gauss_R[[1]][i]-mean_B[1], gauss_R[[2]][i]-mean_B[2])
	vector_trans = t(vector_a)
	dens_b = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	
	
	
	
	temp = (dens_r)*Pc_r/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_R[[i,1]] = temp;	
	max_o = temp;
	pos=1;		

	temp = (dens_g)*Pc_g/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_R[[i,2]] = temp;	
	if(temp>max_o){
	 	pos=2
		max_o =temp
	}

	temp = (dens_b*Pc_b)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_R[[i,3]] = temp;	
	if(temp>max_o){
	 	pos=3
		max_o =temp
	}

	classification_r[[i,1]] = label_class[pos];

}


#calculating posterior and classification training set 2
vector_a = c(0,0);
max_o = 0;
temp=0;
classification_g= matrix(nrow=3000,ncol=1)
post_G = matrix(nrow=3000,ncol=3) #posterior for G dataset

for (i in 1:Nk[2]){			

	vector_a =  c(gauss_G[[1]][i]-mean_R[1], gauss_G[[2]][i]-mean_R[2])
	vector_trans = t(vector_a)
	vector_trans= t(vector_a)
	dens_r = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
		
	vector_a =  c(gauss_G[[1]][i]-mean_G[1], gauss_G[[2]][i]-mean_G[2])
	vector_trans = t(vector_a)
	dens_g = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
				
	vector_a =  c(gauss_G[[1]][i]-mean_B[1], gauss_G[[2]][i]-mean_B[2])
	vector_trans = t(vector_a)
	dens_b = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
		
	temp = (dens_r*Pc_r)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_G[[i,1]] = temp;	
	max_o = temp;
	pos=1;		

	temp = (dens_g*Pc_g)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_G[[i,2]] = temp;	
	if(temp>max_o){
	 	pos=2
		max_o =temp
	}

	temp = (dens_b*Pc_b)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_G[[i,3]] = temp;	
	if(temp>max_o){
	 	pos=3
		max_o =temp
	}

	classification_g[[i,1]] = label_class[pos];

}


#calculating posterior and classification training set 3
vector_a = c(0,0);
max_o = 0;
temp=0;
classification_b= matrix(nrow=3000,ncol=1)
post_B = matrix(nrow=3000,ncol=3) #posterior for B dataset

for (i in 1:Nk[2]){			

	vector_a =  c(gauss_B[[1]][i]-mean_R[1], gauss_B[[2]][i]-mean_R[2])
	vector_trans = t(vector_a)
	vector_trans= t(vector_a)
	dens_r = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
		
	vector_a =  c(gauss_B[[1]][i]-mean_G[1], gauss_B[[2]][i]-mean_G[2])
	vector_trans = t(vector_a)
	dens_g = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
				
	vector_a =  c(gauss_B[[1]][i]-mean_B[1], gauss_B[[2]][i]-mean_B[2])
	vector_trans = t(vector_a)
	dens_b = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
		
	temp = (dens_r*Pc_r)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_B[[i,1]] = temp;	
	max_o = temp;
	pos=1;		
	
	temp = (dens_g*Pc_g)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_B[[i,2]] = temp;	
	if(temp>max_o){
	 	pos=2
		max_o =temp
	}

	temp = (dens_b*Pc_b)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_B[[i,3]] = temp;	
	if(temp>max_o){
	 	pos=3
		max_o =temp
	}

	classification_b[[i,1]] = label_class[pos];

}

#all data together
#final_class shows class of each spot
#final is the matrix with the corresponding spots

final_class=c(0,0);
final = matrix(nrow=N,ncol=2)
for (i in 1:Nk[1]){	

final[i,1]=gauss_R[[1]][i]
final[i,2]=gauss_R[[2]][i]
final_class[i]=classification_r[[i,1]]

final[i+Nk[1],1]=gauss_G[[1]][i]
final[i+Nk[1],2]=gauss_G[[2]][i]
final_class[i+Nk[1]]=classification_g[[i,1]]

final[i+Nk[1]+Nk[2],1]=gauss_B[[1]][i]
final[i+Nk[1]+Nk[2],2]=gauss_B[[2]][i]
final_class[i+Nk[1]+Nk[2]]=classification_b[[i,1]]

}


#data with their max posteriors
post = matrix(nrow=9000,ncol=1)
for (i in 1:Nk[1]){
post[i,1]=max(post_R[[i,1]],post_R[[i,2]],post_R[[i,3]])
post[i+Nk[1],1]=max(post_G[[i,1]],post_G[[i,2]],post_G[[i,3]])
post[i+Nk[1]+Nk[2],1]=max(post_B[[i,1]],post_B[[i,2]],post_B[[i,3]])	
}




#classifying the whole data from beginning to the 3 classes
#1->x,2->y,3->probability
a=0;
a=sum(final_class == "r") #sum of elements in r class
class_r=matrix(nrow=a,ncol=3)


a=0;
a=sum(final_class == "g") #sum of elements in g class
class_g=matrix(nrow=a,ncol=3)


a=0;
a=sum(final_class == "b") #sum of elements in b class
class_b=matrix(nrow=a,ncol=3)



counter1=1;
counter2=1;
counter3=1;
for (i in 1:N){
 	if(final_class[i]=="r"){
		class_r[counter1,1] = final[i,1];
		class_r[counter1,2] = final[i,2];
		class_r[counter1,3] = post[i,1];
		counter1=counter1+1;
	}
	if(final_class[i]=="g"){
		class_g[counter2,1] = final[i,1];
		class_g[counter2,2] = final[i,2];
		class_g[counter2,3] = post[i,1];
		counter2=counter2+1;
	}
	if(final_class[i]=="b"){
		class_b[counter3,1] = final[i,1];
		class_b[counter3,2] = final[i,2];
		class_b[counter3,3] = post[i,1];
		counter3=counter3+1;
	}

}

#calculation to find the outliers
e = which(class_r[,3]<0.4);
tr1=matrix(nrow=length(e),ncol=2)
counter1=1;
#spots in the boundaries
for (i in e){
tr1[counter1,1]=class_r[i,1]
tr1[counter1,2]=class_r[i,2]
counter1= counter1+1
}

#calculation to find the outliers
e = which((class_g[,3])<0.4);
tr2=matrix(nrow=length(e),ncol=2)
counter1=1;
#spots in the boundaries
for (i in e){
tr2[counter1,1]=class_g[i,1]
tr2[counter1,2]=class_g[i,2]
counter1= counter1+1
}

#calculation to find the outliers
e = which(class_b[,3]<0.4);
tr3=matrix(nrow=length(e),ncol=2)
counter1=1;
#spots in the boundaries
for (i in e){
tr3[counter1,1]=class_b[i,1]
tr3[counter1,2]=class_b[i,2]
counter1= counter1+1
}


plot(class_b[,1],class_b[,2],col="blue",xlim=c(-1.6,1.6),ylim=c(-1.7,1.7),pch=1,panel.first = grid(),xlab="x",ylab="y")
points(class_g[,1],class_g[,2],col="green",pch=1)
points(class_r[,1],class_r[,2],col="red",pch=1)

t=( colMeans(tr1) + colMeans(tr2) + colMeans(tr3) )/3
t2=( colMeans(tr1) + colMeans(tr2) )/2
lines( c(-1.7,t[1]),c(t[2],t[2]) ,lwd=1)

t=( colMeans(tr1) + colMeans(tr2) + colMeans(tr3) )/3
t2=(colMeans(tr2) + colMeans(tr3) )/2;
f= (t[2]-t2[2])/(t[1]-t2[1]) ;
x=f*10;
y=f*x
lines( c(t[1],x),c(t[2],y+2) ,lwd=1)

t=( colMeans(tr1) + colMeans(tr2) + colMeans(tr3) )/3
t2=(colMeans(tr1) + colMeans(tr3) )/2;
f= (-t[2]+t2[2])/(t[1]-t2[1]) ;
x=f*10;
y=f*x
lines( c(t[1],x),c(t[2],-y-8.5) ,lwd=1)






#------------------------------------------------------

#install.packages("plot3D") #package for plotting 2d and 3d
library("plot3D")


par(mfrow=c(2,2))
x = class_r[,1]
y = class_r[,2]
z = class_r[,3]
points3D(x,y,z)


x = class_g[,1]
y = class_g[,2]
z = class_g[,3]
points3D(x,y,z)


x = class_b[,1]
y = class_b[,2]
z = class_b[,3]
points3D(x,y,z)



e=which(final_class=="r");
x1= matrix(nrow=length(e),ncol=1) 
y1= matrix(nrow=length(e),ncol=1) 
z1= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x1[counter1,1] = final[i,1]
y1[counter1,1]  = final[i,2]
z1[counter1,1]  = post[i,1];
counter1= counter1+1
}
points3D(x1,y1,z1, col="red",xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),zlim=c(0,1))


e=which(final_class=="g");
x2= matrix(nrow=length(e),ncol=1) 
y2= matrix(nrow=length(e),ncol=1) 
z2= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x2[counter1,1] = final[i,1]
y2[counter1,1]  = final[i,2]
z2[counter1,1]  = post[i,1];
counter1= counter1+1
}
points3D(x2,y2,z2, col="green",add = TRUE)


e=which(final_class=="b");
x1= matrix(nrow=length(e),ncol=1) 
y1= matrix(nrow=length(e),ncol=1) 
z1= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x1[counter1,1] = final[i,1]
y1[counter1,1]  = final[i,2]
z1[counter1,1]  = post[i,1];
counter1= counter1+1
}
points3D(x1,y1,z1, col="blue",add = TRUE)





par(mfrow=c(1,1))
plot(class_b[,1],class_b[,2],col="blue",xlim=c(-1.6,1.6),ylim=c(-1.7,1.7),pch=20,panel.first = grid())
points(class_g[,1],class_g[,2],col="green",pch=20)
points(class_r[,1],class_r[,2],col="red",pch=20)







#==================================================================
#task3

#read csv file in the working directory
gauss_T = read.csv(file="/home/stelios/test_set.csv", header=TRUE, sep=";") 

t=length(gauss_T[,1])





#calculating posterior and classification of test set
vector_a = c(0,0);
max_o = 0;
temp=0;
classification_t= matrix(nrow=t,ncol=1) #classification of the first dataset
post_T = matrix(nrow=t,ncol=1) #Max posterior for test dataset
post_T_all = matrix(nrow=t,ncol=3) #Max posterior for test dataset


for (i in 1:t){			

	vector_a =  c(gauss_T[[1]][i]-mean_R[1], gauss_T[[2]][i]-mean_R[2])
	vector_trans = t(vector_a)
	vector_trans= t(vector_a)
	dens_r = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	

	
	
	
	vector_a =  c(gauss_T[[1]][i]-mean_G[1], gauss_T[[2]][i]-mean_G[2])
	vector_trans = t(vector_a)
	dens_g = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	
			
	vector_a =  c(gauss_T[[1]][i]-mean_B[1], gauss_T[[2]][i]-mean_B[2])
	vector_trans = t(vector_a)
	dens_b = ( 1/( (2*pi)^(3/2) ) ) * (1/((det(S))^(1/2)) ) * exp(-0.5*vector_trans%*%(solve(S))%*%vector_a)
	
	
	
	
	temp = (dens_r)*Pc_r/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_T_all[[i,1]] = temp;	
	max_o = temp;
	pos=1;		

	temp = (dens_g)*Pc_g/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_T_all[[i,2]] = temp;	

	if(temp>max_o){
	 	pos=2
		max_o =temp
	}

	temp = (dens_b*Pc_b)/(dens_r*Pc_r+dens_g*Pc_g+dens_b*Pc_b);
	post_T_all[[i,3]] = temp;

	if(temp>max_o){
	 	pos=3
		max_o =temp
	}
	
	post_T[[i,1]] = max_o ;	
	classification_t[[i,1]] = label_class[pos];

}



#---------------plotting classified training data ---------------------



par(mfrow=c(1,1))
plot(c(0,0),c(0,0),col="white",xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),pch=20,panel.first = grid(),xlab="x",ylab="y")


for (i in 1:t){
	if(classification_t[[i,1]]=="r"){
		points(gauss_T[[1]][i],gauss_T[[2]][i],col="red",pch=1)
	}
	if(classification_t[[i,1]]=="g"){
		points(gauss_T[[1]][i],gauss_T[[2]][i],col="green",pch=1)
	}
	if(classification_t[[i,1]]=="b"){
		points(gauss_T[[1]][i],gauss_T[[2]][i],col="blue",pch=1)
	}
}









plot(c(0,0),c(0,0),col="white",xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),pch=20,panel.first = grid(),xlab="x",ylab="y")

for (i in 1:t){
if(classification_t[[i,1]]=="r"){
points(gauss_T[[1]][i],gauss_T[[2]][i],col="coral",pch=20)
}
if(classification_t[[i,1]]=="g"){
points(gauss_T[[1]][i],gauss_T[[2]][i],col="chartreuse1",pch=20)
}
if(classification_t[[i,1]]=="b"){
points(gauss_T[[1]][i],gauss_T[[2]][i],col="cadetblue2",pch=20)
}
}






#calculation to find the outliers
e = which((post_T)<0.5);
tr1=matrix(nrow=length(e),ncol=2)
counter1=1;
#spots in the boundaries
for (i in e){
tr1[counter1,1]=gauss_T[[1]][i];
tr1[counter1,2]=gauss_T[[2]][i];
counter1= counter1+1
}

points(tr1,col="black",pch=5)




#======================plotting test data=================



e=which(classification_t[,1]=="r");
x1= matrix(nrow=length(e),ncol=1) 
y1= matrix(nrow=length(e),ncol=1) 
z1= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x1[counter1,1] = gauss_T[i,1]
y1[counter1,1]  = gauss_T[i,2]
z1[counter1,1]  = post_T[i,1];
counter1= counter1+1
}
points3D(x1,y1,z1, col="red",xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),zlim=c(0,1))


e=which(classification_t[,1]=="g");
x2= matrix(nrow=length(e),ncol=1) 
y2= matrix(nrow=length(e),ncol=1) 
z2= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x2[counter1,1] = gauss_T[i,1]
y2[counter1,1]  = gauss_T[i,2]
z2[counter1,1]  = post_T[i,1];
counter1= counter1+1
}
points3D(x2,y2,z2, col="green",add = TRUE)


e=which(classification_t[,1]=="b");
x1= matrix(nrow=length(e),ncol=1) 
y1= matrix(nrow=length(e),ncol=1) 
z1= matrix(nrow=length(e),ncol=1) 
counter1=1;
for (i in e){
x1[counter1,1] = gauss_T[i,1]
y1[counter1,1]  = gauss_T[i,2]
z1[counter1,1]  = post_T[i,1];
counter1= counter1+1
}
points3D(x1,y1,z1, col="blue",add = TRUE)



#==========================================================


#calculation to find the outliers
e = which((post_T)<0.5);
tr1=matrix(nrow=length(e),ncol=3)
counter1=1;
#spots in the boundaries
for (i in e){
tr1[counter1,1]=gauss_T[[1]][i];
tr1[counter1,2]=gauss_T[[2]][i];
tr1[counter1,3]=post_T[i,1];
counter1= counter1+1
}

points3D(tr1[,1],tr1[,2],tr1[,3],pch=5,xlim=c(-1,1),ylim=c(-1,1),zlim=c(0,1))




dev.off()





#outlier The resulting decision boundaries, corresponding to the minimum misclassification rate, will occur when two of the posterior probabilities (the two largest) are equal, and so will be defined by linear functions of x, and so again we have a generalized linear model.


