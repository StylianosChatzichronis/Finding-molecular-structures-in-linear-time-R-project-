pdf('standardizing_distant.pdf')
#install robustHD library for doing standardizing
#install.packages("robustHD")
library("robustHD")





#read csv file in the working directory
set1 = read.csv(file="/home/stelios/training_set_1.csv", header=TRUE, sep=";") 
set2 = read.csv(file="/home/stelios/training_set_2.csv", header=TRUE, sep=";")  
set3 = read.csv(file="/home/stelios/training_set_3.csv", header=TRUE, sep=";")  



#-------------------------------initializing-----------------------------


#we know that there are 3 clusters. 
k=3;
#merge all dataset into one list



#case for taking means from each dataset




#all sets into one dataset
dataset = matrix(nrow=length(set1[,1])+length(set2[,1])+length(set3[,1]),ncol=2);

for (i in 1:length(set1[,1])){
	dataset[[i,1]] = set1[[1]][i];
	dataset[[i,2]] = set1[[2]][i];
}

for (i in (length(set1[,1])+1):(length(set1[,1])+length(set2[,1])) ){
	dataset[[i,1]] = set2[[1]][(i-length(set1[,1]))];
	dataset[[i,2]] = set2[[2]][(i-length(set1[,1]))];
}

for (i in ( length(set1[,1])+length(set2[,1])+1):(  length(set1[,1])+length(set2[,1])+length(set3[,1])  )){
	dataset[[i,1]] = set3[[1]][  (i-length(set1[,1]) - length(set2[,1]) ) ];
	dataset[[i,2]] = set3[[2]][  (i-length(set1[,1]) - length(set2[,1]) ) ];
}


#standardizing data mean=0,sd=1
dataset=standardize(dataset);



#counting time for algorithm
ptm <- proc.time()




#kind of initializing mean of each class
mean_Ar = matrix(nrow=k,ncol=2);

#kind of initialization
mean_Ar[[1,1]] = 2.5
mean_Ar[[1,2]] = 2.5

mean_Ar[[2,1]] = -2.5
mean_Ar[[2,2]] = -2.5

mean_Ar[[3,1]] = -2.5
mean_Ar[[3,2]] = 2.5


#rnk matrix show where each data point belongs.
#rnk_matrix[1]=k--> k class
rnk_matrix = matrix(nrow=length(dataset[,1]),ncol=1);



par(mfrow=c(2,2))
plot(dataset[,1],dataset[,2],col="green",xlim=c(-3.2,3.2),ylim=c(-3.2,3.2),pch=20,panel.first = grid(),xlab="x",ylab="y")
points(mean_Ar[,1],mean_Ar[,2],col="blue",pch=4)




#-------------------------------K-means-----------------------------
J=abs(sum(dataset));
previous_J=J+10;


#loop with convergence
while(abs(J-previous_J)>1){
previous_J=J;

#calculating where each data point belongs, rnk criterion

for (i in 1:length(dataset[,1])){
	
	rnk_matrix[i] = 1;	 	
	min = dist(rbind(dataset[i,], mean_Ar[1,]))	

	for (j in 2:k){
		if(  dist(rbind(dataset[i,], mean_Ar[j,])) <min  ){
			min = dist(rbind(dataset[i,], mean_Ar[j,]))				
			rnk_matrix[i] = j;
		}	
	}
}


#calculating new mk
for (i in 1:length(dataset[,1])){
	
	rnk_matrix[i] = 1;	 	
	min = dist(rbind(dataset[i,], mean_Ar[1,]))	

	for (j in 2:k){
		if(  dist(rbind(dataset[i,], mean_Ar[j,])) <min  ){
			min = dist(rbind(dataset[i,], mean_Ar[j,]))				
			rnk_matrix[i] = j;
		}	
	}
}


#for k class
for(t in 1:k){

	#2 for calculating the mk of the k class 
	for (i in 1:length(dataset[,1])){
	
		#select from all the data only those belonging to the k class
		a= which(rnk_matrix==t)	

		sum1=matrix(c(0,0),nrow=2,ncol=1);
		for (j in a){
			sum1 = sum1 + dataset[j,]
		}
			mean_Ar[t,] = sum1/length(a)
	}

}

size_class = c(0,0,0)
#---------------------plot first class-----------------------------

a= which(rnk_matrix==1)	
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}
size_class[1] =counter1;

plot(tr1[,1],tr1[,2],col="red",xlim=c(-3.2,3.2),ylim=c(-3.2,3.2),panel.first = grid(),xlab="x",ylab="y",main=toString(J))

#---------------------plot second class-----------------------------
a= which(rnk_matrix==2)	
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}
size_class[2] =counter1;
points(tr1[,1],tr1[,2],col="green")



#---------------------plot third class-----------------------------
a= which(rnk_matrix==3)	
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}
size_class[3] =counter1;
points(tr1[,1],tr1[,2],col="dodgerblue")
points(mean_Ar[,1],mean_Ar[,2],col="black",pch=4)




J=0;
#J formula
for(g in 1:length(dataset[,1])){	
		  	a=rnk_matrix[i] 
			sum2 = dist(rbind(dataset[g,], mean_Ar[a,]))				
			J=J+sum2	

}

print(J)

}
#--------------------------------------------while loop----------------------------








#--------------------------------------k-means----package-----------------------
#algorithm k-means-->Hartigan-Wong

ptm=proc.time() - ptm
#needs library stats
library("stats")

ptm2 <- proc.time()
rop=kmeans(dataset,3)
ptm2=proc.time() - ptm2


#---------------------plot first class-----------------------------

a= which(rop[1][[1]]==1)
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}


plot(tr1[,1],tr1[,2],col="red",xlim=c(-3.2,3.2),ylim=c(-3.2,3.2),panel.first = grid(),xlab="x",ylab="y",main="using package k-means")

#---------------------plot second class-----------------------------
a= which(rop[1][[1]]==2)	
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}

points(tr1[,1],tr1[,2],col="green")



#---------------------plot third class-----------------------------
a= which(rop[1][[1]]==3)	
tr1=matrix(nrow=length(a),ncol=2)
counter1=1;
for (i in a){
tr1[counter1,1]=dataset[i,1]
tr1[counter1,2]=dataset[i,2]
counter1= counter1+1
}

points(tr1[,1],tr1[,2],col="dodgerblue")
points(rop[[2]][,1],rop[[2]][,2],col="black",pch=4)


#size of each class
t(c(length(which(rnk_matrix==1)),length(which(rnk_matrix==2)),length(which(rnk_matrix==3))))
rop[7]

#centers
mean_Ar
rop[2]

#iterations of kmeans-package
rop[8]


#time
#user kmeans
ptm
#package time
ptm2


dev.off()



