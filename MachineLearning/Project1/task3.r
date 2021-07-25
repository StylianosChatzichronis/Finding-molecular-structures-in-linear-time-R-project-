set.seed(102)


#Saving plots
pdf(file = "task3.pdf")

#prior sd is 4 prior
#hyperparameters h_mean=0 h_sd=2
h_mean = 4
h_sd = 2
 

par(mfrow=c(3,3),xpd=FALSE)
N=1;
x=rnorm(n=N,m=7,sd=4)
#finding observed values
	x<-rnorm(n=N,m=7,sd=4) #creation of observed data
	sampleMean=mean(x);	

	x1 <- seq(-30, 30, 0.1)

	#Finding posterior	
	posterior_mean = (N*(h_sd^2)*sampleMean + 16*0 )  /  ( N*(h_sd^2) + 16)
	posterior_var = (h_sd^2*16) / N*( h_sd^2 + 16 )  
	post.norm<-rnorm(n=N,m=posterior_mean,sd=sqrt(posterior_var)) #creation of observed data



	#finding posterior, observed, prior
	plot(x1, dnorm(x1, 0, 2), type = "l",xlim=c(-20, 20), xlab=N,col="blue")
	lines(x1,dnorm(x1, posterior_mean, sqrt(posterior_var)),col="gold")
	lines(x,dnorm(x1, x, sqrt(0)),col="red")

	legend("topleft", 
       inset=.05, 
       cex = 1,  
       c("data","post","prior"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("red","gold","blue"), 
       box.lty=0)



for (N in c(5,10,20,50,100,1000)){
	print(paste("The number of observed values is", N))
	
	#finding observed values
	x<-rnorm(n=N,m=7,sd=4) #creation of observed data
	sampleMean=mean(x);	

	x1 <- seq(-30, 30, 0.1)

	#Finding posterior	
	posterior_mean = (N*(h_sd^2)*sampleMean + 16*0 )  /  ( N*(h_sd^2) + 16)
	posterior_var = (h_sd^2*16) / N*( h_sd^2 + 16 )  
	post.norm<-rnorm(n=N,m=posterior_mean,sd=sqrt(posterior_var)) #creation of observed data



	#finding posterior, observed, prior
	plot(x1, dnorm(x1, 0, 2), type = "l",xlim=c(-20, 20), xlab=N,col="blue",ylim=c(0, 0.35))
	lines(x1,dnorm(x1, posterior_mean, sqrt(posterior_var)),col="darkgoldenrod1")
	lines(density(x),col="red")

	
  legend("topleft", 
       inset=.05, 
       cex = 1,  
       c("data","post","prior"), 
       horiz=FALSE, 
       lty=c(1,1), 
       lwd=c(2,2), 
       col=c("red","gold","blue"), 
       box.lty=0)
}

dev.off()

#Παρατηρήσεις: Η εκ των προτέρων (prior) κατανομή(με εξαιρεση Ν=1) παραμένει σταθερή,
#η κατανομή των παρατηρήσεων παραμένει σχεδόν σταθερή με αβεβαιότητα
#λόγω της μικρής κορυφής και η posterior, όσο αυξάνονται τα δείγματα
# μικραίνει η διακύμανση και αποκτά μεγάλη κορυφή. Λογικό, όπως
#αναφέρεται στη βιβλιογραφία όσο το Ν αυξάνεται η αβεβαιότητα στη posterior
# μειώνεται αφού μειώνεται η επίδραση της σ.

