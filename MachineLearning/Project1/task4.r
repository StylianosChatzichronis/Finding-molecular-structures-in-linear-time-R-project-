

#saving plots
pdf('task4.pdf')


par(mfrow=c(3,2),xpd=FALSE)
N<-1;
set.seed(32)



#loop for the cases of 2,3,4,5,9
for (N in c(2,3,4,5,9)){
#Creating a sinus in R
q <- 1:10000
q <- q*0.0001
y <- sin(2*pi*q)

plot(q,y,type="l",col="firebrick1",xlab="q",xlim=c(-0.1, 1.1), ylim=c(-6.5, 6.5) )


#sampling the points

#libraries used
library(RMThreshold) #necessary functions
library(hydroGOF);
library(Metrics v0.1.2);

x2 = vector(,10);
l<-length(q)
step=l/10
for (i in 1:10){
 	x2[i] <- q[i*step];
}

#creating the noisy signal
q <- x2;
y <- sin(2*pi*q);
noise<-rnorm(length(q),0,1)
noisy.y <- y +noise
points(q,noisy.y,col='deepskyblue4',lwd=1)

#computing the model with lm function for N=2,3,4,5,9
model <- lm(noisy.y ~ poly(q,N))

#printing the coefficient values from regression model
a2Pval <- summary(model)$coefficients[2, 4]
print ("w is the following")
print (a2Pval)

#finding the residual standard error
summary(model)
print("Residual standard error")
print (summary(model))

 
predicted.intervals <- predict(model,data.frame(x=q),interval='confidence', level=0.95)
#finding RMSE ,actual and predicted
error<-RMSE(y,predicted.intervals[1,])
print(error)

#displaying the fitting curves
predicted.intervals <- predict(model,data.frame(x=q),interval='confidence', level=0.95)
lines(q,predicted.intervals[,1],col='green',lwd=2)
lines(q,predicted.intervals[,2],col='black',lwd=1)
lines(q,predicted.intervals[,3],col='black',lwd=1)



}

dev.off()



