library(readxl)
census1=as.data.frame(read_xlxs("GreekCensus2011Data.xlsx",sheet=1))
census2=as.data.frame(read_xlxs("GreekCensus2011Data.xlsx",sheet=2))
census3=as.data.frame(read_xlxs("GreekCensus2011Data.xlsx",sheet=3))

#answer1 is objective 1
answer1=sum( rowSums ( census1[  2:(dim(census1))[1],2:(dim(census1))[2]  ]) ) #sum is10817485

#answer2 is objective 2
answer2= col ])Sums ( census1[  2:(dim(census1))[1],2:(dim(census1))[2]  ]) 

#barplot is objective 3
barplot(answer2,main="Age Distribution"),xlab="Age",ylab="Absolute Values")

arrayBarplot = answer2/answer1
arrayBarplot= arrayBarplot*100
barplot(arrayBarplot,main="Age Distribution"),xlab="Age",ylab="Percentage")


#creating arrays with FALSE,TRUE
#arrayTrue είναι στην "Ευρώπη"==TRUE,
#arrayTrue3 είναι στην "Ναι"==TRUE,
arrayTrue = census2[,2]
arrayTrue3 = census3[,2]
arrayTrue = arrayTrue == "Ευρώπη"
arrayTrue3 = arrayTrue3 == "Ναι"


#if(arrayTrue[i] == TRUE){
#επιλέγουμε από την πρώτη καρτέλα τις χώρες που είναι στην Ευρώπη
#if(arrayTrue3[i] == TRUE){
#επιλέγουμε από την πρώτη καρτέλα τις χώρες που είναι στην ΕΕ και παίρνουμε
#παίρνουμε άθροισμα των επιθυμητών γραμμών
#ομοίως στο else για τους εκτός ΕΕ

sumEu = 0
sumNonEu=0;
j=1;
for(i in 1:length(arrayTrue)){
	if(arrayTrue[i] == TRUE){
		if(arrayTrue3[j] == TRUE){
			sumEu = sumEu + rowSums(census1[1,2:dim(census1)[2]])
		}
		else
		sumNonEu = sumNonEu +rowSums(census1[i,2:dim(census1)[2]])
		}
	j=j+1
	}
}



#sumEu 10058669
#sumNonEu 495658


slices = c(sumNonEu,sumEeu)
lbl = c("outside European Union","In European Union")
pie(slices,labels=lbl,main="Europeans in the country")



arrayNoEurope = census2[,2]

sumAfrica =0;
sumAmerica =0;
sumAsia =0;
sumOther =0;

#επιλέγουμε για κάθε ήπειρο και αθροίζουμε στο αντίστοιχο άθροισμα

for(i in 1:length(arrayNoEurope)){

if(census2[i,2] == "Αφρική"){
sumAfrica =sumfrica+ rowSums(census1[1,2:dim(census1)[2]]);
}
else{
if(census2[i,2] == "Αμερική" ){
sumAmerica = sumAmerica + rowSums(census1[1,2:dim(census1)[2]]);
}
else{
if((census2[i,2] == "Ασία"){
sumAsia = sumAsia  + rowSums(census1[1,2:dim(census1)[2]]);
}else{
 sumOther = sumOther + rowSums(census1[1,2:dim(census1)[2]]);
}
}
}
}

lbl2 = c("Africa","America","Asia","Other/Unknown")
slices2 = c(sumAfrica,sumAmerica,sumAsia,sumOther)
pie(slices2,labels=lbl2,main="Non Europeans living in the Country")