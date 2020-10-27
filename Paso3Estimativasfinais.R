btamix=rbind( c(1.20   ,3.5,  3.8,  2.50),
              c( -0.5, -5.5, -5.5, 2.2),
          


theta=read.table("theta") 
plot(theta[,1], type="l")
lines(theta[,7])
lines(theta[,13])


DB=read.table("Db") 
plot(DB[,1], type="l")
lines(theta[,7])
lines(theta[,13])

#Elaboration of summary

TT=length(theta[,1])
jump=seq(1,TT,20)
R=DB[jump,]
R=R[-100,]
Est=numeric()
for(i in 1:9){
vetor=pp<-matrix(R[,i])
vetor <- mcmc(vetor)
hpd=HPDinterval(vetor,prob=0.95)
int=hpd[1:2]
est=c(mean(R[,i]),int)
Est=rbind(Est,est)
}






  
Sm=read.table("S") 
TT=length(Sm[,1]) 
jump=TT/2
T=length(Sm[ddi,1])
ddi=jump:TT

IID=numeric()
for(t in 1:m){
  A=Sm[ddi,t]
  IID[t]=as.numeric(names(table(A))[which.max(table(A))])
}


  Year=c(1994,1998,2002,2006,2010)
par(mfrow=c(3,1))
y_sim=t(matrix(y,Tp,m))
plot(Year,y_sim[1,],col="white", type="l",ylim=c(0,1),main="K 3 Grau 3",xaxt='n')
axis(1, at=Year,labels=ano, col.axis="black", las=1)

for(i in 1:m){
if(IID[i]==1)
lines(Year,y_sim[i,],col=IID[i],type="o")
}
 

plot(Year,y_sim[1,],col="white", type="l",ylim=c(0,1),main="K 3 Grau 3",xaxt='n')
axis(1, at=Year,labels=ano, col.axis="black", las=1)

for(i in 1:m){
if(IID[i]==2)
lines(Year,y_sim[i,],col=IID[i],type="o")
}



plot(Year,y_sim[1,],col="white", type="l",ylim=c(0,1),main="K 3 Grau 3",xaxt='n')
axis(1, at=Year,labels=ano, col.axis="black", las=1)

for(i in 1:m){
if(IID[i]==3)
lines(Year,y_sim[i,],col=IID[i],type="o")
}























