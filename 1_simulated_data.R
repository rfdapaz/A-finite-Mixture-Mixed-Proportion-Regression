
           
if(!require(truncnorm)) install.packages("truncnorm");require(truncnorm) 
if(!require(MASS)) install.packages("MASS");require(MASS) 
if(!require(llogistic)) install.packages("llogistic");require(llogistic) 
if(!require(MCMCpack)) install.packages("MCMCpack");require(MCMCpack) 
if(!require(mvtnorm)) install.packages("mvtnorm");require(mvtnorm)

if(!require(LaplacesDemon)) install.packages("LaplacesDemon");require(LaplacesDemon)
setwd("G:\\Meu Drive\\0 UFC\\0 PROJETOS\\Artigo mixture misto ll\\Programas utilizado na correcao do artigo_2020")


##############################################################################
##########      DADOS SIMULADOS
###################################### ########################################
#
m=100    #n?mero de individuos
Tp=5     #n?mero de replica??es no tempo
n=m*Tp     # total de medi??es
k=2       #n?mero de grupos 
grau=2 #number of fixed effects
pcov=grau+1
seed=2005    
w=c(0.4,0.6) 
# 

################## VARIAVEL EXPLICATIVA DO EFEITO FIXO #######
#com polinomio ortogonal na matriz do modelo
Vt=rep(1:Tp,m)
X=model.matrix(~poly(Vt, degree=grau))
#X=model.matrix(~Vt)
cc=ncol(X)
X=matrix(X,n,cc)
#
################## VARIAVEL EXPLICATIVA DO EFEITO ALEAT?RIO #######
#
#####forma de blocos para modelo misto      (nxm)
Z1=X[1:Tp,]
Z=diag(m)%x%Z1 
Zi=Xi=Z1

#
################# DEFINE O GRUPO E CALCULA O PREDITOR LINEAR ##########
#
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 val}
#
set.seed(seed)
pgrupos<-w
Sverd<-numeric(m)
indiv<-NULL
Y<-NULL
for (i in 1:m){
	Sverd[i]<-rDiscreta(pgrupos)
}
#

       

##################  VARIAVEL EXPLICATIVA DO EFEITO FIXO #######
t=rep(seq(1:Tp),m)
           

set.seed(seed)
btamix=rbind( c(1.20,3.5, 3.8,  2.50),
              c( -0.5, -6.5, -5.5, 2.2),
              c( -2.50, 2.0, 6.5,  3.2))  
              
              
 betamix=cbind(btamix[,1:cc],btamix[,length(btamix[1,])]         )  
 thetamix=cbind(btamix[1:k,],w)         
 s=c(t(matrix(rep(Sverd,Tp),m,Tp)) )        
 betamix=betamix[s,]
 ind=c(t(matrix(rep(seq(1:m),Tp),m,Tp)) )
 
 ##################  EFEITOS ALEAT?RIOS b #################
npbi=ncol(Z1)
mu=rep(0,npbi) 
set.seed(seed)
Sigmab=solve(diag(length(mu))+2*runif(npbi))
b=mvrnorm(m, mu, Sigmab)
b=c(t(b))

   
#########################################################
################ generated mixture  #####################
#########################################################


y_sim=numeric()
for(i in 1:n){
eta=X[i,]%*%betamix[i,1:pcov]+Z[i,]%*%as.matrix(b)    ### preditor linear do modelo misto para o individuo i
med=exp(eta)/(1+exp(eta))
ph=exp(betamix[i,(pcov+1)])
y=rllogistic(1 , med, ph)     #### amsotra gerada sem mistura
y_sim=rbind(y_sim,y)
}   

par(mfrow=c(1,2))
y=y_sim
y_sim=t(matrix(y,Tp,m))
 ano=c(1,2,3,4,5)
y_sim=t(matrix(y,Tp,m))
plot(ano,y_sim[1,], type="l",ylim=c(0,1),main="Real Classification", ylab=" Proportion",
 xlab="Point at time",xaxt='n')
for(i in 2:m){
lines(ano,y_sim[i,],col=Sverd[i], type="o",  lty=Sverd[i])

}
axis(1, at=ano,labels=ano, col.axis="black", las=1)

 



#########################################################################3
#################################################3
#####################################################################

                  [,1]       [,2]       [,3]
[1,]  0.8437310 -0.1562690 -0.1562690
[2,] -0.2776464  0.7223536 -0.2776464
[3,] -0.2149233 -0.2149233  0.7850767



   

####### Gerando valores correlacionados 
### a partir de uma matriz de variancias 
## e covariancias


##Matriz de Cholesky

J=ncol(Sigmab)

Db=Sigmab      ### matriz de variancias e covariancias
D=diag(diag(Db))
C=(Db-D)+diag(rep(1,J)) ##corelation matrix
tau=diag(sqrt(diag(Db))) ## standart deviation matrix

L=chol(C)  ### matriz de cholesky para correlação
A=(tau)%*%(L)

z=rnorm(J)

x=(A)%*%(z)   ###valores correlacionados


