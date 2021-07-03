
###################### MIXTURE ###############################
### size of groups
##############################################################
update_ns = function(s)
{
  ns = rep(0,k)
  for(i in 1:k) ns[i] = sum(s==i)
  return(ns)
}


###############################
# Updating indicator variable ## 
###############################

update_s <- function(y,X,Z,b,thetamix)
{
  s.new <- rep(0,m)
  prob=numeric()
  for(i in 1:m){
    idd=which(ind==i)
    for(j in 1:k){
      et=X[idd,]%*%thetamix[j,1:p]+as.numeric(Z[idd,]%*%b) 
      med=exp(et)/(1+exp(et))
      sigm=exp(thetamix[j,np])
      prob[j]=exp(log(thetamix[j,h])+sum(dllogistic( y[ind==i],med,sigm,log=TRUE)) ) 
    }
    
    probs=prob/sum(prob)
    s.new[i]=sample(1:k,1,prob=probs)
    if(max(s.new)!=k){id=sample(1:m,1,prob=rep(1/m,m)); s.new[id]=k}
  }
  return(s.new)
}




######################################################################
##########################  A-Likelihood function                  
######################################################################


Loglink = function(y,X,Z,b,bta){
  eta=X%*%bta[1:p]+Z%*%b 
  med=c(exp(eta)/(1+exp(eta)))
  phi=exp(bta[np])
  som=sum(dllogistic(y,med,phi, log= TRUE))
  return(som)
}


###################################################################                                 
############# POSTERIOR for coeficients ############################
###################################################################

post_betas <- function(y,X,Z,b,bta){  
  Loglink(y,X,Z,b,bta)+ sum(dnorm(bta[1:np],0,100,log=TRUE)) 
}

######################################################################         
#######METROPOLIS hastings for coeficients ###########
###################################################################### 

update_betas=function(y,X,Z,b,bta){
  
  for(l in 1:np){
    btan=bta
    btan[l]=rnorm(1, mean = bta[l],Wb)
    fnew=  post_betas(y,X,Z,b,btan)
    fold=  post_betas(y,X,Z,b,bta)
    wp=min(1,exp(fnew - fold));
    u=runif(1)   
    if(u< wp){bta[l]=btan[l]}
  } 
  
  
  return(bta)
}



###############################################################
########## Posterior for cov-matrix of b #######################
###############################################################         

post_Db <- function(b){   
  if(npbi==1){
    v=1
    V=1
    MM=rgamma(1, (m+v)/2,rate=(sum(b^2)+V)/2)
    post=sqrt(1/MM)
    return(post)  }else{ 
      bl=t(matrix(b,npbi,m))
      del=200
      grl=m+del
      Sb=solve(diag(npbi)*rep(100,npbi)+ t(bl)%*%(bl))
      MM=   rWishart( 1,grl, Sb) 
      MM= matrix(MM,npbi,npbi)
      post=(MM)
      return(post) 
    }
}



######################################################
##### Posterior for random-effects bi's ###############
######################################################


post_bi <- function(yi,Xi,Zi,bta,bi,Db){
  eta=Xi%*%bta[1:p]+Zi%*%bi 
  med=c(exp(eta)/(1+exp(eta)))
  phi=exp(bta[np])
  L=sum(dllogistic(yi,med,phi, log= TRUE))
  mu=rep(0,npbi) 
  if(length(mu)==1){
    Pbi=sum(dnorm(bi, mu, Db, log=TRUE) )
  }else{Pbi=dmvnorm(bi, mu, Db, log=TRUE) } 
  post=L+Pbi 
  return(post)    
}


######################################################
#######   Metropolis-Hastings for updating bi's #######
######################################################                      




update_bis=function(y,X,Zi,bi,theta,Db){
  mi=bi
  Sigmab=diag(Vbi,length(mi),length(mi))
  mi=mvrnorm(1, bi, Sigmab) 
  fnew= post_bi(y,X,Zi,theta,mi,Db)
  fold=post_bi(y,X,Zi,theta,bi,Db)
  wp=min(1,exp(fnew - fold))
  u=runif(1)
  if(u< wp){bi=mi}
  return(bi)
}






##============================================================
############### Declarations for algorithm
#####################################################
ind=c(t(matrix(rep(seq(1:m),Tp),m,Tp)) )
p=length(X[1,])   # number of the coefficients in each component
np=p+1            # total of the coefficients in each component
h=np+1
###=============================================================


#############################################
########## fixed hyper-parameters
#############################################
nu=rep(1,k)

#############################################
########## Initial values
#############################################
mu_zeros=rep(0,npbi) 
Db=diag(length(mu_zeros))
set.seed(seed)
b=mvrnorm(m, mu_zeros, Db)
b=c(t(b))

set.seed(seed)
thetamix=cbind(matrix(rep(runif(1),np*k),k,np),rep(1/k,k)) 

S_init=numeric() 
for (i in 1:m){
  S_init[i] = rDiscreta(rep(1/k,k))
}
s = S_init

ns = update_ns(s)

#############################################
########## Declarations ######################
#############################################


Samp=list(beta=NULL,delta=NULL,omega=NULL,B=NULL,S=NULL,D=NULL,Loglink=NULL )


Vbi=0.3 #variance of proposal function for updating random-effects
Wb=0.2  #variance of proposal function for updating fixed-effects

nit=300000
r = 0



#############################################
########## Algorithm
#############################################


while(r <= nit){  
  r = r + 1
  s=update_s(y,X,Z, b, thetamix)
  ns=update_ns(s)
  thetamix[,h]=rdirichlet(1,nu+ ns)
  
  s_l=c(t(matrix(rep(s,Tp),m,Tp)) )
  
  for(j in 1:k){
    id=which(s_l==j)
    bta= thetamix[j,1:np]
    yj=y[id]
    Xj=X[id,]
    Vj=Z[id,]
    thetamix[j,1:np]=update_betas(yj,Xj,Vj, b,bta)
  }
  
  Db=post_Db(b)
  b_w=t(matrix(b,npbi,m))
  
  for(i in 1:m) {
    bi=b_w[i,]
    ym=y[ind==i]
    b_w[i,]=update_bis(ym,Xi,Zi,bi,  thetamix[ s[i],1:np], Db)
  }
  b=c(t(b_w))
  
  
  if(r==100){Wb=0.01;Vbi=0.01}     
  print(  thetamix)
  print( Db)
  
  L = 0
  for(i in 1:m){
    dk = 0
    for(j in 1:k){ 
      bta=thetamix[j,] 
      me=Xi%*%c(bta[1:p]) + Zi%*%c(b_w[i,]) 
      med_r=exp(me)/(1+exp(me))
      phi_r=exp(bta[np])
      dencomp=log(bta[h]) + sum(dllogistic(y_w[i,],med_r,phi_r,  log= TRUE))
      dk=dk + dencomp
    }
    L = L + dk
  }
  
  
  
  Samp[["beta"]] = rbind(Samp[["beta"]],c(thetamix[,1:pcov]))
  Samp[["delta"]] = rbind(Samp[["delta"]],c(thetamix[,np]))
  Samp[["omega"]] = rbind(Samp[["omega"]],c(thetamix[,h]))
  Samp[["B"]] = rbind(Samp[["B"]],c(b_w))
  Samp[["S"]] = rbind(Samp[["S"]],c(s))
  Samp[["D"]] = rbind(Samp[["D"]],c(Db))
  Samp[["Loglink"]] = c(Samp[["Loglink"]],c(L))
  
  
  flush.console()}











