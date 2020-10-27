  
     
                ###################### MIXTURE ###############################
### size of groups
update.nz = function(z)
{
  nz = rep(0,k)
  for(i in 1:k) nz[i] = sum(z==i)
  return(nz)
}

# Updating indicator variable

update.z <- function(y,X,Z,b,thetamix)
{
  z.new <- rep(0,m)
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
    z.new[i]=sample(1:k,1,prob=probs)
    if(max(z.new)!=k){id=sample(1:m,1,prob=rep(1/m,m)); z.new[id]=k}
  }
    return(z.new)
}



#####################################################################
######################################################################
##########################   Likelihood fuction                  
######################################################################
######################################################################
                  
            
   Loglink = function(y,X,Z,b,bta){
      eta=X%*%bta[1:p]+Z%*%b 
      med=c(exp(eta)/(1+exp(eta)))
      phi=exp(bta[np])
      som=sum(dllogistic(y,med,phi, log= TRUE))
   return(som)
   }

                 
                                 
############# POSTERIOR for coeficients ############################
 
   post.betas <- function(y,X,Z,b,bta){  
   Loglink(y,X,Z,b,bta)+ sum(dnorm(bta[1:np],0,100,log=TRUE)) 
   }
         
#######METROPOLIS hastings for coeficients ###########
 
   update.betas=function(y,X,Z,b,bta){
       
       for(l in 1:np){
       btan=bta
       btan[l]=rnorm(1, mean = bta[l],Wb)
       fnew=  post.betas(y,X,Z,b,btan)
       fold=  post.betas(y,X,Z,b,bta)
       wp=min(1,exp(fnew - fold));
       u=runif(1)   
       if(u< wp){bta[l]=btan[l]}
        } 
       
                      
   return(bta)
   }
  ######################################################################
######################################################################
############# POSTERIOR Db ############################
         
           
 post.Db <- function(b){   
##### hyperparameters
 if(npbi==1){
 v=1
 V=1
 MM=rgamma(1, (m+v)/2,rate=(sum(b^2)+V)/2)
 post=sqrt(1/MM)
 return(post)  }else{ 
 b=t(matrix(b,npbi,m))
 del=3
 grl=m+del
 Sb=solve(diag(npbi)*rep(0.0001,npbi)+ t(b)%*%(b))
 MM=   rWishart( 1,grl, Sb) 
 MM= matrix(MM,npbi,npbi)
 post=(MM)
 return(post) 
 }
 }
       
   
     
       
       
######################################################
#####
######################################################




post.bi <- function(yi,Xi,Zi,bta,bi,Db){
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
    
    

    
                             
   #######METROPOLIS B'IS' ###########
                      
                            
                     
                            

update.bis=function(y,X,Zi,bi,theta,Db){
            mi=bi
            Sigmab=diag(Vbi,length(mi),length(mi))
            mi=mvrnorm(1, bi, Sigmab) 
            fnew= post.bi(y,X,Zi,theta,mi,Db)
            fold=post.bi(y,X,Zi,theta,bi,Db)
            wp=min(1,exp(fnew - fold))
            u=runif(1)
            if(u< wp){bi=mi}
             return(bi)
}





    
##============================================================
############### Components
#####################################################
ind=c(t(matrix(rep(seq(1:m),Tp),m,Tp)) )
p=length(X[1,])     # numero de coeficientes fixos por componente
np=p+1            # total de coeficientes por componente
## parametros
npbi=length(Z1[1,])
h=np+1
###=============================================================
########## huperparameter fixed
 delta=1
 
#############################################
########## initials values
mu=rep(0,npbi) 
Sigmab=diag(length(mu))
b=mvrnorm(m, mu, Sigmab)
b=c(t(b))



thetamix=cbind(matrix(rep(runif(1),np*k),k,np),rep(1/k,k)) 
Sv=numeric() 
for (i in 1:m){
	Sv[i]<-rDiscreta(rep(1/k,k))
}
z=Sv

nz=update.nz(z)
Db=Sigmab
v=list(z=z,thetamix=thetamix,b=b,nz=nz,Db=Db)
 

Vbi=0.3# variance of proposal function
Wb=0.2 # variance of proposal function
nit=20000
 

cat( file="theta", append=FALSE);
cat( file="B", append=FALSE);
cat( file="S", append=FALSE);
cat( file="Db", append=FALSE);
for(r in 1:nit){ for(l in 1:10){  
  v$z=update.z(y,X,Z,v$b,v$thetamix)
  v$nz=update.nz(v$z)
  v$thetamix[,h]=rdirichlet(1,delta+v$nz)
  
  zm=c(t(matrix(rep(v$z,Tp),m,Tp)) )
  
  for(j in 1:k){
    id=which(zm==j)
    bta=v$thetamix[j,1:np]
    yj=y[id]
    Xj=X[id,]
    Vj=Z[id,]
    v$thetamix[j,1:np]=update.betas(yj,Xj,Vj,v$b,bta)
  }
  
  v$Db=post.Db(v$b)
  bm=t(matrix(v$b,npbi,m))
  
  for(i in 1:m) {
    bi=bm[i,]
    ym=y[ind==i]
    bm[i,]=update.bis(ym,Xi,Zi,bi, v$thetamix[v$z[i],1:np],v$Db)
  }
  v$b=c(t(bm))
}
  
  if(r==100){Wb=0.01;Vbi=0.01}     
  print( v$thetamix)
  print(v$Db)
  cat(c(t( v$thetamix)),"\n",file="theta", append=TRUE)
  cat(c(v$z),"\n",file="S", append=TRUE)
  cat(c(v$Db),"\n",file="Db", append=TRUE) 
  cat(c(v$b),"\n",file="B", append=TRUE)                                      
  flush.console()}



   







