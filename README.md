## A finite Mixture Mixed Proportion Regression Model for Classification Problems in Longitudinal Voting data


The aim of this work is to present the scripts for implement the methodology presented in the manuscript: A finite Mixture Mixed Proportion Regression Model for Classification Problems in Longitudinal Voting data.
This work is the result of a partnership between the Federal University of Ceará ( Laboratory of Innovative Technologies from portuguese, LTI, <https://dmontier.pro.br/lti/>), Instituto de Ciências Matemáticas e de Computação, Universidade de São Paulo, both in Brazil, and University of Connecticut, Storrs, USA.


Autors: Rosineide da Paz, Jorge Luis Bazán, Victor Hugo Lachos and Dipak Dey


The methodology is implemented in real and simulated data set. Here we present the code for simulate the data used in the manuscript.


### Packages requerid for simulation



```{r include=FALSE}
if(!require(truncnorm)) install.packages("truncnorm");require(truncnorm) 
if(!require(MASS)) install.packages("MASS");require(MASS) 
if(!require(llogistic)) install.packages("llogistic");require(llogistic) 
if(!require(MCMCpack)) install.packages("MCMCpack");require(MCMCpack) 
if(!require(mvtnorm)) install.packages("mvtnorm");require(mvtnorm)
if(!require(LaplacesDemon)) install.packages("LaplacesDemon");require(LaplacesDemon)

```


