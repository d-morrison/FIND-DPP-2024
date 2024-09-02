

library(tidyverse)

load('LCMdata.RData')

n=10

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

test.matrix = matrix((NA),nrow = 2^n, ncol=n)
for (i in 0:(2^n -1)){
  test.matrix[i+1,] = number2binary(i,n)
}  

## need to generate the result vector
trim.data = data.merged %>%
  dplyr::select(case,SD.IgM,Typhidot.IgM,Widal,Enterocheck,TestIt,CTK.IgG,Tubex,Spectrum.IgM,DPP.pos.90sens)

result.vector= rep(0,dim(test.matrix)[1])
for (i in 1:dim(test.matrix)[1]){
  for (j in 1:dim(trim.data)[1]){
    if (sum(as.numeric(as.numeric(trim.data[j,])== as.numeric(test.matrix[i,])))==dim(test.matrix)[2]){
      result.vector[i] = result.vector[i] + 1
    }
  }
}


N= sum(result.vector)

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
estBetaParams(0.6,0.05**2)

S.alpha= c(57,1,1,1,1,1,1,1,1,1)
S.beta= c(38,1,1,1,1,1,1,1,1,1)
C.alpha=c(99999,1,1,1,1,1,1,1,1,1)
C.beta=c(1,1,1,1,1,1,1,1,1,1)
pi.alpha=1
pi.beta=1


S.vec = rep(NA,n)
C.vec= rep(NA,n)
Y.vec = rep(NA,2^n)
for (i in 1:n){
  S.vec[i]=rbeta(1,S.alpha[i],S.beta[i])
  C.vec[i]=rbeta(1,C.alpha[i],C.beta[i])
}
pi=0.5

sims= 100000

pi.vec= rep(NA,sims)
S.mat= matrix((NA),nrow=n,ncol=sims)
C.mat= matrix((NA),nrow=n,ncol=sims)


for (i in 1:sims){ #loop over simulations
  if (i %% 1000 == 0){print(i)}
  for (k in 1:(2^n)){ #loop over possible test result combinations
    tp = pi
    fn= 1-pi
    for (m in 1:n){ #loop over tests
      tp = tp*(S.vec[m]**test.matrix[k,m])*((1-S.vec[m])**(1-test.matrix[k,m]))
      fn= fn*((1-C.vec[m])**(test.matrix[k,m]))*(C.vec[m]**(1-test.matrix[k,m]))
    }
    Y.vec[k] = rbinom(1,result.vector[k],tp/(tp+fn))
  }
  pi= rbeta(1,sum(Y.vec)+pi.alpha,N-sum(Y.vec)+pi.beta)
  for (j in 1:n){
    S.vec[j]= rbeta(1,sum(test.matrix[,j]*Y.vec) + S.alpha[j],sum((1-test.matrix[,j])*Y.vec) + S.beta[j]) #multiple vector of test results by the test matrix
    C.vec[j]= rbeta(1,sum((1-test.matrix[,j])*result.vector) - sum((1-test.matrix[,j])*Y.vec) + C.alpha[j], 
                    sum(test.matrix[,j]*result.vector) - sum(test.matrix[,j]*Y.vec) + C.beta[j]) #first number is all negatives minus false negatives; second number is all positives minus true positives 
  }
  pi.vec[i]=pi
  S.mat[,i]=S.vec
  C.mat[,i]=C.vec
  }


#blood.culture,SD.IgM,Typhidot.IgM,Widal,Enterocheck,TestIt,CTK.IgG,Tubex,Spectrum.IgM
subsample = seq(50001,100000,100)

df = data.frame(
  pi=pi.vec,
  S.bloodculture= S.mat[1,],
  S.SD.IgM= S.mat[2,],
  S.Typhidot.IgM = S.mat[3,],
  S.Widal = S.mat[4,],
  S.Enterocheck = S.mat[5,],
  S.TestIt = S.mat[6,],
  S.CTK.IgG = S.mat[7,],
  S.Tubex = S.mat[8,],
  S.Spectrum.IgM = S.mat[9,],
  S.DPP_comb = S.mat[10,],
  C.bloodculture= C.mat[1,],
  C.SD.IgM= C.mat[2,],
  C.Typhidot.IgM = C.mat[3,],
  C.Widal = C.mat[4,],
  C.Enterocheck = C.mat[5,],
  C.TestIt = C.mat[6,],
  C.CTK.IgG = C.mat[7,],
  C.Tubex = C.mat[8,],
  C.Spectrum.IgM = C.mat[9,],
  C.DPP_comb = C.mat[10,],
  sim= 1:sims
)


df.long = pivot_longer(df,cols=1:21)

results2=df.long %>%
  group_by(name) %>%
  summarise(med= median(value),q2.5= quantile(value,0.025),q97.5=quantile(value,0.975))

results2$param = c(rep("Specificity",10),rep("Sensitivity",10),"Prevalence")
results2$name = c(rep(c("CTK IgG","DPP","Enterocheck","SD IgM","Spectrum IgM","TestIt","Tubex","Typhidot IgM",
                      "Widal","Blood culture"),2),"Prevalence")


