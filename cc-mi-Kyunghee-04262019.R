# rm(list=ls(all.names=TRUE))

mac <- 0

if (mac==1) {
  # Mac
  indMC <- 5
  m <- 5 # 1, 2, 3, 4, 5, 6
  
} else {
  # Linux
  argFor <- commandArgs(TRUE)
  m <- as.numeric(argFor[1])
  indMC <- as.numeric(argFor[2])
}

library(mitools)
library(survey)
library(gam)

# outer parameter (to be designated)
# N <- 10000; optn <- 'sub'; m <- 1
N <- 10000
gam1 <- -5
xi <- 1.8

Delta <- seq(0,3.5,length.out=6) # log(seq(1,exp(3.5),length.out=6))

# global parameter
n <- c()
M <- 100
B <- 125


logPn <- function(y,x,gam1,gam2,delta,xi) {

  N <- length(x)
  
  tmp <- rep(0,N)
  tmp[which(x>xi)] <- 1
  logitP <- gam1 + gam2*x + delta*(x-xi)*tmp
  # logitQ <- alpha + beta*x
  
  logP <- y*logitP - log(1+exp(logitP))
  # logQ <- y*logitQ - log(1+exp(logitQ))
  
  return(sum(logP[which(logP!=-Inf)]))
}

logQn <- function(y,x,alpha,beta) {
  
  N <- length(x)
  
  # tmp <- rep(0,N)
  # tmp[which(x>xi)] <- 1
  # logitP <- gam1 + gam2*x + delta*(x-xi)*tmp
  logitQ <- alpha + beta*x
  
  # logP <- y*logitP - log(1+exp(logitP))
  logQ <- y*logitQ - log(1+exp(logitQ))
  
  return(sum(logQ[which(logQ!=-Inf)]))
}

mu <- function(x,gam1,gam2,delta,xi) {
  N <- length(x)
  tmp <- rep(0,N)
  
  tmp[which(x>xi)] <- 1
  pr <- exp(gam1+gam2*x+delta*(x-xi)*tmp)/(1+exp(gam1+gam2*x+delta*(x-xi)*tmp))
  
  return(pr)
}


yGen <- function(x,gam1,gam2,delta,xi) {
  N <- length(x)
  tmp <- rep(0,N)
  
  tmp[which(x>xi)] <- 1
  pr <- exp(gam1+gam2*x+delta*(x-xi)*tmp)/(1+exp(gam1+gam2*x+delta*(x-xi)*tmp))
  
  y <- rbinom(rep(1,N),rep(1,N),pr)
  return(y)
}


time1 <- Sys.time()
asdf <- asdf1 <- list()
gam2All <- betaTarget <- alphaTarget <- matrix(nrow=B,ncol=length(Delta))

print(m)

delta <- Delta[m]

lrt <- rep(0,B)
lrStat <- lrStatTrue <- rep(0,B)
betaAll <- interAll <- c()
scoreIpwAtTarget <- c()
ksTest <- c()

for (k in 1:B) {
  
  
  time01 <- Sys.time()
  
  print(c(m,k))
  
  if (mac!=1) {
    write.table(k,paste('/home/khhan/project/multi-impute/simulation/loop-first-04262019-indMC',indMC,'-N',N,'-B',B,'-m',m,'.txt',sep=''))
  }
  
  gam2Grid <- seq(0,1.5,length.out=500)
  
  eps <- 1
  iter <- 0
  
  # target parameter - beta, gam2
  betaTest <- epsTest <- c()
  while (eps > 0.001 && iter < length(gam2Grid)) {
    
    iter <- iter+1
    
    set.seed(B*(indMC-1)+k)
    
    x <- rnorm(N)
    y <- yGen(x,gam1,gam2Grid[iter],delta,xi)
    # x0 <- rep(gam1,length(x))
    
    model0<-glm(y~x,family=binomial)
    betaTest <- c(betaTest,model0$coefficients[2])
    eps <- abs(model0$coefficients[2]-1)
    
    epsTest <- c(epsTest,eps)
    
  }
  
  gam2 <- gam2Grid[which(epsTest==min(epsTest))[1]]
  beta <- betaTest[which(epsTest==min(epsTest))[1]]
  
  gam2All[k,m] <- gam2
  betaTarget[k,m] <- beta
  
  # target parameter - alpha
  alphaTest <- gam1Test <- c()
  alphaTest0 <- c()
  options(warn = -1) 
  for (i in 1:500) {
    
    set.seed(i*5000)
    
    x <- rnorm(N)
    y <- yGen(x,gam1,gam2,delta,xi)
    
    tmp <- rep(0,N)
    tmp[which(x>xi)] <- 1
    z <- (x-xi)*tmp
    df <- data.frame(y=y,x=x,z=z)#,x0=x0)
    
    insample<-c(which(df$y==1),sample(which(df$y==0),sum(df$y)))
    sdf<-df[insample,]
    
    
    model0<-glm(y~x,family=binomial,data=df)
    modelm <- glm(y~x,family=binomial,data=sdf)
    
    modelTmp <- glm(y~x+z,family=binomial,data=sdf)
    
    alphaTest0 <- c(alphaTest0,model0$coefficients[1])
    alphaTest <- c(alphaTest,modelm$coefficients[1])
    
    gam1Test <- c(gam1Test,modelTmp$coefficients[1])
    
  }
  
  alphaFull <- mean(alphaTest0)
  alphaSub <- mean(alphaTest)
  gam1Sub <- mean(gam1Test)
  
  
  
  ### data generation
  
  set.seed(B*(indMC-1)+k)
  
  x <- rnorm(N)
  y <- yGen(x,gam1,gam2,delta,xi)
  # x0 <- rep(gam1,length(x))
  
  tmp <- rep(0,N)
  tmp[which(x>xi)] <- 1
  z <- (x-xi)*tmp
  df <- data.frame(y=y,x=x,z=z)#,x0=x0)
  
  
  
  ### Case-control sampling
  insample<-c(which(df$y==1),sample(which(df$y==0),sum(df$y)))
  n[k] <- length(which(df$y==1))
  sdf<-df[insample,]
  fdf<-df
  fdf$x[-insample]<-NA
  Nmiss<-nrow(df)-length(insample)
  w<-sum(df$y==0)/sum(sdf$y==0)
  sdf$wt<-with(sdf,ifelse(y==1,1,w))
  
  
  # fit
  model0<-glm(y~x,family=binomial,data=df)
  modelm<-glm(y~x,family=binomial,data=sdf)
  models<-glm(y~x,family=binomial,data=sdf,weight=wt)
  
  muTrue <- mu(x,gam1,gam2,delta,xi)
  scoreIpwAtTarget[k] <- mean(((y-muTrue)*x)[insample]*models$data$wt)

  
  # LRT
  modelTrue <- glm(y~x+z,family=binomial,data=sdf)
  modelNearlyTrue <- glm(y~x,family=binomial,data=sdf)
  
  gam1Hat <- modelTrue$coefficients[1]
  gam2Hat <- modelTrue$coefficients[2]
  deltaHat <- modelTrue$coefficients[3]
  
  alphaHat <- modelNearlyTrue$coefficients[1]
  betaHat <- modelNearlyTrue$coefficients[2]
  
  lrStatTrue[k] <- -2*(logQn(y[insample],x[insample],alphaSub,beta) - logPn(y[insample],x[insample],gam1Sub,gam2,delta,xi))
  lrStat[k] <- -2*(logQn(y[insample],x[insample],alphaHat,betaHat) - logPn(y[insample],x[insample],gam1Hat,gam2Hat,deltaHat,xi))
  

  if (m==1) {
    if (lrStat[k] > qchisq(0.95,1)) {
      lrt[k] <- 1
    }
  } else {
    if (deltaHat > 0 & sqrt(lrStat[k]) > qnorm(0.95)) {
      lrt[k] <- 1
    }
  }

  
  
    
    
  
  ### lineartiy test
  linLogitHat <- glm(y~x,family=binomial,data=sdf)
  gamLogitHat <- gam(y~lo(x),family=binomial,data=sdf)
  
  linProb <- linLogitHat$fitted.values
  gamProb <- gamLogitHat$fitted.values
  
  linGof <- mean((sdf$y - linProb)^2/(linProb*(1-linProb)))
  gamGof <- mean((sdf$y - gamProb)^2/(gamProb*(1-gamProb)))
  
  gofHat <- linGof - gamGof
  gofHat
  
  # plot(x,prob,ylim=c(0,1))
  # points(sdf$x,linLogitHat$fitted.values,col=2)
  # points(sdf$x,gamLogitHat$fitted.values,col=3)
  
  gofBoot <- c()
  for(tt in 1:250){
    
    # bootstrap
    
    # indBoot <- sample(1:nrow(sdf),nrow(sdf),replace=TRUE)
    # 
    # # indBoot1 <- sample(which(sdf$y==1),length(which(sdf$y==1)),replace=TRUE)
    # # indBoot0 <- sample(which(sdf$y==0),length(which(sdf$y==0)),replace=TRUE)
    # # indBoot <- c(indBoot1,indBoot0)
    # 
    # xBoot <- x[indBoot]
    # yBoot <- y[indBoot]
    
    # yBoot <- rbinom(length(indBoot),1,prob=linProb[indBoot])
    
    yBoot <- rbinom(nrow(sdf),1,prob=linProb)
    
    linLogitBoot <- glm(yBoot~sdf$x,family=binomial)
    gamLogitBoot <- gam(yBoot~lo(sdf$x),family=binomial)
    
    linProbBoot <- linLogitBoot$fitted.values
    gamProbBoot <- gamLogitBoot$fitted.values
    
    linGofBoot <- mean((yBoot - linProbBoot)^2/(linProbBoot*(1-linProbBoot)))
    gamGofBoot <- mean((yBoot - gamProbBoot)^2/(gamProbBoot*(1-gamProbBoot)))
    
    gofBoot[tt] <- linGofBoot - gamGofBoot
    
  }
  
  # plot(sdf$x,sdf$y)
  # points(sdf$x,linLogitHat$fitted.values,col=2)
  # points(sdf$x,gamLogitHat$fitted.values,col=3)
  # points(sdf$x,linLogitBoot$fitted.values,col=4)
  # points(sdf$x,gamLogitBoot$fitted.values,col=5)
  # 
  # mean(gofBoot>gofHat)
  # hist(gofBoot,breaks=50,main=round(gofHat,3))
  # abline(v=gofHat,col=2)
  
  ksTest[k] <- ifelse(mean(gofBoot>gofHat)<0.05,1,0)
  
  
  
  
    
    
  
  ## Resampling imputation
  p3<-vector("list",M)
  for(i in 1:M){
    p3[[i]]<-fdf
    p3[[i]]$x[-insample]<-with(sdf,sample(x[y==0],Nmiss,replace=TRUE))
  }
  p3list<-imputationList(p3)


  ## Parametric imputation
  p4<-vector("list",M)

  # m0<-with(sdf, mean(x-y*coef(modelm)[2]))
  # s0<-with(sdf, sd(x-y*coef(modelm)[2]))
  
  modelmTmp<-lm(x~y,data=sdf)
  
  m0<-with(sdf, mean(x-y*coef(modelmTmp)[2]))
  s0<-with(sdf, sd(modelmTmp$residuals))
  
  for(i in 1:M){
    p4[[i]]<-fdf
    p4[[i]]$x[-insample]<-rnorm(Nmiss,m0,s0)
  }
  p4list<-imputationList(p4)
  
  
  # fit
  modelmi1<-with(p3list, glm(y~x,family=binomial))
  modelmi2<-with(p4list, glm(y~x,family=binomial))
  
  betaVec<-c(coef(model0)[2],
             coef((modelm))[2],
             coef((models))[2],
             coef(MIcombine(modelmi1))[2],
             coef(MIcombine(modelmi2))[2])
  
 
  
  betaAll <- rbind(betaAll,betaVec)
  
  
  
  time02 <- Sys.time()
  
  if (mac!=1) {
    write.table(time02-time01,paste('/home/khhan/project/multi-impute/simulation/time-first-04262019-indMC',indMC,'-N',N,'-B',B,'-m',m,'.txt',sep=''))
  }
  
}

colnames(betaAll) <- c('mle-full','mle-sub','ipw','mi-boot','mi-para')

Dn <- betaAll[,2]-betaAll[,3]

# betaAll
betaMean <- apply(betaAll,2,'mean')
tmp <- matrix(rep(betaMean,rep(B,length(betaMean))),nrow=B,ncol=length(betaMean))

# MSE and bias-variance decomposition
betaMse <- apply((betaAll-1)^2,2,'mean')
betaBias <- betaMean-1
betaVar <- (apply((betaAll-tmp)^2,2,'mean'))

# dist
tmp1 <- matrix(rep(betaAll[,1],(ncol(betaAll)-1)),nrow=nrow(betaAll),ncol=(ncol(betaAll)-1))
betaDist <- sqrt(apply((betaAll[,-1]-tmp1)^2,2,'mean'))

round(betaMse,3)
round(betaBias,3)
# round(betaVar,3)
round(betaMse,3)-round(betaBias^2,3)

round(mean(n)/N*100,3)

qwer <- list()
qwer[[1]] <- round(sqrt(betaMse),4)
qwer[[2]] <- round(betaBias,4)
qwer[[3]] <- round(sqrt(betaVar),4)
qwer[[4]] <- round(betaDist,4)
qwer[[5]] <- round(mean(n)/N*100,4)
qwer[[6]] <- round(c(mean(gam2All[,m]),delta),4)
qwer[[7]] <- round(c(mean(scoreIpwAtTarget),sd(scoreIpwAtTarget)),4)
qwer[[8]] <- round(cor(Dn,lrStatTrue),4)
qwer[[9]] <- round(c(mean(Dn),sd(Dn)),4)
qwer[[10]] <- round(c(mean(lrStat),sd(lrStat)),4)
qwer[[11]] <- round(c(mean(lrt),sd(lrt)),4)
qwer[[12]] <- round(c(mean(ksTest),sd(ksTest)),4)

names(qwer) <- c('sqrtMse','bias','sqrtVar','dist','caseToCohort',
                 'gam2Delta','avgScoreIpwAtTarget','corrBiaslogLR','Dn','logLR','NPtest','ksTest')

asdf[[m]] <- qwer
  
time2 <- Sys.time()

asdf <- asdf[[m]]
asdf[c(6,1,2,3,4,7,8,11,12)]
# asdf1
time2-time1

# round(1-apply(betaTarget,2,'mean'),4)


if (mac==1) {
  save.image(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/simulation/first-04262019-indMC',indMC,'-N',N,'-B',B,'-m',m,'.RData',sep=''))
} else {
  save.image(paste('/home/khhan/project/multi-impute/simulation/first-04262019-indMC',indMC,'-N',N,'-B',B,'-m',m,'.RData',sep=''))
}













# global parameter
M <- 100
B <- 125
N <- 10000
numB <- 1000/B

Delta <- seq(0,3.5,length.out=6) # log(seq(1,exp(3.5),length.out=6))

betaAllMC <- matrix(nrow=numB*B,ncol=5)
colnames(betaAllMC) <- c('mle-full','mle-sub','ipw','mi-boot','mi-para')

mseAll <- biasAll <- varAll <- matrix(nrow=length(Delta),ncol=5)
powerMC <- ksPowerMC <- c()

for (m in 1:6) {

  DnMC <- c()
  lrStatTrueMC <- c()
  lrtMC <- c()
  ksTestMC <- c()

  for (indMC in 1:numB) {
    print(indMC)

    # print(paste('Delta = ',round(Delta[m],3),sep=''))
    load(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/simulation/first-04262019-indMC',indMC,'-N',N,'-B',B,'-m',m,'.RData',sep=''))
    # print(asdf[c(6,1,2,3,4,7,8,11,12)])

    betaAllMC[(B*(indMC-1)+1):(B*indMC),] <- as.matrix(betaAll)

    DnMC <- c(DnMC,Dn)
    lrStatTrueMC <- c(lrStatTrueMC,lrStatTrue)

    lrtMC <- c(lrtMC,lrt)
    ksTestMC <- c(ksTestMC,ksTest)

    print(round(cor(Dn,lrStatTrue),4))

    # hist(lrStat)

  }

  DnMC <- betaAllMC[,2]-betaAllMC[,3]

  # betaAll
  betaMean <- apply(betaAllMC,2,'mean')
  tmp <- matrix(rep(betaMean,rep(numB*B,length(betaMean))),nrow=numB*B,ncol=length(betaMean))

  # MSE and bias-variance decomposition
  betaMse <- apply((betaAllMC-1)^2,2,'mean')
  betaBias <- betaMean-1
  betaVar <- (apply((betaAllMC-tmp)^2,2,'mean'))

  mseAll[m,] <- sqrt(betaMse)
  biasAll[m,] <- betaBias
  varAll[m,] <- sqrt(betaVar)

  powerMC[m] <- mean(lrtMC)
  ksPowerMC[m] <- mean(ksTestMC)

  qwer <- list()
  qwer[[1]] <- round(sqrt(betaMse),3)
  qwer[[2]] <- round(betaBias,3)
  qwer[[3]] <- round(sqrt(betaVar),3)
  qwer[[4]] <- round(cor(DnMC,lrStatTrueMC),3)
  qwer[[5]] <- round(c(mean(DnMC),sd(DnMC)),3)
  qwer[[6]] <- round(c(mean(lrtMC),sd(lrtMC)),3)
  qwer[[7]] <- round(c(mean(ksTestMC),sd(ksTestMC)),3)

  names(qwer) <- c('sqrtMse','bias','sqrtVar',
                   'corrBiaslogLR','Dn','NPtest','ksTest')

  print(round(c(mean(gam2All[,m]),delta),3))
  print(qwer)

}


mseAll <- mseAll[,c(2,3,5,4)]
biasAll <- biasAll[,c(2,3,5,4)]
varAll <- varAll[,c(2,3,5,4)]
powerMC[1] <- powerMC[1]/2
ksPowerMC <- sort(ksPowerMC)



# par(mar=c(9,6,1,10)+0.1)
# pdf('~/Dropbox/research/Combining multiple imputation with ranking of weights/manuscript/figures/table1.pdf',width=10,height=7)
par(mar=c(9,6,1,10)+0.1)

plot(Delta,mseAll[,1],type='b',lwd=5,pch=1,col=1,
     ylim=c(0.1,0.4),xlim=c(0,Delta[length(Delta)]*1.1),
     axes=FALSE,cex=2,cex.lab=2,
     xlab='',ylab=expression(sqrt('MSE')))
axis(2,round(seq(0.1,0.4,0.05),3),line=0,
     cex.axis=1.75)
axis(1,round(Delta,3),
     cex.axis=1.75)
mtext(expression(delta),1,at=Delta[length(Delta)]*1.15,line=0,
      cex=2)

for (k in 2:ncol(mseAll)) {
  if (k==3 | k==4) {
    points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,lty=4,
           cex=2)
  } else {
    points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,
           cex=2)
  }
}

axis(1,round(Delta,3),labels=round(powerMC,3),line=3,
     cex.axis=1.75)
mtext(' Power (MP)',1,at=Delta[length(Delta)]*1.25,line=3,
      cex=1.75)

axis(1,round(Delta,3),labels=round(ksPowerMC,3),line=6,
     cex.axis=1.75)
mtext('(linearity test)',1,at=Delta[length(Delta)]*1.275,line=6,
      cex=1.75)

abline(v=round(Delta,3),h=round(seq(0.1,0.4,0.05),3),col=8,lty=3)
legend('topright',c('MLE','IPW','MI-P','MI-B'),
       lwd=4,col=1:7,pch=1:7,bty='n',cex=1.5,lty=c(1,1,5,5),
       xpd=TRUE,inset=c(-0.2,0))
# dev.off()














