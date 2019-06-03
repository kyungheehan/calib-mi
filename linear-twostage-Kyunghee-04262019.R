# rm(list=ls(all.names=TRUE))

mac <- 0

if (mac==1) {
  # Mac
  indMC <- 1
  # symm <- 0; m <- 1 # 1, 2, 3, 4, 5, 6

  # source('~/Dropbox/research/Combining multiple imputation with ranking of weights/calib-MI/linear-twostage-Kyunghee-04262019.R')

} else {
  # Linux
  argFor <- commandArgs(TRUE)
  symm <- as.numeric(argFor[1])
  m <- as.numeric(argFor[2])
  indMC <- as.numeric(argFor[3])
}


library(missreg3)
library(survey)
library(mitools)
library(invgamma)
library(MASS)
library(np)

# library(splines)
# library(hdrcde)

B <- 250
M <- 100

N <- 5000
n <- N*0.05



estfun.glm<-function (model)
{
  xmat <- model.matrix(model)
  residuals(model, "working") * model$weights * xmat
}


logPn <- function(y,x,z,gam1,gam2,delta) {
  
  if (symm==0) {
    zstrat <- cut(z, c(-Inf, -2.33,2.33,Inf))
  } else {
    zstrat <- cut(z, c(-Inf, -1.80,1.80,Inf))
  }
  zmidx <- ifelse(as.numeric(zstrat)==2,x,0)
  
  muP <- gam1 + gam2*x + delta*zmidx
  
  logP <- -(y-muP)^2/2 - log(2*pi)/2
  
  return(sum(logP[which(logP!=-Inf)]))
}

logQn <- function(y,x,alpha,beta) {
  
  muQ <- alpha + beta*x
  
  logQ <- -(y-muQ)^2/2 - log(2*pi)/2
  
  return(sum(logQ[which(logQ!=-Inf)]))
}


yGen <- function(x,z,gam1,gam2,delta) {
  N <- length(x)
  
  
  if (symm==0) {
    zstrat <- cut(z, c(-Inf, -2.33,2.33,Inf))
  } else {
    zstrat <- cut(z, c(-Inf, -1.80,1.80,Inf))
  }
  zmidx <- ifelse(as.numeric(zstrat)==2,x,0)
  y <- gam1 + gam2*x + delta*zmidx + rnorm(N)
  
  return(y)
}


if (symm==0) {
  Delta <- log(seq(1,exp(0.3),length.out=6))
} else {
  Delta <- -log(seq(1,exp(0.3),length.out=6))
}


time1 <- Sys.time()
betaAll <- gam2All <- betaTarget <- c()
gam1Hat <- gam2Hat <- deltaHat <- alphaHat <- betaHat <- c() 
lrStatTrue <- lrStat <- c()
lrt <- rep(0,B)
deltaTest <- x2Test <- x3Test <- ksTest <- rep(0,B)
for (k in 1:B) {
  
  
  if (mac!=1) {
    write.table(k,paste('/home/khhan/project/multi-impute/simulation/loop-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.txt',sep=''))
  }
  
  time01 <- Sys.time()
  
  print(k)
  
  delta <- Delta[m]
  gam2Grid <- seq(-0.5,1.5,length.out=250)
  
  eps <- 1
  iter <- 0
  
  # target parameter - beta
  betaTest <- epsTest <- c()
  while (iter < length(gam2Grid)) {
    
    iter <- iter+1
    
    set.seed(B*(indMC-1)+k)
    
    x <- rnorm(N)
    
    if (symm==0) {
      z <- x+rnorm(N)
    } else {
      z <- x*rgamma(N,4,4)
    }
    # z <- z - mean(z)
    
    
    if (symm==0) {
      zstrat <- cut(z, c(-Inf, -2.33,2.33,Inf))
    } else {
      zstrat <- cut(z, c(-Inf, -1.80,1.80,Inf))
    }
    zmidx <- ifelse(as.numeric(zstrat)==2,x,0)
    
    y <- yGen(x,z,0,gam2Grid[iter],delta)
    
    fit <- glm(y~x)
    
    eps <- abs(fit$coefficients[2]-1)
    epsTest <- c(epsTest,eps)
    betaTest <- c(betaTest,fit$coefficients[2])
    
  }
  
  # plot(gam2Grid,epsTest)
  
  gam2 <- gam2Grid[which(epsTest==min(epsTest))]
  beta <- betaTest[which(epsTest==min(epsTest))]
  
  gam2All[k] <- gam2
  betaTarget[k] <- beta
  
  
  # target parameter - alpha
  alphaTest0 <- alphaTest <- gam1Test <- c()
  options(warn = -1) 
  for (i in 1:500) {
    
    set.seed(i*5000)
    
    x <- rnorm(N)
    
    if (symm==0) {
      z <- x+rnorm(N)
    } else {
      z <- x*rgamma(N,4,4)
    }
    # z <- z - mean(z)
    
    
    if (symm==0) {
      zstrat <- cut(z, c(-Inf, -2.33,2.33,Inf))
    } else {
      zstrat <- cut(z, c(-Inf, -1.80,1.80,Inf))
    }
    zmidx <- ifelse(as.numeric(zstrat)==2,x,0)
    
    y <- yGen(x,z,0,gam2,delta)
    
    id <- 1:N
    insample <- id %in% c(which(zstrat==levels(zstrat)[1]), which(zstrat==levels(zstrat)[3]),sample(which(zstrat==levels(zstrat)[2]),n))
    ind <- which(insample==TRUE)
    
    model0 <- glm(y~x)
    modelm <- glm(y[ind]~x[ind])
    modelTmp <- glm(y[ind]~x[ind]+zmidx[ind])
    
    alphaTest0 <- c(alphaTest0,model0$coefficients[1])
    alphaTest <- c(alphaTest,modelm$coefficients[1])
    
    gam1Test <- c(gam1Test,modelTmp$coefficients[1])
    
  }
  
  alphaFull <- mean(alphaTest0)
  alphaSub <- mean(alphaTest)
  gam1Sub <- mean(gam1Test)
  
  
  
  
  ### data generation
  set.seed(B*(indMC-1)+k)
  
  xFull <- x <- rnorm(N)
  
  if (symm==0) {
    z <- x+rnorm(N)
  } else {
    z <- x*rgamma(N,4,4)
  }
  # z <- z - mean(z)
  
  
  if (symm==0) {
    zstrat <- cut(z, c(-Inf, -2.33,2.33,Inf))
  } else {
    zstrat <- cut(z, c(-Inf, -1.80,1.80,Inf))
  }
  zmidx <- ifelse(as.numeric(zstrat)==2,x,0)
  
  y <- yGen(x,z,0,gam2,delta)
  
  id <- 1:N
  insample <- id %in% c(which(zstrat==levels(zstrat)[1]), which(zstrat==levels(zstrat)[3]),sample(which(zstrat==levels(zstrat)[2]),n))
  
  x[which(insample==FALSE)] <- NA
  
  df <- data.frame(y=y,x=x,xFull=xFull,z=z,zstrat=zstrat,zmidx=zmidx,id=id,insample=insample)
  sdf <- df[insample,]
  
  
  fitTrue <- lm(y~x+zmidx,data=sdf)
  deltaTest[k] <- ifelse(summary(fitTrue)$coef[3,4]<0.05,1,0)
  
  fitX2 <- lm(y~x+x2,data=data.frame(sdf,x2=sdf$x^2))
  x2Test[k] <- ifelse(summary(fitX2)$coef[3,4]<0.05,1,0)
  
  fitX3 <- lm(y~x+x2+x3,data=data.frame(sdf,x2=sdf$x^2,x3=sdf$x^3))
  x3Test[k] <- ifelse(summary(fitX3)$coef[4,4]<0.05,1,0)
  
  fitLm <- lm(y~x,data=sdf)
  bwsNp <- npregbw(y~x,data=sdf)
  fitNp <- npreg(bwsNp)
  
  # plot(sdf$y,fitLm$fitted.values)
  # points(sdf$y,fitNp$mean,col=2)
  # 
  # plot(fitNp)
  
  tHat <- mean(fitLm$residuals^2) - fitNp$MSE
  
  seqWB <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
  probWB <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  tRe <- tBoot <- c()
  time1Test <- Sys.time()
  for (tt in 1:250) {
    
    # xFullRe <- xRe <- rnorm(N)
    # 
    # if (symm==0) {
    #   zRe <- xRe+rnorm(N)
    # } else {
    #   zRe <- xRe*rgamma(N,4,4)
    # }
    # # z <- z - mean(z)
    # 
    # 
    # if (symm==0) {
    #   zstratRe <- cut(zRe, c(-Inf, -2.33,2.33,Inf))
    # } else {
    #   zstratRe <- cut(zRe, c(-Inf, -1.80,1.80,Inf))
    # }
    # zmidxRe <- ifelse(as.numeric(zstratRe)==2,xRe,0)
    # 
    # yRe <- yGen(xRe,zRe,0,beta,0)
    # 
    # idRe <- 1:N
    # insampleRe <- idRe %in% c(which(zstratRe==levels(zstratRe)[1]), which(zstratRe==levels(zstratRe)[3]),sample(which(zstratRe==levels(zstratRe)[2]),n))
    # 
    # xRe[which(insampleRe==FALSE)] <- NA
    # 
    # dfRe <- data.frame(y=yRe,x=xRe,xFull=xFullRe,z=zRe,zstrat=zstratRe,zmidx=zmidxRe,id=idRe,insample=insampleRe)
    # sdfRe <- dfRe[insampleRe,]
    # 
    # fitLmRe <- lm(y~x,data=sdfRe)
    # fitNpRe <- npreg(y~x,data=sdfRe,bws=bwsNp$bw)
    # tRe[tt] <- mean(fitLmRe$residuals^2) - fitNpRe$MSE
    
    
    s2Boot <- sum(fitLm$residuals^2)/(nrow(sdf)-2)
    yBoot <- with(sdf, fitLm$fitted.values + fitLm$residuals*sample(seqWB,nrow(sdf),prob=probWB,replace=TRUE))
    
    # yBoot <- rnorm(nrow(sdf),yBoot,sqrt(s2Boot))
    
    fitLmBoot <- lm(yBoot~sdf$x)
    fitNpBoot <- npreg(yBoot~sdf$x,bws=bwsNp$bw)
    tBoot[tt] <- mean(fitLmBoot$residuals^2) - fitNpBoot$MSE
    
  }
  time2Test <- Sys.time()
  
  time2Test - time1Test
  
  # par(mfrow=c(1,2))
  # hist(tRe)
  # abline(v=tHat,col=2)
  # mean(tRe>tHat)
  # 
  # hist(tBoot)
  # abline(v=tHat,col=2)
  # mean(tBoot>tHat)
  
  
  # ksTest[k] <- ifelse(mean(tRe>tHat)<0.05,1,0)
  ksTest[k] <- ifelse(mean(tBoot>tHat)<0.05,1,0)
  
  
  
  ## semi-parametric mle
  df$obstype.name<-ifelse(insample,"retro","strata")
  yCuts<-c(qnorm(0.025),-1,0,1,qnorm(0.975))
  
  options(warn=-1)
  mle0<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="ycutmeth",obstype.name="obstype.name",yCuts=yCuts,errdistn="normal", print.progress=0)
  mle<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="direct",obstype.name="obstype.name",start=mle0$coefficients,errdistn="normal", print.progress=0)
  options(warn=0)
  
  
  
  
  ## imputation model 
  # impmodel1<-glm(x~z,data=sdf)
  impmodel2<-glm(x~y+z,data=sdf)
  
  
  
  
  ## conventional calibration - single imputation
  calmodel <- glm(x~z,data=sdf)
  
  betaCal <- calmodel$coefficients
  fullZ <- cbind(rep(1,nrow(df)),df$z)
  
  xCal <- c(fullZ%*%betaCal)
  xCal <- ifelse(is.na(df$x), xCal,df$x)
  
  modelCal<- glm(y~xCal)
  
  
  
  
  ## calibrated - raking
  
  # calibration model
  predx<-predict(impmodel2,newdata=df)
  df1 <- cbind(df,predx)
  calmodel2<-glm(y~predx,data=df1)
  
  # plot(df1$predx1,df1$xFull,col=df1$insample+1)
  # abline(coef=c(0,1))
  
  # influence function
  eif<-as.data.frame(estfun.glm(calmodel2))
  names(eif)<-c("eif1","eif2")
  
  df2<-cbind(df1,eif)
  
  # calibrated design for raking
  des<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df2, subset=~insample)
  cdes2<-calibrate(des,formula=~(eif1+eif2)*zstrat,calfun="raking",epsilon=1e-6)
  
  
  
  
  ## mi-raking: wild bootstrapping & Bayesian resampling
  p<-vector("list",M)
  
  seqWB <- c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
  probWB <- c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5)))
  
  # seqWB <- c(-1,1)
  # probWB <- c(1/2,1/2)
  
  xHat <- impmodel2$fitted.values
  eHat <- impmodel2$residuals
  
  subYZ <- cbind(rep(1,length(xHat)),sdf$y,sdf$z)
  fullYZ <- cbind(rep(1,nrow(df)),df$y,df$z)
  
  fitTmp <- data.frame(xHat=xHat,eHat=eHat)
  sdfTmp <- cbind(sdf,fitTmp)
  for(i in 1:M){
    p[[i]]<-df
    
    set.seed(i*M)
    
    # wild bootstrap    
    xNew <- with(sdfTmp, xHat + eHat*sample(seqWB,nrow(sdfTmp),prob=probWB,replace=TRUE))
    sdfTmp$xNew <- xNew 
    
    impmodel2WBoot<-glm(xNew~y+z,data=sdfTmp)
    
    s2WBoot <- sum(impmodel2WBoot$residuals^2)/(nrow(sdfTmp)-3)
    betaWBoot <- impmodel2WBoot$coefficients
    
    xWBoot <- c(fullYZ%*%betaWBoot)
    resamxWBoot <- rnorm(nrow(df),xWBoot,sqrt(s2WBoot))
    
    # plot(df$xFull,resamxWBoot)
    # points(df$xFull,predict(impmodel2,newdata=df),col=2)
    
    p[[i]]$resamxWBoot <- resamxWBoot
    p[[i]]$impdxWBoot <- with(p[[i]], ifelse(is.na(x), resamxWBoot,x))
    
    
    set.seed(i*M)
    
    # Bayesian
    s2Bayes <- sum(eHat^2)*rinvchisq(1,df=nrow(sdfTmp)-3)
    betaBayes <- mvrnorm(1,impmodel2$coefficients,s2Bayes*solve(t(subYZ)%*%subYZ))
    
    xBayes <- c(fullYZ%*%betaBayes)
    resamxBayes <- rnorm(nrow(df),xBayes,sqrt(s2Bayes))
    
    # plot(df$xFull,resamxBayes)
    # points(df$xFull,resamxWBoot,col=3)
    # points(df$xFull,predict(impmodel2,newdata=df),col=2)
    
    # plot(density(resamxBayes))
    # lines(density(resamxWBoot),col=3)
    # lines(density(predict(impmodel2,newdata=df)),col=2)
    
    p[[i]]$resamxBayes <- resamxBayes
    p[[i]]$impdxBayes <- with(p[[i]], ifelse(is.na(x), resamxBayes,x))
    
  }
  plist<-imputationList(p)
  
  modelmiWBoot<-with(plist, glm(y~impdxWBoot))
  modelmiBayes<-with(plist, glm(y~impdxBayes))
  
  # modelmiWBoot<-with(plist, glm(y~resamxWBoot))
  # modelmiBayes<-with(plist, glm(y~resamxBayes))
  
  
  # wild bootstrap
  # mean influence function
  modelmiRkWBoot<-with(plist, glm(y~resamxWBoot))
  eifWBoot<-as.data.frame(Reduce("+",lapply(modelmiRkWBoot, estfun.glm))/M)
  names(eifWBoot)<-c("mieif1","mieif2")
  df3<-cbind(df,eifWBoot)
  
  # calibrated design for mi-para
  desmiRkWBoot<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df3, subset=~insample)
  cdesmiRkWBoot<-calibrate(desmiRkWBoot,formula=~(mieif1+mieif2)*zstrat,calfun="raking",epsilon=1e-6)
  
  
  # Bayesian
  # mean influence function
  modelmiRkBayes<-with(plist, glm(y~resamxBayes))
  eifBayes<-as.data.frame(Reduce("+",lapply(modelmiRkBayes, estfun.glm))/M)
  names(eifBayes)<-c("mieif1","mieif2")
  df4<-cbind(df,eifBayes)
  
  # calibrated design for mi-para
  desmiRkBayes<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df4, subset=~insample)
  cdesmiRkBayes<-calibrate(desmiRkBayes,formula=~(mieif1+mieif2)*zstrat,calfun="raking",epsilon=1e-6)
  
  
  
  
  ## estimates - beta
  betaVec <- data.frame(mlf=coef(glm(y~xFull,data=df))[2],
                        mls=summary(mle)$coef.table2[2,1],
                        rk=coef(svyglm(y~x,design=cdes2))[2],
                        miRkWBoot=coef(svyglm(y~x,design=cdesmiRkWBoot))[2],
                        miRkBayes=coef(svyglm(y~x,design=cdesmiRkBayes))[2],
                        miWBoot=coef(MIcombine(modelmiWBoot))[2],
                        miBayes=coef(MIcombine(modelmiBayes))[2],
                        impCal=coef(modelCal)[2])
  # mib=mean(coefBoot),
  # mibKde=coef(svyglm(y~x,design=cdesmiBoot))[2])
  
  round(betaVec-beta,4)
  betaAll <- rbind(betaAll,betaVec)
  

  
  
  ## correlation - difference and log-likelihood ratio
  modelTrue <- glm(y~x+zmidx,data=sdf)
  modelNearlyTrue <- glm(y~x,data=sdf)
  
  gam1Hat[k] <- modelTrue$coefficients[1]
  gam2Hat[k] <- modelTrue$coefficients[2]
  deltaHat[k] <- modelTrue$coefficients[3]
  
  alphaHat[k] <- modelNearlyTrue$coefficients[1]
  betaHat[k] <- modelNearlyTrue$coefficients[2]
  
  lrStatTrue[k] <- -2*(logQn(y[insample],x[insample],alphaSub,beta) - logPn(y[insample],x[insample],z[insample],gam1Sub,gam2,delta))
  lrStat[k] <- -2*(logQn(y[insample],x[insample],alphaHat[k],betaHat[k]) - logPn(y[insample],x[insample],z[insample],gam1Hat[k],gam2Hat[k],deltaHat[k]))
  
  if (m==1) {
    
    if (lrStat[k] > 0 & deltaHat[k]>0)  {
      
      if (sqrt(lrStat[k]) > qnorm(0.95)) {
        
        lrt[k] <- 1
        
      }
    }
    
  } else {
    
    if (lrStat[k] > qchisq(0.95,df=1)) {
      
      lrt[k] <- 1
    }
    
  }
  
  Dn <- betaAll[,2]-betaAll[,3]
  
  
  
  time02 <- Sys.time()
  print(time02-time01)
  
  if (mac!=1) {
    write.table(time02-time01,paste('/home/khhan/project/multi-impute/simulation/time-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.txt',sep=''))
  }
}

time2 <- Sys.time()


colnames(betaAll) <- c('mle-full','mle-sub','raking','miRk-WBoot','miRk-Bayes','mi-WBoot','mi-Bayes','impCal')

betaMean <- apply(betaAll,2,'mean')
tmp <- matrix(rep(betaMean,rep(B,length(betaMean))),nrow=B,ncol=length(betaMean))

betaMse <- apply((betaAll-1)^2,2,'mean')
betaBias <- betaMean-1
betaVar <- (apply((betaAll-tmp)^2,2,'mean'))

Dn <- betaAll[,2]-betaAll[,3]
DnWBoot <- betaAll[,6]-betaAll[,4]
DnBayes <- betaAll[,7]-betaAll[,5]

qwer <- list()
qwer[[1]] <- round(c(mean(gam2All),delta),4)
qwer[[2]] <- round(sqrt(betaMse),4)
qwer[[3]] <- round(betaBias,4)
qwer[[4]] <- round(sqrt(betaVar),4)
qwer[[5]] <- round(cor(Dn,lrStatTrue),4)
qwer[[6]] <- round(cor(DnWBoot,lrStatTrue),4)
qwer[[7]] <- round(cor(DnBayes,lrStatTrue),4)
qwer[[8]] <- round(c(mean(lrt),sd(lrt)),4)
qwer[[9]] <- round(c(mean(deltaTest),sd(deltaTest)),4)
qwer[[10]] <- round(c(mean(ksTest),sd(ksTest)),4)
qwer[[11]] <- round(c(mean(x2Test),sd(x2Test)),4)
qwer[[12]] <- round(c(mean(x3Test),sd(x3Test)),4)

names(qwer) <- c('gam2Delta',
                 'sqrtMse','bias','sqrtVar',
                 'corrBiaslogLR','corrWBoot','corrBayes',
                 'NPtest','deltaTest','ksTest','x2Teset','x3Test')

print(qwer)

print(time2-time1)



if (mac==1) {
    save.image(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/simulation/test/second-04262019-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.RData',sep=''))
} else {
    save.image(paste('/home/khhan/project/multi-impute/simulation/second-04262019-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.RData',sep=''))
}








symm <- 1

# if (symm==1) {
#   numB <- 1
#   B <- 100
# } else {
  B <- 250
  numB <- 1000/B
# }
M <- 100

N <- 5000
n <- N*0.05

if (symm==0) {
  Delta <- log(seq(1,exp(0.3),length.out=6))
} else {
  Delta <- -log(seq(1,exp(0.3),length.out=6))
}


mseAll <- biasAll <- varAll <- matrix(nrow=length(Delta),ncol=8)
power <- c()
ksPower <- c()
x3Power <- c()
for (m in 1:length(Delta)) {

  print(paste('Delta = ',round(Delta[m],3),sep=''))

  betaAllMC <- matrix(nrow=numB*B,ncol=8)
  colnames(betaAllMC) <- c('mle-full','mle-sub','raking','miRk-WBoot','miRk-Bayes','mi-WBoot','mi-Bayes','impCal')
  
  DnMC <- DnWBootMC <- DnBayesMC <- c()
  lrStatTrueMC <- c()
  lrtMC <- c()
  deltaTestMC <- ksTestMC <- x2TestMC <- x3TestMC <- c()


  for (indMC in 1:numB) {
    print(indMC)

    # if (symm==0) {
      load(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/simulation/second-04262019-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.RData',sep=''))
    # } else {
    #   load(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/simulation/test/second-04262019-indMC',indMC,'-N',N,'-n',n,'-B',B,'-m',m,'-symm',symm,'.RData',sep=''))
    # }

    betaAllMC[(B*(indMC-1)+1):(B*indMC),] <- as.matrix(betaAll)

    DnMC <- c(DnMC,Dn)
    DnWBootMC <- c(DnWBootMC,DnWBoot)
    DnBayesMC <- c(DnBayesMC,DnBayes)

    lrStatTrueMC <- c(lrStatTrueMC,lrStatTrue)

    lrtMC <- c(lrtMC,lrt)
    deltaTestMC <- c(deltaTestMC,deltaTest)
    ksTestMC <- c(ksTestMC,ksTest)
    x2TestMC <- c(x2TestMC,x2Test)
    x3TestMC <- c(x3TestMC,x3Test)
    
  }

  betaMean <- apply(betaAllMC,2,'mean')
  tmp <- matrix(rep(betaMean,rep(numB*B,length(betaMean))),nrow=numB*B,ncol=length(betaMean))

  betaMse <- apply((betaAllMC-1)^2,2,'mean')
  betaBias <- betaMean-1
  betaVar <- (apply((betaAllMC-tmp)^2,2,'mean'))

  qwer <- list()
  qwer[[1]] <- round(c(mean(gam2All),delta),3)
  qwer[[2]] <- round(sqrt(betaMse),3)
  qwer[[3]] <- round(betaBias,3)
  qwer[[4]] <- round(sqrt(betaVar),3)
  qwer[[5]] <- round(cor(DnMC,lrStatTrueMC),3)
  qwer[[6]] <- round(cor(DnWBootMC,lrStatTrueMC),3)
  qwer[[7]] <- round(cor(DnBayesMC,lrStatTrueMC),3)
  qwer[[8]] <- round(c(mean(lrtMC),sd(lrtMC)),3)
  qwer[[9]] <- round(c(mean(deltaTestMC),sd(deltaTestMC)),3)
  qwer[[10]] <- round(c(mean(ksTestMC),sd(ksTestMC)),3)
  qwer[[11]] <- round(c(mean(x2TestMC),sd(x2TestMC)),3)
  qwer[[12]] <- round(c(mean(x3TestMC),sd(x3TestMC)),3)

  names(qwer) <- c('gam2Delta',
                   'sqrtMse','bias','sqrtVar',
                   'corrBiaslogLR','corrWBoot','corrBayes',
                   'NPtest','deltaTest','ksTest','x2Teset','x3Test')

  mseAll[m,] <- qwer$sqrtMse
  biasAll[m,] <- qwer$bias
  varAll[m,] <- qwer$sqrtVar

  power[m] <- qwer$NPtest[1]
  ksPower[m] <- qwer$ksTest[1]
  x3Power[m] <- qwer$x3Test[1]

  print(qwer)
}



mseAll <- mseAll[,-1]
biasAll <- biasAll[,-1]
varAll <- varAll[,-1]




if (symm==0) {

  # par(mar=c(9,6,1,10)+0.1)
  # pdf('~/Dropbox/research/Combining multiple imputation with ranking of weights/manuscript/figures/table2.pdf',width=10,height=7)
  par(mar=c(9,6,1,10)+0.1)

  plot(Delta,mseAll[,1],type='b',lwd=5,pch=1,col=1,
       ylim=c(0,0.15),xlim=c(0,Delta[length(Delta)]*1.1),
       axes=FALSE,cex=2,cex.lab=2,
       xlab='',ylab=expression(sqrt('MSE')))
  axis(2,round(seq(0,0.15,0.025),3),line=0,
       cex.axis=1.75)
  axis(1,round(Delta,3),
       cex.axis=1.75)
  mtext(expression(delta),1,at=Delta[length(Delta)]*1.15,line=0,
        cex=2)

  for (k in 2:6) {
    if (k==4 | k==6) {
      points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,lty=4,
             cex=2)
    } else {
      points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,
             cex=2)
    }
  }
  points(Delta,mseAll[,7],type='b',lwd=5,pch=7,col='darkgoldenrod3',lty=1,
         cex=2)
  axis(1,round(Delta,3),labels=round(power,3),line=3,
       cex.axis=1.75)
  mtext(' Power (MP)',1,at=Delta[length(Delta)]*1.25,line=3,
        cex=1.75)

  # axis(1,round(Delta,3),labels=round(x3Power,3),line=6,
  #      cex.axis=1.75)
  # mtext(' (parametric)',1,at=Delta[length(Delta)]*1.25,line=6,
  #       cex=1.75)
  #
  # axis(1,round(Delta,3),labels=round(ksPower,3),line=9,
  #      cex.axis=1.75)
  # mtext('(nonparametric)',1,at=Delta[length(Delta)]*1.3,line=9,
  #       cex=1.75)

  axis(1,round(Delta,3),labels=round(ksPower,3),line=6,
       cex.axis=1.75)
  mtext('(linearity test)',1,at=Delta[length(Delta)]*1.275,line=6,
        cex=1.75)

  abline(v=round(Delta,3),h=round(seq(0,0.15,0.025),3),col=8,lty=3)
  legend('topright',c('MLE','Raking','MIR-Boot','MIR-Bayes','MI-Boot','MI-Bayes','Calib'),
         lwd=4,col=c(1:6,'darkgoldenrod3'),pch=1:7,bty='n',cex=1.5,lty=c(1,1,1,4,1,4,1),
         xpd=TRUE,inset=c(-0.3,0))
  # dev.off()

} else {

  # par(mar=c(12,6,1,10)+0.1)
  # pdf('~/Dropbox/research/Combining multiple imputation with ranking of weights/manuscript/figures/table3.pdf',width=10,height=7)
  par(mar=c(9,6,1,10)+0.1)

  plot(Delta,mseAll[,1],type='b',lwd=5,pch=1,col=1,
       ylim=c(0,0.3),xlim=c(0,Delta[length(Delta)]*1.1),
       axes=FALSE,cex=2,cex.lab=2,
       xlab='',ylab=expression(sqrt('MSE')))
  axis(2,round(seq(0,0.3,0.05),3),line=0,
       cex.axis=1.75)
  axis(1,round(Delta,3),
       cex.axis=1.75)
  mtext(expression(delta),1,at=Delta[length(Delta)]*1.25,line=0,
        cex=2)

  for (k in 2:6) {
    if (k==4 | k==6) {
      points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,lty=4,
             cex=2)
    } else {
      points(Delta,mseAll[,k],type='b',lwd=5,pch=k,col=k,
             cex=2)
    }
  }
  points(Delta,mseAll[,7],type='b',lwd=5,pch=7,col='darkgoldenrod3',lty=1,
         cex=2)
  axis(1,round(Delta,3),labels=round(power,3),line=3,
       cex.axis=1.75)
  mtext('Power (MP)',1,at=Delta[length(Delta)]*1.25,line=3,
        cex=1.75)

  # axis(1,round(Delta,3),labels=round(x3Power,3),line=6,
  #      cex.axis=1.75)
  # mtext(' (parametric)',1,at=Delta[length(Delta)]*1.25,line=6,
  #       cex=1.75)
  #
  # axis(1,round(Delta,3),labels=round(ksPower,3),line=9,
  #      cex.axis=1.75)
  # mtext('(nonparametric)',1,at=Delta[length(Delta)]*1.3,line=9,
  #       cex=1.75)

  axis(1,round(Delta,3),labels=round(ksPower,3),line=6,
       cex.axis=1.75)
  mtext('(linearity test)',1,at=Delta[length(Delta)]*1.275,line=6,
        cex=1.75)

  abline(v=round(Delta,3),h=round(seq(0,0.3,0.05),3),col=8,lty=3)
  legend('topright',c('MLE','Raking','MIR-Boot','MIR-Bayes','MI-Boot','MI-Bayes','Calib'),
         lwd=4,col=c(1:6,'darkgoldenrod3'),pch=1:7,bty='n',cex=1.5,lty=c(1,1,1,4,1,4,1),
         xpd=TRUE,inset=c(-0.3,0))
  # dev.off()

}

