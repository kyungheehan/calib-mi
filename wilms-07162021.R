# rm(list=ls(all.names=TRUE))

# source('~/Dropbox/research/Combining multiple imputation with ranking of weights/data analysis/wilms-07162021.R')


library(survey)
library(mitools)
library(MASS)
library(stats)

invLogit <- function(t) exp(t)/(1+exp(t))

estfun.glm<-function (model)
{
  xmat <- model.matrix(model)
  residuals(model, "working") * model$weights * xmat
}

nwts <- read.table('~/Dropbox/research/Combining multiple imputation with ranking of weights/data analysis/nwts-expanded.txt',header=TRUE)
# head(nwts)

stage34 <- ifelse(nwts$stage%in%c(1,2),0,1)
nwts <- data.frame(nwts,stage34=stage34)

fitFull <- glm(relaps~histol*stage34+age+tumdiam,family=binomial,data=nwts)

impModelFull <- glm(histol~relaps*instit*stage34+age+tumdiam,family=binomial,data=nwts)
xMatFull <- model.matrix(impModelFull)

summary(fitFull)

fitForm <- fitFull$formula
impForm <- impModelFull$formula

B <- 1000
M <- 100
betaAll <- list()
mse <- avg <- matrix(0,nrow=6,ncol=5)
time1MC <- Sys.time()
m <- 0
B0 <- 1
while (B0 <= B) {
  
  try(
    {
      m <- m+1
      time1 <- Sys.time()
      print(c(m,B0))
        
      set.seed(B*m)
      
      # some index sets
      indRe1 <- with(nwts,which(relaps==1))
      indRe0 <- with(nwts,which(relaps==0))
      
      indUH1 <- with(nwts,which(instit==1))
      indUH0 <- with(nwts,which(instit==0))
      
      
      # # relapse sampling
      # subSample <- c(indRe1,sample(indRe0,length(indRe1)))
      # indSample <- (1:nrow(nwts)) %in% subSample
      
      
      # # # instit sampling
      # subSample <- c(indUH1,sample(indUH0,length(indUH1)))
      # indSample <- (1:nrow(nwts)) %in% subSample
      
      
      # # instit-relaps sampling
      # ind00 <- with(nwts,which(instit==0 & relaps==0))
      # ind01 <- with(nwts,which(instit==0 & relaps==1))
      # ind10 <- with(nwts,which(instit==1 & relaps==0))
      # ind11 <- with(nwts,which(instit==1 & relaps==1))
      # 
      # subSample <- c(indUH1,ind01,with(nwts,sample(ind00,length(indRe1)-length(ind10))))
      # indSample <- (1:nrow(nwts)) %in% subSample
      
      
      # stage-instit-relaps sampling
      ind000 <- with(nwts,which(stage34==0 & instit==0 & relaps==0))
      ind001 <- with(nwts,which(stage34==0 & instit==0 & relaps==1))
      ind010 <- with(nwts,which(stage34==0 & instit==1 & relaps==0))
      ind011 <- with(nwts,which(stage34==0 & instit==1 & relaps==1))

      ind100 <- with(nwts,which(stage34==1 & instit==0 & relaps==0))
      ind101 <- with(nwts,which(stage34==1 & instit==0 & relaps==1))
      ind110 <- with(nwts,which(stage34==1 & instit==1 & relaps==0))
      ind111 <- with(nwts,which(stage34==1 & instit==1 & relaps==1))

      subSample <- c(indRe1,ind010,ind110,
                     with(nwts,sample(ind000,length(ind001)+length(ind011)-length(ind010))),
                     with(nwts,sample(ind100,length(ind101)+length(ind111)-length(ind110))))
      indSample <- (1:nrow(nwts)) %in% subSample
      
      
      # data frame
      df <- data.frame(nwts,indSample=indSample)
      df$histol[which(indSample==FALSE)] <- NA
      
      sdf <- df[subSample,]
      # table(data.frame(sdf$instit,sdf$histol))
      # table(data.frame(sdf$instit,sdf$relaps))
      
      
      ## complete case
      fitSub <- glm(fitForm,family=binomial,data=sdf)
      
      
      ## imputation model
      impModel <- glm(impForm,family=binomial,data=sdf)
      
      xMat <- model.matrix(impModel)
      betaHat <- impModel$coefficients
      covHat <- vcov(impModel)
      
      
      ## raking
      
      # calibration model
      histolPred <- ifelse(invLogit(predict(impModel,newdata=df))>0.5,1,0)
      
      # hist(invLogit(predict(impModel,newdata=df)))
      # table(data.frame(histolPred,df$histol))
      # table(data.frame(sdf$instit,sdf$histol))
      # table(data.frame(sdf$relaps,sdf$histol))
      # table(data.frame(df$instit,histolPred))
      
      df1 <- data.frame(df,histolPred=histolPred)
      calibModel <- glm(relaps~histolPred*stage34+age+tumdiam,family=binomial,data=df1)
      
      # influence function
      eif <- as.data.frame(estfun.glm(calibModel))
      names(eif) <- c('eif1','eif2','eif3','eif4','eif5','eif6')
      
      df2 <- data.frame(df1,eif)
      
      # calibrated design for raking
      des <- twophase(id=list(~1,~1), subset=~indSample, strata=list(NULL,~instit), data=df2)
      cdes <- calibrate(des,formula=~(eif1+eif2+eif3+eif4+eif5+eif6)*instit,calfun="raking",epsilon=1e-4)
      
      
      
      ## mi and mi+raking
      p<-vector("list",M)
      
      for(i in 1:M){
        p[[i]]<-df
        
        set.seed(i*B)
        
        # # asymptotic distribution
        # betaRe <- mvrnorm(1,mu=betaHat,Sigma=covHat)
        
        # bootstrap
        indResam <- sample(1:nrow(sdf),nrow(sdf),replace=TRUE)
        
        impModelResam <- glm(impForm,family=binomial,data=sdf[indResam,])
        betaRe <- impModelResam$coefficients
        
        
        probRe <- invLogit(xMatFull%*%betaRe)
        
        # plot(probRe,impModelFull$fitted.values)
        # abline(coef=c(0,1),col=2)
        
        histolRe <- rbinom(nrow(df),1,prob=probRe)
        
        p[[i]]$histolRe <- histolRe
        p[[i]]$histolImp <- with(p[[i]], ifelse(indSample, histol, histolRe))
        
      }
      plist<-imputationList(p)
      
      # mi
      modelmi <- with(plist, glm(relaps~histolImp*stage34+age+tumdiam,family=binomial))
      # modelmi <- with(plist, glm(relaps~histolRe*stage34+age+tumdiam,family=binomial))
      
      # mean influence function
      modelMIR <- with(plist, glm(relaps~histolRe*stage34+age+tumdiam,family=binomial))
      eifMIR <- as.data.frame(Reduce("+",lapply(modelMIR, estfun.glm)))
      names(eifMIR) <- c('eif1','eif2','eif3','eif4','eif5','eif6')
      df3 <- data.frame(df,eifMIR)
      
      # calibrated design for mi-para
      desMIR <- twophase(id=list(~1,~1), subset=~indSample, strata=list(NULL,~instit), data=df3)
      cdesMIR <- calibrate(desMIR,formula=~(eif1+eif2+eif3+eif4+eif5+eif6)*instit,calfun="raking",epsilon=1e-4)
      
      
      
      ## conventional calibration - single imputation
      calmodel <- glm(histol~instit,family=binomial,data=sdf)
      
      betaCal <- calmodel$coefficients
      fullZ <- cbind(rep(1,nrow(df)),df$instit)
      
      xCal <- c(fullZ%*%betaCal)
      xCal <- ifelse(exp(xCal)/(1+exp(xCal))>0.5, 1, 0)
      xCal <- ifelse(is.na(df$histol), xCal,df$histol)
      # table(data.frame(xCal,df$histol))
      
      dfCal <- data.frame(df,histolCal=xCal)
      
      modelCal<- glm(relaps~histolCal*stage34+age+tumdiam,family=binomial,data=dfCal)
      
      
      ## estimates - beta
      betaVec <- data.frame(mlf=coef(fitFull),
                            mls=coef(fitSub),
                            rk=coef(svyglm(fitForm,family=binomial,design=cdes)),
                            mi=coef(MIcombine(modelmi)),
                            mir=coef(svyglm(fitForm,family=binomial,design=cdesMIR)),
                            impCal=coef(modelCal))
      
      betaFull <- matrix(rep(betaVec$mlf,ncol(betaVec)),nrow=length(betaVec$mlf),ncol=ncol(betaVec))
      mse <- mse + ((betaVec-betaFull)^2)[,-1]
      avg <- avg + betaVec[,-1]
      
      betaAll[[m]] <- betaVec
      
      time2 <- Sys.time()
      print(time2-time1)
      

      B0 <- m-sum(unlist(lapply(betaAll,'is.null')))+1
      
    }, silent=TRUE
  )
    
} 
time2MC <- Sys.time()
print(time2MC-time1MC)

err <- rep(1,m)
for (k in 1:m) {
  if (is.null(betaAll[[k]])) {
    err[k] <- 0
  }
}
print(sum(err))

mse <- mse/B
avg <- avg/B

# mse
round(sqrt(mse),3)
round(apply(mse[-1,],2,'sum'),3)

# bias
bias <- avg-betaFull
round(bias,3)
round(apply(bias[-1,]^2,2,'sum'),3)

# var
# var <- round(rbind(mse,sum=apply(mse,2,'sum')),4) - round(rbind(bias^2,sum=apply(bias^2,2,'sum')),4)
round(sqrt(mse-bias^2),3)
round(apply(mse[-1,],2,'sum'),3) - round(apply(bias[-1,]^2,2,'sum'),3)


save.image(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/data analysis/wilms-07162021-B',B,'-M',M,'.RData',sep=''))




### load

B <- 1000 # 500
M <- 100

load(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/data analysis/wilms-07162021-B',B,'-M',M,'.RData',sep=''))

summary(fitFull)

# mse
round(sqrt(mse),3)
round(apply(mse[-1,],2,'sum'),3)

# bias
bias <- avg-betaFull
round(bias,3)
round(apply(bias[-1,]^2,2,'sum'),3)

# var
# var <- round(rbind(mse,sum=apply(mse,2,'sum')),4) - round(rbind(bias^2,sum=apply(bias^2,2,'sum')),4)
round(sqrt(mse-bias^2),3)
round(apply(mse[-1,],2,'sum'),3) - round(apply(bias[-1,]^2,2,'sum'),3)




# varName0 <- c('Intercept','Histology','Stage','Age','Tumor Diameter','Histology*Stage')
# varName <- rep(varName0, each=4*B)
# method <- rep(c('MLE','Raking','MI','MIR'), each=B)
# 
# indFit <- which(err==1)
# fitDev <- matrix(nrow=4*B,ncol=6)
# for (j in 1:6) {
#   tmpDev <-matrix(nrow=B,ncol=4)
#   for (m in 1:B) {
#     tmpDev[m,] <- as.numeric(betaAll[[indFit[m]]][j,-1] - betaFull[j,1])
#   }
#   fitDev[,j] <- c(tmpDev)
# }
# fitDevVec <- c(fitDev)
# 
# data <- data.frame(varName,method,fitDevVec)
# 
# 
# 
# for (j in 1:6) {
# 
#   # par(mar=c(3,5,1,1))
#   # pdf(paste('~/Dropbox/research/Combining multiple imputation with ranking of weights/manuscript/figures/regressor',j,'.pdf',sep=''))
#   par(mar=c(3,5,1,1))
#   fitHist <- data.frame(matrix(fitDev[,j],nrow=B,ncol=4))
#   names(fitHist) <- c('MLE','Raking','MI','MIR')
#   boxplot(fitHist,#col=c(1,2,5,4),
#           ylab=varName0[j],
#           lwd=3,lty=1,
#           col=c('gray','orange','seagreen1','skyblue'),
#           border=rep(c('black','brown','seagreen','navyblue'),3),
#           cex.lab=2.25,cex.axis=2.25)#,
#           # ylim=c(-1.5,1.5))
#   if (j %in% c(1,2)){
#     abline(h=seq(-2,2,0.5),v=1:4,col=8,lty=3)
#   }
#   if (j==3 | j==6) {
#     abline(h=seq(-2,2,0.2),v=1:4,col=8,lty=3)
#   } 
#   if (j==4 | j==5) {
#     abline(h=seq(-2,2,0.01),v=1:4,col=8,lty=3)
#   }
#   boxplot(fitHist,#col=c(1,2,5,4),
#           ylab=varName0[j],
#           lwd=3,lty=1,
#           col=c('gray','orange','seagreen1','skyblue'),
#           border=rep(c('black','brown','seagreen','navyblue'),3),
#           cex.lab=2.25,cex.axis=2.25,add=TRUE)#,
#           # ylim=c(-1.5,1.5))
#   abline(h=0,col=2,lwd=3,lty=2)
#   # dev.off()
# }




