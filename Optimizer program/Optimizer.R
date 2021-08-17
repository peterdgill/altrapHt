##########################################
###Logistic regressions Peak Height vs log10LR subsource fits
##This program fits RFU logistic regressions to log10 LR subsource
##data in order to find closest RFU fit
#run following lines########################
library(rstanarm)
library(parallel)
#set your working directory to find direct transfer file
dat <- read.csv("DiretTransfer.csv")
bayesianFit = TRUE# FALSE  #decide whether stan_glm or ordinary glm should be used

#Ymatrix.csv is prepared - these are subsource
#logistic regressions for log10 1..12 across 0-24h calcualated from
#direct transfer.csv file
#Ymatrix<-read.csv("mydat2.csv")
Ymatrix = read.csv("Ymatrix.csv")#,sep=" ")
Times = Ymatrix[,1] #outcome of time
PrLRModel = Ymatrix[,2:13]
colnames(PrLRModel) = 1:ncol(PrLRModel)

dat$POIcontrRfu #=  #dat$MxC1_Hp*dat$Phexp_Hp
#plot(PGfile$MxC1_Hp,PGfile$Log10LR)
#plot(dat$POIcontrRfu,dat$Log10LR)

#Defined optimise function
optimF<-function(rfu,returnPr=FALSE){#where Gtrx is logistic regression threshold;mydat2 is the logistic data for log10 LR
  dat$ProposeResponse <- ifelse(dat$POIcontrRfu>=rfu,1,0)  #obtain binary response depending on criterion

  if(bayesianFit) {
   set.seed(101)
    #Model=stan_glm( ProposeResponse ~ Time, family=binomial(link="logit"),data=dat, cores=detectCores(),refresh=0)#logistic regression create model
    Model=stan_glm( ProposeResponse ~ Time, family=binomial(link="logit"),data=dat,refresh=0)#logistic regression create model
  } else {
    Model=glm( ProposeResponse ~ Time, family=binomial(link="logit"),data=dat)#logistic regression create model
  }
  prPropose=predict(Model,newdata=data.frame(Time=Times),type="resp")#obtain probabilities

#  plot(Times,prPropose)
  if(returnPr) return(prPropose)
  score = mean((prPropose-prTarget)^2) # MSE
  return(score)
}

#Perform optimization for each LR threshold
# and write to a pdf file
RFUmin = rep(NA,ncol(PrLRModel))
pdf(paste0("Comparison_",bayesianFit,".pdf"),width=12,height=24)
par(mfrow=c(ncol(PrLRModel)/2,2)) #Create 6x2 panel
par(mar=c(3,3,2,2))  #Reduce margin to minimum (no labels)
for (i in 1:ncol(PrLRModel)) {
  print(i)
  prTarget = PrLRModel[,i] #probabilities to compare against

  #Optimize over wide range
  gridSize = 5
  grid = seq(30,500,l=gridSize)
  gridVal = sapply(grid,optimF ) #calculate start
  #plot(grid,gridVal)

  #Optimize over narrow range
  minval = which.min(gridVal)
  lim = c(max(minval-1,1), min(minval+1,gridSize)) #interal to search
  opt = optimize(optimF,interval = grid[lim])
  RFUmin[i] <- rfuLim <- opt$min #this is optimal choice

  #Check curve
  prPropose = optimF(rfuLim,returnPr=TRUE) #obtain fitted
  matplot(cbind(Times,Times),cbind(prTarget,prPropose),ty="l",xlab="",ylab="")#, main="Comparison (LR vs RFU)")
  mtext(paste0("MSE=",signif(opt$objective,3)))
  leg1=paste0("log10LR Threshold=",i)
  leg2=paste0("RFU Threshold=",round(rfuLim,0))
  legend("topright",c(leg1,leg2 ),lty=1:2,col=1:2)
}
dev.off()


