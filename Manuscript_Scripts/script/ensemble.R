#####################################################
# R commands to run ensemble machine learning methods
#####################################################

# Load packages
library(pROC)
library(randomForest)
library(lda)
library(MASS)
library(e1071)
library(glmnet)

# Read data Table S2, relabel and remove missing data
d<-read.csv("TableS2.csv",sep="\t",skip=1)
colnames(d)[grep("Training",colnames(d))]<-"Training"
d<-d[which(!is.na(d$Training)),]
d$Training<-as.factor(d$Training)
miss<-apply(d[,c(18,19,14,13,15,17,16)],1,function(x) length(which(is.na(x))))
d<-d[which(miss==0),]
d<-d[,c(12,18,19,14,13,15,17,16)]

# Model formula
model<-formula("Training ~ LRT.logistic+LRTm.logistic+PolyPhen2+Provean+SiftScore+Gerp..+MAPP")

# Logistic regression
fit<-glm(model,data=d,family=binomial())
pred<-predict(fit,d, type="response")
roc(Training~pred,data=cbind(d,pred),plot=F)$auc

# Support vector machine
fit<-svm(model,data=d)
pred<-predict(fit,d, decision.value=TRUE)
roc(Training~pred,data=cbind(d,pred=as.numeric(attr(pred,"decision.value"))),plot=F)$auc

# Random forest
fit<-randomForest(model,data=d,importance=TRUE)
roc(Training~pred,data=cbind(d,pred=fit$votes[,2]),plot=F)$auc

# Linear discriminant analysis
fit<-MASS::lda(model, data = d)
pred<-predict(fit,d, type="prob")
roc(Training~pred,data=cbind(d,pred=pred$posterior[,1]),plot=F)$auc

# Generalized linear model with penalized maximum likelihood
x<-as.matrix(d[,-1])
y<-as.numeric(as.character(d$Training))
fit = glmnet(x,y,family = "binomial")
pred<-predict(fit,type="response", newx = x,s = min(fit$lambda))
roc(Training~pred,data=cbind(d,pred=as.numeric(pred)),plot=F)$auc

