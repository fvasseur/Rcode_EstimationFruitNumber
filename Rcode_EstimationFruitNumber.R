
library(glmnet)


#################################################################################################
#
# Evaluate models for estimating fruit number from inflorescence skeleton descriptors
#
#################################################################################################


# Import phenotypic data extracted from ImageJ analysis of inflorescence images
#==============================================================================

setwd("/home/ubuntu/Desktop/output_macro_InflorescenceSkeletons") # corresponding to the output folder of the RAPAmacro_InflorescenceSkeleton.txt
lsfile <- dir(, ".xls$")

# Assemble datasets
dataOut <- read.table(lsfile[1], skip=1)
dataOut <- dataOut[,c(2:11)]
names(dataOut) <- c("file", "nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                    "Triplepoints", "QuadPoints", "MaxBranchL")
for(i in c(1)) {dataOut[,i] <- as.factor(as.character(dataOut[,i]))}
for(i in c(2:10)) {dataOut[,i] <- as.numeric(as.character(dataOut[,i]))}  
dataOut <- aggregate(dataOut[,c(2:10)], by=list(dataOut$file), FUN=sum, na.rm=T)
names(dataOut)[1] <- "file"
dataOut$idPot <- strsplit(x=as.character(strsplit(x=as.character(dataOut$file[1]), split="_")[[1]][1]), split="pot")[[1]][2]          
dataOut <- dataOut[,c(11,2:10)]
dataOut$idPot <- rep(NA, length(dataOut$idPot))
dataOut <- na.omit(dataOut)
dataOut0 <- dataOut

for(i in lsfile) {
  dataOut <- read.table(i, skip=1)
  dataOut <- dataOut[,c(2:11)]
  names(dataOut) <- c("file", "nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                      "Triplepoints", "QuadPoints", "MaxBranchL")
  for(j in c(1)) {dataOut[,j] <- as.factor(as.character(dataOut[,j]))}
  for(k in c(2:10)) {dataOut[,k] <- as.numeric(as.character(dataOut[,k]))}  
  dataOut <- aggregate(dataOut[,c(2:10)], by=list(dataOut$file), FUN=sum, na.rm=T)
  names(dataOut)[1] <- "file"
  dataOut$idPot <- strsplit(x=as.character(strsplit(x=as.character(dataOut$file[1]), split="_")[[1]][1]), split="pot")[[1]][2]          
  dataOut <- dataOut[,c(11,2:10)]
  dataOut0 <- rbind(dataOut0, dataOut)
  print(i)
}

aut <- dataOut0
aut$idPot <- as.factor(as.character(aut$idPot))
aut <- aggregate(aut[,c(2:10)], by=list(aut$idPot), FUN=mean, na.rm=T)
names(aut) <- c("idPot", "nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                "Triplepoints", "QuadPoints", "MaxBranchL")
aut <- droplevels(aut)



# Import phenotypic data on fruit number manually counted in a subset of the population
#=========================================================================================
setwd("/home/ubuntu/Desktop")
man <- read.csv("Siliques_man_GT05.csv", header=T, sep=",", dec=".")
man <- man[,c(2,5)]
man <- aggregate(man[,2], by=list(man[,1]), FUN=mean, na.rm=T)
names(man) <- c("idPot", "Count")
man <- man[man$Count>0,]
man$idPot <- as.factor(as.character(man$idPot))
man <- droplevels(man)


# Merge data
##############
counttot <- merge(man, aut, by="idPot", all.y=T)
counttot$idPot <- as.factor(as.character(counttot$idPot))
counttot <- droplevels(counttot)
str(counttot)



# Perform and compare models
#--------------------------------------------------
temp0 <- data.frame(nit=NA,ntrain=NA,lm_accuracy=NA,quad_accuracy=NA,LASSO_accuracy=NA,RIDGE_accuracy=NA)

for(k in seq(10,120,10))
{
  lm_accuracy <- NA
  quad_accuracy <- NA
  LASSO_accuracy <- NA
  RIDGE_accuracy <- NA
  
  for(n in 1:100)
  {
    counttot$set <- NA
    trainsample <- sample(x=levels(man$idPot), k)
    
    for(i in 1:dim(counttot)[1])
    {
      if(counttot[i,"idPot"] %in% trainsample) {
        counttot[i,"set"] <- "train"
      } else{
        counttot[i,"set"] <- "test"
      }
    }
    
    testsample <- sample(levels(as.factor(as.character(counttot[counttot$set=="test" & counttot$Count>0, "idPot"]))), 100)
    
    test <- droplevels(counttot[counttot$idPot %in% testsample, ])
    train <- droplevels(counttot[counttot$idPot %in% trainsample, ])
    train <- na.omit(train)
    test <- na.omit(test)
    test <- droplevels(test)
    train <- droplevels(train)
    
    # lm
    model0 <- lm(Count ~ nbBranches+JunctionsVox+nbJunctions+EndPointsVox+SlabVox+AvBranchL+Triplepoints+QuadPoints+MaxBranchL, data=train)
    test$predSil0 <- predict(model0, newdata=test)
    lm_accuracy[n] <- as.numeric(((cor.test(test$Count, test$predSil0)$estimate)^2)) 
    
    # quadratic
    model1 <- lm(Count ~ poly(nbBranches,2)+poly(JunctionsVox,2)+poly(nbJunctions,2)+poly(EndPointsVox,2)+
                   poly(SlabVox,2)+poly(AvBranchL,2)+poly(Triplepoints,2)+poly(QuadPoints,2)+poly(MaxBranchL,2), data=train)
    test$predSil1 <- predict(model1, newdata=test)
    poly_accuracy[n] <- as.numeric(((cor.test(test$Count, test$predSil1)$estimate)^2)) 
    
    # LASSO model
    model3 <- glmnet(x=as.matrix(train[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                          "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                     y=as.matrix(train[,"Count"]),
                     family="gaussian",
                     alpha = 1) # for lasso, alpha=1, for ridge alpha=0)
    cv.out1 <- cv.glmnet(x=as.matrix(train[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                              "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                         y=as.matrix(train[,"Count"]), alpha = 1)
    bestlam1 <- cv.out1$lambda.min
    test$predSil3 <- predict(model3,
                             newx=as.matrix(test[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                                    "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                             s = bestlam1)
    LASSO_accuracy[n] <- as.numeric(((cor.test(test$Count, test$predSil3)$estimate)^2)) 
    
    # RIDGE model
    model4 <- glmnet(x=as.matrix(train[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                          "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                     y=as.matrix(train[,"Count"]),
                     family="gaussian",
                     alpha = 0) # for lasso, alpha=1, for ridge alpha=0)
    cv.out2 <- cv.glmnet(x=as.matrix(train[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                              "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                         y=as.matrix(train[,"Count"]), alpha = 0)
    bestlam2 <- cv.out2$lambda.min
    test$predSil4 <- predict(model4,
                             newx=as.matrix(test[,c("nbBranches", "nbJunctions", "EndPointsVox", "JunctionsVox", "SlabVox", "AvBranchL",
                                                    "Triplepoints", "QuadPoints", "MaxBranchL")]), 
                             s = bestlam2)
    RIDGE_accuracy[n] <- as.numeric(((cor.test(test$Count, test$predSil4)$estimate)^2)) 
    
    if(n %in% seq(10,100,10)){print(n)}
  }
  nit <- c(1:100)
  ntrain <- rep(k, 100)
  tot_accuracy <- cbind(nit, ntrain)
  tot_accuracy <- cbind(tot_accuracy, lm_accuracy)
  tot_accuracy <- cbind(tot_accuracy, quad_accuracy)
  tot_accuracy <- cbind(tot_accuracy, LASSO_accuracy)
  tot_accuracy <- cbind(tot_accuracy, RIDGE_accuracy)
  
  temp0 <- rbind(temp0, tot_accuracy)
  print(k)
  
  write.table(temp0, "Results_models_FruitPrediction.csv", dec=".", sep=",", row.names = F)
  
}
