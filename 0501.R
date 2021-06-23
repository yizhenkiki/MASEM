require(metaSEM)
require(dplyr)
require(magrittr)

data <- read.csv("meta_dataset0503.csv")
data

dataset <- data[5:14]
dataset

nvar <- 5 
varnames <- c("OUT","BEH","CL","SE","PER")
labels <- list(varnames,varnames)

convert2matrix <- list()

for(i in 1:length(dataset)){
  convert2matrix[[i]] <- vec2symMat(as.matrix(dataset[i,1:10]),
                                 diag = FALSE)
  dimnames(convert2matrix[[i]]) <- labels
}
convert2matrix

#put NA on diagonal if variable is missing
for(i in 1:length(convert2matrix)){
 for(j in 1:nrow(convert2matrix[[i]])){
   if(sum(is.na(convert2matrix[[i]][j,]))==nvar-1)
     {convert2matrix[[i]][j,j] <- NA}
 }
}


for(i in 1:length(convert2matrix)){
  for(j in 1:nrow(convert2matrix[[i]])){
    for(k in 1:nvar){
      if(is.na(convert2matrix[[i]][j,k]==TRUE)
         &is.na(convert2matrix[[i]][j,j])!=TRUE
         &is.na(convert2matrix[[i]][k,k])!=TRUE){

        if(sum(is.na(convert2matrix[[i]])[j,])>sum(is.na(convert2matrix[[i]])[k,]))
               {convert2matrix[[i]][k,k] <- NA}
        if(sum(is.na(convert2matrix[[i]])[j,])<=sum(is.na(convert2matrix[[i]])[k,]))
               {convert2matrix[[i]][j.j] <- NA}
      }
    }
  }
}
convert2matrix

#Stage 1 Random
stage1random <- tssem1(convert2matrix,data$N,method = "REM",RE.type = "Diag",acov="individual")
summary(stage1random)

coef(stage1random, select="random")

A <- create.mxMatrix(
  c(0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    "0.1*b51","0.1*b52","0.1*b53","0.1*b54",0),
  type="Full",
  nrow=5,
  ncol=5,
  byrow=TRUE,
  name="A",
  dimnames=list(varnames,varnames))
A

S <- create.mxMatrix(
  c("1*p11",
    "0.1*p21","1*p22",
    "0.1*p31","0.1*p32","1*p33",
    "0.1*p41","0.1*p42","0.1*p43","1*p44",
    0,0,0,0,"1*p55"),
  type="Symm",
  byrow=TRUE,
  name="S",
  dimnames=list(varnames,varnames))
S


#Stage 2 Random
stage2 <- tssem2(stage1random, Amatrix=A, Smatrix=S,diag.constraints=FALSE, intervals="LB")
stage2

stage2 <- tssem2(stage1random, Amatrix=A, Smatrix=S)
summary(stage2)


#plot
install.packages("semPlot")
require("semPlot")

my.plot <- meta2semPlot(stage2)
semPaths(my.plot,whatLabels="est",layout="tree2")

#Subgroup Analysis

# Data for studies with majority low SES
data_low <- convert2matrix[Roorda11$SES>50]
n_low <- data$N[Roorda11$SES>50]
# Data for studies with majority high SES
data_high <- convert2matrix[Roorda11$SES<=50]
n_high <- data$N[Roorda11$SES<=50]

#Fitting a random-effects Stage 1 model in both subgroups
# Stage 1 analysis per subgroup (random-effects analysis)
stage1_low.fit <- tssem1(my.df = data_low, n = n_low, method = "REM", RE.type = "Diag")
stage1_high.fit <- tssem1(my.df = data_high, n = n_high, method = "REM", RE.type = "Diag")
summary(stage1_low.fit)
#Fitting the Stage 2 model in both subgroups
# Stage 2 analysis per subgroup (random-effect analysis)
stage2_low.fit <- tssem2(stage1_low.fit, Amatrix=A, Smatrix=S)
stage2_high.fit <- tssem2(stage1_high.fit, Amatrix=A, Smatrix=S)
summary(stage2_low.fit)
