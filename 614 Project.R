data = read.csv(file = "project.csv")
normal=function(x){return((x-min(x))/(max(x)-min(x)))}
normdat=as.data.frame(lapply(data,normal))
normdat
library(tidyverse)
library(factoextra)
library(qcc)
results = prcomp(normdat,scale = FALSE,center = TRUE)
results$rotation = -1*results$rotation
results$rotation
names(results)

results$x = -1*results$x
(results$sdev^2 / sum(results$sdev^2))*100
head(results$x)

biplot(results,scale=0)
results$sdev
var_explained = (results$sdev^2 / sum(results$sdev^2))*100
var_explained
#create scree plot
plot(var_explained,
     xlab="Principal Component",
     ylab="Variance Explained",
     ylim=c(0,50),
     main="Scree Plot" )


s=cov(normdat)
s
sum(diag(s))
eigen= eigen(s)
eigen
for (s in eigen$values) {
  print(s / sum(eigen$values))
}
results2 = prcomp(data,scale = TRUE)
fviz_eig(results2,ncp=25,addlabels = TRUE)
fviz_eig(results,ncp=25,addlabels = TRUE)

covdata = cov(normdat)
results3 = prcomp(covdata, scale = TRUE)
fviz_eig(results3, addlabels= TRUE)

covdata2 = cov(data)
results4 = prcomp(covdata2, scale = TRUE)
fviz_eig(results4, addlabels= TRUE)

covdata2
covdata


cordata = cor(data)
results5 = prcomp(cordata,scale = TRUE)
fviz_eig(results5, addlabels= TRUE)

cordata2 = cor(normdat)
results6 = prcomp(cordata2,scale = TRUE)
fviz_eig(results6, addlabels= TRUE)

cordata
cordata2

V = matrix(,nrow = 209, ncol = 209)

for (i in 1:209) {
  for (j in 1:209) {
    if (i == j) {
      V[i,j] = 1/sqrt(covdata2[i,j])
      print(j)
      print(V[i,j])
    } else
      V[i,j] = 0
  }
}

V

mu = colMeans(data)
W= read.csv(file = "projectx_mu.csv")
x_mu = data.matrix(W)

eigenv = eigen(cordata)
e1 = data.matrix(eigenv$vectors)
e2 = as.matrix(e1[,1])
e3 = as.matrix(e1[,2])
e4 = as.matrix(e1[,3])

Z1 = t(e2)%*%V
Z2 = t(e3)%*%V
Z3 = t(e4)%*%V

PC1 = Z1%*%t(x_mu)
PC2 = Z2%*%t(x_mu)
PC3 = Z3%*%t(x_mu)

Q1 = as.vector(PC1)
Q2 = as.vector(PC2)
Q3 = as.vector(PC3)
PC = rbind(PC1,PC2,PC3)
PC
write.csv(PC,file = "PCs.csv")
attributes(PC)

P1 = qcc(Q1,type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
P11 = qcc(Q1,type = "xbar.one", confidence.level = 0.9991, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")

P2 = qcc(Q2,type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC2")
P21 = qcc(Q2,type = "xbar.one", confidence.level = 0.9991, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC2")

P3 = qcc(Q3, type = "xbar.one",nsigmas = 3,  label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC3")
P31 = qcc(Q3, type = "xbar.one",confidence.level = 0.9991,  label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC3")
#
#
#
#
#PHASE II
#Begins here...

X1 = P1$limits
a = c()
for (i in 1:552){
  if (Q1[i]<X1[1]){
    a = c(a,i)
  }
}

for (i in 1:552){
  if (Q1[i]>X1[2]){
    a = c(a,i)
  }
}

X2 = P2$limits
b = c()
for (i in 1:552){
  if (PC2[1,i]<X2[1]){
    b = c(b,i)
  }
}

for (i in 1:552){
  if (PC2[1,i]>X2[2]){
    b = c(b,i)
  }
}

X3 = P3$limits
d = c()
for (i in 1:552){
  if (PC3[1,i]<X3[1]){
    d = c(d,i)
  }
}

for (i in 1:552){
  if (PC3[1,i]>X3[2]){
    d = c(d,i)
  }
}


a
b
d

P11 = qcc(Q1[-a],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
attributes(P11)
a2 = P11$violations$beyond.limits
a1 = c(a,a2)
P12 = qcc(Q1[-a1],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a2 = P12$violations$beyond.limits
a3 = c(a1,a2)
P12 = qcc(Q1[-a3],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a4 = c(a3, P12$violations$beyond.limits)
P14 = qcc(Q1[-a4],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a5 = c(a4, P14$violations$beyond.limits)
P15 = qcc(Q1[-a5],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a6 = c(a5,P15$violations$beyond.limits)
P16 = qcc(Q1[-a6],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a7 = c(a6,P16$violations$beyond.limits)
P17 = qcc(Q1[-a7],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a8 = c(a7, P17$violations$beyond.limits)
a8
P18 = qcc(Q1[-a8],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a9 = c(a8, P18$violations$beyond.limits)
a9
P19 = qcc(Q1[-a9],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a10 = c(a9, P19$violations$beyond.limits)
a10
P110 = qcc(Q1[-a10],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a11 = c(a10, P110$violations$beyond.limits)
a11
P111 = qcc(Q1[-a11],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a12 = c(a11, P111$violations$beyond.limits)
a12
P112 = qcc(Q1[-a12],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a13 = c(a12, P112$violations$beyond.limits)
a13
P113 = qcc(Q1[-a13],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a14 = c(a13, P113$violations$beyond.limits)
a14
P114 = qcc(Q1[-a14],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a15 = c(a14, P114$violations$beyond.limits)
P115 = qcc(Q1[-a15],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a16 = c(a15, P115$violations$beyond.limits)
a16
P116 = qcc(Q1[-a16],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a17 = c(a16, P116$violations$beyond.limits)
a17
P117 = qcc(Q1[-a17],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
a18 = c(a17, P117$violations$beyond.limits)
a18
P118 = qcc(Q1[-a18],type = "xbar.one",nsigmas = 3, label.limits = c("LCL","UCL"),xlab = "Index", ylab = "PC1")
