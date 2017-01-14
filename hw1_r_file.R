### Machine Learning Hw1
### Mengchuan Fu

library(MASS)   # use for the inverse matrix function "ginv"

### Problem 1
polynomial <- function(x,y,d){
  ytrain<-y[trainset]
  xtrain<-x[trainset]
  ytest<-y[-trainset]
  xtest<-x[-trainset]
  n <-nrow(data)*0.5
  msetrain=0
  msetest=0
  x<-rep(1,n)
  for(n in 1:d){
    x<-cbind(x,xtrain^n)
  }
  A<-ginv(t(x)%*%x)%*%t(x)%*%ytrain
  for(i in 1:n){
    msetrain <- msetrain + (ytrain[i]-(A[1]+A[2]*xtrain[i]))^2
    msetest <- msetest + (ytest[i]-(A[1]+A[2]*xtest[i]))^2
  }
  return(list(c(msetrain,msetest),A))
}

mse_d <- function(){
  mse_train<-rep(NA,20)
  mse_test<-rep(NA,20)
  for(i in 1:20){
    mse_train[i] <-polynomial(data[,1],data[,2],i)[[1]][1]
    mse_test[i] <- polynomial(data[,1],data[,2],i)[[1]][2]
  }
  return(list(mse_train,mse_test))
}

# read data
data <- read.csv("hw1a.csv") 
# use cross-validation by randomly splitting the data in to two halves
trainset <- sample(nrow(data),nrow(data)*0.5)
d <-which.min(mse_test)

par(mfrow=c(3,1))
plot(1:20,mse_d()[[1]],type="l",main="Training MSE",xlab="Polynomial order D",ylab="MSE")
plot(1:20,mse_d()[[2]],type="l",main="Testing MSE",xlab="Polynomial order D",ylab="MSE")
plot(data[,1],data[,2],add=T,main="Function overlaid on the data",xlab="x",ylab="y")
abline(a=polynomial(data[,1],data[,2],1)[[2]][1],b=polynomial(data[,1],data[,2],1)[[2]][2],col="red")



### Problem 2
data1 <- read.csv("hw1b.csv")

ridge <- function(lambda,x,y,trainset){
  xtrain <- x[trainset,]
  xtest <- x[-trainset,]
  ytrain <- y[trainset]
  ytest <- y[-trainset]
  
  beta <- solve( t(xtrain) %*% xtrain +  diag(lambda,nrow=ncol(x)), toll=0 ) %*% t(xtrain) %*% ytrain
  error <- y - x %*% beta
  losstrain <- sum(error[trainset]^2)
  losstest <- sum(error[-trainset]^2)
  return ( list(losstrain,losstest))
}

result <- sapply(0:1000, ridge,x = as.matrix(data1[,2:101]),y = as.matrix(data1[,1]),
                 trainset = sample(nrow(data1),nrow(data1)*0.5))
trainloss <- result[1,]
testloss <- result[2,]
par(mfrow=c(2,1))
plot(0:1000,trainloss,main="Training Error",xlab="lambda",ylab="Error")
plot(0:1000,testloss,main="Testing Error",xlab="lambda",ylab="Error")
min(unlist(testloss))
which.min(unlist(testloss))
points(719,24814.25,pch=8,col="red",lwd="10")


### Problem 4
loss_function <- function(x,y,theta){
  n <-length(y)
  h <- 1/(1+exp(-x%*%theta))
  loss<-(1/n)*((t(y)-1)%*%log(1-h)-t(y)%*%log(h))
  return(loss)
}
  
gradient <- function(x,y,theta){
  h <- 1/(1+exp(-x%*%theta))
  gradient <- t(x)%*%(h-y)/length(y)
  return(gradient)
}
  
gra_des <- function(x,y,epsilon,alpha){
  n <-nrow(data)
  theta = as.matrix(rep(0,3))
  old_theta = theta+2*epsilon
  loss_set<-c()
  while (norm(theta-old_theta,"2") > epsilon){
    old_theta<- theta
    theta <- theta - alpha*gradient(x,y,theta)
    loss_set <- c(loss_set,loss_function(x,y,theta))
    print(norm(gradient(x,y,theta)),"2")
  }
  return(list(theta,loss_set))
}

my_plot <- function(data,result){
  par(mfrow=c(2,1))
  x<-data[,1:2]
  pass<-which(data[,4]==1)
  plot(x[pass,1],x[pass,2],type ="p",col="blue",pch=13,xlab="X1",ylab="X2",ylim=c(-1,1),xlim=c(0,1),
       main = "Plot of Decision Boundary")
  lines(x[-pass,1],x[-pass,2],type ="p",col="green",pch=17,xlab="X1",ylab="X2",ylim=c(-1,1),xlim=c(0,1))
  abline(-result[[1]][3]/result[[1]][2],-result[[1]][1]/result[[1]][2],col="red")
  loss<-result[[2]]
  plot(1:length(loss),loss,xlab="Iteration times",ylab="Loss",main="Rate of Convergence")
}
  
data <- read.csv("hw1c.csv")  
trainset <- sample(nrow(data),nrow(data)*0.5) 
x<-as.matrix(data[,1:3])   
y<-as.matrix(data[,4])
  
result <- gra_des(as.matrix(data[-4])[trainset,],as.matrix(data[4])[trainset],epsilon=0.01,alpha=0.5)
result <- gra_des(as.matrix(data[-4])[trainset,],as.matrix(data[4])[trainset],epsilon=0.001,alpha=0.5)
my_plot(data,result)
print(result[1])



### Problem 5
gradient <- function(x,y,theta){
  h <- 1/(1+exp(-x%*%theta))
  gradient <- t(x)%*%(h-y)/length(y)
  return(gradient)
}

hessian <- function(x,y,theta){
  h <- 1/(1+exp(-x%*%theta))
  hes <- t(x)%*% diag(as.numeric(h)) %*% diag(as.numeric(t(1-h)))%*% x
  return(hes)
}

loss_function <- function(x,y,theta){
  n <-length(y)
  h <- 1/(1+exp(-x%*%theta))
  loss<-(1/n)*((t(y)-1)%*%log(1-h)-t(y)%*%log(h))
  return(loss)
}

newton <- function(x,y,epsilon){
  n <-nrow(data)
  theta = as.matrix(rep(0,3))
  old_theta = theta+2*epsilon
  loss_set<-c()
  while (norm(gradient(x,y,theta),"2") > epsilon){
    old_theta<- theta
    theta <- theta - ginv(hessian(x,y,theta))%*%gradient(x,y,theta)
    loss_set <- c(loss_set,loss_function(x,y,theta))
    print(norm(gradient(x,y,theta)),"2")
  }
  return(list(theta,loss_set))
}

my_plot <- function(data,result){
  par(mfrow=c(2,1))
  x<-data[,1:2]
  pass<-which(data[,4]==1)
  plot(x[pass,1],x[pass,2],type ="p",col="blue",pch=13,xlab="X1",ylab="X2",ylim=c(-1,1),xlim=c(0,1),
       main = "Plot of Decision Boundary")
  lines(x[-pass,1],x[-pass,2],type ="p",col="green",pch=17,xlab="X1",ylab="X2",ylim=c(-1,1),xlim=c(0,1))
  abline(-result[[1]][3]/result[[1]][2],-result[[1]][1]/result[[1]][2],col="red")
  loss<-result[[2]]
  plot(1:length(loss),loss,xlab="Iteration times",ylab="Loss",main="Rate of Convergence")
}

data <- read.csv("hw1c.csv")  
trainset <- sample(nrow(data),nrow(data)*0.5) 
x<-as.matrix(data[,1:3])   
y<-as.matrix(data[,4])

result <- newton(as.matrix(data[-4])[trainset,],as.matrix(data[4])[trainset],epsilon=0.01)
result <- newton(as.matrix(data[-4])[trainset,],as.matrix(data[4])[trainset],epsilon=0.001)
my_plot(data,result)
print(result[1])
