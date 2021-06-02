rm(list=ls(all=TRUE))

#########################################################
########################### Connecting to MySQL database
install.packages("RMySQL")
library(RMySQL)

mysqlconnection = dbConnect(MySQL(), user = 'root', password = '*******', dbname='Ali1',
                            host = '127.0.0.1', port=3306)
class(result)

# List the tables available in this database.
dbListTables(mysqlconnection)
###########Querying the Tables
# Query the "students" tables to get all the rows.
result = dbSendQuery(mysqlconnection, "select * from Q1")

# Store the result in a R data frame object. n = 5 is used to fetch first 5 rows.
mydata.frame = fetch(result, n = 5)
print(mydata.frame)
#########Query with Filter Clause
result = dbSendQuery(mysqlconnection, "SELECT Name, History from students WHERE History>19;")

mydata.frame = fetch(result, n = 5)
print(mydata.frame)

# Use the R data frame "mtcars" to create the table in MySql.
# All the rows of mtcars are taken inot MySql.
dbWriteTable(mysqlconnection, "mtcars", mtcars[, ], overwrite = TRUE)
result = dbSendQuery(mysqlconnection, "select * from mtcars")

###########Dropping Tables in MySql
dbSendQuery(mysqlconnection, 'drop table if exists mtcars')


####################################################################
############### Portfolio Optimization
install.packages("quantmod")
library(quantmod)
getSymbols("JPM") #NYSE JP Morgan & Chase from yahoo 
getSymbols("MSFT")
getSymbols("F")
getSymbols("BAC")
getSymbols("AAPL")
getSymbols("PFE")
RateAnCoeff<-252

dat.ret=matrix(, nrow = 2584, ncol = 6)
dat.ret[,1]= diff(log(coredata(JPM$JPM.Close)))
dat.ret[,2]= diff(log(coredata(MSFT$MSFT.Close)))
dat.ret[,3]= diff(log(coredata(F$F.Close)))
dat.ret[,4]= diff(log(coredata(BAC$BAC.Close)))
dat.ret[,5]= diff(log(coredata(AAPL$AAPL.Close)))
dat.ret[,6]= diff(log(coredata(PFE$PFE.Close)))

risk.param <- 0.5
Dmat <- cov(dat.ret)
dvec <- colMeans(dat.ret) * risk.param
# Constraints 1: sum(x_i) = 1
Amat <- matrix(1, nrow=nrow(Dmat))
bvec <- 1
meq <- 1

install.packages('quadprog')
library(quadprog)
#solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
qp$solution
w=qp$solution
mio=w%*%colMeans(dat.ret)*RateAnCoeff
mySD<- sqrt(t(w)%*%(Dmat%*%(w))*RateAnCoeff)

dat.ret2=matrix(, nrow = 2579, ncol = 6)
dat.ret2[,1]= diff(log(coredata(JPM$JPM.Close)))
dat.ret2[,2]= diff(log(coredata(MSFT$MSFT.Close)))
dat.ret2[,3]= diff(log(coredata(F$F.Close)))
dat.ret2[,4]= diff(log(coredata(BAC$BAC.Close)))
dat.ret2[,5]= diff(log(coredata(AAPL$AAPL.Close)))
dat.ret2[,6]= diff(log(coredata(PFE$PFE.Close)))

Dmat <- cov(dat.ret2)
Amat <- matrix(1, nrow=nrow(Dmat))
bvec <- 1
meq <- 1

risk.param<-.2
dvec <- colMeans(dat.ret2) * risk.param
qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
qp$solution
w=qp$solution
mio=w%*%colMeans(dat.ret2)*RateAnCoeff
mySD<- sqrt(t(w)%*%(Dmat%*%(w))*RateAnCoeff)
par(new=TRUE)
plot(mySD, mio, main="Efficient Frontier", 
     xlab="portfolio variance", ylab="portfolio mean",
     xlim=c(0.0, .4), ylim=c(-0.3, .25),col = "red", pch =16)

allmio= c()
allmySD<-c()
for (risk.param in seq(from=0, to=.7, by=.01)){
  print(risk.param)
  dvec <- colMeans(dat.ret2) * risk.param
  qp <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  w=qp$solution
  mio=w%*%colMeans(dat.ret2)*RateAnCoeff
  mySD<- sqrt(t(w)%*%(Dmat%*%(w))*RateAnCoeff)
  allmio <- c(allmio, mio)
  allmySD<- c(allmySD, mySD)
}

par(new=TRUE)

#plot(allmySD, allmio, main="Efficient Frontier", 
#      xlab="portfolio variance", ylab="portfolio mean", col="green")

plot(allmySD, allmio, main="Efficient Frontier", 
     xlab="portfolio variance", ylab="portfolio mean",
     xlim=c(0.0, .4), ylim=c(-0.3, .25), col="blue",pch =16)
lines(seq(0.0, .4, by=.1),rep(0, 5),  col = "red")

plot(allmySD, allmio, main="Efficient Frontier", 
     xlab="portfolio variance", ylab="portfolio mean",
     xlim=c(0.0, .4), ylim=c(-0.3, .25))

# Constraints 2: sum(x_i) = 1 & x_i >= 0
Amat <- cbind(1, diag(nrow(Dmat)))
bvec <- c(1, rep(0, nrow(Dmat)))
meq <- 1

# Constraints: sum(x_i) = 1 & x_i >= 0 & x_i <= 0.15
Amat <- cbind(1, diag(nrow(Dmat)), -1*diag(nrow(Dmat)))
bvec <- c(1, rep(0, nrow(Dmat)), rep(-0.15, nrow(Dmat)))
meq <- 1

###############################################
############ VaR calculation and other tests
#VaR: 1) Historical approach
inv= 10000
miosamples=as.vector(w%*%t(dat.ret))
TomorrowPortV= inv*(miosamples+1)
L<- length(TomorrowPortV)
hist(TomorrowPortV)
quantile(TomorrowPortV,.99)
quantile(TomorrowPortV,.01)

mean(TomorrowPortV)

#VaR: 2) Model Building approach
mymio=(w%*%colMeans(dat.ret))
mySD<- sqrt(t(w)%*%(Dmat%*%(w)))
qnorm(.01, mean= (1+mymio)*inv, sd= mySD*inv)

#############
##Using a built-in function (VaR)
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
myVaR= (VaR(miosamples, p=.01, method="gaussian")+1)*inv
hist((miosamples+1)*inv)
plot(miosamples)
lines(miosamples)
myCVaR= (CVaR(miosamples, p=.01, method="gaussian")+1)*inv

##############
#Another way:
port.price<- na.omit(merge(JPM$JPM.Close, MSFT$MSFT.Close, F$F.Close, BAC$BAC.Close, AAPL$AAPL.Close, PFE$PFE.Close))
port.return<- ROC(port.price, type="continuous")[-1,]
VaR(port.return , p=.01,weights= w, portfolio_method="component", method="gaussian")

##################
#####T-test
T1<- read.csv('/Users/administrator/Desktop/QF/QF104-032717/QFpanel.csv',header = FALSE)
T2<- read.csv('/Users/administrator/Desktop/QF/QF104-032717/QFpanel2.csv',header = FALSE)
#TT1<-T1[-1,]
TT1=as.matrix(T1)
TT2=as.matrix(T2)
t.test(TT1,TT2)
t.test(TT1,TT2,paired=TRUE)
t.test(TT1,mu=0) # Ho: mu=0
