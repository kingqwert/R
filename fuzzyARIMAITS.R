library(astsa)
library(forecast)
library(dplyr)
library(zoo)
library(xts)
library(ggpmisc)
library(scales)
library(openxlsx)
library(tseries)
library(tsbox)
library(truncnorm)
library(tidyr)

create.data.normal = function(intervention, data, l1,l2, mean=0, sd=1){
	Ft = dtruncnorm(c(l1:l2), a=l1, b=l2, mean = mean, sd = sd)
	tmp = data.frame(date = seq(as.Date(intervention)+l1,as.Date(intervention)+l2,by='day'), Ft)	
	res = merge(data,tmp,by="date",all.x=T)
	res$Ft[is.na(res$Ft)] = 0
	
	# to create F_{t-1} - F_{t-l1-l2}
	for(i in 1:(-l1+l2)){
		tmptmp = data.frame(date = seq(as.Date(intervention)+l1+i,as.Date(intervention)+l2,by='day'),
							F_tmp = Ft[1:(-l1+l2+1-i)]
				)
		names(tmptmp) = c('date',paste('F_',i,sep=""))
		res = merge(res,tmptmp,by="date",all.x=T)
	}
	res[is.na(res)] = 0
	return(res)
}

create.data.unif = function(intervention, data, l1,l2){
	Ft = rep(1/(l2-l1+1),l2-l1+1)
	tmp = data.frame(date = seq(as.Date(intervention)+l1,as.Date(intervention)+l2,by='day'), Ft)	
	res = merge(data,tmp,by="date",all.x=T)
	res$Ft[is.na(res$Ft)] = 0
	
	# to create F_{t-1} - F_{t-l1-l2}
	for(i in 1:(-l1+l2)){
		tmptmp = data.frame(date = seq(as.Date(intervention)+l1+i,as.Date(intervention)+l2,by='day'),
							F_tmp = Ft[1:(-l1+l2+1-i)]
				)
		names(tmptmp) = c('date',paste('F_',i,sep=""))
		res = merge(res,tmptmp,by="date",all.x=T)
	}
	res[is.na(res)] = 0
	return(res)
}

df3.norm = create.data.normal('2020-04-07',df2,l1=-6,l2=6,mean=0,sd=1)


# create C matrix
createCmat_poly = function(d, l1, l2){
	intercept =rep(1,(-l1+l2))

	for(i in 1:d){
		intercept = cbind(
			intercept,
			c(1:(-l1+l2))^i
			)
	}
	rbind(
		c(1,rep(0,d)),
		intercept
	)
}

createCmat_MA_L2 = function(l1, l2){
	mat = matrix(c(1,0,0,1),ncol=2)
	for(i in 1:(-l1+l2-1)){
		mat = rbind(
			mat,
			c((2^i-(-1)^i)/3/2^i,1-(2^i-(-1)^i)/3/2^i)
			)
	}
	return(mat)
}

Cmat1 = createCmat_poly(3,l1=-6,l2=6)
X.norm = as.matrix(df3.norm[,c(3:15)]) %*% Cmat1

# Use automated algorithm to identify p/q parameters
# Specify first difference = 0 and seasonal difference = 0

X2.norm = data.frame(X.norm, 
	positive = test[c(6:95),"testedPositive"]/test[c(6:95),"peopleTested"],
	deaths = test[c(6:95),"deaths"],
	avetemp = temp[temp$date>'2020-02-19' & temp$date<'2020-05-20','avetemp'])

X2.norm[is.na(X2.norm$deaths),]$deaths = 0


auto.arima(ts.data.diff1, seasonal = FALSE,
	xreg = as.matrix(X2.norm), 
	start.p=1, max.p = 1,
	max.d=0, max.D=0, method="ML",
	stepwise=FALSE, trace=TRUE, ic="aic"
	)

res.norm = Arima(ts.data.diff1, xreg=as.matrix(X2.norm), order=c(1,0,0), method="ML")

# Check residuals
checkresiduals(res.norm)
Box.test(res.norm$residuals, lag = 6, type = "Ljung-Box")


# counterfactual
X3.norm = X2.norm
X3.norm[,c(1:4)] = 0
res.fc.norm = Arima(ts.data.diff1, xreg=as.matrix(X3.norm),model=res.norm)
res.fc.norm.ts = ts(as.numeric(res.fc.norm$fitted[42:90]), start=c(2020,91), frequency=365)

