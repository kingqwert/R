library(meta)
library(ggplot2)
library(tidyr)
mean =1 # c(0.5,1)
tau2 = 0.7 # c(0.01,0.10,0.50,0.70)

#あとはmean1、tau2=0.5と0.7

result=list()
theta = list()
for(r in 1:3){
	eps=c(3,5,10)[r]
	result_each = matrix(NA,ncol=length(2:50),nrow=1000)
	result_theta = matrix(NA,ncol=length(2:50),nrow=1000)
	for(S in c(2:50)){
		K = S*eps
		if(K %%1 ==0){
			for(j in 1:1000){
				si = 0.25*rchisq(K,1)
				si = ifelse(si <0.009,0.009,ifelse(si>0.6,0.6,si))
				
				yi = numeric(K)
				for(i in 1:K){
					yi[i] = rnorm(1,mean=mean,sd=(si[i]+tau2))
				}
				
				res = metagen(yi,sqrt(si),comb.fixed=F,comb.random=T)
				all_p = pnorm(-res$TE.random*sqrt(sum(res$w.random)))
				#all_p = pnorm(-res$TE.fixed*sqrt(sum(res$w.fixed)))
				tau = res$tau
				wi = 1/(si+tau)
				pi = pnorm(-sqrt(wi)*yi)
				
				sub_si = si[rank(-pi)<=S]
				sub_yi = yi[rank(-pi)<=S]
				
				sub_res = metagen(sub_yi,sqrt(sub_si),comb.fixed=F,comb.random=T)
				sub_p = pnorm(-sub_res$TE.random*sqrt(sum(sub_res$w.random)))
					
				theta_min = -sub_res$TE.fixed-qnorm(0.05)/sqrt(sum(sub_res$w.fixed))
				
				result_each[j,S-1]=ifelse(all_p > 0.05,NA, 
										ifelse( all_p < 0.05 & sub_p>0.05,1,0))
				result_theta[j,S-1]=theta_min
			}
		}
	}
	result[[r]]=result_each
	theta[[r]] = result_theta
}


std = function(x) sqrt(sum(x,na.rm=T)/sum(!is.na(x))*(1-sum(x,na.rm=T)/sum(!is.na(x))))
se_prop = function(x)  sqrt(sum(x,na.rm=T)/sum(!is.na(x))*(1-sum(x,na.rm=T)/sum(!is.na(x))))/sqrt(sum(!is.na(result[[1]][,1])))


x1 = rbind(apply(result[[1]],2,mean,na.rm=T),apply(result[[1]],2,se_prop))
x2 = rbind(apply(result[[2]],2,mean,na.rm=T),apply(result[[2]],2,se_prop))
x3 = rbind(apply(result[[3]],2,mean,na.rm=T),apply(result[[3]],2,se_prop))

df1 = cbind(gather(key = Group, value = Mean, as.data.frame(x1[1,])),gather(key = Group, value = SD, as.data.frame(x1[2,])))
df2 = cbind(gather(key = Group, value = Mean, as.data.frame(x2[1,])),gather(key = Group, value = SD, as.data.frame(x2[2,])))
df3 = cbind(gather(key = Group, value = Mean, as.data.frame(x3[1,])),gather(key = Group, value = SD, as.data.frame(x3[2,])))

df1$S = c(2:50)
df1$r = "r=1/3"
df2$S = c(2:50)
df2$r = "r=1/5"
df3$S = c(2:50)
df3$r = "r=1/10"

df = rbind(df1,df2,df3)
df=df[,-c(1,3)]
df = df[!is.na(df$Mean),]

d = df

d=d[d$S<=30,]
pd <- position_dodge(width = 0.2)
file = paste("/Users/daisuke/Dropbox/論文/Cherry picking/figs/Case2/Case2(theta",mean,"tau",tau2,").pdf",sep="")

g = ggplot(d, aes(x = S, y = Mean, color=r,group=r))
g <- g + geom_line(position = pd,size=1)+geom_point(size = 1, shape = 18, position = pd) 
#g = g + geom_line()+geom_point(size = 2, shape = 16)+ylim(c(0,1))
#g <- g + geom_errorbar(aes(ymin = ifelse(Mean - 1.96*SD<0,0,Mean - 1.96*SD), ymax = ifelse(Mean + 1.96*SD>1,1,Mean + SD), width = 0.3),position = pd)+ scale_x_continuous(expand = c(0, 0),limits=c(0,16),breaks=c(0,2,5,10,15))
g <- g + geom_errorbar(aes(ymin = ifelse(Mean - SD<0,0,Mean - SD), ymax = ifelse(Mean + SD>1,1,Mean + SD), width = 0.3),position = pd)+ scale_x_continuous(expand = c(0, 0),limits=c(0,32),breaks=c(0,2,5,10,15,20,25,30))+ylim(c(0,1))
g = g+theme_classic()+theme(text=element_text(size=30,  family="serif"),legend.position="none")+ labs(x="",y="") 
g
ggsave(file,units="cm",width=28,height=19)




theta_copy = theta
theta[[1]] = ifelse(result[[1]]==0,NA,theta[[1]])
theta[[2]] = ifelse(result[[2]]==0,NA,theta[[2]])
theta[[3]] = ifelse(result[[3]]==0,NA,theta[[3]])

x_theta1 = rbind(apply(theta[[1]],2,mean,na.rm=T),apply(theta[[1]],2,function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))))
x_theta2 = rbind(apply(theta[[2]],2,mean,na.rm=T),apply(theta[[2]],2,function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))))
x_theta3 = rbind(apply(theta[[3]],2,mean,na.rm=T),apply(theta[[3]],2,function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))))

df_theta1 = cbind(gather(key = Group, value = Mean, as.data.frame(x_theta1[1,])),gather(key = Group, value = SD, as.data.frame(x_theta1[2,])))
df_theta2 = cbind(gather(key = Group, value = Mean, as.data.frame(x_theta2[1,])),gather(key = Group, value = SD, as.data.frame(x_theta2[2,])))
df_theta3 = cbind(gather(key = Group, value = Mean, as.data.frame(x_theta3[1,])),gather(key = Group, value = SD, as.data.frame(x_theta3[2,])))

df_theta1$S = c(2:50)
df_theta1$r = "r=1/3"
df_theta2$S = c(2:50)
df_theta2$r = "r=1/5"
df_theta3$S = c(2:50)
df_theta3$r = "r=1/10"

df_theta = rbind(df_theta1,df_theta2,df_theta3)
df_theta=df_theta[,-c(1,3)]
df_theta = df_theta[!is.na(df_theta$Mean),]


df_theta=df_theta[df_theta$S<=30,]
df_theta$Mean = ifelse(df_theta$Mean<0,0,df_theta$Mean)
pd <- position_dodge(width = 0.2)

file = paste("/Users/daisuke/Dropbox/論文/Cherry picking/figs/Case2/Case2_thetamin(theta",mean,"tau",tau2,").pdf",sep="")
g = ggplot(df_theta, aes(x = S, y = Mean, color=r))
g <- g + geom_line(position = pd,size=1) +geom_point(size = 1, shape = 18, position = pd) 
#g = g + geom_line()+geom_point(size = 2, shape = 16)+ylim(c(0,1))
#g <- g + geom_errorbar(aes(ymin = ifelse(Mean - 1.96*SD<0,0,Mean - 1.96*SD), ymax = ifelse(Mean + 1.96*SD>1,1,Mean + SD), width = 0.3),position = pd)+ scale_x_continuous(expand = c(0, 0),limits=c(0,16),breaks=c(0,2,5,10,15))
g <- g + geom_errorbar(aes(ymin =ifelse(Mean - SD<0,0,Mean-SD), ymax = Mean + SD, width = 0.3),position = pd)+ scale_x_continuous(expand = c(0, 0),limits=c(0,32),breaks=c(0,2,5,10,15,20,25,30))+ scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))
g = g+theme_classic()+theme(text=element_text(size=30,  family="serif"),legend.position="none")+ labs(x="",y="")
g

ggsave(file,units="cm",width=28,height=19)

res = list(result = result, theta = theta_copy)
file = paste("/Users/daisuke/Dropbox/論文/Cherry picking/results/Case2_res(theta",mean,"tau",tau2,").RData",sep="")
save(res,file=file)

