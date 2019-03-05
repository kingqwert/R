library(MASS)
library(mvtnorm)
library(speedglm)
library(proxy)

source("dist.cal.gene.R")
source("sw.move.R")
source("bdk.R")
source("bdl.R")
source("change.R")
source("change.model.R")
source("create.matrix.R")
source("logL2.R")
source("cost.R")
source("accept.R")
source("make.init.trees.R")



################################################################################################################################################
# make simulation dataset
################################################################################################################################################
K.true = 3
true.logic = list()

# make design matrix
P = 2000 # number of covs
n = 10000 # number of samples
p.cov = runif(P) # parameter for binomial dist of covs
design = rbinom(n, 1, p.cov[1])

for(i in 2:P){design = cbind(design, rbinom(n, 1, p.cov[i]))} # design matrix

# make true distant matrix
e.true = 0
f.true = 1
g.true = 1000
margin = sample(1:100, P-1, replace = T)
snploc = c(1,cumsum(margin))

# make true gene-gene-interaction
num.genes = 10
gene.snp.set = sort(sample(1:num.genes, size = P, replace=T)) # make SNP-Gene relationship table
gene.gene.inte = matrix(NA,num.genes,num.genes)
gene.gene.inte[upper.tri(gene.gene.inte,diag=T)]=c(rep(0,8),1,rep(0,8),1,1,rep(0,36)) 
gene.gene.inte[lower.tri(gene.gene.inte,diag=T)]=c(rep(0,20),1,0,1,rep(0,6),1,rep(0,25)) 

# create true distance matrix
dist.true = dist.cal.gene(P,snploc, gene.snp.set, gene.gene.inte, e.true, f.true, g.true, lambda=0) # true distance matrix. 0 means there are no correlation between genes

# make true logics
# pick up main SNPs
mainsnp.gene3 = sample(which(gene.snp.set==3)[1]:tail(which(gene.snp.set==3),1), 1) # decide main snp randomly from gene3
mainsnp.gene4 = sample(which(gene.snp.set==4)[1]:tail(which(gene.snp.set==4),1), 1) # decide main snp randomly from gene4
mainsnp.gene6 = sample(which(gene.snp.set==6)[1]:tail(which(gene.snp.set==6),1), 1) # decide main snp randomly from gene4

# pick up SNPs located in neiborhood
cov.set3 = sample(which(gene.snp.set==3)[1]:tail(which(gene.snp.set==3),1), 
     4, 
     prob=dist.true[mainsnp.gene3,which(gene.snp.set==3)[1]:tail(which(gene.snp.set==3),1)])
cov.set4 = sample(which(gene.snp.set==4)[1]:tail(which(gene.snp.set==4),1), 
     3, 
     prob=dist.true[mainsnp.gene4,which(gene.snp.set==4)[1]:tail(which(gene.snp.set==4),1)])
cov.set6 = sample(which(gene.snp.set==6)[1]:tail(which(gene.snp.set==6),1), 
     3, 
     prob=dist.true[mainsnp.gene6,which(gene.snp.set==6)[1]:tail(which(gene.snp.set==6),1)])

# sample operators
op.set3 = sample(c("&","|"), length(cov.set3)-1, replace=T)
comp.set3 = sample(c("","c"), length(cov.set3), replace=T)
true.logic[[1]] = list(substr(paste(paste0(paste0(paste0("X",cov.set3),comp.set3),op.set3),collapse=""), 1, nchar(paste(paste0(paste0(paste0("X",cov.set3),comp.set3),op.set3),collapse=""))-1),
     lk = length(cov.set3),
     covs = cov.set3,
     ops = op.set3,
     comp = comp.set3)
op.set4 = sample(c("&","|"), length(cov.set4)-1, replace=T)
comp.set4 = sample(c("","c"), length(cov.set4), replace=T)
true.logic[[2]] = list(substr(paste(paste0(paste0(paste0("X",cov.set4),comp.set4),op.set4),collapse=""), 1, nchar(paste(paste0(paste0(paste0("X",cov.set4),comp.set4),op.set4),collapse=""))-1),
     lk = length(cov.set4),
     covs = cov.set4,
     ops = op.set4,
     comp = comp.set4)
op.set6 = sample(c("&","|"), length(cov.set6)-1, replace=T)
comp.set6 = sample(c("","c"), length(cov.set6), replace=T)
true.logic[[3]] = list(substr(paste(paste0(paste0(paste0("X",cov.set6),comp.set6),op.set6),collapse=""), 1, nchar(paste(paste0(paste0(paste0("X",cov.set6),comp.set6),op.set6),collapse=""))-1),
     lk = length(cov.set6),
     covs = cov.set6,
     ops = op.set6,
     comp = comp.set6)


# make true parameter b and outcome y
b = c(-1,1,1,1)
Li.true = create.matrix(design=design, logic=true.logic)
y = rbinom(n, 1, prob=1/(1+exp(-Li.true%*%b)))




#################################################
# Simulations
#################################################

###################
# True sets of covs and distance matrix
###################
colnames(design) = paste("V", c(1:P),sep="")
dist.tmp=dist.true
diag(dist.tmp)=1
dist.tmp = dist.tmp+0.0000000000001

# set the number of iteration, number of trees, update timing and tempterture of simulated annealing
niter = 10000
num.trees = 5
update = 500
temper = 10^seq(to=-1,from=1,by=-1/25)

# set tuning parameters
lambda1=1
lambda2=1
lambda3=1

# sample 100 from 10000 population 
id = sample(1:10000, 100, replace=F)
design.sample = design[id,]
y.sample = y[id] 

####################################################################################
# Main part
# proposed method by using spatial logic
####################################################################################
# make initial trees
init.trees = make.init.trees(P=P,K.init=3,num.trees=num.trees,min.lk.init=1,max.lk.init=5)

acceptprob = list()
trees = init.trees
update.thres = update * num.trees

for (k in 1:length(temper)){
     show(k)
     update.count = 0

     for(s in 1:niter){
          show(s)
          for(t in 1:num.trees){
               if(k ==1){
                    logic.old = init.trees[[t]]
               }else{logic.old = trees[[t]]}

               model = sw.move(P=P, logic.old=logic.old,max.num.tree=5,max.num.leaf=5,dist=dist.true)
               logic.new = change.model(model=model,logic=logic.old)
               other.regs = trees[-t]
     
               acceptprob = accept(y.sample,P,dist=dist.tmp,design.sample,model,logic.old,logic.new,other.regs,lambda1=lambda1,lambda2=lambda2,alpha=lambda3,temp=temper[k])

               if(runif(1)<acceptprob){
                    update.count = update.count+1
                    trees[[t]] = logic.new
               }else{
                    trees[[t]] = logic.old   
               }
     
          }
          if(update.count > update.thres) break
     }  
}

# select one optimal tree that has highest likelihood 
opt.tree.num = which.max(sapply(trees,function(x){logL2(y.sample,L=create.matrix(design.sample,x))$logL}))
opt.tree = trees[[opt.tree.num]]

# get regression coefficients
reg = glm(y.sample~-1+create.matrix(design.sample, opt.tree),family=binomial(link=logit))


