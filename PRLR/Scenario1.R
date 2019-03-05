#########################################
#This code is for Spatial logic regression with L2 penarty.
#Author: Daisuke Yoneoka (blue.sky.sea.dy@gmail.com)
#First modification: July 13, 2016
#Last modification: July 13, 2016
#########################################

#Model
#Logit(P(Y=1)) = b0+b1*L_1+b2*L_2+...+bK*L_K

#Definitions
#i = samples (i = 1,...n)
#j = covariates (j = 1,...p)
#K = number of trees (k = 1,...K)
#l_k = number of leaves within tree k (l = 1,...lk)
#L_k = tree with boolean logics (= (PO_(1,k),PO_(2,k),...PO_(l_k,k)))
#PO_(l,k) = Position within one tree (= (X_(l,k),O_(l,k)))
#X_(l,k) = covariate at position l within tree k
#O(l,k) = operator at position l within tree k
#bk = coefficients

#Parameters 
#bk ~ N(b,v^2)
#b ~ Laplace (a,b)  # By incorporating laplace dist, we can get sparse solution (= lasso L2). Hyper parameters a and b should be fixed previously.
#K ~ geometric(c)  #Hyper parameter c should be fixed previously
#l_k ~ U(1,d) #Hyper parameter e should be fixed previously. In general, d allows very big number (i.e. allows many leaves)
#T={t1: Birth or Death of one tree | t2: Birth or Death of one leaf | t3: Change cov or operator}
#kt = number of trees at iteration t
#lkt = number of leaves within tree k at iteration t
#s1 =P(BDK = Birth| T=t1) : probability of birth of one tree when T=t1.
#e,f,g ~ Uniform(min=e1(or f1,g1),max=e2(or f2,g2))


library(proxy)
library(smoothmest) 
library(MASS)
library(mvtnorm)
library(VGAM)
library(LogicReg)
library(glmnet)
library(sac)
library(scrime)
library(speedglm)
library(ROCR)

#library(foreach)
#library(doParallel)

#registerDoParallel(detectCores())
#cl <- makeCluster(10) 
#registerDoParallel(cl)


################################################################################################################################################
# make functions
################################################################################################################################################
# make a distant matrix with gene sets
dist.cal.gene = function(P,snploc,gene.snp.set,gene.gene.inte,e,f,g,lambda){
     dist.mat = as.matrix(dist(snploc, function(x,y){
          e+f*exp(-(x-y)^2/g^2)
          }))

     gene.snp.inte = matrix(NA,P,P)
     
     for(snp1 in 1:P){
          for(snp2 in 1:P){
               gene1 = gene.snp.set[snp1]
               gene2 = gene.snp.set[snp2] 
               gene.snp.inte[snp1,snp2] = gene.gene.inte[gene1,gene2]
          } 
     }

     for(s1 in 1:P){
          for(s2 in 1:P){
               g1 = gene.snp.set[s1]
               g2 = gene.snp.set[s2]
               if(g1!=g2){
                    dist.mat[s1,s2] = 0
               }
          }
     }

     dist.mat+lambda*gene.snp.inte
}



# switch move type
sw.move = function(P, logic.old,max.num.tree,max.num.leaf,dist){ 
     
     u = runif(1)
     if (u < 1/3) 
         move.type = 1
     else if (u < 2/3) 
         move.type = 2
     else move.type = 3
     
     switch(move.type, 
          bdk(P,logic.old,max.num.tree,max.num.leaf,dist), 
          bdl(P,logic.old,max.num.leaf,dist), 
          change(logic.old)
     )  
}


# functions for different moves in trees
bdk = function(P, logic.old,max.num.tree,max.num.leaf,dist){

     u = runif(1)
     if (u < 0.5) 
        bdk.move.type = 1
     else bdk.move.type = 2
     
     lk = sample(1:max.num.leaf, 1)
     comp = sample(c("c",""), lk, replace =T)
     
     if(lk>1){
          op = sample(c("&","|"), lk-1, replace =T)    
     }else{
          op = c("&") 
     }

     new.x.main = sample(1:P, 1)
     new.x = c(new.x.main, sample((1:P)[-new.x.main], lk-1,prob=dist[new.x.main,-new.x.main]))

     if(bdk.move.type==2 & length(logic.old)==1){bdk.move.type=1}
     if(bdk.move.type==1 & length(logic.old)==max.num.tree){bdk.move.type=2}
     
     switch(bdk.move.type, 
          return(list(ftype=1, new.x=new.x, new.comp=comp, new.op=op)), #this is f1
          return(list(ftype=2, k.minus=sample(1:length(logic.old),1))) #this is f2
     )
}

# functions for different moves in leaves
bdl = function(P,logic.old,max.num.leaf,dist){

     u = runif(1)
     if (u < 0.5){bdl.move.type = 1}else{bdl.move.type = 2}
          
     k1 = sample(1:length(logic.old),1) #which tree should be modified
     covs.in.leaf = logic.old[[k1]]$covs #extract the information about covariates included in tree k1  
     
     if(length(covs.in.leaf)==1){
          xj = sample((1:P)[-covs.in.leaf], 1, prob=dist[-covs.in.leaf,covs.in.leaf])
     }else{
          xj = sample((1:P)[-covs.in.leaf], 1, prob=apply(dist[-covs.in.leaf,covs.in.leaf],1,sum))  #pick up one covariate from design matrix without covariates, which already included in tree k1          
     }

     lbirth = sample(0:logic.old[[k1]]$lk,1) #pick up one position within L_k = (PO_(1,k),PO_(2,k),...PO_(l_k,k)) 
     ldeath = sample(1:logic.old[[k1]]$lk, 1) #pick up one position within L_k = (PO_(1,k),PO_(2,k),...PO_(l_k,k)) 
     
     comp = sample(c("c",""), 1, replace =T)
     op = sample(c("&","|"), 1, replace =T)
     
     if(bdl.move.type==2 & length(covs.in.leaf)==1 & length(logic.old)==1){bdl.move.type = 1}
     if(bdl.move.type==1 & length(covs.in.leaf)>=max.num.leaf){bdl.move.type=2}

     switch(bdl.move.type, 
          return(list(ftype=3, k.bdl=k1, new.x=xj, position=lbirth, new.comp=comp, new.op=op)), #this is f3
          return(list(ftype=4, k.bdl=k1, position=ldeath)) #this is f4
     )    
}

# functions for different moves in change
change = function(logic.old){
     u = runif(1)
     if (u < 0.5) 
          change.move.type = 1
     else change.move.type = 2
     
     k1 = sample(1:length(logic.old),1) #which tree should be modified
     l.cov = sample(1:logic.old[[k1]]$lk,1) #pick up one covariate in tree k1 to change
     l.op = sample(1:(logic.old[[k1]]$lk-1),1) #pick up one operator in tree k1 to change

     if(l.op ==0) {l.op=1}
     
     switch(change.move.type, 
          return(list(ftype=5, change='cov', k.change=k1, position=l.cov)), #this is f5
          return(list(ftype=6, change='operator', k.change=k1, position=l.op)) #this is f6
     )  
}


# update model representation
change.model = function(model,logic) {
     if(model$ftype==1){
     
          new.tree = list(logic = substr(paste(paste0(paste0(paste0("X",model$new.x),model$new.comp),model$new.op),collapse=""),
               1,
               nchar(paste(paste0(paste0(paste0("X",model$new.x),model$new.comp),model$new.op),collapse=""))-1),
               lk = length(model$new.x),
               covs = model$new.x,
               ops = model$new.op,
               comp = model$new.comp, 
               new = list(ftype=1,new.x=model$new.x,new.comp=model$new.comp,new.op=model$new.op)
          )
          logic = c(logic,list(new.tree))

     }else if(model$ftype==2){
     
          logic[[model$k.minus]] = NULL

     }else if(model$ftype==3){
     
          covs = append(logic[[model$k.bdl]]$covs, model$new.x, after = model$position)
          comp = append(logic[[model$k.bdl]]$comp, model$new.comp, after = model$position)
          
          if(is.na(logic[[model$k.bdl]]$ops[1])){
               ops = model$new.op
          }else{
               ops = append(logic[[model$k.bdl]]$ops, model$new.op, after = model$position)
          }
          
          logic[[model$k.bdl]] = list(logic = substr(paste(paste0(paste0(paste0("X",covs),comp),ops),collapse=""), 
              1, 
              nchar(paste(paste0(paste0(paste0("X",covs),comp),ops),collapse=""))-1),
              lk = logic[[model$k.bdl]]$lk+1,
              covs = covs,
              ops = ops,
              comp = comp,
              new = list(ftype=3,k.bdl=model$k.bdl, nex.x=model$new.x, position=model$position)
          )    

     }else if(model$ftype==4){
     
          if(logic[[model$k.bdl]]$lk==1){
               logic[[model$k.bdl]] = NULL  # when the tree with one leaf, remove the tree               
          }else{
               covs = logic[[model$k.bdl]]$covs[-model$position]
               if(length(logic[[model$k.bdl]]$ops)==1){
                ops = NA
               }else{
                ops = logic[[model$k.bdl]]$ops[-model$position]
          }
               comp = logic[[model$k.bdl]]$comp[-model$position]
          
               logic[[model$k.bdl]] = list(logic = substr(paste(paste0(paste0(paste0("X",covs),comp),ops),collapse=""), 
                  1, 
                  nchar(paste(paste0(paste0(paste0("X",covs),comp),ops),collapse=""))-1),
                  lk = logic[[model$k.bdl]]$lk-1,
                  covs = covs,
                  ops = ops,
                  comp = comp,
                  new = list(ftype=4,k.dbl=model$k.bdl, position=model$position)
               )
          }

     }else if(model$ftype==5){
     
          if(logic[[model$k.change]]$comp[model$position]=="c"){
               logic[[model$k.change]]$comp[model$position] = ""
          }else{
               logic[[model$k.change]]$comp[model$position] = "c"
          }
          comp = logic[[model$k.change]]$comp
     
          logic[[model$k.change]] = list(logic = substr(paste(paste0(paste0(paste0("X",logic[[model$k.change]]$covs),comp),logic[[model$k.change]]$ops),collapse=""), 
               1,
               nchar(paste(paste0(paste0(paste0("X",logic[[model$k.change]]$covs),comp),logic[[model$k.change]]$ops),collapse=""))-1),
               lk = logic[[model$k.change]]$lk,
               covs = logic[[model$k.change]]$covs,
               ops = logic[[model$k.change]]$ops,
               comp = comp,
               new = list(ftype=5,k.change=model$k.change, position=model$position)
               )
     
     }else if(model$ftype==6){
     
          if(is.na(logic[[model$k.change]]$ops[1])){
            logic[[model$k.change]]$ops =NA
          }else if(logic[[model$k.change]]$ops[model$position]=="&"){
               logic[[model$k.change]]$ops[model$position] = "|"
          }else{
               logic[[model$k.change]]$ops[model$position] = "&"
          }
          ops = logic[[model$k.change]]$ops
     
          logic[[model$k.change]] = list(logic = substr(paste(paste0(paste0(paste0("X",logic[[model$k.change]]$covs),logic[[model$k.change]]$comp),ops),collapse=""), 
               1, 
               nchar(paste(paste0(paste0(paste0("X",logic[[model$k.change]]$covs),logic[[model$k.change]]$comp),ops),collapse=""))-1),
               lk = logic[[model$k.change]]$lk,
               covs = logic[[model$k.change]]$covs,
               ops = ops,
               comp = logic[[model$k.change]]$comp,
               new = list(ftype=6,k.change=model$k.change, position=model$position)
               )

    }else{print("move type error!")}

    return(logic)
}

# update model matrix
create.matrix = function(design,logic){
     Li = matrix(NA, ncol=length(logic), nrow=dim(design)[1])

     for(l in 1:length(logic)){
          L = design[,logic[[l]]$covs]
          
          if(is.vector(L)){
               Li[,l] = L
          }else{
               for (ll in 1:dim(L)[2]){
                    if(logic[[l]]$comp[ll]=="c"){L[,ll]=(1-L[,ll])} # transform design matrix according to the result of moves for compliments.
               }
     
               Li[,l] = L[,1]
               for (ll in 2:dim(L)[2]){
                    if(logic[[l]]$ops[ll-1]=="&"){  # main part of calculation of logic
                         Li[,l] = as.numeric(Li[,l]&L[,ll]) 
                    }else if(logic[[l]]$ops[ll-1]=="|"){
                         Li[,l] = as.numeric(Li[,l]|L[,ll])
                    }else{print("operator error! check create.matrix()")}  
               }    
          }         
     }
     return(cbind(1,Li))
}

#################
# calculate log likelihood
#################
#logL = function(y,L){
#     reg = glm(y~-1+L,family=binomial)
#     
#     return(list(coefs= coef(reg),logL = logLik(reg)[1])) #this log-likelihood is 
#}

logL2 = function(y,L){
      reg = try(speedglm(y~-1+L,family=binomial(link=logit)))
      if(inherits(reg,'try-error')){reg = try(glm(y~-1+L,family=binomial(link=logit)))}

      return(list(coefs= coef(reg),logL = logLik(reg)[1])) #this log-likelihood is 
}

#################
# calculate cost function
#################
cost = function(y,design,logic,dist,other.regs,lambda1,lambda2,lambda3){
  this.covs = lapply(logic,function(y){y$covs})

  -logL2(y, L=create.matrix(design, logic))$logL+lambda1*sim.A(dist=dist,logic=logic,other.regs=other.regs)+lambda2*sim.B(dist=dist,logic=logic,this.covs=this.covs)-lambda3*sim.C(dist=dist,logic=logic,this.covs=this.covs)
}

cost2 = function(y,design,logic,dist,other.regs,lambda1,lambda2){
     this.covs = lapply(logic,function(y){y$covs})

     -logL2(y, L=create.matrix(design, logic))$logL+lambda1*log(sim.A(dist=dist,logic=logic,other.regs=other.regs))+lambda2*log(sim.D(dist=dist,logic=logic,this.covs=this.covs))
}

cost3 = function(y,design,logic,dist,other.regs,lambda1,lambda2,alpha){
     this.covs = lapply(logic,function(y){y$covs})

     -logL2(y, L=create.matrix(design, logic))$logL+lambda1*log(sim.A(dist=dist,logic=logic,other.regs=other.regs))+lambda2*log(sim.B(dist=dist,logic=logic,this.covs=this.covs))+alpha*mean(sapply(this.covs,length))
}

sim.A = function(dist,logic,other.regs){
  all.covs = sort(unique(unlist(lapply(logic,function(y){y$covs}))))
  covs.list=lapply(other.regs,function(x){sort(unique(unlist(lapply(x,function(y){y$covs}))))}) 
  
  mean(sapply(1:length(covs.list),
    function(x){
      mean(dist[all.covs,covs.list[[x]]])
      }))
}

sim.B = function(dist,logic,this.covs){
  mean(apply(expand.grid(1:length(this.covs),1:length(this.covs)),1,
    function(x){
      mean(dist[this.covs[[x[1]]],this.covs[[x[2]]]])
      }))
}

sim.C = function(dist,logic,this.covs){
  mean(sapply(1:length(this.covs),
    function(x){
      mean(dist[this.covs[[x]],this.covs[[x]]])
      }))
}

sim.D = function(dist,logic,this.covs){
     mean(apply(expand.grid(1:length(this.covs),1:length(this.covs)),1,
     function(x){
          mean(dist[this.covs[[x[1]]],this.covs[[x[2]]]])
          }))/mean(sapply(1:length(this.covs),
          function(x){
               mean(dist[this.covs[[x]],this.covs[[x]]])
               }))
}


#################
# Acceptance rate with temperature
#################
accept = function(y,P,dist,design,model,logic.old,logic.new,other.regs,lambda1,lambda2,lambda3,temp){ # acceptance prob for change in length of theta (and beta) and lk (i.e. length of L_k)
     cost.old = cost(y,design,logic=logic.old,dist,other.regs,lambda1,lambda2,lambda3)
     cost.new = cost(y,design,logic=logic.new,dist,other.regs,lambda1,lambda2,lambda3)

     return(min(1,exp(-(cost.new-cost.old)/temp)))
}

accept2 = function(y,P,dist,design,model,logic.old,logic.new,other.regs,lambda1,lambda2,temp){ # acceptance prob for change in length of theta (and beta) and lk (i.e. length of L_k)
     cost.old = cost2(y,design,logic=logic.old,dist,other.regs,lambda1,lambda2)
     cost.new = cost2(y,design,logic=logic.new,dist,other.regs,lambda1,lambda2)

     return(min(1,exp(-(cost.new-cost.old)/temp)))
}

accept3 = function(y,P,dist,design,model,logic.old,logic.new,other.regs,lambda1,lambda2,alpha,temp){ # acceptance prob for change in length of theta (and beta) and lk (i.e. length of L_k)
     cost.old = cost3(y,design,logic=logic.old,dist,other.regs,lambda1,lambda2,alpha)
     cost.new = cost3(y,design,logic=logic.new,dist,other.regs,lambda1,lambda2,alpha)

     return(min(1,exp(-(cost.new-cost.old)/temp)))
}

################
# prepare initial trees and logics
################
make.init.trees = function(P,K.init,num.trees,min.lk.init,max.lk.init){
  init.trees = list()

  for(i in 1:num.trees){
    lk.init = sample(min.lk.init:max.lk.init, K.init,replace=T) #number of leaves in each trees

      mainsnp.init = sample(1:P, K.init) # decide main snp randomly
      logic.init = list()
      
      for (t in 1:(K.init)){
          cov.set = c(mainsnp.init[t],sample(1:P, lk.init[t]-1)) # pick up some snps with correlation
          op.set = sample(c("&","|"), lk.init[t]-1, replace=T)
          if(lk.init[t]==1){op.set = "&"}
          comp.set = sample(c("","c"), lk.init[t], replace=T)
          logic.init[[t]] = list(substr(paste(paste0(paste0(paste0("X",cov.set),comp.set),op.set),collapse=""), 1, nchar(paste(paste0(paste0(paste0("X",cov.set),comp.set),op.set),collapse=""))-1),
              lk = lk.init[t],
              covs = cov.set,
              ops = op.set,
              comp = comp.set,
              new = NULL)
      }

      init.trees[[i]] = logic.init
  }
  return(init.trees)
}

# performance measure for predictions
perf = function(cut, fit, y){
   yhat = (fit>cut)
   w = which(y==1)
   sen = mean( yhat[w] == 1 ) 
   spe = mean( yhat[-w] == 0 ) 
   return(c(sen,spe))
}





################################################################################################################################################
# make simulation dataset
################################################################################################################################################
set.seed(12345)

# set ups for true parameters of trees and leaves
#K.true = rposgeom(n=1, prob=0.3) #number of trees K=1 in this dataset
#K.true = 2

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

num.genes = 10
gene.snp.set = sort(sample(1:num.genes, size = P, replace=T)) # make SNP-Gene relationship table
gene.gene.inte = matrix(NA,num.genes,num.genes)
gene.gene.inte[upper.tri(gene.gene.inte,diag=T)]=c(rep(0,8),1,rep(0,9),0,rep(0,7),0,rep(0,28)) 
gene.gene.inte[lower.tri(gene.gene.inte,diag=T)]=c(rep(0,20),1,rep(0,34)) 


dist.true = dist.cal.gene(P,snploc, gene.snp.set, gene.gene.inte, e.true, f.true, g.true, lambda=0) # true distance matrix. 0 means there are no correlation between genes
#dist.true = dist.true/apply(dist.true,1,sum)





# make true parameter b and outcome y
b = c(0)
y = rbinom(n, 1, prob=1/(1+exp(b)))



#################################################
# Simulations
#################################################

###################
# True sets of covs
###################
#true = sort(unique(unlist(lapply(true.logic,function(x){x$covs}))))
#true.comb = unlist(lapply(true.logic,FUN=function(x){outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=",")[lower.tri(outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=","))]}))
#true.comb.forFBLR = unlist(lapply(true.logic,FUN=function(x){outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=",")[upper.tri(outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=","))]}))
#true.logL = logL2(y,L=create.matrix(design,true.logic))$logL

colnames(design) = paste("V", c(1:P),sep="")

# store results 
selected.lasso = list()
selected.logic = list()
selected.comb = list()
selected.pro = list()
selected.FBLR = list()

# pars related to lasso or logic reg
log.pars = list()
logic.fit = list()

myanneal = logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 500)


# set the number of iteration and burn-in in mcmc
niter = 10000
dist.tmp=dist.true
diag(dist.tmp)=1
dist.tmp = dist.tmp+0.0000000000001

res = matrix(NA,ncol=7,nrow=10)

args = commandArgs(T)
seed = as.numeric(args[1])
samples = as.numeric(args[2])
set.seed(seed)

for(m in 1:10){
 # sample 500 persons from 1000 population 
     id = sample(1:10000, samples, replace=F)
     design.sample = design[id,]
     y.sample = y[id] 
     
     ####################################################################################
     #proposed method by using spatial logic
     num.trees =5
     update = 500
     init.trees = make.init.trees(P=P,K.init=3,num.trees=num.trees,min.lk.init=1,max.lk.init=5)
     
     acceptprob = list()
     trees = init.trees
     temper = 10^seq(to=-1,from=1,by=-1/25)
     update.thres = update * num.trees
     
     for (k in 1:length(temper)){
          show(k)
          update.count = 0

          for(s in 1:niter){
               for(t in 1:num.trees){
                    if(k ==1){
                         logic.old = init.trees[[t]]
                    }else{logic.old = trees[[t]]}
     
                    model = sw.move(P=P, logic.old=logic.old,max.num.tree=3,max.num.leaf=5,dist=dist.true)
                    logic.new = change.model(model=model,logic=logic.old)
                    other.regs = trees[-t]
          
                    acceptprob = accept3(y.sample,P,dist=dist.tmp,design.sample,model,logic.old,logic.new,other.regs,lambda1=1,lambda2=0.01,alpha=10,temp=temper[k])

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
     
     opt.tree.num = which.max(sapply(trees,function(x){logL2(y.sample,L=create.matrix(design.sample,x))$logL}))
     opt.tree = trees[[opt.tree.num]]
     selected.pro = sort(unique(unlist(lapply(opt.tree,function(x){x$covs}))))
     selected.comb = unlist(lapply(opt.tree,FUN=function(x){outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=",")[lower.tri(outer(as.character(sort(unlist(x$covs))),as.character(sort(unlist(x$covs))),FUN=paste ,sep=","))]}))
     reg_tmp = glm(y.sample~-1+create.matrix(design.sample, opt.tree),family=binomial(link=logit))
     auc.pro = as.numeric(performance(prediction(predict(reg_tmp,type=c("response")), y.sample),"auc")@y.values)
     senspe.pro = perf(0.5,predict(reg_tmp,type=c("response")),y.sample)

     # Performance evaluation     
     proportion1.pro = length(selected.pro)
     proportion1.pro = ifelse(is.na(proportion1.pro),0,proportion1.pro)
          
     proportion2.pro = length(selected.pro)
     proportion2.pro = ifelse(is.na(proportion2.pro),0,proportion2.pro)

     

     # two-SNP interactions
     # Combination w/o ops
     
     res1.pro.comb = length(unlist(selected.comb))
     res1.pro.comb = ifelse(is.na(res1.pro.comb),0,res1.pro.comb)
     res2.pro.comb = length(unlist(selected.comb))
     res2.pro.comb = ifelse(is.na(res2.pro.comb),0,res2.pro.comb)


     res[m,] = c(proportion1.pro,
        proportion2.pro,
        res1.pro.comb,
        res2.pro.comb,
        auc.pro,
        senspe.pro)
    
     file = paste("~/result/parallel/Add_result_parallel_Scenario1_",samples,"_",seed,".csv",sep="")
     write.csv(res, file = file, append=T)

}

