# calculate cost function
cost = function(y,design,logic,dist,other.regs,lambda1,lambda2,alpha){
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