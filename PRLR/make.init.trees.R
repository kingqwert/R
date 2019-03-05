# prepare initial trees and logics
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
