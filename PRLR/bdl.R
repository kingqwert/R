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