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
