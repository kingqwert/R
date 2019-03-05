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