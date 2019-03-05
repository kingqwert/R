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
