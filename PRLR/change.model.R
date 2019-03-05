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