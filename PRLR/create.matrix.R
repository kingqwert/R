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