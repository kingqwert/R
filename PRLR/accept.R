# Acceptance rate with temperature
accept = function(y,P,dist,design,model,logic.old,logic.new,other.regs,lambda1,lambda2,alpha,temp){ # acceptance prob for change in length of theta (and beta) and lk (i.e. length of L_k)
     cost.old = cost(y,design,logic=logic.old,dist,other.regs,lambda1,lambda2,alpha)
     cost.new = cost(y,design,logic=logic.new,dist,other.regs,lambda1,lambda2,alpha)

     return(min(1,exp(-(cost.new-cost.old)/temp)))
}
