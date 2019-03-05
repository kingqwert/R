# calculate log likelihood
logL2 = function(y,L){
      reg = try(speedglm(y~-1+L,family=binomial(link=logit)))
      if(inherits(reg,'try-error')){reg = try(glm(y~-1+L,family=binomial(link=logit)))}

      return(list(coefs= coef(reg),logL = logLik(reg)[1])) #this log-likelihood is 
}