# make a distant matrix with gene sets
dist.cal.gene = function(P,snploc,gene.snp.set,gene.gene.inte,e,f,g,lambda){
     dist.mat = as.matrix(dist(snploc, function(x,y){
          e+f*exp(-(x-y)^2/g^2)
          }))

     gene.snp.inte = matrix(NA,P,P)
     
     for(snp1 in 1:P){
          for(snp2 in 1:P){
               gene1 = gene.snp.set[snp1]
               gene2 = gene.snp.set[snp2] 
               gene.snp.inte[snp1,snp2] = gene.gene.inte[gene1,gene2]
          } 
     }

     for(s1 in 1:P){
          for(s2 in 1:P){
               g1 = gene.snp.set[s1]
               g2 = gene.snp.set[s2]
               if(g1!=g2){
                    dist.mat[s1,s2] = 0
               }
          }
     }

     dist.mat+lambda*gene.snp.inte
}
