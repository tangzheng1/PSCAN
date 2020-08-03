library(MASS)
library(SKAT) 
library(igraph)
library(rgl)


# type: "mean" for PSCAN-M, "variance" for PSCAN-V
# p.comb: use "Cauchy" or "minP" method to combine p-values across windows and obtain global (gene-level) p-value 
# U: a vector of score statistics with names as SNP ID in chr:position format
# V: a covariance matrix
# MAC: a vector of minor allele counts
# weight: a vector contains the weight for each SNP (default = NULL, flat weight)
# N.MC: number of permutation for signal detection
# FWER: family-wise error rate in signal detection procedure (default = 0.05). If set to NULL, will not perform signal detection step. 
# f: a cutoff for the fraction of overlap in signal detection procedure (default f=0, assume signal windows are not overlapping with each other) 
# details: by default, only gene-level p-value and the SNP ID of the signal variants are output. If details=TRUE, then the p-value for each signal regions and the coordinates of the SNPs on the protein space will also be output
PSCAN <- function(type="mean", p.comb="Cauchy", U, V, MAC, weight=NULL, N.MC=1000, FWER=0.05, f=0, details=FALSE, plot3D=FALSE){

  if(!is.vector(U)){
    warning("score statistic is not a vector")
    U = as.vector(U) 
  }
  
  m = length(U)
  
  SNPID = names(U)
  if(is.null(SNPID)){
    stop("chr:position SNP ID should be the name of the score statistic vector")  
  }
  
  tmp = strsplit(SNPID,":")
  chr = as.numeric( tmp[[1]][1] )
  pos = sapply(tmp, function(x) as.numeric(x[2])) 
  
  if( sum(is.na(chr)) + sum(is.na(pos)) > 0  ){
    stop("Invalid SNP ID") 
  }
  
  
  if(is.null(weight)){
    weight = rep(1,m)
  }
  
  if(!is.vector(weight)){
    warning("weight is not a vector")
    weight = as.vector(weight) 
  }
  
  if(!is.matrix(V)){
    warning("Covariance is not a matrix")
    V = as.matrix(V)  
  }
  
  if(!isSymmetric.matrix(V)){
    stop("Covariance matrix is not symmetric") 
  }
  # if( det(V)<=0 ){
  #   stop("Covariance matrix is not positive-definite") 
  # }
  
  if(m != length(weight) | m != nrow(V)){
    stop("Dimensions among U, V, and/or weight do not agree")
  }
  
  ####################################
  #
  #  Map SNPs to protein space. First, use experimental structure, then use computational structure
  #
  ####################################
  structure.type = NA
  flag.mapro = 0
  index.pos = NULL
  
  for(i in 1:m){
    
    pos.i = pos[i]
    index.tmp = which.min( abs(PSCAN:::exp.start[[chr]] - pos.i) )
    diff = pos.i-PSCAN:::exp.start[[chr]][index.tmp] 
    
    #print(i)
    if( abs(diff)<3 ){
      
      if( diff>=0 ){
        
        index.pos = c(index.pos, index.tmp)
        
      }else{
        
        if(index.tmp==1){
          
          index.pos = c(index.pos, NA)
          
        }else{
          
          diff2 = pos.i-PSCAN:::exp.start[[chr]][index.tmp-1] 
          if(diff2>=0 & diff2<3){
            index.pos = c(index.pos, index.tmp-1)
          }else{
            index.pos = c(index.pos, NA)
          }
          
        }
        
      }
      
      #print(index.tmp)
    }else{
      index.pos = c(index.pos, NA)
    }
    
    
  }
  
  flag.mapro = as.numeric(sum(!is.na(index.pos))>0)
  
  if( flag.mapro ){
    
    structure.type = "PDB"
    idx = which(!is.na(index.pos))
    SNPID.mapro = SNPID[idx]
    
    idx.mapro = index.pos[idx]
    COORD = PSCAN:::exp.coord[[chr]][idx.mapro,,drop=FALSE]
    rownames(COORD) = SNPID.mapro
    pdbid.mapro = PSCAN:::exp.pdbid[[chr]][idx.mapro]
    pdbid.mapro.list = unique(pdbid.mapro)
    n.pdbid = length(pdbid.mapro.list)
    
    
  }else{
    
    flag.mapro.mod = 0
    
    index.pos = NULL
    
    for(i in 1:m){
      
      pos.i = pos[i]
      index.tmp = which.min( abs(PSCAN:::mod.start[[chr]] - pos.i) )
      diff = pos.i-PSCAN:::mod.start[[chr]][index.tmp] 
      
      #print(i)
      if( abs(diff)<3 ){
        
        if( diff>=0 ){
          
          index.pos = c(index.pos, index.tmp)
          
        }else{
          
          if(index.tmp==1){
            
            index.pos = c(index.pos, NA)
            
          }else{
            
            diff2 = pos.i-PSCAN:::mod.start[[chr]][index.tmp-1] 
            if(diff2>=0 & diff2<3){
              index.pos = c(index.pos, index.tmp-1)
            }else{
              index.pos = c(index.pos, NA)
            }
            
          }
          
        }
        
        #print(index.tmp)
      }else{
        index.pos = c(index.pos, NA)
      }
      
      
    }
    
    flag.mapro.mod = as.numeric(sum(!is.na(index.pos))>0)
    
    if( flag.mapro.mod ){
      
      structure.type = "Modbase"
      idx = which(!is.na(index.pos))
      SNPID.mapro = SNPID[idx]
      
      idx.mapro = index.pos[idx]
      COORD = PSCAN:::mod.coord[[chr]][idx.mapro,,drop=FALSE]
      rownames(COORD) = SNPID.mapro
      pdbid.mapro = PSCAN:::mod.pdbid[[chr]][idx.mapro]
      pdbid.mapro.list = unique(pdbid.mapro)
      n.pdbid = length(pdbid.mapro.list)
      
      flag.mapro = 1
      
    }
    
  }  
  
  
  ####################################
  #
  #  Perform Tests
  #
  ####################################   
  PSCAN.PVAL = NA
  PSCAN.SIGNAL = NA
  SIGNAL.PVAL = NA
  n.signal.region = NULL
  COORD.ls = NULL
  
  if( flag.mapro ){
    
    if(length(idx) == length(U)){
      SNPID.imapro = NULL 
    }else{
      SNPID.imapro = SNPID[-idx]
    }
    
    SNPID.ls = list(); COORD.ls = list();
    for(k in 1:n.pdbid){
      
      PDBID = pdbid.mapro.list[k]
      
      idx.pdb = which(pdbid.mapro == PDBID)
      
      SNPID.ls = c(SNPID.ls, list(SNPID.mapro[idx.pdb]))
      COORD.ls = c(COORD.ls, list(COORD[idx.pdb, , drop=FALSE]))
      
    }
    
    if(!is.null(SNPID.imapro)){
      SNPID.ls = c(SNPID.ls, list(SNPID.imapro) )
    }
    
    SNPID.all = unlist(SNPID.ls)
    reord.idx = match( SNPID.all, SNPID )
    U.ls = U[reord.idx]
    V.ls = V[reord.idx, reord.idx, drop=FALSE]
    MAC.ls = MAC[reord.idx]  
    weight.ls = weight[reord.idx]
    

    
    #### PSCAN ####
    if(type=="mean"){
      
      tmp = try( PSCAN.M(p.comb=p.comb, U.ls, V.ls, MAC.ls, weight=weight.ls, COORD.ls, N.MC=N.MC, FWER=FWER, f=f) )
      
      if(!inherits(tmp, 'try-error') ){
      
        PSCAN.PVAL = tmp$Q.scan2.pval
        PSCAN.SIGNAL = list()
        n.signal.region = length(tmp$signal.region)
        if(n.signal.region>0){
          for(i in 1:n.signal.region){
            PSCAN.SIGNAL[[i]] = SNPID.all[tmp$signal.region[[i]]]
          }
        }
        SIGNAL.PVAL = tmp$signal.pval
        
      }
      

    }
    
    if(type=="variance"){
      
      tmp = try( PSCAN.V(p.comb=p.comb, U.ls, V.ls, MAC.ls, weight=weight.ls, COORD.ls, N.MC=N.MC, FWER=FWER, f=f) )
      
      if(!inherits(tmp, 'try-error') ){
        
        PSCAN.PVAL = tmp$Q.scan3.pval
        PSCAN.SIGNAL = list()
        n.signal.region = length(tmp$signal.region)
        if(n.signal.region>0){
          for(i in 1:n.signal.region){
            PSCAN.SIGNAL[[i]] = SNPID.all[tmp$signal.region[[i]]]
          }
        }
        SIGNAL.PVAL = tmp$signal.pval
        
      }
      
    
    }
    
    
  }
  else{# perform burden/SKAT test if protein structure is not availiable
    
    if(type=="mean"){
      tmp.U = weight %*% U
      tmp.V = weight %*% V %*% weight  
      PSCAN.PVAL = pchisq(tmp.U*tmp.U/tmp.V,1, lower.tail=F)
    }
    
    if(type=="variance"){
      PSCAN.PVAL = Q.pval(U, V, weight=weight, r=0)
    }
    
    
    PSCAN.SIGNAL = NA
    SIGNAL.PVAL = NA
    
  }
  
  ####################################
  #
  #  3D plot for signal regions in protein space
  #
  ####################################  
  if(n.signal.region!=0 & !is.null(n.signal.region)){
    
    j.dom = which(sapply(SNPID.ls, function(x) sum(!is.na(match(unlist(PSCAN.SIGNAL), x)))) > 0)[1]  # if signal region spread over multiple structures, then the first structure will be plot
    if (j.dom > n.pdbid) {
      
      structure.type=NA  #Signal region(s) only contain unmapped SNPs
    }
    
  }

  
  if (plot3D & !is.null(n.signal.region)) {
    
    
    if (n.signal.region > 0) {
      
      if (j.dom > n.pdbid) {
        print("Signal region(s) only contain unmapped SNPs. 3D plot not generated.")
      }
      else {
        COORD = COORD.ls[[j.dom]]
        m = nrow(COORD)
        start = 0
        if (j.dom > 1) {
          for (j in 1:(j.dom - 1)) {
            start = start + nrow(COORD.ls[[j]])
          }
        }
        idx = (start + 1):(start + m)
        signal.idx = match(unlist(PSCAN.SIGNAL), rownames(COORD))
        signal.idx = signal.idx[!is.na(signal.idx)]
        col.list = rep("gray", m)
        col.list[signal.idx] = "deeppink"
        
        open3d()
        rgl.bringtotop()
        plot3d(COORD[, 1], COORD[, 2], COORD[, 3], col = col.list, 
               type = "s", radius = 1, xlab = "", 
               ylab = "", zlab = "", box = TRUE, axes = FALSE)
        rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = "white")
      }
    }
    else {
      print("No signal region detected. 3D plot not generated.")
    }
  }
  
  if(length(COORD.ls)>0){
    names(COORD.ls) = pdbid.mapro.list
  }
  
  if(is.null(FWER)){
    return(list(pscan.pval = PSCAN.PVAL ))
    
  }else{
    if(details){
      return(list(pscan.pval = PSCAN.PVAL, signal = PSCAN.SIGNAL, signal.pval = SIGNAL.PVAL, structure.type=structure.type, structure.coord = COORD.ls ))
      
    }else{
      return(list(pscan.pval = PSCAN.PVAL, signal = PSCAN.SIGNAL ))
      
    }
  }
  
}


ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    Pvals[which(Pvals==1)] = 0.99  # v1.1 07/30/2020 update. Cauchy combined p-value becomes 1 if any p-value in Pvals is exactly 1 
    #warning("There are p-values that are exactly 1!")
    #return(1)
  }
  
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  
  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}

# U: a vector of score statistics
# V: a covariance matrix
# MAC: a vector of minor allele counts
# weight: a vector contains the weight for each SNP (default = NULL, flat weight)
# COORD.ls: a lists with each component being a matrix contins variant x-, y-, and z- coordinates in a fragment of protein structure
# N.MC: number of simulation for signal detection (default = 1000)
# FWER: family-wise error rate in signal detection procedue (default = 0.05). If set to NULL, will not perform signal detection step. 
# f: a cutoff for the fraction of overlap in signal detection procedure (default f=0, assume signal windows are not overlapping with each other) 

## Components in U, V, weight are already ordered according to COORD.ls, (mapped SNPs appears first, then the unmapped ones)

## test for mean
PSCAN.M <- function(p.comb="Cauchy", U, V, MAC, weight=NULL, COORD.ls, N.MC=1000, FWER=0.05, f=0){
  
  MAC.LB = 10
  
  if(sum(MAC) < MAC.LB){
    
    Q.max = NA; pval= NA; signal.region = NA; signal.pval = NA; C = NA;
    
  }else{
    
    m = length(U)
    n.dom = length(COORD.ls)
    m.ls = sapply(COORD.ls, nrow)
    m.mapro = sum(m.ls)
    m.imapro = m - m.mapro
    
    flag.search = as.numeric(m.ls>1)
    
    
    cluster.list = NULL
    
    m.idx = 0
    for(d in 1:n.dom){
      
      if(flag.search[d]){
        
        cords.dist = as.matrix(dist(COORD.ls[[d]]))
        dist.uniq = sort(unique(  as.numeric(cords.dist[lower.tri(cords.dist)])  ))
        n.poss = length(dist.uniq)
        cluster.list.d = as.list(c(1:m.ls[d])) 
        cluster.pre = cluster.list.d
        cluster.num.pre = m.ls[d]
        
        for(j in 1:n.poss){
          tmp = clusters(graph_from_adjacency_matrix(0+(cords.dist<=dist.uniq[j]), mode="undirected"))
          cluster.num.cur = tmp$no
          
          if(cluster.num.cur == cluster.num.pre){
            next
          }else{
            cluster.cur = split(1:m.ls[d], tmp$membership)
            cluster.add = cluster.cur[!(cluster.cur %in% cluster.pre)]
            cluster.list.d = c(cluster.list.d, cluster.add)
            
            cluster.pre = cluster.cur
            cluster.num.pre = cluster.num.cur
          }
          if(cluster.num.cur ==1){
            break
          }
        }
        
        for(j in 1:length(cluster.list.d)){
          cluster.list.d[[j]] = m.idx+cluster.list.d[[j]]
        }
        
        cluster.list = c(cluster.list, cluster.list.d)
        
      }else{
        
        cluster.list = c(cluster.list, list(m.idx + 1))
      }
      
      m.idx = m.idx + m.ls[d]
      
    }
    
    
    if( m.imapro ){
      
      cluster.list = c(cluster.list, list( (1:m)[-c(1:m.mapro)] ) )
      
    }
    
    
    n.scan = length(cluster.list)
    C = matrix(0, nrow=n.scan, ncol=m)
    for(j in 1:n.scan){
      C[j, cluster.list[[j]]] = 1
    }
    
    if( length( which(rowSums(C)==m) )==0 ){
      C = rbind(C, rep(1, m))
      n.scan = n.scan + 1
      cluster.list = c(cluster.list, list(rep(1:m)) )
    }
    
    # remove window if its MAC<10
    C.MAC = C %*% MAC
    rm.idx = which(C.MAC < MAC.LB)
    if(length(rm.idx)>0){
      C = C[-rm.idx,,drop=FALSE]
      n.scan = n.scan - length(rm.idx)
      cluster.list = cluster.list[-rm.idx]
    }
    
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    UC = C %*% U
    VC = C %*% V %*% t(C)
    Q = UC*UC / diag(VC)
    Q.max =  max(Q)
    
    Qp = pchisq(Q, 1, lower.tail=F)
    
    if(p.comb=="Cauchy"){
      pval = ACAT( Qp )
    }
    

    signal.region = NULL
    signal.pval = NULL
    if(!is.null(FWER) | p.comb=="minP"){
      
      set.seed(11)
      if(is.null(N.MC)){
        N.MC = 50/FWER
      }
      Q.max.perm = NULL
      for(l in 1:N.MC){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        tmp = C %*% U.perm
        Q.max.perm = c(Q.max.perm, max( tmp*tmp/diag(VC) ) )
      }
      
      if(p.comb=="minP"){
        pval = sum(Q.max.perm >= Q.max)/N.MC
      }
      
      if(!is.null(FWER)){
        
        signal.thresh = quantile(Q.max.perm, 1-FWER)
        cand.idx = which(Q >= signal.thresh)
        
        if( length(cand.idx)>0 ){
          
          Q.max.idx = which(Q==Q.max)
          if( length(Q.max.idx)>1 ){
            Q.max.idx = Q.max.idx[which.min(sapply(cluster.list[Q.max.idx], length))]
          }
          signal.region[[1]] = cluster.list[[Q.max.idx]]
          signal.pval[1] = pchisq(Q[Q.max.idx],1, lower.tail=F)
          signal.cand = cluster.list[cand.idx] 
          Q.cand = Q[cand.idx]
          signal.cand.max.region =  cluster.list[[Q.max.idx]]
          #nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
          if(f==0){
            lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
          }else{
            lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) )/length(x) < f) )
          }
          
          while( length(lowoverlap.idx)>0  ){
            signal.cand = signal.cand[lowoverlap.idx] 
            Q.cand = Q.cand[lowoverlap.idx]
            Q.max.idx = which(Q.cand==max(Q.cand))
            if( length(Q.max.idx)>1 ){
              Q.max.idx = Q.max.idx[which.min(sapply(signal.cand[Q.max.idx], length))]
            }
            signal.region = c(signal.region , signal.cand[Q.max.idx])
            signal.pval = c(signal.pval, pchisq(Q.cand[Q.max.idx],1, lower.tail=F) )
            signal.cand.max.region = signal.cand[[Q.max.idx]]
            #nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
            if(f==0){
              lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) ) == 0) )
            }else{
              lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.max.region) )/length(x) < f) )
            }
            
          }
          
        }        
      }

    }
    
    
  }

  
  return(list(Q.scan2.stat = Q.max, Q.scan2.pval=pval, signal.region = signal.region, signal.pval = signal.pval, C = C ))
  
}



Q.pval <- function(U, V, weight=NULL, r = 0){
  
  U =  U * weight
  V =  t(t(V * weight) * weight)
  pval = Met_SKAT_Get_Pvalue(U, V, r, method="davies")$p.value
  return(pval)
  
}


## test for variance
PSCAN.V <- function(p.comb="Cauchy", U, V, MAC, weight=NULL, COORD.ls, N.MC=1000, FWER = 0.05, f=0){
  
  MAC.LB = 10
  if(sum(MAC) < MAC.LB){  
    
    Qp.min = NA; pval= NA; signal.region = NA; signal.pval = NA; C = NA;
   
  }else{
    
    m = length(U)
    n.dom = length(COORD.ls)
    m.ls = sapply(COORD.ls, nrow)
    m.mapro = sum(m.ls)
    m.imapro = m - m.mapro
    
    flag.search = as.numeric(m.ls>1)
    
    
    cluster.list = NULL
    
    m.idx = 0
    for(d in 1:n.dom){
      
      if(flag.search[d]){
        
        cords.dist = as.matrix(dist(COORD.ls[[d]]))
        dist.uniq = sort(unique(  as.numeric(cords.dist[lower.tri(cords.dist)])  ))
        n.poss = length(dist.uniq)
        cluster.list.d = as.list(c(1:m.ls[d])) 
        cluster.pre = cluster.list.d
        cluster.num.pre = m.ls[d]
        
        for(j in 1:n.poss){
          tmp = clusters(graph_from_adjacency_matrix(0+(cords.dist<=dist.uniq[j]), mode="undirected"))
          cluster.num.cur = tmp$no
          
          if(cluster.num.cur == cluster.num.pre){
            next
          }else{
            cluster.cur = split(1:m.ls[d], tmp$membership)
            cluster.add = cluster.cur[!(cluster.cur %in% cluster.pre)]
            cluster.list.d = c(cluster.list.d, cluster.add)
            
            cluster.pre = cluster.cur
            cluster.num.pre = cluster.num.cur
          }
          if(cluster.num.cur ==1){
            break
          }
        }
        
        for(j in 1:length(cluster.list.d)){
          cluster.list.d[[j]] = m.idx+cluster.list.d[[j]]
        }
        
        cluster.list = c(cluster.list, cluster.list.d)
        
      }else{
        
        cluster.list = c(cluster.list, list(m.idx + 1))
      }
      
      m.idx = m.idx + m.ls[d]
      
    }
    
    
    if( m.imapro ){
      
      cluster.list = c(cluster.list, list( (1:m)[-c(1:m.mapro)] ) )
      
    }
    
    n.scan = length(cluster.list)
    C = matrix(0, nrow=n.scan, ncol=m)
    for(j in 1:n.scan){
      C[j, cluster.list[[j]]] = 1
    }
    
    if( length( which(rowSums(C)==m) )==0 ){
      C = rbind(C, rep(1, m))
      n.scan = n.scan + 1
      cluster.list = c(cluster.list, list(rep(1:m)) )
    }
    
    
    #  remove window if its MAC<10
    C.MAC = C %*% MAC
    rm.idx = which(C.MAC < MAC.LB)
    if(length(rm.idx)>0){
      C = C[-rm.idx,,drop=FALSE]
      n.scan = n.scan - length(rm.idx)
      cluster.list = cluster.list[-rm.idx]
    }
    
    if(!is.null(weight)){
      
      C = t(t(C) * weight)
      
    }
    
    # observed stat
    Qp = NULL
    for(j in 1:n.scan){
      idx = which(C[j,]!=0)
      Qp = c(Qp, Q.pval(U[idx], V[idx,idx], weight=C[j,idx], r = 0) )
    }
    Qp[Qp<0] = 0  
    Qp.min = min(Qp)
    
    if(p.comb=="Cauchy"){
      
      pval = ACAT( Qp )
    
    }
    
    
    
    signal.region = NULL
    signal.pval = NULL
    if(!is.null(FWER) | p.comb=="minP"){
      
      set.seed(11)
      if(is.null(N.MC)){
        N.MC = 50/FWER
      }
      Qp.min.perm = NULL
      for(l in 1:N.MC){
        U.perm = mvrnorm(n = 1, rep(0, m), V, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        
        Qp.perm = NULL
        for(j in 1:n.scan){
          idx = which(C[j,]!=0)
          Qp.perm = c(Qp.perm, Q.pval(U.perm[idx], V[idx,idx], weight=C[j,idx], r = 0) )
        }
        Qp.perm[Qp.perm<0] = 0  
        Qp.min.perm = c(Qp.min.perm, min(Qp.perm) )
      }
      
      if(p.comb=="minP"){
        pval = sum(Qp.min.perm <= Qp.min)/N.MC
      }
      
      if(!is.null(FWER)){
        
        signal.thresh = quantile(Qp.min.perm, FWER)
        cand.idx = which(Qp <= signal.thresh)
        
        if( length(cand.idx)>0 ){
          
          Qp.min.idx = which(Qp==Qp.min)
          if( length(Qp.min.idx)>1 ){
            Qp.min.idx = Qp.min.idx[which.min(sapply(cluster.list[Qp.min.idx], length))]
          }
          signal.region[[1]] = cluster.list[[Qp.min.idx]]
          signal.pval[1] = Qp[Qp.min.idx]
          
          signal.cand = cluster.list[cand.idx] 
          Qp.cand = Qp[cand.idx]
          signal.cand.min.region =  cluster.list[[Qp.min.idx]]
          #nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
          if(f==0){
            lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
          }else{
            lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) )/length(x) < f) )
          }
          
          while( length(lowoverlap.idx)>0  ){
            signal.cand = signal.cand[lowoverlap.idx] 
            Qp.cand = Qp.cand[lowoverlap.idx]
            Qp.min.idx = which(Qp.cand==min(Qp.cand))        
            if( length(Qp.min.idx)>1 ){
              Qp.min.idx = Qp.min.idx[which.min(sapply(signal.cand[Qp.min.idx], length))]
            }
            signal.region = c(signal.region , signal.cand[Qp.min.idx])
            signal.pval = c(signal.pval, Qp.cand[Qp.min.idx])
            signal.cand.min.region = signal.cand[[Qp.min.idx]]
            #nonoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
            if(f==0){
              lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) ) == 0) )
            }else{
              lowoverlap.idx = which( sapply(signal.cand, function(x) length( intersect(x, signal.cand.min.region) )/length(x) < f) )
            }
          }
        }        
      }

      
    }
    
    
  }

  
  return(list(Q.scan3.stat = Qp.min, Q.scan3.pval=pval, signal.region = signal.region, signal.pval = signal.pval, C=C))
  
}



Beta.Weights<-function(MAF,weights.beta, Cutoff=1, Is.MAF=TRUE){
  # no change
  n<-length(MAF)
  weights<-rep(0,n)
  Sign<-rep(1,n)
  #print(MAF)
  
  IDX1<-which(MAF > 0.5)
  if(length(IDX1) > 0){
    Sign[IDX1]<--1
    MAF[IDX1]<-1-MAF[IDX1]
  }
  
  IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
  if(length(IDX_0) == n){
    #stop("No polymorphic SNPs")
    weights<-rep(0,n)
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  
  weights = weights * Sign	
  return(weights)
  
}
SKAT_META_Optimal_Get_Q<-function(Score, r.all){
  # no change
  n.r<-length(r.all)
  Q.r<-rep(0,n.r)
  
  for(i in 1:n.r){
    r.corr<-r.all[i]
    Q.r[i]<-(1-r.corr) * sum(Score^2) + r.corr * sum(Score)^2
  }
  Q.r = Q.r /2
  
  re<-list(Q.r=Q.r)
  return(re)
}
SKAT_META_Optimal_Param<-function(Phi,r.all){
  # no change
  p.m<-dim(Phi)[2]
  r.n<-length(r.all)
  
  # ZMZ
  Z.item1.1<- Phi %*% rep(1,p.m)
  ZZ<-Phi
  ZMZ<- Z.item1.1 %*% t(Z.item1.1) / sum(ZZ)
  
  # W3.2 Term : mixture chisq
  W3.2.t<-ZZ - ZMZ
  lambda<-SKAT:::Get_Lambda(W3.2.t)
  
  # W3.3 Term : variance of remaining ...
  W3.3.item<-sum(ZMZ *(ZZ-ZMZ)) * 4
  
  # tau term 
  z_mean_2<- sum(ZZ)/p.m^2
  tau1<- sum(ZZ %*% ZZ) / p.m^2 / z_mean_2
  
  # Mixture Parameters
  MuQ<-sum(lambda)
  VarQ<-sum(lambda^2) *2 + W3.3.item
  KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
  Df<-12/KerQ
  
  # W3.1 Term : tau1 * chisq_1
  tau<-rep(0,r.n)
  for(i in 1:r.n){
    r.corr<-r.all[i]
    term1<-p.m^2*r.corr * z_mean_2 + tau1 * (1-r.corr)
    tau[i]<-term1
  }
  
  out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau,
            z_mean_2=z_mean_2, p.m=p.m,
            tau.1 = tau1,
            tau.2= p.m*z_mean_2 )
  
  #param2<<-out
  return(out)
}
SKAT_META_Optimal_Get_Pvalue<-function(Q.all, Phi, r.all, method){
  # no change
  n.r<-length(r.all)
  n.q<-dim(Q.all)[1]
  p.m<-dim(Phi)[2]
  
  lambda.all<-list()
  for(i in 1:n.r){
    r.corr<-r.all[i]
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi_rho<- L %*% (Phi %*% t(L))
    lambda.all[[i]]<-SKAT:::Get_Lambda(Phi_rho)
  }
  
  # Get Mixture param 
  param.m <- SKAT_META_Optimal_Param(Phi,r.all)
  Each_Info <- SKAT:::SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
  pmin.q<-Each_Info$pmin.q
  pval <- rep(0,n.q)
  
  # added
  pmin<-Each_Info$pmin
  
  for(i in 1:n.q){
    pval[i]<-SKAT:::SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])
  }
  
  # Check the pval 
  # Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
  # To correct conservatively, we use min(p-values) * 3
  
  multi<-3
  if(length(r.all) < 3){
    multi<-2
  }
  
  for(i in 1:n.q){
    pval.each<-Each_Info$pval[i,]
    IDX<-which(pval.each > 0)
    
    pval1<-min(pval.each) * multi
    if(pval[i] <= 0 || length(IDX) < length(r.all)){
      pval[i]<-pval1
    }
    # if pval==0, use nonzero min each.pval as p-value
    if(pval[i] == 0){
      if(length(IDX) > 0){
        pval[i] = min(pval.each[IDX])
      }
    }
    
  }
  
  return(list(p.value=pval,p.val.each=Each_Info$pval))
}
SKAT_META_Optimal <- function(Score, Phi, r.all, method="davies"){
  # no change
  # if r.all >=0.999 ,then r.all = 0.999
  IDX<-which(r.all >= 0.999)
  if(length(IDX) > 0){
    r.all[IDX]<-0.999	
  }
  
  p.m<-dim(Phi)[2]
  n.r<-length(r.all)
  
  ###########################################
  # Compute Q.r and Q.r.res
  ##########################################
  out.Q <- SKAT_META_Optimal_Get_Q(Score, r.all)
  Q.res=NULL
  Q.all<-rbind(out.Q$Q.r, Q.res) 
  
  ##################################################
  # Compute P-values 
  #################################################
  
  out<-SKAT_META_Optimal_Get_Pvalue(Q.all, Phi/2, r.all, method)
  
  param<-list(p.val.each=NULL, q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)
  
  id_temp<-which(param$p.val.each == min(param$p.val.each))
  id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
  if(length(id_temp1) > 0){
    param$rho[id_temp1] = 1
  }
  
  param$rho_est<-param$rho[id_temp]
  p.value<-out$p.value[1]
  re<-list(p.value = p.value, param=param)
  return(re)	
}


# library(CompQuadForm)
pchisqsum2 <- function (Q, lambda, delta = rep(0, length(lambda)), method = c("saddlepoint", "integration", "liu"), acc = 1e-07) 
{
  method <- match.arg(method)
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
  if (method == "saddlepoint") {
    p = saddle(Q, lambda, delta)
    if (is.na(p)) {
      method <- "integration"
    }
    else {
      return(list(p = p, errflag = 0))
    }
  }
  if (method == "integration") {
    tmp <- suppressWarnings( CompQuadForm::davies(q = Q, lambda = lambda, 
                                delta = delta, acc = acc) )
    if (tmp$ifault > 0) {
      lambda <- zapsmall(lambda, digits = 2)
      delta <- delta[lambda > 0]
      lambda <- lambda[lambda > 0]
      tmp <- suppressWarnings( CompQuadForm::farebrother(q = Q, lambda = lambda, 
                                       delta = delta) )
    }
    Qq <- if ("Qq" %in% names(tmp)) 
      tmp$Qq
    else tmp$res
    return(list(p = Qq, errflag = 0))
  }
  if (method == "liu") {
    tmp <- suppressWarnings( CompQuadForm::liu(Q, lambda = lambda, delta = delta) )
    return(list(p = tmp, errflag = 0))
  }
}

Met_SKAT_Get_Pvalue<-function(Score, Phi, r.corr, method){
  # change SKAT 
  Q.res = NULL
  p.m<-nrow(Phi)
  # if Phi==0
  if(sum(abs(Phi)) == 0){
    warning("No polymorphic SNPs!",call.=FALSE)
    return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
  }
  
  if(length(Phi) <=1){
    r.corr=0
  } else{
    if(ncol(Phi) <=10){
      if(qr(Phi)$rank <= 1){
        r.corr=0
      }
    }
  }
  
  if(length(r.corr) > 1){
    re = SKAT_META_Optimal(Score, Phi, r.corr, method=method)
    return(re)
  } 
  
  if (r.corr == 0){
    Q<-sum(Score^2)/2
  } else if (r.corr==1){
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    a<- as.matrix(sum(Phi))
    re <- SKAT:::Get_Liu_PVal(Q, a, Q.res)
    return(re)
  } else {
    # like r.corr = 0.1 or 0.2
    Q <- SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Phi<- L %*% (Phi %*% t(L))
  }
  lambda <- SKAT:::Get_Lambda(Phi/2)
  re1 <- pchisqsum2(Q, lambda, method = "integration", acc=1e-20)
  re <- list(p.value = re1$p, errflag = re1$errflag)
  # re<-SKAT:::Get_Davies_PVal(Q, Phi)
  return(re)
}
