
getmode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


RemoveDuplicate = function(X, thresh = 0.95, listreturn = T){
  
  S = cor(X, use = "pairwise.complete.obs" )
  
  diag(S) = 0
  
  index.arr.over.thresh = which(abs(S) > thresh , arr.ind = T)
  
  index = vector()
  
  while(length(index.arr.over.thresh) != 0){
    
    H = getmode(as.vector(index.arr.over.thresh))
    
    index = c(index, H)
    
    index.arr.over.thresh =
      index.arr.over.thresh[ (index.arr.over.thresh[,1] !=  H) &
                               (index.arr.over.thresh[,2] !=  H),]
    
  }
  
  if(length(index) == 0){
    
    indlist =  colnames(X)
    
  }else{
    
    indlist =  colnames(X)[-index]
    
  }
  
  if(listreturn){
    
    return(indlist)
    
  }else{
    
    return(X[,indlist])
  }
  
}

na.process <- function(X){
  apply(X,2,function(x)
  {
    x[is.na(x)] <- mean(x,na.rm=TRUE);
    return(x);
  });
} # Imputation of NA values

Z.to.Mg <- function(x, n){
  
  return(x /sqrt(x^2 + n - 2))
  
}
