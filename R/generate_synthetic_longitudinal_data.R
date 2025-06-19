generate_synthetic_longitudinal_data <- function(nsamps=1000){
  
  # training set
  Tx = cbind(rnorm(nsamps),rnorm(nsamps)) # T1, T2
  colnames(Tx) = c(1,2)
  X = matrix(rnorm(nsamps*2),ncol=2)
  Y = matrix(rnorm(nsamps*15),ncol=15)
  colnames(Y) = c(rep(3,5),rep(4,5),rep(5,5))
  
  # test set
  Txte = cbind(rnorm(nsamps),rnorm(nsamps)) # T1, T2
  colnames(Txte) = c(1,2)
  Xte = matrix(rnorm(nsamps*2),ncol=2)
  Yte = matrix(rnorm(nsamps*15),ncol=15)
  colnames(Yte) = c(rep(3,5),rep(4,5),rep(5,5))

  return(list(Tx=Tx,X=X,Y=Y,
              Txte=Txte,Xte=Xte,Yte=Yte))
}
