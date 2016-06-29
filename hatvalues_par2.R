
hatvalues.ivreg_par <-function (model, cl=NULL, ...) 
{
  #xz <- model.matrix(model, component = "projected") 
  #x <- model.matrix(model, component = "regressors") 
  #z <- model.matrix(model, component = "instruments") 
  #solve_qr <- function(x) chol2inv(qr.R(qr(x)))
  #sqr_xz <- solve_qr(xz)
  #sqr_z <- solve_qr(z)


  solve_qr <- function(x) chol2inv(qr.R(qr(x)))
  
  l_ma <- list(model.matrix(model, component = "regressors"), 
               t(model.matrix(model, component = "regressors")),
                 solve_qr(model.matrix(model, component = "instruments")) ) 
  l_mb <- list(solve_qr(model.matrix(model, component = "projected")),
               model.matrix(model, component = "instruments"),
               t(model.matrix(model, component = "instruments")))             
  
  f<-foreach(i = 1:3) %dopar%  {
    print(paste('inicio',i))
    fx<-l_ma[[i]] %*% l_mb[[i]]
    print(paste('fim',i))
    fx
  }
  
  diag(f[[1]]%*%f[[2]]%*%f[[3]])
}

