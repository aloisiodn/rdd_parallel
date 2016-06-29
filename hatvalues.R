hatvalues.ivreg <-function (model, ...) 
{
  xz <- model.matrix(model, component = "projected")
  x <- model.matrix(model, component = "regressors")
  z <- model.matrix(model, component = "instruments")
  solve_qr <- function(x) chol2inv(qr.R(qr(x)))
  diag(x %*% solve_qr(xz) %*% t(x) %*% z %*% solve_qr(z) %*% t(z))
}
