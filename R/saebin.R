#' @export
saebebin <- function(y,  n){
  result <- list(theta_i_hat_EB2 = NA, Parameter = list(alpha = NA, beta = NA), MSE.theta_i_hat_EB2
                 = list(method = "Jackknife", mse = NA), direct = list(est = NA, mse = NA))

  if (any(is.na(y)))
    stop("Argument y=", response, " contains NA values.")
  if (any(is.na(n)))
    stop("Argument n=", sampel, " contains NA values.")

  m <- length(y) #jumlah area
  nT <- sum(n) #jumlah seluruh sampel
  w <- n/nT #bobot

  #menghitung dugaan langsung proporsi dan variansnya
  theta_i <- y/n #dugaan langsung proporsi
  mse.theta_i <- theta_i * (1-theta_i) / n #varians dugaan langsung
  result$direct$est <- theta_i
  result$direct$mse <- mse.theta_i

  #menghitung proporsi dan ragam proporsi
  theta_ib <- w * theta_i
  theta_hat <- sum(theta_ib) #rataan terboboti
  s_theta2 <- w * (theta_i-theta_hat)^2
  sum_s_theta2 <- sum(s_theta2) #ragam terboboti

  #menduga parameter sebaran beta-binomial alpha dan beta dengan metode momen Kleinman
  n2nT <- (n^2) / nT
  sum_n2nT <- sum(n2nT)
  k11 <- (nT*sum_s_theta2) - theta_hat*(1-theta_hat) * (m-1)
  k12 <- theta_hat*(1-theta_hat) * (nT-sum_n2nT-(m-1))
  k1 <- k11/k12

  #nilai alpha beta
  alpha <- theta_hat*(1-k1)/k1
  beta <- alpha*(1/theta_hat-1)

  #pendugaan bayes dan ragam posterior bagi theta_i
  k21 <- (y+alpha)*(n-y+beta)
  k22 <- (n+alpha+beta+1)*(n+alpha+beta)^2
  theta_i_hat_EB1 <- (y+alpha)/(n+alpha+beta) #penduga bayes
  var_theta_i_hat_EB <- k21/k22 #ragam posterior

  #pendugaan empirical bayes bagi theta_i
  gamma <- n/(n+alpha+beta) #nilai gamma
  theta_i_hat_EB2 = gamma*theta_i+(1-gamma)*theta_hat #penduga empirical bayes
  result$Parameter$alpha <- alpha
  result$Parameter$beta <- beta
  result$theta_i_hat_EB2 <- theta_i_hat_EB2

  #Menghitung MSE penduga empirical bayes
  jackknife <- function(y, n, l){
    theta_i_jk <- 0
    n.jk <- 0
    theta_hat_jk <- 0
    theta_i_jk <- y[-l]/n[-l] #dugaan langsung proporsi
    nT.jk <- sum(n[-l])
    theta_ib_jk <- n[-l]/nT.jk * theta_i_jk
    theta_hat_jk <- sum(theta_ib_jk) #rataan terboboti
    s_theta2_jk <- n[-l]/nT.jk * (theta_i_jk-theta_hat_jk)^2
    sum_s_theta2_jk <- sum(s_theta2_jk) #ragam terboboti
    n2nT.jk <- (n[-l]^2) / nT.jk
    sum_n2nT.jk <- sum(n2nT.jk)
    k11.jk <- (nT.jk*sum_s_theta2_jk) - theta_hat_jk*(1-theta_hat_jk) * (m-1)
    k12.jk <- theta_hat_jk*(1-theta_hat_jk) * (nT.jk-sum_n2nT.jk-(m-1))
    k1.jk <- k11.jk/k12.jk
    alpha <- theta_hat_jk*(1-k1.jk)/k1.jk #alpha
    beta <- alpha*(1/theta_hat_jk-1) #beta
    k21 <- (y+alpha)*(n-y+beta)
    k22 <- (n+alpha+beta+1)*(n+alpha+beta)^2
    theta_i_hat_EB1 <- (y+alpha)/(n+alpha+beta) #penduga bayes
    var_theta_i_hat_EB <- k21/k22 #ragam posterior
    #pendugaan empirical bayes bagi theta_i
    gamma <- n/(n+alpha+beta) #nilai gamma
    theta_i_hat_EB2 = gamma*theta_i+(1-gamma)*theta_hat #penduga empirical bayes
    result <- list(gamma = gamma, theta_i_hat_EB2 = theta_i_hat_EB2, var_theta_i_hat_EB = var_theta_i_hat_EB)
    return(result)
  }
  jk <- lapply(1:m, function(l) jackknife(y, n, l))
  M1 <- sapply(1:m, function(i){
    m1 <- var_theta_i_hat_EB[i]-(m-1)/m*sum(sapply(1:m, function(l){
      return(jk[[l]]$var_theta_i_hat_EB[i]-var_theta_i_hat_EB[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:m, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(l){
      return((jk[[l]]$theta_i_hat_EB2[i]-theta_i_hat_EB2[i])^2)
    }))
    return(m2)
  })
  mse <- M1+M2
  result$MSE.theta_i_hat_EB2$mse <- mse
  return(result)
}

#' @export
saelognom <- function(y, n){
  result <- list(pi_eb = NA, Parameter = list(miu = NA, sigma = NA), MSE.pi_eb
                 = list(method = "Jackknife", mse = NA), direct = list(est = NA, mse = NA))

  if (any(is.na(y)))
    stop("Argument y=", response, " contains NA values.")
  if (any(is.na(n)))
    stop("Argument n=", sampel, " contains NA values.")

  m <- length(y) #jumlah area

  #menghitung dugaan langsung proporsi dan variansnya
  pi = y/n #dugaan langsung proporsi
  mse.pi <- pi * (1-pi) / n #varians dugaan langsung
  result$direct$est <- pi
  result$direct$mse <- mse.pi


  #menduga parameter sebaran logit-normal miu dan sigma
  logit = log(pi/(1-pi))

  #nilai miu dan sigma
  miu=mean(logit)#perhitungan rata-rata logit
  sigma=sd(logit) #perhitungan standard deviasi logit


  #ragam posterior bagi pi
  g = ((y+miu)*(n-y+sigma)) / ((miu+n+sigma+1)*(miu+n+sigma)^2)

  #pendugaan empirical bayes bagi pi
  zi <- (logit-miu)/sigma
  par = miu + (sigma * zi)
  h1 = exp(par) / (1+exp(par))
  h2 = par*yi-ni*log(1+exp(par))
  nz=(1/(2*3.14159))*exp(-1/2*((zi^2)))
  a = h1*exp(h2)*nz
  b = exp(h2)*nz
  pi_eb = a/b

  result$Parameter$miu <- miu
  result$Parameter$sigma <- sigma
  result$pi_eb <- pi_eb

  #Menghitung MSE penduga empirical bayes
  jackknife <- function(y, n, l){
    #pi_jk <- 0
    #n.jk <- 0
    pi_jk <- y[-l]/n[-l] #dugaan langsung proporsi
    mse.pi_jk <- pi_jk * (1-pi_jk) / n[-l] #varians dugaan langsung
    logit_jk = log(pi_jk/(1-pi_jk))

    #nilai miu dan sigma
    miu_jk=mean(logit_jk)#perhitungan rata-rata logit
    sigma_jk=sd(logit_jk) #perhitungan standard deviasi logit

    #pendugaan empirical bayes bagi pi
    zi_jk <- (logit_jk-miu_jk)/sigma_jk
    par_jk = miu_jk + (sigma_jk * zi_jk)
    h1_jk = exp(par_jk) / (1+exp(par_jk))
    h2_jk = par_jk*y[-l]-n[-l]*log(1+exp(par_jk))
    nz_jk=(1/(2*3.14159))*exp(-1/2*((zi_jk^2)))
    a_jk = h1_jk*exp(h2_jk)*nz_jk
    b_jk = exp(h2_jk)*nz_jk
    pi_eb_jk = a_jk/b_jk #penduga bayes
    g_jk = ((y+miu_jk)*(n-y+sigma_jk)) / ((miu_jk+n+sigma_jk+1)*((miu_jk+n+sigma_jk)^2)) #ragam posterior
    result <- list(pi_eb = pi_eb_jk, g = g_jk)
    return(result)
  }
  jk <- lapply(1:m, function(l) jackknife(y, n, l))
  M1 <- sapply(1:m, function(i){
    m1 <- g[i]-(m-1)/m*sum(sapply(1:m, function(l){
      return(jk[[l]]$g[i]-g[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:m, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(l){
      return((jk[[l]]$pi_eb[i]-pi_eb[i])^2)
    }))
    return(m2)
  })
  mse <- M1
  result$MSE.pi_eb$mse <- mse

  return(result)
}
