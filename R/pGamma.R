
lambda1 = function(p1) {
  unlist(lapply(p1, function(x) log(x) - log(1 - x)))
}


lambda2 = function(p0) {
  unlist(lapply(p0, function(x) log(x) - log(1 - x)))
}


expGamma = function(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta, S) {
  sum = lambda1g * gamma_g %*% gamma_g_1
  sum = sum + lambda2g * gamma_g %*% (1 - gamma_g_1)
  sum = sum + gamma_g %*% (Theta * S) %*% gamma_g
  return(exp(sum))
}


expGamma1 = function(gamma_g, lambda1g, Theta, S) {
  sum = lambda1g * gamma_g %*% gamma_g
  sum = sum + gamma_g %*% (Theta * S) %*% gamma_g
  return(exp(sum))
}


Zloop = function(vecs, gamma_g_1, lambda1g, lambda2g, Theta, S) {
  sum = apply(vecs, 2, function(x) expGamma(x, gamma_g_1, lambda1g, lambda2g, Theta, S))
  return(sum(sum))
}


Zloop1 = function(vecs, lambda1g, Theta, S) {
  sum = apply(vecs, 2, function(x) expGamma1(as.numeric(x), lambda1g, Theta, S))
  return(sum(sum))
}


Z = function(K, gamma_g_1, lambda1g, lambda2g, Theta, S) {
  vecs = t(as.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE))))
  sum = Zloop(vecs, gamma_g_1, lambda1g, lambda2g, Theta, S)
  return(sum)
}


pGamma_g = function(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta, S) {
  K = length(gamma_g)
  ratio = expGamma(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta, S) /
    Z(K, gamma_g_1, lambda1g, lambda2g, Theta, S)
  return(ratio)
}


Z1 = function(K, lambda1g, Theta, S) {
  vecs = t(as.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE))))
  sum = Zloop1(vecs, lambda1g, Theta, S)
  return(sum)
}


pGamma_1 = function(gamma_g, lambda1g, Theta, S) {
  K = length(gamma_g)
  ratio = expGamma1(gamma_g, lambda1g, Theta, S) / Z1(K, lambda1g, Theta, S)
  return(ratio)
}


Zprod = function(GAMMA, l1, l2, Theta, S) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Zvec = c(Z1(K, l1[1], Theta, S),
           unlist(lapply(2:G, function(g) Z(K, GAMMA[, g - 1], l1[g], l2[g], Theta, S) )))
  return(Zvec)
}


Zprod2 = function(GAMMA, l1, Theta, S) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Zvec = unlist(lapply(1:G, function(g) Z(K, GAMMA[, g], l1[g], 0, Theta, S) ))
  return(Zvec)
}


GammaGprod = function(Zvec, GAMMA, l1, l2, Theta, S) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Gammavec = c(expGamma1(GAMMA[, 1], l1[1], Theta, S),
               unlist(lapply(2:G, function(g) expGamma(GAMMA[, g], GAMMA[, g - 1], l1[g], l2[g], Theta, S) )))
  Gammavec = exp(log(Gammavec) - log(Zvec))
  return(Gammavec)
}


GammaGprod2 = function(Zvec, GAMMA, l1, Theta, S) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Gammavec = unlist(lapply(1:G, function(g) expGamma(GAMMA[, g], GAMMA[, g], l1[g], 0, Theta, S) ))
  Gammavec = exp(log(Gammavec) - log(Zvec))
  return(Gammavec)
}


pGamma_ratio = function(GAMMAnew, GAMMAold, k, Xk, Yk, tauk, sigma0, nu0, Gammavec_old,
                        l1, l2, Theta, S) {
  Gammavec_new = Gammavec_old
  c = GAMMAnew[k, ] - GAMMAold[k, ]
  c = which(c != 0)
  for (i in 1:length(c)) {
    if (c[i] == 1) {
      Gammavec_new[1] = pGamma_1(GAMMAnew[, 1], l1[1], Theta, S)
    } else if (c[i] < length(l1)) {
      Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i] - 1],
                                    l1[c[i]], l2[c[i]], Theta, S)
      Gammavec_new[c[i] + 1] = pGamma_g(GAMMAnew[, c[i] + 1], GAMMAnew[, c[i]],
                                        l1[c[i] + 1], l2[c[i] + 1], Theta, S)
    } else {
      Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i] - 1],
                                    l1[c[i]], l2[c[i]], Theta, S)
    }
  }

  gamma_new = which(GAMMAnew[k, ] == 1)
  gamma_old = which(GAMMAold[k, ] == 1)
  Sigma_gammanew = tauk * Xk[, gamma_new] %*% t(Xk[, gamma_new])
  Sigma_gammaold = tauk * Xk[, gamma_old] %*% t(Xk[, gamma_old])
  diag(Sigma_gammanew) = diag(Sigma_gammanew) + 1
  diag(Sigma_gammaold) = diag(Sigma_gammaold) + 1
  Sigma_gammanew = sigma0 * Sigma_gammanew
  Sigma_gammaold = sigma0 * Sigma_gammaold
  mulTnew = dmvt(as.numeric(Yk), sigma = Sigma_gammanew, df = 2 * nu0, log = TRUE)
  mulTold = dmvt(as.numeric(Yk), sigma = Sigma_gammaold, df = 2 * nu0, log = TRUE)

  ratio = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * exp(mulTnew - mulTold)
  if(is.na(ratio)) ratio = 0
  return(list(ratio = ratio,
              Gammavec = Gammavec_new))
}


pGamma_ratio2 = function(GAMMAnew, GAMMAold, k, Xk, Yk, tauk, sigma0, nu0, Gammavec_old,
                         l1, Theta, S) {
  Gammavec_new = Gammavec_old
  c = GAMMAnew[k, ] - GAMMAold[k, ]
  c = which(c != 0)
  for (i in 1:length(c)) {
    Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i]],
                                  l1[c[i]], 0, Theta, S)
  }

  gamma_new = which(GAMMAnew[k, ] == 1)
  gamma_old = which(GAMMAold[k, ] == 1)
  Sigma_gammanew = tauk * Xk[, gamma_new] %*% t(Xk[, gamma_new])
  Sigma_gammaold = tauk * Xk[, gamma_old] %*% t(Xk[, gamma_old])
  diag(Sigma_gammanew) = diag(Sigma_gammanew) + 1
  diag(Sigma_gammaold) = diag(Sigma_gammaold) + 1
  Sigma_gammanew = sigma0 * Sigma_gammanew
  Sigma_gammaold = sigma0 * Sigma_gammaold
  mulTnew = dmvt(as.numeric(Yk), sigma = Sigma_gammanew, df = 2 * nu0, log = TRUE)
  mulTold = dmvt(as.numeric(Yk), sigma = Sigma_gammaold, df = 2 * nu0, log = TRUE)

  ratio = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * exp(mulTnew - mulTold)
  if(is.na(ratio)) ratio = 0
  return(list(ratio = ratio,
              Gammavec = Gammavec_new))
}


pTheta_ratio = function(Thetanew, Thetaold, GAMMA, k, kp, alpha1, alpha2, Gammavec_new, Gammavec_old) {
  thetanew = Thetanew[k, kp]
  thetaold = Thetaold[k, kp]
  res = dgamma(thetanew, shape = alpha1, rate = alpha2) /
    dgamma(thetaold, shape = alpha1, rate = alpha2)
  res = res * exp(sum(log(Gammavec_new) - log(Gammavec_old)))
  if(is.na(res)) res = 0
  return(res)
}


pL1_ratio = function(l1_new, l1_old, l1o, sigmal, Gammavec_new, Gammavec_old) {
  if (length(unique(l1_new)) == 1) {
    dnew = dnorm(l1_new[1], l1o[1], sigmal)
    dold = dnorm(l1_old[1], l1o[1], sigmal)
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * dnew / dold
    if(is.na(res)) res = 0
    return(res)
  } else {
    dnew = dtmvnorm(l1_new, l1o, sigmal * diag(length(l1o)), log = TRUE)
    dold = dtmvnorm(l1_old, l1o, sigmal * diag(length(l1o)), log = TRUE)
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old)) + dnew - dold)
    if(is.na(res)) res = 0
    return(res)
  }
}


pL2_ratio = function(l2_new, l2_old, l2o, sigmal, Gammavec_new, Gammavec_old) {
  if (length(unique(l2_new)) == 1) {
    dnew = dnorm(l2_new[1], l2o[1], sigmal)
    dold = dnorm(l2_old[1], l2o[1], sigmal)
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * dnew / dold
    if(is.na(res)) res = 0
    return(res)
  } else {
    dnew = dtmvnorm(l2_new, l2o, sigmal * diag(length(l2o)), log = TRUE)
    dold = dtmvnorm(l2_old, l2o, sigmal * diag(length(l2o)), log = TRUE)
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old)) + dnew - dold)
    if(is.na(res)) res = 0
    return(res)
  }
}


Gamma_prop = function(GAMMA, k) {
  Gamma = GAMMA[k, ]
  sors = sample(c(0, 1), 1)
  if (sors == 0 & length(which(Gamma == 1)) > 1) {
    ones = which(GAMMA[k, ] == 1)
    zeros = which(GAMMA[k, ] == 0)
    Gamma[sample(ones, 1)] = 0
    Gamma[sample(zeros, 1)] = 1
  } else if (sors == 1 & length(which(Gamma == 1)) > 1) {
    c1 = sample(1:length(Gamma), 1)
    Gamma[c1] = 1 - Gamma[c1]
  } else if (length(which(Gamma == 1)) <= 1) {
    zeros = which(GAMMA[k, ] == 0)
    Gamma[sample(zeros, 1)] = 1
  }

  GAMMA[k, ] = Gamma
  return(GAMMA)
}


Theta_prop = function(Theta, alpha1, k, kp) {
  theta = Theta[k, kp]
  theta = rgamma(1, shape = alpha1, rate = alpha1 / theta)
  Theta[k, kp] = theta
  Theta[kp, k] = theta
  return(Theta)
}


L1_prop = function(l1, sigmal) {
  if (length(unique(l1)) == 1) {
    l1new = as.numeric(rtmvnorm(1, l1[1], sigmal, upper =  0))
    return(rep(l1new, length(l1)))
  } else {
    l1new = rep(0, length(l1))
    for (i in 1:length(l1new)) {
      l1new[i] = as.numeric(rtmvnorm(1, l1[i], sigmal, upper =  0, algorithm = "gibbs"))
    }
    return(l1new)
  }
}


L2_prop = function(l2, sigmal) {
  if (length(unique(l2)) == 1) {
    l2new = as.numeric(rtmvnorm(1, l2[1], sigmal, upper =  0))
    return(rep(l2new, length(l2)))
  } else {
    l2new = rep(0, length(l2))
    for (i in 1:length(l2new)) {
      l2new[i] = as.numeric(rtmvnorm(1, l2[i], sigmal, upper =  0, algorithm = "gibbs"))
    }
    return(l2new)
  }
}


sigma_prop = function(sigma, GAMMA, X, Y, tau, nu0, sigma0) {
  shape = (sum(unlist(lapply(Y, length))) + 2 * nu0) / 2
  scale = 0
  for (k in 1:dim(GAMMA)[1]) {
    Xk = X[[k]]
    Yk = Y[[k]]
    Gamma = which(GAMMA[k, ] == 1)
    SigmaY = tau[k] * Xk[, Gamma] %*% t(Xk[, Gamma]) + diag(length(Yk))
    scale = scale + t(Yk) %*% solve(SigmaY) %*% Yk
  }
  scale = scale / 2 + nu0 * sigma0
  return(MCMCpack::rinvgamma(1, shape = shape, scale = scale))
}


Yk_prop = function(Y, Yobs, delta, X, GAMMA, tau, nu0, sigma0, k) {
  cs = which(delta[[k]] == 0)
  uncs = which(delta[[k]] == 1)
  if (length(cs) == 0) {
    Y[[k]] = Yobs[[k]]
  } else {
    Xk = X[[k]]
    Yk = Y[[k]]
    gamma = which(GAMMA[k, ] == 1)
    Sigma_gamma = tau[k] * Xk[, gamma] %*% t(Xk[, gamma])
    diag(Sigma_gamma) = diag(Sigma_gamma) + 1
    Sigma_gamma = sigma0 * Sigma_gamma
    Sigma11 = Sigma_gamma[-cs, ]
    Sigma11 = Sigma11[, -cs]
    Sigma21 = Sigma_gamma[cs, ]
    Sigma21 = Sigma21[, -cs]
    Sigma12 = Sigma_gamma[, cs]
    Sigma12 = Sigma12[-cs, ]
    invSigma11 = solve(Sigma11)
    mu21 = Sigma21 %*% invSigma11 %*% Yk[-cs]
    Sigma221 = Sigma_gamma[cs, cs] - Sigma21 %*% invSigma11 %*% Sigma12
    ll = Yobs[[k]][cs]
    p1 = length(uncs)
    d1 = as.numeric(Yk[-cs] %*% invSigma11 %*% Yk[-cs])
    impute = rtmvt(1, mean = as.numeric(mu21), sigma = (2 * nu0 + d1) / (2 * nu0 + p1) * Sigma221,
                   df = 2 * nu0 + p1, lower = ll, algorithm = "gibbs")
    impute = as.numeric(impute)
    Y[[k]][cs] = impute
  }
  return(Y)
}

