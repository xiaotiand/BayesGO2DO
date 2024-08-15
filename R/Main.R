
iter = function(Yobs, delta, X, K, nchain,
                tau, nu0, sigma0, alpha1, alpha2, l1o, l2o, sigmal,
                GAMMA_vec, l1_vec, l2_vec, Theta_vec, BETA_vec, sigma_vec, Y_vec,
                cur, update_BETA, dependency, strength, ETA, S) {
  for (ich in 1:nchain) {
    Y = Y_vec[[ich]][[cur]]
    l1 = l1_vec[[ich]][[cur]]
    l2 = l2_vec[[ich]][[cur]]
    GAMMA = GAMMA_vec[[ich]][[cur]]
    Theta = Theta_vec[[ich]][[cur]]
    sigma = sigma_vec[[ich]][cur]
    BETA = BETA_vec[[ich]][[length(BETA_vec[[ich]])]]

    Gammavec = GammaGprod(Zprod(GAMMA, l1, l2, Theta, S),
                          GAMMA, l1, l2, Theta, S)
    for (k in 1:K) {
      locur = 1
      repeat{
        GAMMAnew = Gamma_prop(GAMMA, k)
        kres = pGamma_ratio(GAMMAnew, GAMMA, k, X[[k]], Y[[k]],
                            tau[k], sigma0, nu0, Gammavec,
                            l1, l2, Theta, S)
        locur = locur + 1
        alpha_ratio = runif(1, 0, 1)
        if(kres$ratio > alpha_ratio) {
          GAMMA = GAMMAnew
          Gammavec = kres$Gammavec
          break
        } else if(locur > 10) {
          break
        }
      }
    }
    GAMMA_vec[[ich]][[cur + 1]] = GAMMA

    if (dependency) {
      locur = 1
      repeat{
        l1new = L1_prop(l1, sigmal)
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1new, l2, Theta, S), GAMMA, l1new, l2, Theta, S)
        locur = locur + 1
        new_ratio = pL1_ratio(l1new, l1, l1o, sigmal, Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0, 1)
        if(new_ratio > alpha_ratio) {
          l1 = l1new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      l1_vec[[ich]][[cur + 1]] = l1

      locur = 1
      repeat{
        l2new = L2_prop(l2, sigmal)
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2new, Theta, S), GAMMA, l1, l2new, Theta, S)
        locur = locur + 1
        new_ratio = pL2_ratio(l2new, l2, l2o, sigmal, Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0, 1)
        if(new_ratio > alpha_ratio) {
          l2 = l2new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      l2_vec[[ich]][[cur + 1]] = l2
    } else {
      locur = 1
      repeat{
        l2new = L2_prop(l2, sigmal)
        l1new = l2new
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1new, l2new, Theta, S), GAMMA, l1new, l2new, Theta, S)
        locur = locur + 1
        new_ratio = pL2_ratio(l2new, l2, l2o, sigmal, Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0, 1)
        if(new_ratio > alpha_ratio) {
          l2 = l2new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      l2_vec[[ich]][[cur + 1]] = l2
      l1_vec[[ich]][[cur + 1]] = l1
    }

    if (strength) {
      for (k in 1:(K - 1)) {
        for (kp in (k + 1):K) {
          Thetanew = Theta_prop(Theta, alpha1, k, kp)
          if (dependency) {
            Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2, Thetanew, S), GAMMA, l1, l2, Thetanew, S)
          } else {
            Gammavec_new = GammaGprod2(Zprod2(GAMMA, l1, Thetanew, S), GAMMA, l1, Thetanew, S)
          }
          locur = locur + 1
          new_ratio = pTheta_ratio(Thetanew, Theta, GAMMA, k, kp, alpha1, alpha2, Gammavec_new, Gammavec)
          alpha_ratio = runif(1, 0, 1)
          if(new_ratio > alpha_ratio) {
            Theta = Thetanew
            Gammavec = Gammavec_new
          }
        }
      }
    } else {
      Theta = matrix(0, K, K)
    }
    Theta_vec[[ich]][[cur + 1]] = Theta

    sigma = sigma_prop(sigma, GAMMA, X, Y, tau, nu0, sigma0)
    sigma_vec[[ich]] = c(sigma_vec[[ich]], sigma)

    for (k in 1:K) {
      Y = Yk_prop(Y, Yobs, delta, X, GAMMA, tau, nu0, sigma0, k)
    }
    Y_vec[[ich]][[cur + 1]] = Y

    # if (update_BETA) {
    #   BETA = Beta_prop(BETA, GAMMA, X, Y, tau, sigma)
    #   BETA_vec[[ich]][[length(BETA_vec[[ich]]) + 1]] = BETA
    # }
  }

  cur = cur + 1
  result = list(cur = cur,
                GAMMA_vec = GAMMA_vec,
                l1_vec = l1_vec,
                l2_vec = l2_vec,
                Theta_vec = Theta_vec,
                BETA_vec = BETA_vec,
                sigma_vec = sigma_vec,
                Y_vec = Y_vec)
  return(result)
}


#' Main model fitting function
#'
#' This is the main model fitting function for Bayesian framework with GO and DO information.
#' @details
#' This function is used to fit a group of Bayesian regression models.
#' @import GOSemSim
#' @import DOSE
#' @import MCMCpack
#' @import crch
#' @import glmnet
#' @import mvtnorm
#' @import tmvtnorm
#'
#' @param Y A list of response values in regressions
#' @param delta A list of survival status (1 being deceased/uncensored and 0 being alive/censored)
#' @param X A list of design matrices in regressions
#' @param iteration The number of MCMC iterations
#' @param burnin The number of burn-in iterations
#' @param nchain The number of MCMC chains
#' @param dependency Whether to create dependencies among genes/biomarkers; The default is TRUE
#' @param strength Whether to borrow information/strength across regressions; The default is TRUE
#' @param ETA A G-by-G matrix of pairwise gene ontology similarity scores (generated by mgeneSim function in GOSemSim)
#' @param S A K-by-K symmetric matrix of disease ontology similarities (generated by doSim function in DOSE)
#' @return A large list of all the parameter values
#' @examples
#' require(mvtnorm)
#' require(DOSE)
#' require(GOSemSim)
#' K = 5
#' G = 10
#' lower = c(0.55, 0.55, 0.55, 0.55, 0.55)
#' upper = c(0.6, 0.6, 0.6, 0.6, 0.6)
#' Nk = c(10, 20, 30, 40, 50)
#' GAMMA = matrix(0, nrow = K, ncol = G)
#' GAMMA[1, 1:3] = 1
#' GAMMA[2, 1:3] = 1
#' GAMMA[3, 1:3] = 1
#' GAMMA[4, 1:3] = 1
#' GAMMA[5, 1:3] = 1
#' Y = lapply(Nk, function(l) rep(0, l))
#' BETA = matrix(1, nrow = K, ncol = G)
#' sigma = lapply(Nk, function(l) rep(0, l))
#' ck = list(N1 = rep(1, Nk[1]), N2 = rep(2, Nk[2]),
#'           N3 = rep(3, Nk[3]), N4 = rep(4, Nk[4]),
#'           N5 = rep(5, Nk[5]))
#' delta = list()
#' X = list()
#'
#' gene_list = c("TRPM3", "PKD2", "KCNJ16", "KCNJ15", "SLC2A11",
#'               "SLC2A9", "OSCP1", "SLC27A1", "SLC27A2", "SLC25A4")
#' hsGO = godata('org.Hs.eg.db', keytype = "SYMBOL", ont = "MF", computeIC = FALSE)
#' ETA = mgeneSim(gene_list, semData = hsGO, measure = "Wang",
#'                combine = "BMA", verbose = FALSE)
#' S = doSim(c("DOID:4471", "DOID:4467", "DOID:4465", "DOID:3910", "DOID:3907"),
#'           c("DOID:4471", "DOID:4467", "DOID:4465", "DOID:3910", "DOID:3907"), measure = "Wang")
#'
#' for (k in 1:K) {
#'  X[[k]] = matrix(rmvnorm(Nk[k], sigma = diag(G)), nrow = Nk[k], byrow = TRUE)
#'  sigma[[k]] = rnorm(Nk[k], 0, 1)
#'  Y[[k]] = X[[k]] %*% (GAMMA[k, ] * BETA[k, ]) + sigma[[k]]
#'  ck[[k]] = runif(Nk[k], quantile(Y[[k]], lower[k]), quantile(Y[[k]], upper[k]))
#'  delta[[k]] = as.numeric(Y[[k]] <= ck[[k]])
#'  Y[[k]] = pmin(Y[[k]], ck[[k]])
#' }
#'
#' t0 = Sys.time()
#' #res = Main(Y, delta, X, iteration = 500, burnin = 350, nchain = 1, ETA = ETA, S = S)
#' Sys.time() - t0
#'
#' @export

Main = function(Y, delta, X, iteration, burnin, nchain,
                tau = NULL, nu0 = 5, sigma0 = NULL,
                p1 = NULL, p0 = NULL, sigmal = NULL,
                alpha1 = NULL, alpha2 = NULL,
                dependency = TRUE, strength = TRUE,
                ETA, S) {
  require(mvtnorm)
  require(tmvtnorm)
  require(crch)
  require(MCMCpack)

  K = length(Y)
  G = dim(X[[1]])[2]
  if(is.null(tau)) {
    tau = rep(0.5, K)
  }
  if(is.null(sigma0)) {
    sigma0 = sd(unlist(Y)) ^ 2
  }
  if(is.null(sigmal)) {
    sigmal = 1
  }
  if(is.null(alpha1)) {
    alpha1 = 5
  }
  if(is.null(alpha2)) {
    alpha2 = 1
  }
  Xobs = X
  dobs = delta

  eta_vec = unlist(lapply(1:(dim(ETA)[1] - 1), function(i) ETA[i, i + 1]))
  eta_vec = eta_vec - 0.01
  p1 = (eta_vec / mean(eta_vec)) / (1 + eta_vec / mean(eta_vec))
  p1 = c(0.5, p1)
  xi_vec = 1 - eta_vec
  p0 = (xi_vec / mean(xi_vec)) / (1 + xi_vec / mean(xi_vec))
  p0 = c(0.5, p0)
  p1 = rep(0.2, G)
  p0 = rep(0.01, G)
  l1o = lambda1(p1)
  l2o = lambda2(p0)
  l2_vec = list()
  for (ich in 1:nchain) {
    l2_vec[[ich]] = list()
    if (length(unique(p0)) == 1) {
      l2_vec[[ich]][[1]] = rep(l2o[1], G)
    } else {
      l2_vec[[ich]][[1]] = l2o
    }
  }
  if (dependency) {
    l1_vec = list()
    for (ich in 1:nchain) {
      l1_vec[[ich]] = list()
      if (length(unique(p1)) == 1) {
        l1_vec[[ich]][[1]] = rep(l1o[1], G)
      } else {
        l1_vec[[ich]][[1]] = l1o
      }
    }
  } else {
    l1_vec = l2_vec
  }
  Theta_vec = list()
  for (ich in 1:nchain) {
    Theta_vec[[ich]] = list()
    if (strength) {
      Theta_vec[[ich]][[1]] = matrix(runif(1, 0, 1.5), K, K)
      diag(Theta_vec[[ich]][[1]]) = 0
    } else {
      Theta_vec[[ich]][[1]] = matrix(0, K, K)
    }
  }
  GAMMA_vec = list()
  for (ich in 1:nchain) {
    GAMMA_vec[[ich]] = list()
    GAMMA = matrix(0, K, G)
    GAMMA[, 1] = sample(c(0, 1), K, replace = TRUE, prob = c(1 - p1[1], p1[1]))
    for (g in 2:G) {
      for (k in 1:K) {
        if(GAMMA[k, g - 1] == 1) {
          GAMMA[k, g] = sample(c(0, 1), 1, prob = c(1 - p1[g], p1[g]))
        } else {
          GAMMA[k, g] = sample(c(0, 1), 1, prob = c(1 - p0[g], p0[g]))
        }
      }
    }
    GAMMA_vec[[ich]][[1]] = GAMMA
  }
  BETA_vec = list()
  sigma_vec = list()
  for (ich in 1:nchain) {
    BETA_vec[[ich]] = list()
    BETA_vec[[ich]][[1]] = matrix(0, K, G)
    sigma_vec[[ich]] = vector()
    sigma_vec[[ich]][1] = sigma_prop(sigma0, GAMMA_vec[[ich]][[1]], X, Y, tau, nu0, sigma0)
  }
  Y_vec = list()
  for (ich in 1:nchain) {
    Y_vec[[ich]] = list()
    Y_vec[[ich]][[1]] = Y
  }

  cur = 1
  for (index in 1:iteration) {
    if(cur <= burnin) {
      update_BETA = FALSE
    } else if(cur > burnin)  {
      update_BETA = TRUE
    }
    iter_res = iter(Yobs = Y, delta, X, K, nchain,
                    tau, nu0, sigma0, alpha1, alpha2, l1o, l2o, sigmal,
                    GAMMA_vec, l1_vec, l2_vec, Theta_vec, BETA_vec, sigma_vec, Y_vec,
                    cur, update_BETA, dependency, strength, ETA, S)
    cur = iter_res$cur
    GAMMA_vec = iter_res$GAMMA_vec
    l1_vec = iter_res$l1_vec
    l2_vec = iter_res$l2_vec
    Theta_vec = iter_res$Theta_vec
    BETA_vec = iter_res$BETA_vec
    sigma_vec = iter_res$sigma_vec
    Y_vec = iter_res$Y_vec
  }

  return(iter_res)
}

