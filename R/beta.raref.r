#' Rarefaction of true beta diversity (Whittaker's multiplicative beta)
#' @description Calculates true beta diversity rarefied to given number of samples.
#' @author David Zeleny (zeleny.david@@gmail.com). Heavily builds on the code from function \code{specaccum} in package \code{vegan} written by Roeland Kindt and Jari Oksanen to construct sample-based rarefaction curves.
#' @param comm Community matrix (or data frame) for which rarified beta diversity is calculated.
#' @param sites Number of sites for which rarefied true beta should be calculated.
#' @param conditioned Logical; should the result be conditioned on the provided dataset (\code{comm} matrix) or not? Default = \code{TRUE}.
#' @param gamma Estimator of gamma diversity used to estimate number of species in the community (in case that \code{conditioned = FALSE}). \code{gamma = c('chao', 'jack1', 'jack2', 'boot')}.
#' @return
#' The function returns list with the following items:
#' \itemize{
#' \item \code{sites}: number of sites for which rarefaction is calculated;
#' \item \code{gamma}: estimated number of species (gamma diversity) for given number of \code{sites};
#' \item \code{gamma.sd}: standard deviation of \code{gamma} estimate;
#' \item \code{alpha}: mean number of species per sample (alpha diversity);
#' \item \code{beta}: true beta diversity, calculated as \code{gamma / alpha};
#' \item \code{beta.sd}: standard deviation of true beta diversity, calculated as \code{gamma.sd / alpha} (perhaps wrongand better not to use, needs proof!);
#' }
#' @examples 
#' sc <- sample.comm (simul.comm (totS = 100), Np= 100)
#' beta.raref (comm = sc$a.mat, sites = 10)
#' 
#' # True beta diversity calculated on pair of samples (minus one) 
#' # is very close to mean pairwise Sorensen:
#' beta.raref (comm = sc$a.mat, sites = 2)$beta-1
#' library (vegan)
#' mean (vegdist (decostand (sc$a.mat, 'pa')))
#' @importFrom vegan specpool
#' @export
beta.raref <- function (comm, sites, conditioned = TRUE, gamma = 'jack1')
{
  comm <- as.matrix(comm)
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  n <- nrow(comm)
  p <- ncol(comm)
  alpha <- mean (rowSums (comm > 0))
  if (p == 1) {
    comm <- t(comm)
    n <- nrow(comm)
    p <- ncol(comm)
  }
  freq <- colSums(comm > 0)
  freq <- freq[freq > 0]
  f <- length(freq)
  ldiv <- lchoose(n, 1:n)
  result <- ifelse(n - freq < sites, 0, exp(lchoose(n - freq, sites) - ldiv[sites]))
  specaccum <- sum (1 - result)
  if (conditioned) {
    V <- result * (1 - result)
    tmp1 <- cor(comm > 0)
    ind <- lower.tri(tmp1)
    tmp1 <- tmp1[ind]
    tmp1[is.na(tmp1)] <- 0
    tmp2 <- outer(sqrt(V), sqrt(V))[ind]
    cv <- 2 * sum(tmp1 * tmp2)
    V <- sum(V)
    sdaccum <- sqrt(V + cv)
  } else {
    Stot <- vegan::specpool(comm)[, gamma]
    sdaccum1 <- sum((1 - result)^2)
    sdaccum2 <- specaccum^2/Stot
    sdaccum <- sqrt(sdaccum1 - sdaccum2)
  }
  out <- list(sites = sites, gamma = specaccum, gamma.sd = sdaccum, alpha = alpha, beta = specaccum/alpha, beta.sd = sdaccum/alpha)
  out
}