#' Co-occurrence based metric of species habitat specialization
#' @description A set of functions for estimating species niche breadth, based on species composition dataset using co-occurrence-based \emph{theta} metric. Introduced by Fridley et al. (2007).
#' @author David Zeleny (zeleny.david@@gmail.com). Partly based on codes written by Jason Fridley (Fridley et al. 2007) and David Zeleny (Zeleny 2009), extended for other published algorithms and optimised for speed and applicability on large datasets. Function \code{beals.2} is based on function \code{beals} from \code{vegan}, written by Miquel De Caceres and Jari Oksanen.
#' @param comm Community data (\code{matrix} or \code{data.frame}, samples x species). If data are not presence-absence, the matrix will be automatically transformed into presence-absence and warning will be printed.
#' @param species.data Species data (\code{matrix} or \code{data.frame}). If suplied, it should have at least two columns - the first containing species name, the second containing layer. 
#' @param thresh Minimal frequency of species. Habitat specialization will be calculated for species occurring in number of samples equal or higher than minimal frequency threshold. Default = \code{5}.
#' @param psample Size of one random subsample (number of samples) for methods based on subsampling (argument \code{method = c('additive', 'multiplicative', 'multi.sorensen', 'multi.simpson', 'beals', 'rao')}). This value should not be higher than mimal frequency of species (argument \code{thresh}). For default setting (\code{method = 'multiplicative', rarefaction = TRUE} this value number of samples on rarefaction curve on which all the beta diversity calculation is standardized (equivalent to number of subsamples). Default = \code{5}.
#' @param reps Number of random subsamples. Specifies how many times the fixed number of samples (specified by \code{psample}) will be randomly drawn from all samples containing target species. Default = \code{10}.
#' @param method Beta-diversity algorithm used to calculate theta measure. Partial match to \code{'additive'}, \code{'multiplicative'}, \code{'pairwise.jaccard'}, \code{'pairwise.sorensen'}, \code{'pairwise.simpson'}, \code{'multi.sorensen'}, \code{'multi.simpson'}, \code{'rao'}, \code{'beals'}). See Details for available options.
#' @param q Generalization of Whittaker's multiplicative beta diversity for abundance data (only if \code{method = 'multiplicative'}). 
#' @param rarefaction Logical value, which applies for \code{method = 'multiplicative'} and \code{q = 0}: should the Whittaker's multiplicative beta be calculated by rarefaction (\code{rarefaction = TRUE}) or by subsampling (\code{rarefaction = FALSE})? Default = \code{TRUE}.
#' @param beals.file Contains pre-calculated matrix of species co-occurrences. Can be used if \code{method = 'beals'} to speed up repeated calculation.
#' @param pa.transform Should the compositional data be transformed into presence-absence form? This argument applies only if \code{method} is \code{"rao"}, since this method can be meaningfully calculated with either presence absence or abundance data. For \code{method = 'multi'}, the data need to be presence-absence for \code{q = 0} (true beta diversity) and abundance for \code{q = 1} or \code{q = 2} (with presence-absence species values, the \code{q = 1 or 2} returns the same beta diversity as the method with \code{q = 0}). Other methods (Jaccard, Sorensen, Simpson) are calculated on species presences only.
#' @param force.subsample Logical; should the subsampling be forced even for beta diversity metrics which are not influenced by sample size? Default behaviour is \code{force.subsample = FALSE}, meaning that beta diversity metrics dependent on sample size (\code{method = c('additive', 'multiplicative', 'beals', 'multi.sorensen', 'multi.simpson')}) will subsample number of plots equal to \code{psample}, while method not dependent on sample size (\code{method = c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'rao')}) will be calculated using all plots containing target species. If \code{force.subsample = TRUE}, all indices will be calculated using subsampling. 
#' @param parallel Logical; should be the parallel calculation used?
#' @param no.cores Number of cores (if \code{parallel = TRUE}). Note that in case of large datasets the calculation may be limited by RAM of the computer, and increasing number of cores may result in saturation of RAM and calculation collapse.
#' @param remove.out Logical; should be the algorithm removing outliers (sensu Botta-Dukat 2012) applied? 
#' @param out.metric Dissimilarity metric used to calculate outliers which should be removed (\code{out.metric = c('sorensen', 'euclidean', 'binary.euclidean')}). Default value is \code{'sorensen'}, which is compatible with Whittaker's multiplicative beta; if using other \code{method}, consider to change it into \code{'binary.euclidean'} recommended by Botta-Dukat (2012).
#' @param temp.matrix Internal argument; matrix with species composition of plots containing target species.
#' @param sp Internal argument; the order of the species for which the current calculation is done.
#' @param sci.name Internal argument; the name of the species for which the current calculation is done.
#' @param x Internal argument of \code{beals.2} function; input compositional matrix.
#' @param include Internal argument of \code{beals.2} function; include argument from \code{vegan::beals}.
#' @details
#' \describe{
#' Function \code{calculate.theta} calculates theta metric of species habitat specialization using range of proposed beta diversity measures. It uses internal functions \code{calculate.theta.0}, \code{beals.2} (modified from the library \code{vegan} to calculate the sample species pool using Beals smoothing method).
#' The function \code{calculate.theta} offers the following \code{method} argument to calculate beta diversity among samples:
#' \item{\code{additive}}{This is the original algorithm published by Fridley et al. (2007), in which beta diversity among samples containing given species is calculated by additive beta diversity measure.}
#' \item{\code{multiplicative}}{This is the default method, which uses the multiplicative Whittaker's measure of beta diversity instead of the original additive measure, as suggested by Zeleny (2009). Two options are available - using rarefaction of true beta diversity (if \code{rarefaction = TRUE}) to given number of samples (argument \code{psample}), or by repeated subsampling of \code{psample} from the dataset \code{reps}-times (if \code{rarefaction = FALSE}); both methods give comparable results, and the rarefaction one is usually faster. Modification of argument \code{q} calculates multiplicative beta diversity based on number equivalents (or number of effective species, Jost 2007). \code{q = 0} calculates Whittaker's beta, which weights all species equally (meaning that rare species, which are the most susceptible to undersampling, are weighted equally to abundant species); \code{q = 1} calculates number equivalents for Shannon diversity and \code{q = 2} for Simspon diversity. Uses function \code{d} from the packages \code{vegetarian}.}
#' \item{\code{beals}}{Multiplicative beta on species pool. Algorithm suggested by Botta-Dukat (2012), calculating the beta diversity using species pool matrix instead of the original species data matrix. Species pool matrix is calculated using Beals smoothing method (invented by Ewald 2002). While the previous multiplicative beta diversity method gives unbiased results only in case of not-saturated communities, this method should give unbiased results also in case of saturated communities. See Zeleny (2009) and Botta-Dukat (2012) for detail discussion of this saturated/not-saturated communities issue. Argument \code{q} have no effect, since the recalculated species pool data are presence-absence only.}
#' \item{\code{pairwise.jaccard}, \code{pairwise.sorensen}, \code{pairwise.simpson}, \code{multi.sorensen} and \code{multi.simpson}}{Mean pairwise Jaccard, Sorensen and Simpson dissimilarity, and multiple Sorensen and Simpson dissimilarity based on reccomendations of Manthey & Fridley (2009). Authors suggested that neither the original additive algorithm (introduced by Fridley et al. 2007), neither the modified version using the multiplicative beta diversity (Zeleny 2009) is the best solution, and introduced other alternatives, using pairwise or multiple site beta diversity algorithm. Mean pairwise Jaccard dissimilarity (or Sorensen and Simpson, respectively) is based on calculating mean of Jaccard (or Sorensen and Simpson, respectively) dissimilarities among all pairs of samples in each subset, while multiple Sorensen (or Simpson, respectively) is using multiple-site Sorensen (or Simpson, respectively) algorithm introduced by Baselga et al. (2007). Multiple-site Sorensen index is a linear function of Whittaker's beta diversity.}
#' \item{\code{rao}}{Rao index of dissimilarity; this option has been introduced and used by Boulangeat et al. (2012). Advantage of Rao index is a possibility to incorporate known relationships among species using the among-species distance matrix. The formula used here is based on de Bello et al. (2010) \emph{beta.rao = (gamma - mean.alpha)/(1 - mean.alpha)} and is calculated by function \code{RaoRel} from package \code{cati}.}
#' }
#' @return
#' The function \code{calculate.theta} returns data.frame, with species in rows and the following columns:
#' \itemize{
#' \item \code{sci.name}: scientific name of the species;
#' \item \code{local.avgS}: average local species richness (average number of species in plots containing target species);
#' \item \code{occur.freq}: occurrence frequency: number of plots in which species occurs;
#' \item \code{meanco}: mean number of co-occurring species in subset of selected plots;
#' \item \code{meanco.sd}: sd of the number of co-occurring species in subset of selected plots;
#' \item \code{meanco.u, meanco.l}: upper and lower confidence interval of the number of co-occuring species in subset of selected plots;
#' \item \code{theta}: calculated theta value;
#' \item \code{theta.sd}: standard deviation of calculated theta values for individual subsets (not available for metrics which are not calculated by subsampling).
#' }
#' @references
#' Baselga A., Jimenez-Valverde A. & Niccolini G. (2007): A multiple-site similarity measure independent of richness. \emph{Biology Letters}, 3: 642-645.
#' 
#' Baselga A., Orme D., Villeger S., Bortoli J. & Leprieur F. (2013): betapart: Partitioning beta diversity into turnover and nestedness components. R package version 1.3. http://CRAN.R-project.org/package=betapart
#' 
#' Botta-Dukat Z. (2012): Co-occurrence-based measure of species' habitat specialization: robust, unbiased estimation in saturated communities. \emph{Journal of Vegetation Science}, 23: 201-207.
#' 
#' Boulangeat I., Lavergne S., Van Es J., Garraud L. & Thuiller W. (2012): Niche breadth, rarity and ecological characteristics within a regional flora spanning large environmental gradients. \emph{Journal of Biogeography}, 39: 204-214.
#' 
#' De Bello F., Lavergne S., Meynard C.N., Leps J. & Thuiller W. (2010): The partitioning of diversity: showing Theseus the way out of the labyrinth. \emph{Journal of Vegetation Science}, 21: 992-1000.
#' 
#' Fridley J.D., Vandermast D.B., Kuppinger D.M., Manthey M. & Peet R.K. (2007): Co-occurrence based assessment of habitat generalists and specialists: a new approach for the measurement of niche width. \emph{Journal of Ecology}, 95: 707-722.
#' 
#' Jost L. (2007): Partitioning diversity into independent alpha and beta components. \emph{Ecology}, 88: 2427-2439.
#' 
#' Manthey M. & Fridley J.D. (2009): Beta diversity metrics and the estimation of niche width via species co-occurrence data: reply to Zeleny. \emph{Journal of Ecology}, 97: 18-22.
#' 
#' Munzbergova Z. & Herben T. (2004): Identification of suitable unoccupied habitats in metapopulation studies using co-occurrence of species. \emph{Oikos}, 105: 408-414.
#' 
#' Zeleny D. (2009): Co-occurrence based assessment of species habitat specialization is affected by the size of species pool: reply to Fridley et al. (2007). \emph{Journal of Ecology}, 97: 10-17.
#' @examples
#' require (simcom)
#' sc <- sample.comm (simul.comm (S = 100), Np = 100)
#' niches <- sc$args.simcom$niche
#' additive <- calculate.theta (sc$a.mat, method = 'add')
#' multi <- calculate.theta (sc$a.mat, method = 'multiplicative')
#' beals <- calculate.theta (sc$a.mat, method = 'beals')
#' pairs (cbind (niches = niches[names (niches) %in% additive$sci.name], 
#'  additive = additive$theta, multi = multi$theta, beals = beals$theta))
#' @importFrom vegetarian d
#' @import RColorBrewer
#' @import simcom
#' @import ade4
#' @import parallel
#' @import stats
#' @import graphics
#' @import utils
#' @rdname calculate.theta
#' @export
calculate.theta <- function (comm, species.data = NULL, thresh = 5, psample = 5, reps = 10, method = "multiplicative", q = 0, rarefaction = TRUE, beals.file = NULL, pa.transform = FALSE, force.subsample = FALSE, parallel = FALSE, no.cores = 2, remove.out = F, out.metric = 'euclidean') 
{
  METHODS <- c('additive', 'multiplicative', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'rao', 'beals')
  method.n <- pmatch(method, METHODS)
  if (is.na (method.n)) stop ('Invalid name of beta-diversity algorithm')
  if (method.n == -1) stop ('The method values equals to NA') 
  method <- METHODS[method.n]
  if (is.na (reps) || reps < 2) 
    stop ("Number of random subsamples must be integer >= 2")
  if (is.na (thresh) || thresh < 2) 
    stop ("Minimum frequency of species must be integer >= 2")
  if (thresh < psample) 
    stop ("Minimum frequency of species must be >= size of the random subsamples")
  
  
  if (!is.matrix (comm)) comm <- as.matrix (comm)  # if comm is dataframe, changes into matrix
  if (is.null (row.names (comm))) row.names (comm) <- seq (1, nrow (comm))  # if comm has no row.names, these are created as sequence of integers

if (method %in% c('additive', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'beals') || (method %in% 'multiplicative' & q == 0)) pa.transform <- TRUE
if (pa.transform) comm <- ifelse (comm > 0, 1, 0)

  # For which species to calculate theta metric:
  Nplots <- nrow (comm)
  plots.per.spp <- colSums (comm > 0)  # uses only presence-absence data, since it needs count of plots per species, not sum of abundances
  select.spp <- plots.per.spp[plots.per.spp >= thresh]
  Nspp <- length (select.spp)
  
  
# For beals method, it transforms comm into beals smoothed form:  
  if (method == "beals")
  {
    if (is.null (beals.file))
    {
      beals.matrix <- beals.2 (comm, include = T)
      txt.pb <- txtProgressBar (min = 1, max = ncol (comm), initial = 0, style = 3) 
      for (co in seq (1, ncol (comm)))
      {
        if (sum (comm[,co]) > 0)
        {
          beals.temp <- beals.matrix[,co][as.logical (comm[,co])]
          stats.temp <- fivenum (beals.temp)
          iqr <- diff (stats.temp [c(2,4)])
          beals.thresh <- min (beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
          beals.matrix[,co] <- as.numeric (beals.matrix[,co] >= beals.thresh)  
        } else beals.matrix[,co] <- 0
      setTxtProgressBar (txt.pb, co)
      }
      close (txt.pb)
      } else  {
      beals.matrix <- as.matrix (read.delim (file = beals.file, row.names = 1, check.names = F))
      if (!all (dim (beals.matrix) == dim (comm))) stop ('Beals matrix has different size than species matrix!')
      }
    }

  if (!parallel) 
  {
    txt.pb <- txtProgressBar (min = 1, max = Nspp, initial = 0, style = 3) 
    temp.res <- lapply (1:Nspp, FUN = function (sp)
    {
      if (method == 'beals') temp.matrix <- beals.matrix[comm [,colnames (comm) == names (select.spp[sp])]>0,] else 
        temp.matrix <- comm[comm [,colnames (comm) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, out.metric = out.metric, thresh = thresh, psample = psample, reps = reps, method = method, q = q, rarefaction = rarefaction, force.subsample = force.subsample, parallel = parallel)
    })
      
    close (txt.pb)
  }
  
  if (parallel)
  {
    workers <- makeCluster (no.cores)
    if (file.exists ('GS-progress.txt')) file.remove ('GS-progress.txt')
    clusterExport (workers, c('calculate.theta.0', 'comm', 'select.spp', 'remove.out', 'thresh', 'psample', 'reps', 'method', 'parallel'), envir = environment ())
    temp.res <- parLapply (workers, 1:Nspp, fun = function (sp) 
    {
      if (method == 'beals') temp.matrix <- beals.matrix[comm [,colnames (comm) == names (select.spp[sp])]>0,] else 
        temp.matrix <- comm[comm [,colnames (comm) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, out.metric = out.metric, thresh = thresh, psample = psample, reps = reps, method = method, q = q, rarefaction = rarefaction, force.subsample = force.subsample, parallel = parallel) 
    }
    )
    stopCluster (workers)
  }

  theta.out <- do.call (rbind.data.frame, temp.res)
  rownames (theta.out) <- NULL
  
  names (theta.out) <- c ('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'theta', 'theta.sd')
  theta.out$sci.name <- as.character (theta.out$sci.name)  # otherwise this column would be factor, which may cause troubles
  if (!is.null (species.data)) theta.out <- as.data.frame (cbind (sci.name = theta.out[,1], species.data[as.character (theta.out[,'sci.name']),1:2], theta.out[,-1]))
  return (theta.out)
}

#' @name calculate.theta
#' @export
#' 
calculate.theta.0 <- function (temp.matrix, sci.name, sp, remove.out, out.metric, thresh, psample, reps, method, rarefaction, q, pa.transform, force.subsample, parallel)
{
  if (parallel) write (paste (sp, '\n'), file = 'GS-progress.txt', append = T)
  
  #performs outlier analysis sensu Botta-Dukat (2012):  
  if (remove.out)
  {
    if (out.metric == 'euclidean') veg.dist <- as.matrix (dist (temp.matrix))
    if (out.metric == 'sorensen') veg.dist <- as.matrix (vegan::vegdist (temp.matrix > 0))
    if (out.metric == 'binary.euclidean') veg.dist <- as.matrix (dist (temp.matrix > 0))
    diag (veg.dist) <- NA
    distances <- rowMeans (veg.dist, na.rm = T)
    outliers <- distances > (mean (distances) + 2*sd (distances))
    temp.matrix <- temp.matrix[!outliers,]
    temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  }
  
  # first method - use subsampling
  if (method %in% c('additive', 'multiplicative', 'multi.sorensen', 'multi.simpson', 'beals', 'rao') & !(method %in% c('multiplicative', 'beals') & rarefaction) || (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson') & force.subsample))
  {
    if (!nrow (temp.matrix) < thresh)  
    {
    rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = nrow (temp.matrix), byrow = F)
    sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
    
    mc.mat <- array(0,dim=c(psample,ncol (temp.matrix),reps))  
    for(i in 1:reps) mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
    total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
    mean.alpha <- colMeans (apply (mc.mat > 0, c(1,3), sum))
    
    if (method == "multiplicative") if (q == 0) Wbeta.vec <- total.rich/mean.alpha else Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) vegetarian::d (mc.mat[,,i], lev = 'beta', q = q)))  # generalized Whittaker's multiplicative beta - for q = 0 it's classical Whittaker
    if (method == "beals") Wbeta.vec <- total.rich/mean.alpha
    if (method == "additive") Wbeta.vec <- total.rich-mean.alpha 
    if (method == "pairwise.jaccard") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'jaccard')$beta.jac)))
    if (method == "pairwise.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sor)))
    if (method == "pairwise.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sim)))
    if (method == "multi.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SOR))
    if (method == "multi.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SIM))
    if (method == "rao") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) cati::RaoRel (t (mc.mat[,,i]), dfunc = NULL, dphyl = NULL, Jost = TRUE)$TD$Beta_prop))

    theta <- mean(Wbeta.vec)      #mean beta diversity value for all reps (= theta metric)
    theta.sd <- sd(Wbeta.vec)			#s.d. of above
    meanco <- mean(total.rich)			#mean cooccurrences in "psample" plots
    meanco.sd <- sd(total.rich)		#s.d. of above
    
    #sci.name <- sci.name	#scientific name
    local.avgS <- mean(mean.alpha)				#approximate mean local richness
    occur.freq <- nrow (temp.matrix)							#total number of plots
    
#    meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
#    meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
    result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd)
    return (result)
    }
  }
  # second method - not to use subsampling (only for subset of methods which are not dependent on sample size)
  if (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson') & !force.subsample)
  {
    if (!nrow (temp.matrix) < thresh)  
    {
      total.rich <- sum (colSums (temp.matrix) > 0)
      mean.alpha <- mean (rowSums (temp.matrix > 0))
      
      if (method == "pairwise.jaccard") theta <- mean (betapart::beta.pair (temp.matrix, index = 'jaccard')$beta.jac)
      if (method == "pairwise.sorensen") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sor)
      if (method == "pairwise.simpson") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sim)
      meanco <- total.rich			#mean cooccurrences in "psample" plots
      sci.name <- sci.name	#scientific name
      local.avgS <- mean.alpha				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)							#total number of plots
      
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd = 0)
      return (result)
    }
  }
  # third method - only for multiplicative with q = 0 and beals - use beta diversity rarefaction to calculate mean true beta
  if (method %in% c('multiplicative', 'beals') & rarefaction)
  {
    if (!nrow (temp.matrix) < thresh)  
    {
      theta <- beta.raref (comm = temp.matrix, sites = psample)  # contains also sd
      total.rich <- sum (colSums (temp.matrix) > 0)
      meanco <- total.rich			#number of co-occurring species in the subset
      local.avgS <- theta$alpha				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta$beta, theta$beta.sd)
      return (result)
    }
  }
}

#' @name calculate.theta
#' @export
beals.2 <- function (x, include = TRUE) # method of beals from vegan, for only p/a data
{
  x <- as.matrix(x)
  x [x > 0] <- 1
  refX <- x
  incSp <- include
  refX <- as.matrix(refX)
  M <- crossprod(refX, refX)
  C <- diag(M)
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
  S <- rowSums(x)
  b <- x
  for (i in 1:nrow(x))
    b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  SM <- rep(S, ncol(x))
  if (!incSp) SM <- SM - x
  b <- b/replace(SM, SM == 0, 1)
  b
}