#' Calculates species habitat specialization using co-occurrence based \emph{theta} metric.
#' @author David Zeleny (zeleny.david@@gmail.com). Partly based on codes written by Jason Fridley (Fridley et al. 2007) and David Zeleny (Zeleny 2009), extended for other published algorithms and optimised for speed and applicability on large datasets. Function \code{beals.2} is based on function \code{beals} from \code{vegan}, written by Miquel De Caceres and Jari Oksanen.
#' @param input.matrix Community data (\code{matrix} or \code{data.frame}, samples x species). If data are not presence-absence, the matrix will be automatically transformed into presence-absence and warning will be printed.
#' @param species.data Species data (\code{matrix} or \code{data.frame}). If suplied, it should have at least two columns - the first containing species name, the second containing layer. 
#' @param psample Minimal frequency of species. Habitat specialization will be calculated for species occurring in number of samples equal or higher than minimal frequency threshold. Default = \code{5}.
#' @param reps Number of random subsamples. Specifies how many times the fixed number of samples (specified by \code{psample}) will be randomly drawn from all samples containing target species. Default = \code{10}.
#' @param method Beta-diversity algorithm used to calculate theta measure. Currently available are \code{'additive'}, \code{'multiplicative'}, \code{'pairwise.jaccard'}, \code{'pairwise.sorensen'}, \code{'pairwise.simpson'}, \code{'multi.sorensen'}, \code{'multi.simpson'}, \code{'rao'}, \code{'beals'} and \code{'beta.div'}). See Details for available options.
#' @param beta.div.method Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param beta.div.sqrt.D Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param beta.div.samp Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param beals.file Contains pre-calculated matrix of species co-occurrences. Can be used if \code{method = 'beals'} to speed up repeated calculation.
#' @param pa.transform Logical; should the compositional data be transformed into presence-absence form? This choice applies only if \code{method} is \code{"rao"} or \code{"beta.div"}, since the other methods must be calculated on presence-absence data only (for these methods the matrix is automatically transformed into p/a form).
#' @param force.subsample Logical; should the subsampling be forced even for beta diversity metrics which are not influenced by sample size (\code{c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'rao', 'beta.div')})? Default behaviour is \code{force.subsample = FALSE}, meaning that beta diversity metrics dependent on sample size (\code{method = c('additive', 'multiplicative', 'beals', 'multi.sorensen', 'multi.simpson')}) will subsample number of plots equal to \code{psample}, while method not dependent on sample size (\code{method = c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'beta.div')}) will be calculated using all plots containing target species. If \code{force.subsample = TRUE}, even methods not dependent on sample size will be calculated using subsampling. 
#' @param parallel Logical; should be the parallel calculation used?
#' @param no.cores Number of cores (if \code{parallel = TRUE}). Note that in case of large datasets the calculation may be limited by RAM of the computer, and increasing number of cores may result in saturation of RAM and calculation collapse.
#' @param remove.out Logical; should be the algorithm removing outliers (sensu Botta-Dukat 2012) applied? 
#' @param verbal Logical; if \code{TRUE}, tcltk progress bar will popup during the calculation.
#' @param juicer Logical argument specific for launching the function from JUICE software; logical (default = F) - is the function launched from JUICE? If \code{juicer = TRUE}, function is expecting that \code{species.data} have JUICE-specific structure, which enables to import data back to JUICE.
#' @param tcltk Logical argument specific for launching the function from JUICE sofware.
#' @param temp.matrix Internal argument; matrix with species composition of plots containing target species.
#' @param sp Internal argument; the order of the species for which the current calculation is done.
#' @param sci.name Internal argument; the name of the species for which the current calculation is done.
#' @param win.pb Internal argument.
#' @param x Internal argument of \code{beals.2} function; input compositional matrix.
#' @param include Internal argument of \code{beals.2} function; include argument from \code{vegan::beals}.

#' @details
#' Function \code{calculate.theta} calculates theta metric of species habitat specialization using range of proposed beta diversity measures; it uses internal functions \code{calculate.theta.0} and \code{beals.2} (modified from the library \code{vegan}) . Function \code{calculate.theta.tcltk} launches tcltk clickable interface, which enables to select methods and parameters used for calculation; this function is primarily used to be lounched externally, e.g. from JUICE program.
#' 
#' The function \code{calculate.theta} offers the following \code{method} argument to calculate beta diversity among samples:
#' \itemize{
#' \item \code{additive}: This is the original algorithm published by Fridley et al. (2007), in which beta diversity among samples containing given species is calculated by additive beta diversity measure.
#' \item \code{multiplicative}: Uses the multiplicative Whittaker's measure of beta diversity instead of the original additive measure, as suggested by Zeleny (2009).
#' \item \code{beals}: Multiplicative beta on species pool. Algorithm suggested by Botta-Dukat (2012), calculating the beta diversity using species pool matrix instead of the original species data matrix. Species pool matrix is calculated using Beals smoothing method (invented by Ewald 2002). While the previous multiplicative beta diversity method gives unbiased results only in case of not-saturated communities, this method should give unbiased results also in case of saturated communities. See Zeleny (2009) and Botta-Dukat (2012) for detail discussion of this saturated/not-saturated communities issue.
#' \item \code{pairwise.jaccard}, \code{pairwise.sorensen}, \code{pairwise.simpson}, \code{multi.sorensen} and \code{multi.simpson}: Mean pairwise Jaccard, Sorensen and Simpson dissimilarity, and multiple Sorensen and Simpson dissimilarity based on reccomendations of Manthey & Fridley (2009). Authors suggested that neither the original additive algorithm (introduced by Fridley et al. 2007), neither the modified version using the multiplicative beta diversity (Zeleny 2009) is the best solution, and introduced other alternatives, using pairwise or multiple site beta diversity algorithm. Mean pairwise Jaccard dissimilarity (or Sorensen and Simpson, respectively) is based on calculating mean of Jaccard (or Sorensen and Simpson, respectively) dissimilarities among all pairs of samples in each subset, while multiple Sorensen (or Simpson, respectively) is using multiple-site Sorensen (or Simpson, respectively) algorithm introduced by Baselga et al. (2007) 
#' \item \code{rao}: Rao index of dissimilarity; this option has been introduced and used by Boulangeat et al. (2012). Advantage of Rao index is a possibility to incorporate known relationships among species using the among-species distance matrix. If this is not supplied and Rao index is calculated using presence-absence data, it is equal to Gini-Simpson diversity index and results are similar to pairwise Jaccard index.
#' \item \code{beta.div}: Calculating the beta diversity as the variation in community matrix, using the concept introduced by Legendre & De Caceres (2013) and function \code{beta.div} written by Pierre Legendre (see Appendix S4.r of Legendre & De Caceres 2013; the version used here is from Legendre's website http://adn.biol.umontreal.ca/~numericalecology/labo/fonctions_r/beta-diversity.zip; however, I keep using \code{method = "percentagedifference"} instead of "\%difference" from the new version, since I had a problem with roxygenizing the R help files with \% sign). Three additional arguments can be specified if \code{method = "beta.div"}, namely \code{beta.div.method}, \code{beta.div.sqrt.D} and \code{beta.div.samp} (the original arguments in the function \code{beta.div} are \code{method}, \code{sqrt.D} and \code{samp}). \code{beta.div.method} is choosing one of 21 distance metrics (from \code{c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "percentagedifference", "ruzicka", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai")}). Argument \code{beta.div.sqrt.D} (logical) decides whether square root of distance should be used for calculation (important for non-euclidean distances like Bray-Curtis, called \code{"percentagedifference"} in \code{beta.div} function). Argument \code{beta.div.samp} is logical; if \code{beta.div.samp = TRUE}, the abundance-based distances (\code{c("ab.jaccard", "ab.sorensen", "ab.ochiai", "ab.simpson")}) are computed for sample data. If \code{beta.div.samp= FALSE}, they are computed for true population data.
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
#' @examples
#' sc <- sample.comm (simul.comm (totS = 100), Np= 100)
#' niches <- sc$simul.comm$range
#' additive <- calculate.theta (sc$a.mat, method = 'add')
#' multi <- calculate.theta (sc$a.mat, method = 'multiplicative')
#' beals <- calculate.theta (sc$a.mat, method = 'beals')
#' bray <- calculate.theta (sc$a.mat, method = 'beta.div', beta.div.method = 'percentagedifference', beta.div.sqrt.D = TRUE)
#' # Visualize the relationship using function pairs with Spearmann's correlation in the boxes above diagonal (see Examples in ?pairs)
#' panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
#' {
#'   usr <- par("usr"); on.exit(par(usr))
#'   par(usr = c(0, 1, 0, 1))
#'   r <- abs(cor(x, y, method = 'spearman'))
#'   txt <- format(c(r, 0.123456789), digits = digits)[1]
#'   txt <- paste0(prefix, txt)
#'   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#'   text(0.5, 0.5, txt, cex = cex.cor * r)
#' }
#'pairs (cbind (niches = niches[names (niches) %in% additive$sci.name], additive = additive$theta, multi = multi$theta, beals = beals$theta, bray = bray$theta), upper.panel = panel.cor)

#' @rdname calculate.theta
#' @export
calculate.theta <- function (input.matrix, species.data = NULL, psample = 5, reps = 10, method = "multiplicative", beta.div.method = 'hellinger', beta.div.sqrt.D = FALSE, beta.div.samp = TRUE, beals.file = NULL, pa.transform = FALSE, force.subsample = FALSE, parallel = FALSE, no.cores = 2, remove.out = F, verbal = F, juicer = F, tcltk = F) 
{
  require (tcltk)
  METHODS <- c('additive', 'multiplicative', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'rao', 'beals', 'beta.div')
  method.n <- pmatch(method, METHODS)
  if (is.na (method.n)) stop ('invalid method')
  if (method.n == -1) stop ('ambiguous method')
  method <- METHODS[method.n]
  if (!verbal) win.pb <- NULL
  if (is.na (reps) || reps < 2) 
    if (verbal) {tkmessageBox (type = "ok", message = "Number of random subsamples must be integer >= 2"); stop ()} else stop ("Number of random subsamples must be integer >= 2")
  if (is.na (psample) || psample < 2) 
    if (verbal) {tkmessageBox (type = "ok", message = "Minimal frequency of species must be integer >= 2"); stop ()} else stop ("Minimal frequency of species must be integer >= 2")

  if (!is.matrix (input.matrix)) input.matrix <- as.matrix (input.matrix)  # if input.matrix is dataframe, changes into matrix
  if (is.null (row.names (input.matrix))) row.names (input.matrix) <- seq (1, nrow (input.matrix))  # if input.matrix has no row.names, these are created as sequence of integers

if (method %in% c('additive', 'multiplicative', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'beals')) pa.transform <- TRUE
if (pa.transform) input.matrix <- ifelse (input.matrix > 0, 1, 0)

  # For which species to calculate theta metric:
  Nplots <- nrow (input.matrix)
  plots.per.spp <- colSums (input.matrix > 0)  # uses only presence-absence data, since it needs count of plots per species, not sum of abundances
  select.spp <- plots.per.spp[plots.per.spp >= psample]
  Nspp <- length (select.spp)
  
  
# For beals method, it transforms input.matrix into beals smoothed form:  
  if (method == "beals")
  {
    if (is.null (beals.file))
    {
      beals.matrix <- beals.2 (input.matrix, include = T, verbal = verbal)
      if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix), initial = 0, width = 300) 
      for (co in seq (1, ncol (input.matrix)))
      {
        if (verbal) setWinProgressBar (win.pb, co + 1, label = paste ("Prepare beals smoothed table: species ", co))
        if (sum (input.matrix[,co]) > 0)
        {
          beals.temp <- beals.matrix[,co][as.logical (input.matrix[,co])]
          stats.temp <- fivenum (beals.temp)
          iqr <- diff (stats.temp [c(2,4)])
          beals.thresh <- min (beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
          beals.matrix[,co] <- as.numeric (beals.matrix[,co] >= beals.thresh)  
        } else beals.matrix[,co] <- 0
        
      }
      if (verbal) close (win.pb)
      if (tcltk) write.table (beals.matrix, file = 'beals-data.txt', sep = '\t', col.names = TRUE)
    } else  
    {
      if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix)+1, initial = 0, width = 300) 
      if (verbal) setWinProgressBar (win.pb, 1, label = "Reading Beals smoothed table")
      beals.matrix <- as.matrix (read.delim (file = beals.file, row.names = 1, check.names = F))
      if (! all (dim (beals.matrix) == dim (input.matrix))) {tkmessageBox (type = "ok", message = paste ("Selected Beals matrix has different size than species matrix! \nYou need to calculate new beals smoothed species pool data using current species data. Close the JUICE-R application and run it again from JUICE, and calculate the Multiplicative beta on species pool analysis without selecting the Beals smoothing table.")); stop ('Beals matrix has different size than species matrix!')}
#      input.matrix <- beals.matrix
      if (verbal) close (win.pb)
    }
  }

  if (!parallel) 
  {
    if (verbal) win.pb <- winProgressBar (title = "Calculation progress", label = paste ("Species no. ", 1), min = 1, max = Nspp, initial = 0, width = 300) 
    temp.res <- lapply (1:Nspp, FUN = function (sp)
    {
      if (method == 'beals') temp.matrix <- beals.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,] else 
        temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, psample = psample, reps = reps, method = method, beta.div.method = beta.div.method, beta.div.sqrt.D = beta.div.sqrt.D, beta.div.samp = beta.div.samp, force.subsample = force.subsample, parallel = parallel, win.pb = win.pb, verbal = verbal, juicer = juicer)
    })
      
    if (verbal) close (win.pb)
  }
  
  if (parallel)
  {
    require (parallel)
    workers <- makeCluster (no.cores)
    if (verbal) if (file.exists ('GS-progress.txt')) file.remove ('GS-progress.txt')
    clusterExport (workers, c('calculate.theta.0', 'input.matrix', 'select.spp', 'remove.out', 'psample', 'reps', 'method', 'parallel'), envir = environment ())
    temp.res <- parLapply (workers, 1:Nspp, fun = function (sp) 
    {
      if (method == 'beals') temp.matrix <- beals.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,] else 
        temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, psample = psample, reps = reps, method = method, beta.div.method = beta.div.method, beta.div.sqrt.D = beta.div.sqrt.D, beta.div.samp = beta.div.samp, force.subsample = force.subsample, parallel = parallel, win.pb = NULL, verbal = verbal, juicer = juicer) 
    }
    )
    stopCluster (workers)
  }

  theta.out <- do.call (rbind.data.frame, temp.res)
  rownames (theta.out) <- NULL
  
  if (ncol (theta.out) == 9) names (theta.out) <- c ('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'meanco.sd', 'meanco.u', 'meanco.l', 'theta', 'theta.sd') else names (theta.out) <- c ('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'theta')
  theta.out$sci.name <- as.character (theta.out$sci.name)  # otherwise this column would be factor, which may cause troubles
  if (!is.null (species.data)) theta.out <- as.data.frame (cbind (sci.name = theta.out[,1], species.data[as.character (theta.out[,'sci.name']),1:2], theta.out[,-1]))
  if (juicer) write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = T, row.names = F) else return (theta.out)
      
  if (juicer) write.table (file = "theta_import.species.data.via.clipboard.txt", theta.out[,c('full.sci.name', 'layer', 'theta')], quote = F, row.names = F, col.names = F, sep = '\t')
  if (juicer) write.table (file = "clipboard", theta.out[,c('full.sci.name', 'layer', 'theta')], quote = F, row.names = F, col.names = F, sep = '\t')
      
  if (verbal) tkmessageBox (type = "ok", message = paste ("Species theta values have been copied into clipboard - you can import them directly into JUICE (Edit > Paste Clipboard to BLACK species names).\n\nResult files were saved into", getwd (), "\n\nYou can also use the file theta_import.species.data.via.clipboard.txt to import the species theta values to JUICE (Edit > Paste Clipboard to BLACK species names)."))
}

#' @name calculate.theta
#' @export
#' 
calculate.theta.0 <- function (temp.matrix, sci.name, sp, remove.out, psample, reps, method, beta.div.method, beta.div.sqrt.D, beta.div.samp, force.subsample, parallel, win.pb, verbal, juicer)
{
  if (verbal) if (parallel) write (paste (sp, '\n'), file = 'GS-progress.txt', append = T) else setWinProgressBar (win.pb, sp, label = paste ("Species no. ", sp))
  
  #performs outlier analysis sensu Botta-Dukat (2012):  
  if (remove.out)
  {
    veg.dist <- as.matrix (dist (temp.matrix))
    diag (veg.dist) <- NA
    distances <- rowMeans (veg.dist, na.rm = T)
    outliers <- distances > (mean (distances) + 2*sd (distances))
    temp.matrix <- temp.matrix[!outliers,]
    temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  }
  
  # first method - use subsampling
  
  if (method %in% c('additive', 'multiplicative', 'multi.sorensen', 'multi.simpson', 'beals')||(method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'rao', 'beta.div') & force.subsample))
  {
    if (!nrow (temp.matrix) < psample)  
    {
    rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = nrow (temp.matrix), byrow = F)
    sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
    
    mc.mat <- array(0,dim=c(psample,ncol (temp.matrix),reps))  
    for(i in 1:reps) mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
    total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
    mean.alpha <- colMeans (apply (mc.mat > 0, c(1,3), sum))
    
    if (method == "multiplicative" | method == "beals") Wbeta.vec <- total.rich/mean.alpha
    if (method == "additive") Wbeta.vec <- total.rich-mean.alpha 
    if (method == "pairwise.jaccard") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'jaccard')$beta.jac)))
    if (method == "pairwise.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sor)))
    if (method == "pairwise.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sim)))
    if (method == "multi.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SOR))
    if (method == "multi.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SIM))
    if (method == "rao") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (ade4::disc (as.data.frame (t (mc.mat[,,i]))))))
    if (method == "beta.div") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) beta.div (mc.mat[,,i], method = beta.div.method, sqrt.D = beta.div.sqrt.D, nperm = 0)$SStotal_BDtotal[2]))
    
    theta <- mean(Wbeta.vec)      #mean beta diversity value for all reps (= theta metric)
    theta.sd <- sd(Wbeta.vec)			#s.d. of above
    meanco <- mean(total.rich)			#mean cooccurrences in "psample" plots
    meanco.sd <- sd(total.rich)		#s.d. of above
    
    #sci.name <- sci.name	#scientific name
    local.avgS <- mean(mean.alpha)				#approximate mean local richness
    occur.freq <- nrow (temp.matrix)							#total number of plots
    
    meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
    meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
    result <- list(sci.name, local.avgS, occur.freq, meanco, meanco.sd, meanco.u, meanco.l, theta, theta.sd)
    return (result)
    }
  }
  # second method - not to use subsampling (only for subset of methods which are not dependent on sample size)
  if (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'rao', 'beta.div') & !force.subsample)
  {
    if (!nrow (temp.matrix) < psample)  
    {
      total.rich <- sum (colSums (temp.matrix) > 0)
      mean.alpha <- mean (rowSums (temp.matrix > 0))
      
      if (method == "pairwise.jaccard") theta <- mean (betapart::beta.pair (temp.matrix, index = 'jaccard')$beta.jac)
      if (method == "pairwise.sorensen") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sor)
      if (method == "pairwise.simpson") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sim)
      if (method == "rao") theta <- mean (ade4::disc (as.data.frame (t (temp.matrix))))
      if (method == "beta.div") theta <- beta.div (temp.matrix, method = beta.div.method, sqrt.D = beta.div.sqrt.D, nperm = 0)$SStotal_BDtotal[2]
      
      meanco <- total.rich			#mean cooccurrences in "psample" plots
      sci.name <- sci.name	#scientific name
      local.avgS <- mean.alpha				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)							#total number of plots
      
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta)
      return (result)
    }
  }
}

#' @name calculate.theta
#' @export
calculate.theta.tcltk <- function (input.matrix, species.data = NULL, juicer = T)
{
  require (tcltk)
  cancel <- tclVar (0)
  end.end <- F
  beals.file <- NULL
  
  GSmethods <- c ("Additive beta diversity (Fridley et al. 2007)", "Multiplicative beta diversity (Zeleny 2009)", "Multiplicative beta on species pool (Botta-Dukat 2012)", "Pairwise Jaccard dissimilarity (Manthey & Fridley 2009)", "Multiple Sorensen dissimilarity (Manthey & Fridley 2009)", "Multiple Simpson dissimilarity (Manthey & Fridley 2009)", "Rao index of dissimilarity (Boulangeat et al. 2012)")
  base <- tktoplevel()
  tkwm.title(base, "Generalists-specialists")
  
  spec.frm <- tkframe (base, borderwidth=2)
  frame.title <- tkframe (spec.frm, relief = 'groove', borderwidth = 2, background = 'grey')
  frame.a <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.b <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.c <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.d <- tkframe (spec.frm, borderwidth = 2)
  frame.e <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.e.1 <- tkframe (frame.e)
  frame.e.2 <- tkframe (frame.e)
  frame.f <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.f.1 <- tkframe (frame.f)
  frame.g <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  
  
  GSmethod <- tclVar ("additive")
  label.radio <- tklabel (frame.a, text = "Which beta diversity algorithm to use?")
  radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = "additive", variable = GSmethod)
  radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = "multiplicative", variable = GSmethod)
  radio3 <- tkradiobutton (frame.a, text = GSmethods[3], value = "beals", variable = GSmethod)
  radio4 <- tkradiobutton (frame.a, text = GSmethods[4], value = "pairwise.jaccard", variable = GSmethod)
  radio5 <- tkradiobutton (frame.a, text = GSmethods[5], value = "multi.sorensen", variable = GSmethod)
  radio6 <- tkradiobutton (frame.a, text = GSmethods[6], value = 'multi.simpson', variable = GSmethod)
  radio7 <- tkradiobutton (frame.a, text = GSmethods[7], value = 'rao', variable = GSmethod)
  radio8 <- tkradiobutton (frame.a, text = GSmethods[7], value = 'beta.div', variable = GSmethod)
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  parallel <- tclVar (0)
  no.cores <- tclVar (2)
  remove.out <- tclVar (0)
  label.entry1 <- tklabel (frame.b, text = "Minimal frequency of species ")
  entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)
  
  label.entry2 <- tklabel (frame.c, text = "Number of random subsamples ")
  entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = "Calculate", width = 10, height = 2, command = function () calculate.theta (input.matrix = input.matrix, species.data = species.data, psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), method = as.character (tclvalue (GSmethod)), beals.file = beals.file, parallel = as.logical (as.numeric (tclvalue (parallel))), no.cores = as.numeric (tclvalue (no.cores)), remove.out = as.logical (as.numeric (tclvalue (remove.out))), verbal = T, juicer = T, tcltk = T))
  
  
  choose.label <- tklabel (frame.e.2, text = 'Select the file with beals smoothed data')
  choose.button <- tkbutton (frame.e.1, text = 'Select', command = function () assign ('beals.file', choose.files (), inherits = T))
  tkpack (choose.button)
  tkpack (choose.label)
  tkpack (tklabel (frame.e, text = 'Beals smoothing (included in method of Botta-Dukat 2012)'), anchor = 'w')
  tkpack (frame.e.1, frame.e.2, side = 'left',ipady = 5, ipadx = 5, padx = 5, pady = 5)
  
  
  parallel.label <- tklabel (frame.f, text = 'Parallel calculation (enable only if you have more than one core)')
  parallel.no.cores.label <- tklabel (frame.f.1, text = 'number of cores: ')
  parallel.no.cores.entry <- tkentry (frame.f.1, width = 2, textvariable = no.cores)
  parallel.checkbutton <- tkcheckbutton (frame.f.1, text = 'use parallel computing,', variable = parallel)
  
  tkpack (tklabel (frame.g, text = 'Outlier analysis (McCune & Mefford 1999, suggested by Botta-Dukat 2012)'), tkcheckbutton (frame.g, text = 'remove outlier samples (with very different species composition)', variable = remove.out), anchor = 'w')
  
  tkpack (label.radio, radio1, radio2, radio4, radio5, radio6, radio7, radio3, anchor = 'w')
  tkpack (label.entry1, entry1, anchor = 'w', side = 'left')
  tkpack (label.entry2, entry2, anchor = 'w', side = 'left')
  tkpack (button1)
  tkpack (parallel.checkbutton, parallel.no.cores.label, parallel.no.cores.entry, side = 'left')
  tkpack (parallel.label,  frame.f.1, anchor = 'w')
  
  tkpack (tklabel (frame.title, text = paste ('Calculation of generalists and specialists using co-occurrence species data \n Author: David Zeleny (zeleny.david@gmail.com)', if (juicer) '\n JUICE-R application (www.bit.ly/habitat-specialists)', '\n Version of library genspe: ', as.character (packageVersion ('genspe')), '\nNumber of samples: ', dim (input.matrix)[1], ', number of species: ', dim (input.matrix)[2], sep = '')), ipady = 10, ipadx = 10, padx = 10, pady = 10)
  
  tkpack (frame.title, side = 'top', expand = T, fill = 'both')
  tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.e, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.f, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.g, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.b, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10, expand = T, fill = 'both')
  tkpack (frame.d, side = "bottom", pady = 10, padx = 10, expand = T, fill = 'both')
  
  tkpack (spec.frm)
  tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2)  
  
  tkraise (base)
  tkwait.variable (cancel)
}

#' @name calculate.theta
#' @export
beals.2 <- function (x, include = TRUE, verbal = FALSE) # method of beals from vegan, for only p/a data and with progress bar
{
  if (verbal) win.pb2 <- winProgressBar (title = 'Beals', label = 'start', min = 1, max = nrow (x)+2, initial = 0, width = 300)
  x <- as.matrix(x)
  x [x > 0] <- 1
  refX <- x
  incSp <- include
  refX <- as.matrix(refX)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'Crossprod')
  M <- crossprod(refX, refX)
  C <- diag(M)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'First sweep')
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
  S <- rowSums(x)
  b <- x
  for (i in 1:nrow(x)) {
    if (verbal) setWinProgressBar (win.pb2, i+1, label = i)
    b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  }                       
  SM <- rep(S, ncol(x))
  if (!incSp) SM <- SM - x
  b <- b/replace(SM, SM == 0, 1)
  if (verbal) close (win.pb2)
  b
}

beta.div <- function(Y, method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
  #
  # Compute estimates of total beta diversity as the total variance in Y, 
  # for 21 dissimilarity coefficients or analysis of raw data (not recommended). 
  # LCBD indices are tested by permutation within columns of Y.
  # This version includes direct calculation of the Jaccard, Sorensen and Ochiai 
  # coefficients for presence-absence data, done by package ade4.
  #
  # Arguments --
  # 
  # Y : community composition data matrix.
  # method : name of one of the 20 dissimilarity coefficients, or "none" for
#          direct calculation on Y (also the case with method="euclidean").
# sqrt.D : If sqrt.D=TRUE, the distances in matrix D are square-rooted before 
#          computation of SStotal, BDtotal and LCBD. 
# samp : If samp=TRUE, the abundance-based distances (ab.jaccard, ab.sorensen,
#        ab.ochiai, ab.simpson) are computed for sample data. If samp=FALSE, 
#        they are computed for true population data.
# nperm : Number of permutations for test of LCBD.
# save.D : If save.D=TRUE, the distance matrix will appear in the output list.
# clock : If clock=TRUE, the computation time is printed in the R console.
#
# Reference --
#
# Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of 
# community data: dissimilarity coefficients and partitioning. 
# Ecology Letters 16: 951-963. 
#
# License: GPL-2 
# Author:: Pierre Legendre, December 2012, April-May 2013, April 2015
{
  ### Internal functions
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {	One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  mat.cen <- mat %*% D %*% mat
  }
  ###
  BD.group1 <- function(Y, method, save.D, per, n)
  {
    if(method=="profiles") Y = decostand(Y, "total")
    if(method=="hellinger") Y = decostand(Y, "hellinger")
    if(method=="chord") Y = decostand(Y, "norm")
    if(method=="chisquare") Y = decostand(Y, "chi.square")
    #
    s <- scale(Y, center=TRUE, scale=FALSE)^2   # eq. 1
    SStotal <- sum(s)          # eq. 2
    BDtotal <- SStotal/(n-1)   # eq. 3
    if(!per) { SCBD<-apply(s,2,sum)/SStotal }else{ SCBD<-NA }  # eqs. 4a and 4b
    LCBD <- apply(s, 1, sum)/SStotal  # eqs. 5a and 5b
    #
    D <- NA
    if(!per & save.D)   D <- dist(Y)
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), SCBD=SCBD, LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  BD.group2 <- function(Y, method, sqrt.D, n)
  {
    if(method == "divergence") {
      D = D11(Y)		
      
    } else if(any(method == 
                  c("jaccard","sorensen","ochiai"))) 
    {
      if(method=="jaccard") D = dist.binary(Y, method=1) # ade4 takes sqrt(D)
      if(method=="sorensen")  D = dist.binary(Y, method=5) #ade4 takes sqrt(D)
      if(method=="ochiai") D = dist.binary(Y, method=7) # ade4 takes sqrt(D)
      
    } else if(any(method == 
                  c("manhattan","canberra","whittaker","percentagedifference","ruzicka","wishart"))) 
    {
      if(method=="manhattan") D = vegdist(Y, "manhattan")
      if(method=="canberra")  D = vegdist(Y, "canberra")
      if(method=="whittaker") D = vegdist(decostand(Y,"total"), "manhattan")/2
      if(method=="percentagedifference") D = vegdist(Y, "bray")
      if(method=="ruzicka")   D = RuzickaD(Y)
      if(method=="wishart")   D = WishartD(Y)
    } else {
      if(method=="modmeanchardiff") D = D19(Y)
      if(method=="kulczynski")  D = vegdist(Y, "kulczynski")
      if(method=="ab.jaccard")  D = chao(Y, coeff="Jaccard", samp=samp)
      if(method=="ab.sorensen") D = chao(Y, coeff="Sorensen", samp=samp)
      if(method=="ab.ochiai")   D = chao(Y, coeff="Ochiai", samp=samp)
      if(method=="ab.simpson")  D = chao(Y, coeff="Simpson", samp=samp)
    }
    #
    if(sqrt.D) D = sqrt(D)
    SStotal <- sum(D^2)/n      # eq. 8
    BDtotal <- SStotal/(n-1)   # eq. 3
    delta1 <- centre(as.matrix(-0.5*D^2), n)   # eq. 9
    LCBD <- diag(delta1)/SStotal               # eq. 10b
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  ###
  epsilon <- sqrt(.Machine$double.eps)
  method <- match.arg(method, c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "percentagedifference", "ruzicka", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai","none"))
  #
  if(any(method == c("profiles", "hellinger", "chord", "chisquare", "manhattan", "modmeanchardiff", "divergence", "canberra", "whittaker", "percentagedifference", "kulczynski"))) require(vegan)
  if(any(method == c("jaccard","sorensen","ochiai"))) require(ade4)
  #
  if(is.table(Y)) Y <- Y[1:nrow(Y),1:ncol(Y)]    # In case class(Y) is "table"
  n <- nrow(Y)
  if((n==2)&(dist(Y)[1]<epsilon)) stop("Y contains two identical rows, hence BDtotal = 0")
  #
  aa <- system.time({
    if(any(method == 
           c("euclidean", "profiles", "hellinger", "chord", "chisquare","none"))) {
      note <- "Info -- This coefficient is Euclidean"
      res <- BD.group1(Y, method, save.D, per=FALSE, n)
      #
      # Permutation test for LCBD indices, distances group 1
      if(nperm>0) {
        p <- ncol(Y)
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group1(Y.perm, method, save.D, per=TRUE, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, SCBD=res$SCBD, 
                  LCBD=res$LCBD, p.LCBD=p.LCBD, method=method, note=note, D=D)
      
    } else {
      #
      if(method == "divergence") {
        note = "Info -- This coefficient is Euclidean"
      } else if(any(method == c("jaccard","sorensen","ochiai"))) {
        note = c("Info -- This coefficient is Euclidean because dist.binary ",
                 "of ade4 computes it as sqrt(D). Use beta.div with option sqrt.D=FALSE")
      } else if(any(method == 
                    c("manhattan","canberra","whittaker","percentagedifference","ruzicka","wishart"))) {
        if(sqrt.D) {
          note = "Info -- In the form sqrt(D), this coefficient, is Euclidean"
        } else {
          note = c("Info -- For this coefficient, sqrt(D) would be Euclidean", 
                   "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
        }
      } else {
        note = c("Info -- This coefficient is not Euclidean", 
                 "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
      }
      #
      res <- BD.group2(Y, method, sqrt.D, n)
      #
      # Permutation test for LCBD indices, distances group 2
      if(nperm>0) {
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group2(Y.perm, method, sqrt.D, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(sqrt.D) note.sqrt.D<-"sqrt.D=TRUE"  else  note.sqrt.D<-"sqrt.D=FALSE"
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, LCBD=res$LCBD,  
                  p.LCBD=p.LCBD, method=c(method,note.sqrt.D), note=note, D=D)
    }
    #
  })
  aa[3] <- sprintf("%2f",aa[3])
  if(clock) cat("Time for computation =",aa[3]," sec\n")
  #
  class(out) <- "beta.div"
  out
}

RuzickaD <- function(Y)
  #
  # Compute the Ruzicka dissimilarity = (B+C)/(A+B+C) (quantitative form of Jaccard).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2015
{
  n = nrow(Y)
  mat.sq = matrix(0, n, n)
  # A = W = sum of minima in among-site comparisons
  # B = sum_site.1 - W = K[1] - W   # sum of differences for sp(site1) > sp(site2)
  # C = sum_site.2 - W = K[2] - W   # sum of differences for sp(site2) > sp(site1)
  W <- matrix(0,n,n)          # matrix that will receive the sums of minima (A)
  K <- apply(Y,1,sum)         # row sums: (A+B) or (A+C)
  for(i in 2:n) for(j in 1:(i-1)) W[i,j] <- sum(pmin(Y[i,], Y[j,])) # sums of minima (A)
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      mat.sq[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) } # (B+C)/(A+B+C)
  }
  mat = as.dist(mat.sq)
}

D11 <- function(Y, algo=1)
  #
  # Compute Clark's coefficient of divergence. This is
  # coefficient D11 in Legendre and Legendre (2012, eq. 7.51).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = no. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  if(algo==1) {   # Faster algorithm
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        num <- (Y[i,]-Y[j,])
        den <- (Y[i,]+Y[j,])
        sel <- which(den > 0)
        D[i,j] = sqrt(sum((num[sel]/den[sel])^2)/pp[i,j])
      }
    }
    #
  } else {   # Slower algorithm 
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        temp = 0
        for(p2 in 1:p) {
          den = Y[i,p2] + Y[j,p2]
          if(den > 0) {
            temp = temp + ((Y[i,p2] - Y[j,p2])/den)^2
          }
        }
        D[i,j] = sqrt(temp/pp[i,j])
      }
    }
    #
  }	
  DD <- as.dist(D)
}

D19 <- function(Y)
  #
  # Compute the Modified mean character difference. This is
  # coefficient D19 in Legendre and Legendre (2012, eq. 7.46).
  # Division is by pp = number of species present at the two compared sites
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = n. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  D <- vegdist(Y, "manhattan")
  DD <- as.dist(as.matrix(D)/pp)
}

WishartD <- function(Y)
  #
  # Compute dissimilarity = (1 - Wishart similarity ratio) (Wishart 1969).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, August 2012
{
  CP = crossprod(t(Y))
  SS = apply(Y^2,1,sum)
  n = nrow(Y)
  mat.sq = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(i-1)) { mat.sq[i,j] = CP[i,j]/(SS[i] + SS[j] - CP[i,j]) }
  }
  mat = 1 - as.dist(mat.sq)
}

chao <- function(mat, coeff="Jaccard", samp=TRUE)
  #
  # Compute Chao et al. (2006) abundance-based indices.
  #
  # Arguments -
  # mat = data matrix, species abundances
  # coef = "Jaccard" : modified abundance-based Jaccard index
  #        "Sorensen": modified abundance-based Sørensen index
  #        "Ochiai"  : modified abundance-based Ochiai index
  #        "Simpson" : modified abundance-based Simpson index
  # samp=TRUE : Compute dissimilarities for sample data
  #     =FALSE: Compute dissimilarities for true population data
#
# Details -
# For coeff="Jaccard", the output values are identical to those
# produced by vegan's function vegdist(mat, "chao").
#
# Help received from A. Chao and T. C. Hsieh in July 2012 for the computation  
# of dissimilarities for true population data is gratefully acknowledged.
#
# Reference --
# Chao, A., R. L. Chazdon, R. K. Colwell and T. J. Shen. 2006. 
# Abundance-based similarity indices and their estimation when there 
# are unseen species in samples. Biometrics 62: 361–371.
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2012
{
  require(vegan)
  nn = nrow(mat)
  res = matrix(0,nn,nn)
  if(samp) {   # First for sample data
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        #cat("k =",k,"  j =",j,"\n")
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        N.j = sum(v1)   # Sum of abundances in vector 1
        N.k = sum(v2)   # Sum of abundances in vector 2
        shared.sp = v1.pa * v2.pa   # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          C.j = sum(shared.sp * v1)   # Sum of shared sp. abundances in v1
          C.k = sum(shared.sp * v2)   # Sum of shared sp. abundances in v2
          # a1.j = sum(shared.sp * v1.pa)
          # a1.k = sum(shared.sp * v2.pa)
          a1.j = length(which((shared.sp * v2) == 1)) # Singletons in v2
          a1.k = length(which((shared.sp * v1) == 1)) # Singletons in v1
          a2.j = length(which((shared.sp * v2) == 2)) # Doubletons in v2
          if(a2.j == 0) a2.j <- 1
          a2.k = length(which((shared.sp * v1) == 2)) # Doubletons in v1
          if(a2.k == 0) a2.k <- 1
          # S.j = sum(v1[which(v2 == 1)]) # Sum abund. in v1 for singletons in v2
          # S.k = sum(v2[which(v1 == 1)]) # Sum abund. in v2 for singletons in v1
          sel2 = which(v2 == 1)
          sel1 = which(v1 == 1)
          if(length(sel2)>0) S.j = sum(v1[sel2]) else S.j = 0
          if(length(sel1)>0) S.k = sum(v2[sel1]) else S.k = 0
          
          U.j = (C.j/N.j) + ((N.k-1)/N.k) * (a1.j/(2*a2.j)) * (S.j/N.j) # Eq. 11
          if(U.j > 1) U.j <- 1
          U.k = (C.k/N.k) + ((N.j-1)/N.j) * (a1.k/(2*a2.k)) * (S.k/N.k) # Eq. 12
          if(U.k > 1) U.k <- 1
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U.j*U.k/(U.j + U.k - U.j*U.k))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U.j*U.k/(U.j + U.k))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U.j*U.k))
          } else if(coeff == "Simpson") { 
            # Simpson (1943), or Lennon et al. (2001) in Chao et al. (2006)
            res[k,j] = 1 -
              (U.j*U.k/(U.j*U.k+min((U.j-U.j*U.k),(U.k-U.j*U.k))))
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
    
  } else {   # Now for complete population data
    
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        shared.sp = v1.pa * v2.pa    # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          N1 = sum(v1)   # Sum of abundances in vector 1
          N2 = sum(v2)   # Sum of abundances in vector 2
          U = sum(shared.sp * v1)/N1   # Sum of shared sp. abundances in v1
          V = sum(shared.sp * v2)/N2   # Sum of shared sp. abundances in v2
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U*V/(U + V - U*V))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U*V/(U + V))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U*V))
          } else if(coeff == "Simpson") { # "Simpson"
            res[k,j] = 1 - (U*V/(U*V+min((U-U*V),(V-U*V)))) # Eq. ?
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
  }
  res <- as.dist(res)
}

######## End of beta.div function
