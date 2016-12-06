#' Simulation of ecological data structured by single ecological gradient
#' @author David Zeleny (zeleny.david@@gmail.com) - based on the scripts of Jason Fridley (Fridley et al. 2007, Appendix S2) and David Zeleny (Zeleny 2009, Appendix S1)


#' @param totS Total number of species in simulation.
#' @param gr.length Length of the gradient in units.
#' @param niche.type Shape of species response curves ('random', 'normal', 'skewed').
#' @param max.niche.breath Maximum niche breath along the gradient (default = \code{gr.length})
#' @param min.niche.breath Minimum niche breath along the gradient (default = \code{10}). Cannot be more than \code{max.niche.breath}.
#' @param spec.optima Distribution of species optima along simulated gradient. Either one of \code{'random'}, \code{'unimodal'} or \code{'skewed'}, or a vector (of the length equal to \code{totS}) with numbers (ranging between 1 and \code{gr.length}). Default = \code{'random'}.
#' @param prop.random.species Proportion of random species, i.e. species which are not related to the gradient. Should be real number from [0,1]. These species are generated by randomizing occurrence of given proportion of species among samples. Default = 0.
#' @param seed Set random seed for reproducing always the same result.
#' @param plotting Should the species response curves along simulated gradient be drawned? Default = \code{FALSE}.
#' @param highlight.species Vector of species numbers which should be highlighted by color in ploted diagram (if \code{plotting = TRUE}). Default = \code{NULL}, which means that if plot is drawned, two species will be highlighted - the most generalized and the most specialized one.
#' @param simul.comm Object created by \code{simul.comm} function with parameters of individual species response curves.
#' @param Np Number of samples to be generated.
#' @param sample.x Positions of sampling locations along the gradient. Default = \code{"random"}, meaning that random locations are generated. Other options include \code{"equal"} with samples distributed in equal distances, \code{"biased"} with samples accumulated toward higher values of the gradient. Can be also vector of the same length as \code{Np} with exact positions of the samples.
#' @param no.ind Mean number of individuals to be sampled from the species pool (if \code{based.on = 'individuals'}).
#' @param cv.no.ind Coefficient of variation of number of indiduals sampled from the species pool. Increasing it allows to create set of samples with various levels of undersampling.
#' @param k Mean proportion of species from species pool to be drawn into local community (if \code{based.on = 'species'}). If only one number is given, in each sample given proportion of species will be sampled from species pool. If vector of two numbers (range of \code{k} valus) is given, random \code{k} in given range will be generated (with uniform distribution) for each sample; this allows to vary the degree of undersampling among generated samples.
#' @param pa Result table will be generated in presence-absence form (default \code{pa = FALSE}).
#' @param based.on Should the sampling of species from species pool be based on \code{'individuals'} or \code{'species'}? (default = \code{'individuals'})
#' @param standardize.rowsums Standardizes the abundances of species in samples so as the sum of species abundances in sample is 100 (standardizes to rowsums equal to 100). Applies only if \code{pa = FALSE}. Default = \code{TRUE}.
#' @return The function \code{simul.comm} returns \code{list} with 10 items, describing the set of parameters used to simulate community of species response curves:
#' \itemize{
#' \item \code{totS} Total number of species in simulation.
#' \item \code{gr.length} Length of the gradient.
#' \item \code{niche.type} Shape of species response curves.
#' \item \code{Ao} Vector of species amplitudes for the simulated gradient (heights of species response curves, corresponding to maximum probability of species to be selected to community in species optima).
#' \item \code{m} Vector of species optima along the simulated gradient.
#' \item \code{r} Vector of generated species ranges along the simulated gradient (generated niche breaths).
#' \item \code{range} Vector of realised species ranges along the simulated gradient, considering the truncation of species response curves by gradient margins (these differ from \code{r} especially at the gradient margins, where the generated species niche may be wide, but since the margin cuts the species occurrences, realised species niche is narrower.)
#' \item \code{a, g} Vectors of shape parameters for curves (used in beta function to generate the shape of the species response curve)
#' \item \code{A.all} Matrix (dim = gradient length x number of species) with simulated probabilites of individual specie at individual location along the simulated gradient. 
#' }
#' The function \code{sample.comm} returns \code{list} of 5 items with parameters of generated community data; the last item contains also all items returned by \code{simul.comm} function:
#' \itemize{
#' \item \code{a.mat} Matrix (sample x species) of species abundances in samples.
#' \item \code{p.mat} Matrix (sample x species) of species occurrence probabilities in samples.
#' \item \code{sample.x} Vector with positions of samples along the simulated gradient (environmental variable).
#' \item \code{sample.comm} List of 7 items storing initial settings of arguments in \code{sample.comm} (namely arguments \code{Np, based.on, no.ind, k, seed, pa} and \code{standardize.rowsums}).
#' \item \code{simul.comm} List of 10 items returned by function \code{simul.comm}.
#' }




#' @rdname simul.comm
#' @export
simul.comm <- function (totS = 300, gr.length = 5000, niche.type = 'random', max.niche.breath = gr.length, min.niche.breath = 10, spec.optima = 'random', prop.random.species = 0, seed = NULL, plotting = F, highlight.species = NULL)
{
  if (!is.null (seed)) set.seed (seed)
  
  #This is beta function for generating niches
  curve <- function(Ao,m,r,a,g,x)  {
    (Ao*((((x-m)/r)+(a/(a+g)))^a)*((1-(((x-m)/r)+(a/(a+g))))^g))/(((a/(a+g))^a)*((1-(a/(a+g)))^g))
  }
  
  x <- seq(1, gr.length, by = 1)
  S <- totS #number of species
  Ao<-rlnorm(S,2,1)  		#amplitude vector (lognormal distribution)
  if (spec.optima == 'random')  m<-sample(seq(5,max(x)-5),S, replace = T) else
    if (spec.optima == 'unimodal') {
      starting <- seq (1, gr.length, by = gr.length/10)
      biased <- c(20,40,60,80,100,100,80,60,40,20)/600  # keep the same structure as in Zeleny 2009, Appendix S1
      m <- NULL
      for (i in seq (1, 10))
        m <- append (m, sample (seq (starting[i], starting[i]+gr.length/10-1), round (biased[i]*totS)))
      } else
      if (spec.optima == 'skewed'){
        starting <- seq (1, gr.length, by = gr.length/10)
        biased <- c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200)/1100 # keep the same structure as in Zeleny 2009, Appendix S1
        m <- NULL
        for (i in seq (1, 10))
          m <- append (m, sample (seq (starting[i], starting[i]+gr.length/10-1), round (biased[i]*totS)))
      } else m = spec.optima
  r<-runif(S,min=min.niche.breath, max=max.niche.breath)		#range along gradient (niche breadth)
  names (r) <- paste ('spec_', 1:S, sep = '')
  
  # species values for random niches
  if(niche.type=="random") {
    a<-(runif(S,min=.1,max=4))		#shape parameter (alpha)
    g<-(runif(S,min=.1,max=4))		#shape parameter (gamma)
  }
  
  # species values for normal niches
  if(niche.type=="normal") {
    a<-rep(1.99,S)				#shape parameter (alpha)
    g<-rep(1.99,S)				#shape parameter (gamma)
  }
  
  # species values for skewed niches
  if(niche.type=="skewed") {
    #produces skews in either direction, randomly
    a.1<-rep(1.99,S)
    g.1<-rep(.25,S)
    samp<-sample(c(1:S),S/2,replace=FALSE)
    a.1[samp]<-.25
    g.1[samp]<-1.99
    a<-a.1				#shape parameter (alpha)
    g<-g.1				#shape parameter (gamma)
 }
  
  
  A.all <- matrix(0, nrow = gr.length, ncol = S) #response abundances for all points along gradient 1
  
  for(L in 1:S){
    A.all[,L] <- curve(Ao[L],m[L],r[L],a[L],g[L], x)
  }
  
  A.all[is.na (A.all)] <- 0

  if (prop.random.species[1] > 0)
    for (co in sample (1:ncol (A.all), round (prop.random.species[1]*ncol (A.all))))
      A.all[,co] <- sample (A.all[,co])

  range <- colSums (A.all > 0)
  names (range) <- paste ('spec_', 1:S, sep = '')

  if (plotting)
  {
    # Plot species distributions along gradient
    plot(x,A.all[,1],xlim=c(0,max(x)),ylim=c(0,max(Ao)),type="l",xlab="Gradient",ylab="Abundance",cex.lab=1.7,cex.axis=1.5,lwd=2)
    for(L in 2:S) {
      lines(x,A.all[,L])
    }
    if (is.null (highlight.species))
    {
      lines(x,A.all[,which (range == max (range))[1]],lwd=3,col=2)
      lines(x,A.all[,which (range == min (range))[1]],lwd=3,col=3)
    } else 
      for (co in seq (1, length (highlight.species)))
        lines(x,A.all[,highlight.species[co]],lwd=3,col=co+1)
 }  

  
  result <- list (totS = totS, gr.length = gr.length, niche.type = niche.type, Ao = Ao, m = m, r = r, range = range, a = a, g = g, A.all = A.all)
  result
}

#' @rdname simul.comm
#' @export
sample.comm <- function (simul.comm = NULL, Np = 300, sample.x = "random", no.ind = 100, cv.no.ind = 0.01, k = 0.2, seed = NULL, pa = FALSE, based.on = 'individuals', standardize.rowsums = TRUE)
{
  if (!is.null (seed)) set.seed (seed)
  if (is.null (simul.comm)) simul.comm <- simul.comm ()
  sc <- simul.comm
  BASED.ON <- c('individuals', 'species')
  based.on <- match.arg (based.on, BASED.ON)
  
  #Random sample intervals along gradient
  if (sample.x == "random") sample.x <- trunc(sample(c(2:sc$gr.length)-1,Np, replace = T)) else
  #Equal sample intervals along gradient
  if (sample.x == "equal") sample.x <- trunc(seq(2,sc$gr.length-1,length=Np)) else
  #Biased sample intervals along gradient
  if (sample.x == 'biased') {
    exp.x<-sort(rexp(Np,rate=20))
    exp.sample.x<-trunc(exp.x/(max(exp.x))*(sc$gr.length-1))
    while( length(unique(exp.sample.x))!=Np )
      exp.sample.x[duplicated(exp.sample.x)]<-exp.sample.x[duplicated(exp.sample.x)]+1
    sample.x <- rev(sc$gr.length-exp.sample.x)	#switch to other side of gradient
  }

  A <- simul.comm$A.all[sample.x, ]

  p.mat <- A
  p.mat[is.na (p.mat)] <- 0
  
  a.mat <- p.mat*0 # prepared abundance matrix
  
  #output data frames
  samp.out.rand <- matrix(0, nrow = Np, ncol = sc$totS)
  
  #Sampling for random-sample-interval based on individuals
  if (based.on == 'individuals')
  {
    draws.rand <- round(rnorm(Np, mean = no.ind, sd = cv.no.ind * no.ind))
    draws.rand [draws.rand <= 0] <- 1  # if the number of individuals should be zerro or negative, replace by one (to avoid empty samples)
    for(i in 1:Np) {
      samp.prob<-p.mat[i,]		#probabilities of sampling each species in given location (based on rel abundance)
      if (sum (samp.prob) > 0) 
      {
        tab.samp <- table(sample(c(1:sc$totS), size = draws.rand[i], prob = samp.prob, replace = T)) 	#tabulated vector of spp identities after choosing "draws" no. of individuals
        a.mat[i,][as.numeric(names(tab.samp))] <- tab.samp
      } else a.mat[i,] <- 0
    } 
  }
  #Sampling for random-sample-interval based on no of species
  if (based.on == 'species')
    for (i in 1:Np) {
      samp.prob <- p.mat [i,]
      spec.pool.size <- sum (samp.prob > 0)
      #no.spec <- round (rnorm (1, k * spec.pool.size))
      if (length (k) == 1) no.spec <- round (k*spec.pool.size) else
        if (length (k) == 2) no.spec <- round (runif (1, min = k[1], max = k[2]) * spec.pool.size)
      if (no.spec <= 0) no.spec <- 1
      if (no.spec > sum (samp.prob > 0)) no.spec <- sum (samp.prob > 0) # if number of selected species should be higher than number of nonzerro probabilities, it must decrease
      if (sum (samp.prob) > 0) 
      {
        tab.samp <- table(sample(c(1:sc$totS), size=no.spec, prob=samp.prob, replace=FALSE))  #tabulated vector of spp identities after choosing no of species
        a.mat[i,][as.numeric(names(tab.samp))] <- tab.samp*samp.prob[as.numeric(names(tab.samp))]
      } else a.mat[i,] <- 0
    }
  
  colnames (a.mat) <- paste ('spec_', 1:dim (a.mat)[2], sep = '')
  if (standardize.rowsums) a.mat <- vegan::decostand (a.mat, method = 'total')*100
  if (pa == TRUE) a.mat[a.mat>0] <- 1		#presence-absence version
  
  result <- list (a.mat = a.mat, p.mat = p.mat, sample.x = sample.x, sample.comm = list (Np = Np, based.on = based.on, no.ind = if (based.on == 'individuals') no.ind else NULL, k = if (based.on == 'species') k else NULL, seed = seed, pa = pa, standardize.rowsums = standardize.rowsums), simul.comm = simul.comm)
  result
}
