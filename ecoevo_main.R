## To run, either execute within R or enter at the command prompt ($):
## $ Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]

require(deSolve) ## for solving ordinary differential equations (ODEs)
require(tidyverse) ## for manipulating and visualizing data

source("./plotting_functions.R") ## various functions for plotting final data


## ---------------------------- input parameters --------------------------------


clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) { ## command-line arguments
  S <- as.numeric(clargs[1]) ## number of species per trophic level
  vbar <- as.numeric(clargs[2]) ## mean genetic variance
  dbar <- as.numeric(clargs[3]) ## mean dispersal rate
  model <- clargs[4] ## "baseline", "trophic", "Tdep", or "Tdep_trophic"
  replicate <- as.numeric(clargs[5]) ## for seeding random number generator
  outfile <- clargs[6] ## name of file to save data in (w/ path & extension)
} else { ## sample input parameters, if no command line arguments are given
  S <- 50 ## fifty species per trophic level
  vbar <- 1e-1 ## average genetic variance = 0.1 celsius squared
  dbar <- 1e-5 ## average dispersal = 1e-5 (100 meters per year)
  model <- "baseline" ## baseline model (one trophic level, constant comp coeffs)
  replicate <- 1 ## replicate number = 1
  outfile <- "" ## no output file; make plot instead
}


## --------------------------------functions ------------------------------------


## apply smoothed step function to an arbitrary array n: n < 0 are set to 0,
## 0 <= n <= 1 are set to 10*n^3-15*n^4+6*n^5, and n > 1 are set to 1
smoothstep <- function(n) {
  return(ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1))
}

## temperature as a function of space (x) and time (t), where x goes
## from 0 (north pole) to 1 (equator)
Temp <- function(x, t, t0, tE, Cmax, Cmin, Tmax, Tmin) {
  return((Tmax-Tmin)*x+Tmin + ((Cmin-Cmax)*x+Cmax)*smoothstep((t-t0)/tE))
}

## return matrix W[i,j], which is nonzero if consumer i eats resource j;
## SR is the number of resource, SC the number of consumer species
generate_network <- function(SR, SC) {
  w <- matrix(0, SR+SC, SR+SC) ## initialize adjacency matrix
  for (i in 1:SR) { ## determine which resources each consumer eats: it must eat
    indices <- sort(c(i, sample((1:SC)[-i], SR/2-1))) ## the one with
    w[i+SR,indices] <- 1 ## matching trait value, plus a fixed number of
  } ## randomly assigned ones (in this case, to make the connectance 1/2)
  omega <- numeric(0) ## initialize matrix of consumption efforts
  rsum <- rowSums(w) ## omega[i,j] is the proportion of i's consumption rate,
  for (i in 1:(SR+SC)) omega <- cbind(omega, rsum) # targeted at consuming j
  omega[omega!=0] <- 1/omega[omega!=0] ## if not 0, set to proportion
  W <- unname(w*omega) ## only the product of w and omega is used
  return(W)
}

## type II functional response
funcresp <- function(n, Th, arate, W) {
  return(arate*(W%*%diag(n))/as.vector(1+arate*Th*(W%*%(n))))
}

## differential equations for the change in species densities and trait
## means through time and over space
eqs <- function(time, state, pars) {
  ## arrange species densities and traits into SxL matrices
  n <- matrix(state[1:(pars$S*pars$L)], pars$S, pars$L)
  m <- matrix(state[(pars$S*pars$L+1):(2*pars$S*pars$L)], pars$S, pars$L)
  n[n<1e-6] <- 0 ## threshold of extinction
  ## reserve memory for arrays
  dndt <- matrix(0, pars$S, pars$L) ## time derivative of density i in patch k
  dmdt <- matrix(0, pars$S, pars$L) ## time derivative of mean trait i in patch k
  b <- matrix(0, pars$S, pars$L)
  g <- matrix(0, pars$S, pars$L)
  an <- matrix(0, pars$S, pars$L)
  Fn <- matrix(0, pars$S, pars$L)
  bn <- matrix(0, pars$S, pars$L)
  mign <- matrix(0, pars$S, pars$L)
  summig <- matrix(0, pars$S, pars$L)
  mnmu <- matrix(0, pars$S, pars$L)
  ## calculate ingredient functions
  q <- smoothstep(n/pars$nmin) ## reduction of genetic variance at low abundances
  h2 <- q*pars$vmat/(q*pars$vmat+pars$venv) ## heritability (V_genetic/V_total)
  w <- pars$bw-pars$aw*m ## temperature tolerance widths
  T <- Temp(seq(from=0, to=1, l=pars$L), time, pars$t0, pars$tE, pars$Cmax,
            pars$Cmin, pars$Tmax, pars$Tmin) ## temperature in each patch
  for (i in 1:pars$S) {
    mign[i,] <- as.numeric(pars$mig[i,,]%*%n[i,]) ## sum_l mig_ikl*n_il
    summig[i,] <- as.numeric(rep(1, pars$L)%*%pars$mig[i,,]) ## sum_l mig_ilk
    mnmu[i,] <- as.numeric(pars$mig[i,,]%*%
                             (n[i,]*m[i,])) ## sum_l mig_ikl*m_il*n_il
  }
  for (k in 1:pars$L) {
    ef <- exp(-(T[k]-m[,k])^2/(2*((w[,k])^2+pars$s))) ## exp. factor
    b[,k] <- (pars$rho/w[,k])*(w[,k]/sqrt((w[,k])^2+pars$s))*ef - pars$kappa
    g[,k] <- (pars$rho/w[,k])*(w[,k]*pars$s/(
      ((w[,k])^2+pars$s)^(3/2)))*ef*(T[k]-m[,k])
    F <- funcresp(n[,k], pars$Th, pars$arate, pars$W)
    Fn[,k] <- rowSums(pars$eps*F) - as.vector(t(F)%*%n[,k])
  }
  if (model %in% c("Tdep", "Tdep_trophic")) {
    ## sum matrix of phenotypic variances
    sv <- outer(pars$s[1:SR], pars$s[1:SR], FUN="+")
    for (k in 1:pars$L) {
      ## difference matrix of trait means
      dm <- outer(m[1:SR,k], m[1:SR,k], FUN="-")
      an[,k] <- (-exp(-dm^2/(2*sv+pars$eta^2))*pars$eta/
                   sqrt(2*sv+pars$eta^2))%*%n[1:SR,k] ## sum_j alpha_ijk*n_jk
      bn[,k] <- (-2*exp(-dm^2/(2*sv+pars$eta^2))*s[1:SR]*pars$eta*(-dm)/
                   (2*sv+pars$eta^2)^(3/2))%*%n[1:SR,k] ## sum_j beta_ijk*n_jk
    }
  } else {
    for (k in 1:pars$L) an[,k] <- (-pars$a)%*%n[,k] ## sum_j alpha_ijk*n_jk
  }
  ninv <- 1/(n+(1e-10)) ## ninv = 1/n_ik (plus 1e-10 to avoid division by 0)
  ## set up equations
  dndt <- (n*b+n*an+Fn)*smoothstep(n/(1e-6))+mign-n*summig
  dmdt <- h2*(g+bn+ninv*(mnmu-m*mign))
  ## implement periodic boundaries
  dndt[,1] <- dndt[,1] + pars$mig[,1,2]*n[,2] - pars$mig[,2,1]*n[,1]
  dndt[,pars$L] <- dndt[,pars$L] + pars$mig[,pars$L,pars$L-1]*n[,pars$L-1] -
    pars$mig[,pars$L-1,pars$L]*n[,pars$L]
  dmdt[,1] <- dmdt[,1] + h2[,1]*pars$mig[,1,2]*n[,2]*ninv[,1]*(m[,2] - m[,1])
  dmdt[,pars$L] <- dmdt[,pars$L] + h2[,pars$L]*pars$mig[,pars$L,pars$L-1]*
    n[,pars$L-1]*ninv[,pars$L]*(m[,pars$L-1] - m[,pars$L])
  ## return equations by first flattening them back into a single vector
  return(list(c(as.numeric(dndt), as.numeric(dmdt))))
}

## put the results of the numerical integration into a tidy table
organize_data <- function(out, times, pars) {
  dat <- out %>% as.data.frame %>% as_tibble ## convert to tibble
  dat <- dat %>% filter(time %in% times) ## only keep specified time points
  names(dat)[1] <- "time" ## name the first column "time"
  index <- 2 ## keep track of which column we are naming
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { ## name columns showing densities
      names(dat)[index] <- paste0("n_", i, "_", k) ## naming convention:
      index <- index + 1 ## "type_species_patch" - type is either m (trait),
    } ## or n (density)
  }
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { ## name columns showing trait values
      names(dat)[index] <- paste0("m_", i, "_", k) ## same naming convention
      index <- index + 1
    }
  }
  dat %>%
    gather("variable", "v", 2:ncol(dat)) %>% ## tidy up the table
    ## split "variable" into value type (density or trait), species, and patch
    separate(variable, c("type", "species", "patch"), sep="_") %>%
    ## convert species & patch from string ("1","2",...) to integer (1,2,...)
    mutate(species=as.integer(species), patch=as.integer(patch)) %>%
    ## split trait and abundance values into two columns
    spread(type, v) %>%
    ## trophic level (tl): species with index over SR are consumers ("C"),
    ## the rest are resources ("R")
    mutate(tl=ifelse(species>SR, "C", "R")) %>%
    return
}


## ------------------------------- parameters -----------------------------------


## number of species and number of patches
SR <- S ## number of resource species
SC <- 0 ## number of consumer species: 0, unless we have...
if (model %in% c("trophic", "Tdep_trophic")) SC <- S ## ...consumer species
S <- SR + SC ## set S to be the total number of species
L <- 50 ## number of patches

## random- and trophic-dependent quantities
set.seed(1000*replicate+321) ## set random seed for reproducibility
v <- runif(SR, 0.5*vbar, 1.5*vbar) ## resource genetic variances
d <- runif(SR, 0.1*dbar, 10.0*dbar) ## resource dispersal rates
rho <- runif(SR, 0.9, 1.1) ## resource growth-tolerance tradeoff parameter
a <- matrix(0, S, S) ## initialize full competition matrix (resources+consumers)
aP <- matrix(runif(SR*SR, 0.15*0.5, 0.15*1.5), SR, SR) ## resource comp coeffs
diag(aP) <- runif(SR, 0.2*0.5, 0.2*1.5) ## resource intraspecific comp coeffs
a[1:SR,1:SR] <- aP ## top left block: resources
W <- matrix(0, S, S) ## create feeding network: nothing if no consumers
Th <- rep(1, S) ## handling times in type II f.r. (dummy value if no consumers)
arate <- rep(1, S) ## attack rates in type II f.r. (dummy value if no consumers)
if (model %in% c("trophic", "Tdep_trophic")) {
  v <- c(v, runif(SC, 0.5*vbar, 1.5*vbar)) ## add consumer genetic variances
  d <- c(d, runif(SC, 0.1*dbar, 10.0*dbar)) ## add consumer dispersal rates
  rho <- c(rho, runif(SC, 0.9*0.1, 1.1*0.1)) ## add consumer tradeoff parameters
  aH <- matrix(0, SC, SC) ## initialize competition matrix (consumers)
  a[(SR+1):S,(SR+1):S] <- aH ## bottom right: consumers
  W <- generate_network(SR, SC) ## trophic feeding network
  Th[(SR+1):S] <- runif(S-SR, 0.5, 1) ## handling times in type II f.r.
  arate[(SR+1):S] <- runif(S-SR, 1, 10) ## attack rates in type II f.r.
}

## all other parameters
kappa <- 0.1 ## intrinsic mortality parameter
venv <- vbar ## environmental variance
vmat <- matrix(rep(v, L), S, L) ## genetic variances at each patch
s <- v + venv ## species' total phenotypic variances
eta <- 1 ## competition width (centigrade; only for Tdep and Tdep_trophic)
eps <- c(rep(0, SR), rep(0.05, SC)) ## feeding efficiency of consumers
nmin <- 1e-5 ## genetic variance is reduced below this threshold density
aw <- 0.1 ## (negative) slope of trait-dependence of tol. width
bw <- 4 ## intercept of trait-dependence of tol. width
Tmax <- 25.0 ## initial mean temperature at equator
Tmin <- -10.0 ## initial mean temperature at poles
Cmax <- 9.66 ## projected temperature increase at poles
Cmin <- 1.26 ## projected temperature increase at equator
tE <- 300.0 ## saturation time for climate change
t0 <- 1000.0 ## time at which climate change starts

## dispersal matrices
mig <- array(0, c(S, L, L)) ## initialize dispersal tensor
for (i in 1:S) { ## fill migration tensor: each species can only migrate
  for (k in 2:L) mig[i,k-1,k] <- d[i] ## to the two nearest-neighbor patches
  mig[i,,] <- mig[i,,] + t(mig[i,,]) ## can access the first, and vice versa
}

## initial conditions
ninit <- matrix(0, S, L) ## reserve memory for initial densities
## initial trait means
muinit <- matrix(seq(Tmin, Tmax, l=SR), SR, L)
## initial temperatures
Tempinit <- Temp(seq(from=0, to=1, l=L), 0, t0, tE, Cmax, Cmin, Tmax, Tmin)
for (i in 1:SR) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
## initial traits and densities for consumers
if (model %in% c("trophic", "Tdep_trophic")) {
  muinit <- rbind(muinit, matrix(seq(Tmin, Tmax, l=SC), SC, L))
  for (i in (SR+1):S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
}

## coerce parameters into list
pars <- list(SR=SR, SC=SC, S=S, L=L, rho=rho, kappa=kappa, a=a, eta=eta,
               eps=eps, W=W, venv=venv, vmat=vmat, s=s, nmin=nmin, aw=aw, bw=bw,
               Tmax=Tmax, Tmin=Tmin, Th=Th, arate=arate, Cmax=Cmax, Cmin=Cmin,
               tE=tE, t0=t0, mig=mig, model=model)


## --------------------------- integrate ODEs -----------------------------------


## Define sampling points along the time axis
tmax <- 3500 ## time units to simulate for
stepout <- 1 ## spacing of temporal output
time <- seq(0, tmax, by=stepout) ## define time axis
## Solve the system of ODEs
dat <- ode(func=eqs, y=c(ninit,muinit), times=time, parms=pars, method="rk4") %>%
  ## put results in organized tibble
  organize_data(times=seq(from=t0, to=tmax, by=100), pars=pars) %>%
  ## add replicate, genetic var., dispersal rate, and model as new columns
  mutate(replicate=replicate, vbar=vbar, dbar=dbar, model=model) %>%
  ## merge average genetic variance and dispersal into a single column
  mutate(parameterization=paste("V:", vbar, " d:", dbar)) %>%
  ## create regions:
  mutate(region=case_when(
    (patch<=round(max(patch)/3))   ~ "polar", ## top third of patches are "polar"
    (patch>=round(2*max(patch)/3)) ~ "tropical", ## bottom third are "tropical"
    TRUE                           ~ "temperate")) ## the rest are "temperate"


## --------------------------- generate output ----------------------------------


if (outfile!="") { ## if data file to save to was not specified as empty (""):
  write_csv(dat, path=outfile) ## save data to specified file
} else { ## otherwise, create a plot:
  ## replace function below with any function from "plotting_functions.R"
  plot_timeseries(filter(dat, time %in% c(1000,1100,1200,1300,2500,3400,3500)))
}
