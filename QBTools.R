###############################################################################
#                                                                             #
# Quantitative Biodiversity Functions Source Code                             #
#                                                                             #
# Written by: Michael Santel                                                  #
#                                                                             #
# Last update: 2015/03/06                                                     #
#                                                                             #
# Notes: This file contains functions to calculate metrics of taxonomic and   #
# phylogenetic diversity                                                      #
#                                                                             #
###############################################################################

require("vegan")||install.packages("vegan");require("vegan")

# My_Function: returns estimate of beta biodiversity
# Inputs: x, y, z

My_Function <- function(x, y, z){
  ...
}

# sem: Standard Error of the Mean
# Inputs: x

sem <- function(v){
  sd(na.omit(v))/sqrt(length(v))
}

# Alpha Diversity

# S.obs: Calculates observed species richness
# Inputs: x

S.obs <- function(x= ""){
  rowSums( x >  0  ) * 1
}

# C: Calculates Good's Coverage for estimation of how well a site was sampled
# Inputs: x

C <- function(x = ""){
  1 - (sum(x == 1) / rowSums (x)) 
}

# S.chao1: abundance based estimator for single site richness
# Inputs: x

S.chao1 <- function(x = ""){
  S.obs(x) + (sum(x == 1)^2) / (2 * sum(x == 2))
}

# S.chao2: incidence based estimator for richness across multiple sites
# Inputs: x

S.chao2 <- function(site = "", SbyS = ""){
  SbyS = as.data.frame(SbyS)
  x = SbyS[site, ]
  SbyS.pa <- (SbyS > 0) * 1
  Q1 = sum(colSums(SbyS.pa) == 1)
  Q2 = sum(colSums(SbyS.pa) == 2)
  S.chao2 = S.obs(x) + (Q1^2)/(2 * Q2)
  return(S.chao2)
}

# RAC: rank abundance curve
# Inputs: x

RAC <- function(x = ""){
  x = as.vector(x)
  x.ab = x[x > 0]
  x.ab.ranked = x.ab[order(x.ab, decreasing = TRUE)]
  return(x.ab.ranked)
}

# SimpE: Simpson's Evenness
# Inputs: x

SimpE <- function(x = ""){
  x = as.data.frame(x)
  D <- diversity(x, "inv")
  S <- S.obs(x)
  E <- (D)/S
  return(E)
}

# Evar: Smith and Wilson's Evenness Index
# Inputs: x

Evar <- function(x){
  x <- as.vector(x[x>0])
  1 - (2/pi)*atan(var(log(x)))
}

# Beta Diversity

# beta.w: Whittaker's Beta Diversity

beta.w <- function(site1 = "", site2 = ""){
  site1 = subset(site1, select = site1 > 0)
  site2 = subset(site2, select = site2 > 0)
  gamma = union(colnames(site1), colnames(site2))
  s     = length(gamma)
  a.bar = mean(c(specnumber(site1), specnumber(site2)))
  b.w   = round(s/a.bar - 1, 3)
  return(b.w)
}

