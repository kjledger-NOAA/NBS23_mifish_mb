
library(nnet)

set.seed(101)
mlogis = function(invlogis_pars) exp(invlogis_pars) / sum(exp(invlogis_pars))

# Settings
n_obs = 10
n_species = 5

# Define probability for each species
p_c = mlogis( rnorm(n_species) )
names(p_c) = paste0("species_", seq_len(n_species))

# Simulate samples
Y_ic = t(rmultinom( n=5, size=100, prob=p_c ))
dimnames(Y_ic) = list("sample"=1:nrow(Y_ic), "species"=1:length(p_c))

# Binomial GLM
if( length(p_c) == 2 ){
  Glm1 = glm( Y_ic ~ 1, family = binomial )
  phat1_c = c( plogis(Glm1$coef[1]), 1-plogis(Glm1$coef[1]) )
}else{
  phat1_c = rep( NA, length(p_c) )
}

# Poisson GLM
dat = cbind(expand.grid(dimnames(Y_ic)), "count"=as.vector(Y_ic))
Glm2 = glm( count ~ 0 + factor(species) + factor(sample), data=dat, family=poisson )
phat2_c = mlogis(Glm2$coef[seq_along(p_c)])

# Multinom package
# https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/
Glm3 = multinom( species ~ 1, weights=count, data=dat )
phat3_c = mlogis( c(0, coef(Glm3)) )

# Compare estimators
round(cbind( 
  "True" = p_c,
  "stats::glm" = phat1_c,
  "nnet::multinom" = phat3_c,
  "PMT" = phat2_c
), 3)


