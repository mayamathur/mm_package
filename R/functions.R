
################################ FORMAT RESULTS ################################

#' Round with trailing zeroes
#'
#' @param x Value to round
#' @param digits Number of digits to keep
#' @export
#' @examples
my_round = function(x, digits) {
  formatC( round( x, digits ), format='f', digits=digits )
}

#' Format CI as string
#'
#' @param lo Lower limit
#' @param hi Upper limit
#' @param digits Number of digits to keep
#' @export
format_CI = function( lo, hi, digits ) {
  paste( "[", my_round( lo, digits ), ", ", my_round( hi, digits ), "]", sep="" )
}

#' Nicely format p-value
#'
#' @param p The p-value
#' @export
format_pval = function(p) {
  if (p >= 0.01) return( my_round( p, 2 ) )
  if (p < 0.01 & p > 10^-5 ) return( formatC( p, format = "e", digits = 0 ) )
  if ( p < 10^-5 ) return("< 1e-05")
}


################################ META-ANALYSIS ################################

#' Convert Pearsons's r to Fisher's z
#'
#' @param r The correlation
#' @export
r_to_z = Vectorize( function(r) {
  .5 * ( log(1 + r) - log(1 - r) )
}, vectorize.args = "r" )

#' Convert Fisher's z to Pearson's r
#'
#' @param z The Fisher's z
#' @export
z_to_r = Vectorize( function(z) {
  ( exp( 2 * z ) - 1 ) / ( exp( 2 * z ) + 1 )
}, vectorize.args = "z" )

# gives CI for tau from meta-analysis fit in metafor

#' Calculate CI for tau from metafor
#'
#' @param meta Object from metafor
#' @param z.to.r Should we convert from Fisher's z to r?
#' @export
tau_CI = function( meta, z.to.r = FALSE ) {
  t2.lb = meta$tau2 - qnorm(1 - 0.05/2) * meta$se.tau2
  t2.ub = meta$tau2 + qnorm(1 - 0.05/2) * meta$se.tau2

  if ( t2.lb > 0 ) tau.lb = sqrt(t2.lb) else tau.lb = 0
  tau.ub = sqrt(t2.ub)

  if( z.to.r == TRUE ) {
    tau.lb = z_to_r(tau.lb)
    tau.ub = z_to_r(tau.ub)
  }

  if ( tau.lb < 0 ) tau.lb = 0

  return( c(tau.lb, tau.ub))
}




################################ MODELS FOR CORRELATED DATA ################################

#' Calculate Wald p-values for lmer
#'
#' These will agree exactly with confint( , method="Wald" ), unlike
#' lmerTest's too-conservative Satterthwaite approximation.
#' @param .m The merMod object from lmer
#' @export
z_pvals = function( .m ) {
  coefs = data.frame( coef( summary( .m ) ) )
  p = data.frame( 2 * ( 1 - pnorm( abs( coefs$t.value ) ) ) )
  row.names(p) = row.names(coefs)
  return(p)
}




