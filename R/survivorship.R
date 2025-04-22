#' Survivorship function
#' @param x is a vector containing the absolute frequencies of individuals in each age class
#' @return Running the function will give a vector equal in length to x containing survival probabilities for each age class. 
#' @references Price, M., Wolfhagen, J. and Otárola-Castillo, E. (2016). Confidence Intervals in the Analysis of Mortality and Survivorship Curves in Zooarchaeology. Am. Antiq. 81: 157–173.
survivorship = function( x ){
  rel.freq = function( x ) {
    tot = sum( x )
    q = vector()
    for ( i in 1 : length( x ) ) {
      q[ i ] = x[ i ] / tot
      }
    q
  }
  # relative frequencies
  qx = rel.freq( x )
  # cumulative sum
  Qx = cumsum( qx )
  # survivorship probability
  s = 1 - Qx
  s
}
