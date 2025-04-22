#' Correct Counts function used for building mortality profiles. 
#' Function takes as input a string of letters A:I corresponding to Payne's Age classes.
#' Creates a data.frame tallying age classes. 
#' Function corrects counts of age class occurrences by adding fractions to total count based on number of different age classes suggested for an individual.
#' @param a column containing Payne's Age Groups
#' @param probability.correction logical. Whether to correct counts using Vigne and Helmer's (2007) probability correction rule. Default FALSE. If TRUE, combines age groups EF, and HI
#' @returns A two-column dataframe of age class (V) and corrected count (n)
#' 
#'@description
#'Function works as follows... to correct counts of age groups we need to:
#' 1. create a vector of valid groups, v (i.e., A-I)
#' 2. create a table, t to store the number of occurrences, n of each group in v that exist in the input table 
#' 3. for each item, i that has 0 groups listed, skip.
#' 4. for each item i, that has only one group listed, add 1 to n in table t for group v[i]
#' 5. for multiple groups listed, create a vector, g to store the groups listed.
#' 6. calculate the fraction, f that must be added to each group (i.e., 1/length(g))
#' 7. add f to the sum of occurrences for each group in g
#' 8. update t with the new count

correct_counts = function(a, probability.correction = FALSE){
  # vector of groups (A through I)
  v = LETTERS[1:9]
  # table to store counts
  t = data.frame(v=v, n=0)
  
  for(i in 1:length(a)){ # iterate through rows
    if(nchar(a[i])==0){  # if blank, do nothing
      next
    }
    if(nchar(a[i])==1){  # if group is valid add one to count
      t$n[t$v==a[i]] = t$n[t$v==a[i]] + 1
    }
    if(nchar(a[i])>1){   # if multiple groups, parse and add to count for each group proportionally
      g = str_match(a[i], v) # collect age classes
      g = g[,1][!is.na(g)]   # create vector, dropping NA's from match fun
      f = 1/length(g)        # calculate fraction
      for(l in 1:length(g)){ # add fraction to table for each letter represented in age class
        t$n[t$v==g[l]] = t$n[t$v==g[l]] + f
      }
    }
  }
  
  rel.freq = function( x ) {
    tot = sum( x )
    q = vector()
    for ( i in 1 : length( x ) ) {
      q[ i ] = x[ i ] / tot
    }
    q
  }
  
  if( probability.correction ){
    t2 = data.frame(v = t$v[1:4], n = t$n[1:4])
    t2 = rbind.data.frame(t2, data.frame(v = c("EF", "G", "HI"), n = c( t$n[5] + t$n[6], t$n[7], t$n[8] + t$n[9] )))
    
    # relative frequencies
    t2$qx = rel.freq( t2$n )
    
    #-- probability corrections after Vigne and Helmer 2007: 
    #-- for age groups A,   B,   C,  D, EF,   G,   HI
    #-- probability  1/6, 1/3, 1/2, 1,   2,  2,   4
    t2$fd = t2$qx * c(    6,   3,   2, 1, 0.5, 0.5, 0.25 ) # 1/p correction
    # t2$v[5:7] = c("EF", "G", "HI")
    # t2$n[5:7] = c( t$n[5] + t$n[6], t$n[7], t$n[8] + t$n[9] )
    t = t2
  }
  
  return(as.data.frame(t))
}
