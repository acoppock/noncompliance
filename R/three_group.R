
## Three Group Design Functions
# This builds on the replication file from Gerber et al PA 2010

# LL function for three-group design
#' @export
ThreeG <- function( par, v_b, a_b, v_t_c, a_t_c, v_t_u, a_t_u, v_p_c, a_p_c,
                    v_p_u, a_p_u)   {

  alpha <- par[1]
  p_r   <- par[2]
  p_nr  <- par[3]
  tau   <- par[4]

  fvalue <-
    v_b   * log(     alpha * (p_r-tau) + (1-alpha) * p_nr ) +
    a_b   * log(1 -( alpha * (p_r-tau) + (1-alpha) * p_nr )) +
    v_t_c * log( alpha * (   p_r) ) +
    a_t_c * log( alpha * (1-(p_r)) ) +
    v_t_u * log( (1 -  alpha) * (    p_nr) ) +
    a_t_u * log( (1 -  alpha) * (1 - p_nr) ) +
    v_p_c * log( alpha * (   p_r - tau) ) +
    a_p_c * log( alpha * (1-(p_r - tau)) ) +
    v_p_u * log( (1 -  alpha) * (    p_nr) ) +
    a_p_u * log( (1 -  alpha) * (1 - p_nr) )

  if(is.na(fvalue))   { fvalue <- -9999999 }
  if(fvalue == Inf | fvalue == -Inf) { fvalue <- -9999999 }
  return(fvalue)
}


# LL function for baeeline and treatment group comparison
#' @export
ThreeGBT <- function( par, v_b, a_b, v_t_c, a_t_c, v_t_u, a_t_u)    {

  alpha <- par[1]
  p_r   <- par[2]
  p_nr  <- par[3]
  tau   <- par[4]

  fvalue <-
    v_b   * log(     alpha * (p_r-tau) + (1-alpha) * p_nr ) +
    a_b   * log(1- ( alpha * (p_r-tau) + (1-alpha) * p_nr )) +
    v_t_c * log( alpha * (   p_r) ) +
    a_t_c * log( alpha * (1-(p_r)) ) +
    v_t_u * log( (1 -  alpha) * (    p_nr) ) +
    a_t_u * log( (1 -  alpha) * (1 - p_nr) )

  if(is.na(fvalue))   { fvalue <- -9999999 }
  if(fvalue == Inf | fvalue == -Inf)  { fvalue <- -9999999 }
  return(fvalue)
}

# LL function for placebo and treatment group comparison
#' @export
ThreeGPT <- function( par, v_t_c, a_t_c, v_t_u, a_t_u, v_p_c, a_p_c,
                      v_p_u, a_p_u)   {

  alpha <- par[1]
  p_r   <- par[2]
  p_nr  <- par[3]
  tau   <- par[4]

  fvalue <-

    v_t_c * log( alpha * (   p_r) ) +
    a_t_c * log( alpha * (1-(p_r)) ) +
    v_t_u * log( (1 -  alpha) * (    p_nr) ) +
    a_t_u * log( (1 -  alpha) * (1 - p_nr) ) +
    v_p_c * log( alpha * (   p_r - tau) ) +
    a_p_c * log( alpha * (1-(p_r - tau)) ) +
    v_p_u * log( (1 -  alpha) * (    p_nr) ) +
    a_p_u * log( (1 -  alpha) * (1 - p_nr) )

  if(is.na(fvalue))   { fvalue <- -9999999 }
  if(fvalue == Inf | fvalue == -Inf)  { fvalue <- -9999999 }

  return(fvalue)

}
