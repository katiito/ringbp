#' Set up R0 values for between age groups

#' @author Katie Atkins
#' @param rel.infectiousness.c relative infectiousness of children to adults
#' @param rel.susceptibility.c relative susceptibility of children to adults
#' @param r0community R0 within whole population

#' @return list of age by age group R0s



parameter_setup <- function(rel.infectiousness.c = NULL, rel.susceptibility.c = NULL, r0community = NULL){

  
      # relative contact rates of children(0-20y) and adults (21+y)
      rel.contactrate.cc <- 1.0000000
      rel.contactrate.ca <- 0.7938829
      rel.contactrate.ac <- 0.2584824
      rel.contactrate.aa <- 1.0542519
      
      unscaled_ngm <- matrix(
                    c(rel.infectiousness.c * rel.susceptibility.c * rel.contactrate.cc,
                      rel.infectiousness.c * rel.contactrate.ca,
                      rel.susceptibility.c * rel.contactrate.ac, 
                      rel.contactrate.aa), 
                    nrow = 2, byrow = TRUE)
      dom.eval <- max(eigen(unscaled_ngm)$values)
      
      scaling.const <- r0community / dom.eval
      
      ngm <- unscaled_ngm * scaling.const
      
      out <- list(r0community_cc = ngm[1,1],
                      r0community_ca = ngm[1,2],
                      r0community_ac = ngm[2,1],
                      r0community_aa = ngm[2,2])
      return(out)
      
      
}