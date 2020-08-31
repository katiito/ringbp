# source other R files instead of using the ringbp library
source("R/outbreak_setup.R")
source("R/outbreak_step.R")
source("R/outbreak_model.R")
source("R/aux_functions.R")

#' Run a specified number of simulations with identical parameters
#' @author Joel Hellewell / Katie Atkins
#' @param n.sim number of simulations to run
#' @param num.initial.cases Initial number of cases in each initial cluster
#' @param num.initial.clusters Number of initial clusters
#' @param initial.case.adult TRUE for adult, FALSE for child
#' @param prop.ascertain Probability that cases are ascertained by contact tracing
#' @param cap_max_days Maximum number of days to run process for
#' @param cap_cases Maximum number of cases to run process for
#' @param r0isolated basic reproduction number for isolated cases
#' @param r0community basic reproduction number for non-isolated cases
#' @param rel.infectiousness.c relative infectiousness of children to adults
#' @param rel.susceptibility.c relative susceptibility of children to adults
#' @param prop.seq Probability that cases are sequenced
#' @param disp.iso dispersion parameter for negative binomial distribution for isolated cases
#' @param disp.com dispersion parameter for negative binomial distribution for non-isolated cases
#'
#' @importFrom purrr safely
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' res <- scenario_sim(n.sim = 5,
#' num.initial.cases = 5,
#' cap_max_days = 365,
#' cap_cases = 2000,
#' r0isolated = 0,
#' r0community = 2.5,
#' disp.iso = 1,
#' disp.com = 0.16,
#' k = 0.7,
#' delay_shape = 2.5,
#' delay_scale = 5,
#' prop.asym = 0,
#' prop.ascertain = 0)
#' #' }
#'
scenario_sim <- function(n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, cap_cases = NULL,
                         r0isolated = NULL, r0community = NULL, 
                         rel.infectiousness.c = NULL, rel.susceptibility.c = NULL,
                         disp.iso = NULL, disp.com = NULL, k = NULL, initial.case.adult = NULL,
                         num.initial.cases = NULL, prop.asym = NULL,
                         prop.seq = NULL, quarantine = NULL) {
 
  # Run n.sim number of model runs and put them all together in a big data.frame
  
  tt <- purrr::map(.x = 1:n.sim, ~ outbreak_model(num.initial.cases = num.initial.cases,
                                                   prop.ascertain = prop.ascertain,
                                                   cap_max_days = cap_max_days,
                                                   cap_cases = cap_cases,
                                                   r0isolated = r0isolated,
                                                   r0community = r0community,
                                                   rel.infectiousness.c = rel.infectiousness.c,
                                                   rel.susceptibility.c = rel.susceptibility.c,
                                                   initial.case.adult = initial.case.adult,
                                                   disp.iso = disp.iso,
                                                   disp.com = disp.com,
                                                   k = k,
                                                   prop.asym = prop.asym,
                                                   prop.seq = prop.seq,
                                                   quarantine = quarantine))
  
  ## output table/s for FAVITES
  
  
  # extract the infection events and times and fix index case
  inf <- purrr::map(tt, ~
                     select(., infector, caseid, exposure) %>%
                     mutate(., infector = replace(infector, infector==0, "None")))
  rem <- purrr::map(tt, ~
                     select(., infector = caseid, caseid = caseid, exposure = endinfectious) %>%
                     mutate(infector = as.character(infector)))
  din <- purrr::map2(inf, rem, ~bind_rows(.x,.y))
  
  # write to separate files
  list(data = din, sim.num = 1:n.sim) %>%
    purrr::pmap(output_csv)
  
  # res <- purrr::map(.x = 1:n.sim, ~ outbreak_model(num.initial.cases = num.initial.cases,
  #                                            prop.ascertain = prop.ascertain,
  #                                            cap_max_days = cap_max_days,
  #                                            cap_cases = cap_cases,
  #                                            r0isolated = r0isolated,
  #                                            r0community = r0community,
  #                                            disp.iso = disp.iso,
  #                                            disp.com = disp.com,
  #                                            delay_shape = delay_shape,
  #                                            delay_scale = delay_scale,
  #                                            k = k,
  #                                            prop.asym = prop.asym,
  #                                            quarantine = quarantine))


  

 
  # res[, sim := rep(1:n.sim, rep(floor(cap_max_days / 7) + 1, n.sim)), ]
  return(tt)
}
