
#' Run a single instance of the branching process model
#' @author Joel Hellewell
#' @inheritParams outbreak_step
#' @param delay_shape numeric shape parameter of delay distribution
#' @param delay_scale numeric scale parameter of delay distribution
#' @param rel.infectiousness.c relative infectiousness of children to adults
#' @param rel.susceptibility.c relative susceptibility of children to adults
#'
#' @return data.table of cases by week, cumulative cases, and the effective reproduction number of the outreak
#' @export
#'
#' @importFrom data.table rbindlist
#'
#' @examples
#'
#'\dontrun{
#' incfn <- dist_setup(dist_shape = 2.322737,dist_scale = 6.492272)
#' # delay distribution sampling function
#' delayfn <- dist_setup(2, 4)
#' # generate initial cases
#' case_data <- outbreak_setup(num.initial.cases = 5,
#'                             incfn=incfn,
#'                             delayfn = delayfn,
#'                             k=1.95,
#'                             prop.asym=0)
#' # generate next generation of cases
#' case_data <- outbreak_step(case_data = case_data,
#'                            disp.iso = 1,
#'                            disp.com = 0.16,
#'                            r0isolated = 0,
#'                            r0community = 2.5,
#'                            prop.asym = 0,
#'                            incfn = incfn,
#'                            delayfn = delayfn,
#'                            prop.ascertain = 0,
#'                            k = 1.95,
#'                            quarantine = FALSE)
#'}
outbreak_model <- function(num.initial.cases = NULL, initial.case.adult = NULL, prop.ascertain = NULL,
                           cap_max_days = NULL, cap_cases = NULL,
                           r0isolated = NULL, r0community = NULL,
                           rel.infectiousness.c = NULL, rel.susceptibility.c = NULL,
                           disp.iso = NULL, disp.com = NULL,
                           k = NULL, delay_shape = NULL,
                           delay_scale = NULL, sample_shape = NULL, 
                           sample_scale = NULL, prop.asym = NULL, prop.seq = NULL,
                           quarantine = NULL) {

  
  # Set up functions to sample from distributions
  # incubation period sampling function
  incfn <- dist_setup(dist_shape = 2.322737,
                      dist_scale = 6.492272)
  # incfn <- dist_setup(dist_shape = 3.303525,dist_scale = 6.68849) # incubation function for ECDC run
  # onset to isolation delay sampling function
  delayfn <- dist_setup(delay_shape = 1.651524,
                        delay_scale = 4.287786)
  # onset to sampling delay sampling function
  samplefn <- dist_setup(sample_shape,
                         sample_scale)
  
  
  #Set parameter values for age group Rs
  r0params <- parameter_setup(rel.infectiousness.c = rel.infectiousness.c,
                              rel.susceptibility.c = rel.susceptibility.c, 
                              r0community = r0community)

  # Set initial values for loop indices
  total.cases <- num.initial.cases
  latest.onset <- 0
  extinct <- FALSE

  # Initial setup
  case_data <- outbreak_setup(num.initial.cases = num.initial.cases,
                            initial.case.adult = initial.case.adult,
                            incfn = incfn,
                            prop.asym = prop.asym,
                            prop.seq = prop.seq,
                            delayfn = delayfn,
                            samplefn = samplefn,
                            k = k)
  # Preallocate
  effective_r0_vect <- c()
  cases_in_gen_vect <- c()


  # Model loop
  while (latest.onset < cap_max_days & total.cases < cap_cases & !extinct) {
    
    out <- outbreak_step(case_data = case_data,
                             disp.iso = disp.iso,
                             disp.com = disp.com,
                             r0isolated = r0isolated,
                             r0community = r0community,
                             r0community_aa = r0params$r0community_aa,
                             r0community_cc = r0params$r0community_cc,
                             r0community_ac = r0params$r0community_ac,
                             r0community_ca = r0params$r0community_ca,
                             incfn = incfn,
                             delayfn = delayfn,
                             samplefn = samplefn,
                             prop.ascertain = prop.ascertain,
                             prop.seq = prop.seq,
                             k = k,
                             quarantine = quarantine,
                             prop.asym = prop.asym)


    case_data <- out[[1]]
    effective_r0_vect <- c(effective_r0_vect, out[[2]])
    cases_in_gen_vect <- c(cases_in_gen_vect, out[[3]])
    total.cases <- nrow(case_data)
    latest.onset <- max(case_data$onset)
    extinct <- all(case_data$isolated)
  }

  # Prepare output, group into weeks
  weekly_cases <- case_data[, week := floor(onset / 7)
                            ][, .(weekly_cases = .N), by = week
                              ]
  # maximum outbreak week
  max_week <- floor(cap_max_days / 7)
  # weeks with 0 cases in 0:max_week
  missing_weeks <- (0:max_week)[!(0:max_week %in% weekly_cases$week)]

  # add in missing weeks if any are missing
  if (length(missing_weeks > 0)) {
    weekly_cases <- data.table::rbindlist(list(weekly_cases,
                                               data.table(week = missing_weeks,
                                                          weekly_cases = 0)))
  }
  # order and sum up
  weekly_cases <- weekly_cases[order(week)
                               ][, cumulative := cumsum(weekly_cases)]
  # cut at max_week
  weekly_cases <- weekly_cases[week <= max_week]

  # Add effective R0
  weekly_cases <- weekly_cases[, `:=`(effective_r0 = mean(effective_r0_vect,
                                                          na.rm = TRUE),
                                        cases_per_gen = list(cases_in_gen_vect))]
  # return
  # return(weekly_cases)
  return(case_data)
}
