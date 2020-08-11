# Generating transmission trees using a branching process model

This code has been adapted from cmmid/ringbp developed by @jhellewell14

## Usage

### Set up

Set your working directory to the home directory of this project 

```r
library(data.table)
source("R/scenario_sim.R")
```
### Input: Model and parameters

A transmission chain is generated using a branching process model. The following assumptions are made:

- There are two age groups (0-19, child , 20+ adult) who contact each other as per the POLYMOD study (parameterised using socialmixr R package). This contact matrix is C
- The relative infectiousness and susceptibility of a child infection can be set as arguments
- The next generation matrix is calculated by scaling the 


### Generate n transmission chains

Run n.sim simulations.

```r

tt <- scenario_sim(scenario_sim(n.sim = 1,  # number of simulations
                              cap_cases = 100, # minimum number of cases
                              cap_max_days = 100, # minimum number of days (simulation stops after either min cases or min days has been reached in n-1 generation
                              num.initial.cases = 1,  # number of initial infected individuals
                              initial.case.adult = FALSE, # whether initial case is an adult or a child
                              prop.asym = 0,   # probability of asymptomatic infection
                              prop.seq = 1,   # proportion of symptomatic cases that are sequenced
                              r0community = 2.5,  # population R0 
                              disp.com = 0.16,  # dispersal parameter (of neg binom distribution) of community reproduction numbers
                              rel.infectiousness.c = 1, # relative infectiousness of child infections to adult infections
                              rel.susceptibility.c = 1, # relative susceptibility of child infections to adult infections
                              k = 0,  # skew parameter to sample serial interval from incubation period distribution
                              sample_shape = 3,  # param1 for delay from onset to sampling (Weibull distribution)  
                              sample_scale = 2,  # param2 for delay from onset to sampling (Weibull distribution) 
                              ####### final parameters not needed for tranmission chain generation (assuming no isolation / quarantine)
                              r0isolated = 0,  # R0 for isolated cases
                              disp.iso = 1,  # dispersal parameter (of neg binom distribution) of reproduction numbers
                              prop.ascertain = 0, # probability of detecting case for isolation
                              quarantine = FALSE) #whether quarantine is implemented
```

