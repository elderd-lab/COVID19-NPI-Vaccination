# COVID19-NPI-Vaccination
Matlab code associated with the manuscript:
Combined effects of vaccination and non-pharmaceutical interventions on COVID-19 transmission and disease burden.

1. SEAIR_model_vac_NPI_dynamics.m - ODE model of population dynamics for a (S)usceptible-(E)xposed-(A)symptomatic-(I)nfected-(R)ecovered that
    contains an asymptomatic class, vaccination, and non-pharmaceutical interventions (NPI) where one can vary the proportion of the population practicing
    NPI. Parameters in the model are taken from the best estimate of the parameter published in other studies.

2. SEAIR_model_vac_NPI_gridsearch.m - uses the ODE model from SEAIR_model_vac_NPI_dynamics.m where one can varying both vaccination ratio for individuals
    practicing NPI and those not practicing NPI while also varying the percent of the population practicing NPI. Parameters in the model are taken from the best estimate
    of the parameter published in other studies.

3. SEAIR_model_vac_NPI_gridsearch_distributions.m - uses the ODE model from SEAIR_model_vac_NPI_dynamics.m where one can varying both vaccination ratio for individuals
    practicing NPI and those not practicing NPI while also varying the percent of the population practicing NPI.  In comparison to SEAIR_model_vac_NPI_gridsearch.m, this
    version of the model randomly draws model parameters from a distribution. Parameters in the model along with their distributions are taken from published estimates
    in other studies.

4. SEAIR_MCMC_vac_NPI_gridsearch.m - uses a stochastic version of the ODE model from SEAIR_model_vac_NPI_dynamics.m where one can varying both vaccination ratio
    for individuals practicing NPI and those not practicing NPI while also varying the percent of the population practicing NPI.  Parameters in the model are taken from
    the best estimate of the parameter published in other studies.t
