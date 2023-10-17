// replicate the parameter search many times
quietly {
  clear
  matrix final_results = J(1, 7, .)
  forval i = 1/$reps {
    quietly parameter_search params

  	matrix final_results = final_results \ r(RE)
  }

  // rename the columns
  clear
  svmat final_results, names(final_results)
  rename final_results1 n_subj
  rename final_results2 n_pop
  rename final_results3 n_rock
  rename final_results4 beta_1
  rename final_results5 mean_estimate
  rename final_results6 mean_se
  rename final_results7 p_value
  drop in 1

  save "./data/final_results.dta", replace
}
