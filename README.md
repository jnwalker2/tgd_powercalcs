# tgd_powercalcs
Code in support of the manuscript 'Statistical reasons to include transgender and gender diverse sex variables in health research'.

This repository contains two scripts:

power_calc_script: A script which computes sample size under each model (continuous outcome and binary outcome with no misspecification, two-sided misspecification and one-sided misspecification) and power under each model when assuming the baseline models to be true. This code uses the method of simulated power analysis to compute power given a fixed sample size and estimates the minimum sample size necessary to acheive a given power via a bisection method. 

multi_samp_size: A script similar to power_calc_script except that it computes the sample size estimates multiple times, allowing for uncertainty in the estimates to be captured. 
