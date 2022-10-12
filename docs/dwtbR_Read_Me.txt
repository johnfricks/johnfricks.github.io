The R package dwtbR contains R code to implement the stochastic EM algorithm described in Bernstein and Fricks (2015).  
R is an open source statistical computing software that can be downloaded at https://cran.r-project.org/.  
R Studio (https://www.rstudio.com/) is a popular integrated development environment for R that can also be used as an R interface.  

After the package has been downloaded, the package can be installed in R Studio by going to 
Tools -> Install Packages.  The package is loaded in R with the command library(dwtbR).

The package contains five functions:

1) simPath          # Simulates a trajectory of the model
2) particle_filter  # Runs the particle filter
3) update_params    # Implements the M step
4) stochastic_EM    # Runs the stochastic EM algorithm
5) compute_OIM      # Computes the observed information matrix

The command ?foo in R will open documentation and examples for a particular function.  For 
example, to see the documentation and example for simulating a path, type ?simPath.

For questions or comments about the R package, please contact: 
Jason Bernstein (jib5317@psu.edu).
10/27/2015