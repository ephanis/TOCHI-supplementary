
============================================
author: Theophanis Tsandilas, Inria
date: February 2018
email: theophanis.tsandilas@inria.fr

Supplementary material for TOCHI article: Fallacies of Agreement: A Critical Review of Consensus Assessment Methods for Gesture Elicitation

Project page: https://agreement.lri.fr
============================================

Some of these scripts require the installation of additional R libraries. 
You can install the following libraries to avoid any error messages:

- install.packages("coin")
- install.packages("boot")
- install.packages("ggplot2")
- install.packages("zipfR")
- install.packages("partitions")
- install.packages("Rmisc")


===============================================
Agreement Coefficients and Confidence Intervals 
===============================================
+ coefficients/ 
- agreement.coefficients.R : Code for agreement coefficients (Based on implementation by Kilem Gwet)
- agreement.CI.R : Code for calculating confidence intervals
- specific.agreement.CI.R : Code for calculating specific agreement 


==================================
Section 3 - Measuring Agreement
==================================
- section3.bias.distribution.plots.R
Used to produce the bias probability distributions of Figure 3

- section3.bias.monte.carlo.R
Used for Experiments 3.1 and 3.2


==================================
Section 5 - Magnitude of Agreement
==================================
- section5.vatavu2015.partitions.R
Used for Table IV. It derives possible participant partitions according to Vatavu and Wobbrock (2015). It can also calculate their frequencies (number of configurations for each partition) and calculate the mean agreement rate (AR). Finally, it produces the probability distribution presented in Fig. 4.
  
- section5.agreement.distribution.monte.carlo.R
Used for Fig. 5. It produces the null probability distribution of agreement rates (AR) of Fleiss' Kappa through Monte Carlo simulations. 


===========================================================================
Section 6 - Statistical Inference
NOTE: Some evaluation scripts can take quite long to complete (several hours)
===========================================================================
- section6.delta.distributions.R
Generates the distributions of agreement differences in Figure 6.


- section6.within.vatavu2015.Vrd.bailly et al.R
Tests our implementation of the Vrd test on Bailly's dataset. Our implementation accurately replicates Vatavu and Wobbrock (2015)'s reported statistic value (chi-squared) for the overall effect of referent type: Vrd(41) = 1466.818 (p < .001)

- section6.within.vatavu2015.Vrd.evaluation.zipf.R
Code for Experiment 6.1 (Zipf-Mandelbrot distributions)

- section6.within.vatavu2015.Vrd.evaluation.normal.R
Code for Experiment 6.1 (Half-Normal distributions)

- section6.within.vatavu2015.Vrd.evaluation.normal+zipf.R
Code for Experiment A.2 (on linear combinations of Half-Normal and Zipf-Mandelbrot distributions)

- section6.within.vatavu2015.Vrd.evaluation.bailly.R
Code for Experiment 6.2 (on populations generated with Bailly et al.'s dataset)



- section6.between.vatavu2016.Vb.R
Implements Vatavu and Wobbrock's (2016) Vb statistical test for two independent groups

- section6.between.vatavu2016.Vb.evaluation.zipf.R
Code for Experiment 6.3 (Zipf-Mandelbrot distributions)

- section6.between.vatavu2016.Vb.evaluation.normal.R
Code for Experiment 6.3 (Half-Normal distributions)

- section6.between.vatavu2016.Vb.evaluation.normal+zipf.R
Code for Experiment A.5 (on linear combinations of Half-Normal and Zipf-Mandelbrot distributions)

- section6.within.vatavu2015.Vb.evaluation.bailly.R
Code for Experiment 6.4 (on populations generated with Bailly et al.'s dataset)



- section6.within.jackknife.evaluation.zipf.R
Code for Experiment 6.5 (Zipf-Mandelbrot distributions)

- section6.within.jackknife.evaluation.normal.R
Code for Experiment 6.5 (Half-Normal distributions)


- section6.between.bootstrap.R
Implements the bootstrap method for comparing agreement rates of two independent groups

- section6.between.bootstrap.evaluation.zipf.R
Code for Experiment 6.6 (Zipf-Mandelbrot distributions)

- section6.between.bootstrap.evaluation.normal.R
Code for Experiment 6.6 (Half-Normal distributions)



-section6.jackknife.evaluation.coverage.R
Code for Experiment 6.7

- section6.between.bootstrap.kappa.R
Implements the bootstrap method for comparing Fleiss' Kappa of two independent groups

-section6.between.bootstrap.evaluation.bailly.R
Code for Experiment 6.8



- section6.within.bootstrap.evaluation.R
Code for Experiment A.7 (Not reported in the main article)


============================================================
Section 7 - Re-Analysis of Past Gesture Elicitation Studies
NOTE: We only include Bailly et al.'s (2013) data.
I have no permission to share other datasets.
=============================================================
+ data/bailly et al 2013.R - dataset.csv
Dataset of gesture elicitation study by Bailly et al. (2013) 

- section7.bailly et al 2013 - basic analysis.R
Basic analysis of agreement for Bailly et al's dataset. 

- section7.bailly et al 2013 - specific agreement.R
Analysis of specific agreement for Bailly et al's dataset. 

- section7.bailly et al 2013 - women vs men.Vb.R
Using the Vb test for the analysis of agreement differences between women and men for Bailly et al's dataset.

- section7.bailly et al 2013 - women vs men.bootstrap.R
Using the bootstrap method for the analysis of agreement differences between women and men for Bailly et al's dataset. It can be very slow.



================
Acknowledgments
================
Pierre Dragicevic has contributed code for the calculation of agreement indices and their confidence intervals. 
See: Tsandilas, T., and Dragicevic, P. (2016). Accounting for Chance Agreement in Gesture Elicitation Studies. Research Report 1584, LRI - CNRS, University Paris-Sud, Feb 2016.




