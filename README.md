# Comparison Between Lotka-Volterra and Multivariate Autoregressive Models of Ecological Interaction Systems

Abstract:
Lotka-Volterra (LV) and Multivariate Autoregressive (MAR) models are computational frameworks with different mathematical structures that have both been proposed for the same purpose of extracting governing features of dynamic interactions among coexisting populations of different species from observed time series data. We systematically compare the feasibility of the two model structures, using four artificially generated datasets and seven ecological datasets from the literature. We also introduce a new parameter estimation method for LV models, the Algebraic Lotka-Volterra Inference (ALVI) method. We conclude that each of the two approaches is superior in different situations. LV models outperform MAR models in most cases and are generally superior representing cases where the dependent variables deviate greatly from their steady states, when the populations exhibit high abundances with a large signal-to-noise ratio. By contrast, MAR models are better suited for simulations of populations with low abundances and for investigations where the quantification of noise is relevant. 


This repository contains the code used in the paper.
To run it, please check if the working directory is the one you download the files to. You can check this by using the command getwd() and change it using setwd('path to directory').
