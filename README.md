# Inference of Dynamic Interaction Networks:
A Comparison Between Lotka-Volterra and Multivariate 
Autoregressive Models

Abstract:
Networks are ubiquitous throughout biology, spanning the entire range from molecules to food webs 
and global environmental systems. Yet, despite substantial efforts by the scientific community, the 
inference of these networks from data still presents a problem that is unsolved in general. 
One frequent strategy of addressing the structure of networks is the assumption that the interactions 
among molecular or organismal populations are static and correlative. 
While often successful, these methods are no panacea. They usually ignore the asymmetry of 
relationships between two species and inferences become more challenging if the network nodes 
represent dynamically changing quantities. Overcoming these challenges, two very different network 
inference approaches have been proposed in the literature: Lotka-Volterra (LV) models and 
Multivariate Autoregressive (MAR) models.
These models are computational frameworks with different mathematical structures which, 
nevertheless, have both been proposed for the same purpose of inferring the interactions within 
coexisting population networks from observed time-series data. Here, we assess these dynamic 
network inference methods for the first time in a side-by-side comparison, using both synthetically 
generated and ecological datasets.
MAR and LV models are mathematically equivalent at the steady state, but the results of our 
comparison suggest that LV models are generally superior in capturing the dynamics of networks with 
non-linear dynamics, whereas MAR models are better suited for analyses of networks of populations 
with process noise and close-to linear behavior.



This repository contains the code used in the paper.
To run it, please check if the working directory is the one you download the files to. You can check this by using the command getwd() and change it using setwd('path to directory').
