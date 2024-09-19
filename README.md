# BayesGO2DO
A Bayesian genomic selection framework incorporating gene ontology (GO) and disease ontology (DO) information

An R package for a Bayesian genomic selection framework. The framework incorporates GO and DO as prior biological information when performing genomic selections in multiple cancer subtypes. 
The framework can account for functional similarities of genes and ontology similarities of diseases (i.e., cancer subtypes). The response variable in the framework is the (censored) survival time measurement.

The package can be installed by running:

> devtools::install_github("xiaotiand/GO2DO")

Then run 

> library(GO2DO)
> 
> ?Main

to get an example. The example simulate a simple dataset to test the function. Also, see

> ?Predict

for the predicting function.



Citations:

To be added
