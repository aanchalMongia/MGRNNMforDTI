# MGRNNMforDTI
## A matlab tool for prediction of Drug-target interactions


### Brief Description
This repository has been created to hold the source code for the following paper

[Computational prediction of Drug-Disease association based on Graph-regularized one bit Matrix completion](https://www.biorxiv.org/content/10.1101/2020.04.02.020891v1.abstract)

### About 
In this work, we propose a novel matrix completion framework which makes use of the side-information associated with drugs/diseases
for the prediction of drug-disease indications modelled as neighborhood graph: Graph regularized 1-bit matrix compeltion (GR1BMC).
The algorithm is specially designed for binary data and uses parallel proximal algorithm to solve the aforesaid minimization problem
taking into account all the constraints including the neighborhood graph incorporation and restricting predicted scores within the
specified range.


### Implementation
> main2.m

