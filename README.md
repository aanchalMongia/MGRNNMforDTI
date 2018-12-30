# MGRNNMforDTI
## A matlab tool for prediction of Drug-target interactions


### Brief Description
This repository has been created to hold the source code for the following paper

[Drug-Target Interaction prediction using Multi Graph Regularized Nuclear Norm Minimization](https://www.biorxiv.org/content/early/2018/12/28/455642)

### About 
MGRNNM predicts the interactions between drugs and proteins from three inputs: known drug-target interaction network, similarities over drugs and those over targets. The proposed method focuses on finding a low-rank interaction matrix that is structured by the proximities of drugs and targets encoded by graphs. Previous works on Drug Target Interaction (DTI) prediction have shown that incorporating drug and target similarities helps in learning the data manifold better by preserving the local geometries of the original data. But, there is no clear consensus on which kind and what combination of similarities would best assist the prediction task. Hence, we propose to use various multiple drug-drug similarities and target-target similarities as multiple graph Laplacian (over drugs/targets) regularization terms to capture the proximities exhaustively.


### Implementation
> run.m

### References
1. Code for baselines: [Ezzat A, Wu M, Li XL, Kwoh CK. Computational prediction of drug–target interactions using chemogenomic approaches: an empirical survey.Briefings in bioinformatics. 2018](https://github.com/alizat/Chemogenomic-DTI-Prediction-Methods)
2. Dataset: [Yamanishi Y, Araki M, Gutteridge A, Honda W, Kanehisa M. Prediction of drug–target interaction networks from the integration of chemical and genomic spaces Bioinformatics. 2008;24(13):i232–i240.](http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/).

