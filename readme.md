## GWnets
Chowdhury, S. and Mémoli, F. The Gromov-Wasserstein distance between networks and stable network invariants. 2018.
https://arxiv.org/abs/1808.04337

mexEMD code from https://github.com/gpeyre/2016-ICML-gromov-wasserstein

_GWnets is a Matlab/Octave package for computing lower bounds on the Gromov-Wasserstein distance between networks._

_To run a simple example, try this (after loading the directories into Matlab)._

k = 30;  
A = rand(k,k);  
mA= (1/k)*ones(k,1);  

res = emd2RTLB(A,A,mA,mA)

_You should get res = 0._

_Next do the following:_

j = 40;  
B = 10*rand(j,j) + 10;  
mB= (1/j)*ones(j,1);  

res = emd2RTLB(A,B,mA,mB)  
res = emd2RTLB(B,A,mB,mA)  

_You should get the same result in both cases._

### Real Data

To play with some real data, look at the migration dataset in data/worldbankdata/. 
The script buildResMigration.m can be used to compute pairwise distances
between migration networks.
