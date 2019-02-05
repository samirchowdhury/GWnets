## GWnets

_GWnets is a Matlab/Octave package for computing lower bounds on the Gromov-Wasserstein distance between networks._

_To run a simple example, try this._

k = 30;  
A = rand(k,k);  
mA= (1/k)*ones(k,1);  

res = get2RTLB(A,A,mA,mA,200,50)

_You should get res = 0._

_Next do the following:_

j = 40;  
B = rand(j,j);  
mB= (1/j)*ones(j,1);  

res = get2RTLB(A,B,mA,mB,200,50)  
res = get2RTLB(B,A,mB,mA,200,50)  

_You should get the same result in both cases._
