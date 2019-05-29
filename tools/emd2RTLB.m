function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB(A,B,mA,mB)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % mA, mB are the measures
  % No entropic regularization, just EMD
  
  
   
  
  % RTLB consists of an LOP over a cost matrix, each entry of which is given
  % by a formula.
  % For inner LOPs, need to get inner cost matrices.
  % Do this by looping over rows of A and rows of B (ecc-out case)
  
  n     = size(A,1);
  m     = size(B,1);
  
  [r,c] = meshgrid(1:n,1:m);
  
  getOuterCostAnon    = @(ii,jj) getOuterCost(ii,jj,A,B,mA,mB);
  
  [eccout_cst,eccin_cst] = arrayfun(getOuterCostAnon,r,c);
  eccout_cst             = eccout_cst';
  eccin_cst              = eccin_cst';
  
  % using 2-TLB, get squared cost matrix
  eccout_cst             = eccout_cst.^2;
  eccin_cst              = eccin_cst.^2;
  
  
  % EMD
  [dist_out,gamma_out] = mexEMD(mA,mB,eccout_cst);
  [dist_in,gamma_in] = mexEMD(mA,mB,eccin_cst);
  
  %gamma_out
  %gamma_in
  
  % square root
  dist_out = sqrt(dist_out);
  dist_in  = sqrt(dist_in);
  
  % calculate best lower bound
  res_out = dist_out/2;
  res_in  = dist_in/2;
  res = max(res_out,res_in);
    
  
  end






