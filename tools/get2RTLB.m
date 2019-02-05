function [res, eccout, eccin, mu_eccout, mu_eccin] = get2RTLB(A,B,mA,mB,niter,sharpen_by, scaling_bd)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % mA, mB are the measures
  % This Adaptive TLB code determines the regularization parameter
  % lambda from A, B, mA, and mB.   
  % niter is the number of sinkhorn iterations
  % sharpen_by is typically an integer >= 1; represents the orders of magnitude
  % of difference (roughly) between the top two entries of the exponentiated
  % matrix. In the case of TLB(A,A) type problems, want sharpen_by >= 3.
  
  % scaling_bd represents the orders of magnitude up to which the scaling terms
  % can grow without warning
  
  

  
  if nargin < 5 || isempty(niter)
    niter = 50;
  end
  
  
  
  if nargin < 6 || isempty(sharpen_by)
    sharpen_by = 3;
  end
  
  if nargin < 7 || isempty(scaling_bd)
    scaling_bd = 200;
  end
  
  
  res = 0;
  
  
  
  % save original inputs
  origA     = A;
  origB     = B;
  
  % get the 2-sizes of A and B, also largest absolute weights
  sizeOrigA  = getDiam(origA,mA);
  sizeOrigB  = getDiam(origB,mB);
  
  maxPosA   = max(max(abs(origA)));
  maxPosB   = max(max(abs(origB)));
  
  
  % rescale inputs and get new 2-size
  A         = A/(maxPosA);
  B         = B/(maxPosB);
  %sizeA     = getDiam(A,mA);
  %sizeB     = getDiam(B,mB);
  
  
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
  
  
  % outer sinkhorn
  
  lambda_out  = getBestLambda(eccout_cst,mA,mB,sharpen_by,scaling_bd);
  lambda_in   = getBestLambda(eccin_cst,mA,mB,sharpen_by,scaling_bd);
  
  K_outer     = exp(-lambda_out*eccout_cst);
  K_inner     = exp(-lambda_in*eccin_cst);
  mu_eccout   = do_sinkhorn(mA,mB,K_outer,niter);
  mu_eccin    = do_sinkhorn(mA,mB,K_inner,niter);
  
  % calculate TLB cost
  opt_outer_cst   = sqrt(sum(sum(eccout_cst.*mu_eccout)));
  opt_inner_cst   = sqrt(sum(sum(eccin_cst.*mu_eccin)));
  
  % get best lower bound for scaled problem
  % use both in and out eccentricity to get best bound
  scaled_tlb          = 0.5*max(opt_outer_cst,opt_inner_cst);
  res                 = scaled_tlb;
  
  
  
  % finally rescale to get 0.5*tlb for original problem
  
  s   = 0.5*sizeOrigA;
  t   = 0.5*sizeOrigB;
  res = maxPosA*maxPosB*(scaled_tlb^2) - (s^2*maxPosB)/(maxPosA) ...
        - (t^2*maxPosA)/(maxPosB) + s^2 + t^2;
  res = sqrt(max(0,res));
  
  end


function diam = getDiam(A,mA)
  % get the squared 2-size of a measure network
  % writing diam instead of size to avoid confusion with matlab size
  A = abs(A);
  diam = mA'*(A.^2)*mA;
  diam = sqrt(diam);
  end

function [mu,u,v] = do_sinkhorn(mA,mB,K,niter)
  u = ones(size(mA));
  for ii = 1:niter
    v = mB./(K'*u);
    u = mA./(K*v);
  end
  mu = diag(u)*K*diag(v);  
  end


function v = getSmallestEntries(cost)
  % return 2 (unique) smallest entries in each column of cost
  % if 2 unique smallest entries not possible, return just 1 row.
  num_rows  = size(cost,1);
  min_vals  = min(cost);
  dummy     = repmat(min_vals,[num_rows,1]);
  zero_out  = dummy - cost;
  mod_cost  = cost.*logical(zero_out);
  mod_cost(mod_cost == 0) = inf;
  second_min= min(mod_cost);
  
  v         = min_vals;
  if (sum(isinf(second_min)) == 0)
    v         = [min_vals;second_min];
    end
  end

  
function lambda = getBestLambda(cost, mA,mB, sharpen_by, scaling_bd)
  
  col_lambdas = [];
  row_lambdas = [];
  
  % get smallest 2 elements in each col
  col_smallest  = getSmallestEntries(cost);
  row_smallest  = getSmallestEntries(cost');
  
  % get differences (if possible) and candidate lambdas
  if size(col_smallest,1) == 2 
    col_diffs     = col_smallest(2,:) - col_smallest(1,:);
    col_lambdas   = (sharpen_by*log(10)*ones(size(col_diffs)))./(col_diffs);
  end
  if size(row_smallest,1) == 2
    row_diffs     = row_smallest(2,:) - row_smallest(1,:);
    row_lambdas   = (sharpen_by*log(10)*ones(size(row_diffs)))./(row_diffs);
  end
  
  lambda_pool     = [200,col_lambdas,row_lambdas];
  lambda_pool = sort(lambda_pool);
  
  %lambda_found= 0; 
  %ii = 0;
  
  lambda = lambda_pool(1);
  %% Binary search
  L = 1;
  R = length(lambda_pool);
  while L ~= R
    midterm     = ceil((L+R)/2);
    lambda      = lambda_pool(midterm);
    test_cst    = exp(-lambda*cost);
    max_col     = -getSmallestEntries(-test_cst);
    max_row     = -getSmallestEntries(-test_cst');
    if  max(mA'./max_row(1,:)) > 10^scaling_bd || max(mB'./max_col(1,:)) > 10^scaling_bd
      R         = midterm - 1;
    else
      L         = midterm;
    end
  end
  lambda = lambda_pool(L);
  
end


  
function [eccout_cst_ii_jj,eccin_cst_ii_jj] = getOuterCost(ii,jj,...
                                            A,B,mA,mB)
                                    
  vA_out = A(ii,:);
  vB_out = B(jj,:);
  
  vA_in  = A(:,ii);
  vB_in  = B(:,jj);

  % fill in outer cost matrix
  eccout_cst_ii_jj = compareRealDistributions(vA_out,vB_out,mA,mB);
  eccin_cst_ii_jj  = compareRealDistributions(vA_in,vB_in,mA,mB);
end


function dist = compareRealDistributions(vA,vB,mA,mB)
  % computes the 2-TLB linear program over the real line using generalized inverse     
  [sorted_vA,sorted_idA] = sort(vA);
  [sorted_vB,sorted_idB] = sort(vB);
  
  sorted_mA = mA(sorted_idA);
  sorted_mB = mB(sorted_idB);
  
    
  cmA = cumsum(sorted_mA);
  cmB = cumsum(sorted_mB);
  
  % Integrating over the interval [0,1], 
  % functions of interest are the CDFs.
  % Only care about the places where cmA, cmB change
  sorted_all = union(cmA,cmB);
  
  %next stop prevents bugs from rounding errors, needs Matlab to run
  %sorted_all = uniquetol(sorted_all, 10^(-10));
  sorted_all = [0;sorted_all];
  
  summand = 0;
  for ii = 1:length(sorted_all) - 1
    % find first index where cmA exceeds sorted_all(ii)
    % use the index to get the generalized inverse of the distribution function
    idxA = find(cmA > sorted_all(ii) , 1);
    ginvA = sorted_vA(idxA);
    % likewise for B
    idxB = find(cmB > sorted_all(ii) , 1);
    ginvB = sorted_vB(idxB);
    
    riem_int = (abs(ginvA - ginvB))^2; %2-TLB here
    summand  = summand + riem_int*(sorted_all(ii+1) - sorted_all(ii));
  end
  dist = sqrt(summand);
end

