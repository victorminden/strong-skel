function ie_square(n,occ,p,rank_or_tol,method)
% IE_SQUARE  An example usage of strong skeletonization, solving a
%  first-kind integral equation (Laplace single-layer potential) on the 
%  unit square.  Sane defaults are provided for all parameters.
%
%  IE_SQUARE(N,OCC,P,RANK_OR_TOL,SKIP,METHOD) runs the example on a regular
%  N-by-N grid with the other parameters defined as follows:
%  - OCC:         The occupancy parameter, specifying the maximum number of 
%                 points a node in the quadtree can contain before it is
%                 subdivided.  This therefore gives an upper bound on the 
%                 number of points in a leaf node.
%  - P:           The number of proxy points to use to discretize the proxy 
%                 circle used during skeletonization.
%  - RANK_OR_TOL: If a natural number, the maximum number of skeletons to
%                 select during a single box of skeletonization.  If a
%                 float between 0 and 1, an approximate relative tolerance
%                 used to automatically select the number of skeletons.
%  - METHOD:      The type of strong skeletonization to use.  Options are
%                 'srskelf' (standard) or 'srskelf_hybrid' (alternating 
%                 strong and weak skeletonization).

  % Set default parameters to sane values
  if nargin < 1 || isempty(n)
    n = 32;
  end
  if nargin < 2 || isempty(occ)
    occ = 64;
  end
  if nargin < 3 || isempty(p)
    p = 64;
  end
  if nargin < 4 || isempty(rank_or_tol)
    rank_or_tol = 1e-6;
  end
  if nargin < 5 || isempty(method)
    method = 'srskelf';
  end

  % Initialize the grid of points and the proxy circle
  [x1,x2] = ndgrid((1:n)/(n));
  x = [x1(:) x2(:)]';
  N = size(x,2);
  theta = (1:p)*2*pi/p;
  proxy = 1.5*[cos(theta); sin(theta)];
  clear x1 x2

  % Compute the simple adaptive diagonal quadrature term
  h = 1/n;
  intgrl = 4*dblquad(@(x,y)(-1/(2*pi)*log(sqrt(x.^2 + y.^2))),0,h/2,0,h/2);
  
  % Factor the matrix using skeletonization (verbose mode)
  opts = struct('verb',1);
  
  if strcmp(method,'srskelf')
    F = srskelf(@Afun,x,occ,rank_or_tol,@pxyfun,opts);
  elseif strcmp(method,'srskelf_hybrid')
    F = srskelf_hybrid(@Afun,x,occ,rank_or_tol,@pxyfun,opts);
  else
      error('No valid method specified for factorization!');
  end
  w = whos('F');
  fprintf([repmat('-',1,80) '\n'])
  fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)

  % Because the problem is discretized on a regular grid, the forward
  % operator can be applied efficiently using the fast Fourier transform.
  a = reshape(Afun(1:N,1),n,n);
  B = zeros(2*n-1,2*n-1);
  B(  1:n  ,  1:n  ) = a;
  B(  1:n  ,n+1:end) = a( : ,2:n);
  B(n+1:end,  1:n  ) = a(2:n, : );
  B(n+1:end,n+1:end) = a(2:n,2:n);
  B(:,n+1:end) = flipdim(B(:,n+1:end),2);
  B(n+1:end,:) = flipdim(B(n+1:end,:),1);
  G = fft2(B);

  % Test accuracy of forward and inverse operator using randomized power
  % method
  X = rand(N,1);
  X = X/norm(X);

  % For the forward operator, we approximate NORM(A - F)/NORM(A)
  tic
  ntrials = 5;
  for i=1:ntrials
    srskelf_mv_p(F,X);
  end
  t = toc / ntrials;

  [e,niter] = snorm(N,@(x)(mv(x) - srskelf_mv_p(F,x)),[],[],1);
  e = e/snorm(N,@mv,[],[],1);
  fprintf('mv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)

  % For the inverse operator, we approximate the upper bound 
  % NORM(I - A*INV(F)) >= NORM(INV(A) - INV(F))/NORM(INV(A))
  tic
  ntrials = 5;
  for i=1:ntrials
    srskelf_sv_p(F,X);
  end
  t = toc / ntrials;
  [e,niter] = snorm(N,@(x)(x - mv(srskelf_sv_p(F,x))), ...
                      @(x)(x - srskelf_sv_p(F,mv(x))));
  fprintf('sv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)
  
  
  % To evaluate performance as a preconditioner, we first run 
  % unpreconditioned CG to get a baseline iteration count
  Y = rand(N,1);
  Y = Y/norm(Y);
  X = mv(Y);
  [~,~,~,iter] = pcg(@mv,X,1e-12,128);

  % Now we run preconditioned PCG to see the effectiveness of the 
  % preconditioner
  tic
  [Z,~,~,piter] = pcg(@mv,X,1e-12,32,@(x)(srskelf_sv_p(F,x)));
  t = toc;
  e1 = norm(Z - Y)/norm(Z);
  e2 = norm(X - mv(Z))/norm(X);
  fprintf('cg: %10.4e / %10.4e / %4d (%4d) / %10.4e (s)\n',e1,e2, ...
          piter,iter,t)
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Functions used to generate entries of the matrix to be factorized
  function K = Kfun(x,y)
    % KFUN(X,Y) computes the Laplace single layer potential evaluated
    % pairwise between points in X and points in Y (does not handle the
    % singularity)  
    dx = bsxfun(@minus,x(1,:)',y(1,:));
    dy = bsxfun(@minus,x(2,:)',y(2,:));
    K = -1/(2*pi)*log(sqrt(dx.^2 + dy.^2));
  end

  function A = Afun(i,j)
    % AFUN(I,J) computes entries of the matrix A to be factorized at the
    % index sets I and J.  This handles the singularity at the diagonal.
    A = Kfun(x(:,i),x(:,j))/N;
    [I,J] = ndgrid(i,j);
    A(I == J) = intgrl;
  end

  function [Kpxy,nbr] = pxyfun(x,slf,nbr,l,ctr)
    % PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
    % X(:,SLF) and the set of proxy points by scaling the proxy circle to 
    % appropriately contain a box at level L centered at CTR and then
    % calling KFUN
    pxy = bsxfun(@plus,proxy*l,ctr');
    Kpxy = Kfun(pxy,x(:,slf))/N;
  end

  function y = mv(x)
    % MV(X) implements multiplication by the matrix A to be factorized by
    % using the FFT to exploit grid structure
    y = ifft2(G.*fft2(reshape(x,n,n),2*n-1,2*n-1));
    y = reshape(y(1:n,1:n),N,1);
  end
end