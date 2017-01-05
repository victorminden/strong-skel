% First-kind integral equation on the unit square, Laplace single-layer.

function ie_square1_s(n,occ,p,rank_or_tol,skip,method)
  rng(1);
  % set default parameters
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
  if nargin < 5 || isempty(skip)
    skip = 0;
  end
  if nargin < 6 || isempty(method)
    method = 'srskelf';
  end

  % initialize
  [x1,x2] = ndgrid((1:n)/(n));
  x = [x1(:) x2(:)]';
  N = size(x,2);
  theta = (1:p)*2*pi/p;
  proxy = 1.5*[cos(theta); sin(theta)];
  clear x1 x2

  % compute diagonal quadratures
  h = 1/n;
  intgrl = 4*dblquad(@(x,y)(-1/(2*pi)*log(sqrt(x.^2 + y.^2))),0,h/2,0,h/2);
  % factor matrix
  opts = struct('verb',1,'skip',skip);
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

  % set up FFT multiplication
  a = reshape(Afun(1:N,1),n,n);
  B = zeros(2*n-1,2*n-1);
  B(  1:n  ,  1:n  ) = a;
  B(  1:n  ,n+1:end) = a( : ,2:n);
  B(n+1:end,  1:n  ) = a(2:n, : );
  B(n+1:end,n+1:end) = a(2:n,2:n);
  B(:,n+1:end) = flipdim(B(:,n+1:end),2);
  B(n+1:end,:) = flipdim(B(n+1:end,:),1);
  G = fft2(B);

  % test accuracy using randomized power method
  X = rand(N,1);
  X = X/norm(X);

  % NORM(A - F)/NORM(A)
  tic
  ntrials = 5;
  for i=1:ntrials
    srskelf_mv_p(F,X);
  end
  t = toc / ntrials;

  [e,niter] = snorm(N,@(x)(mv(x) - srskelf_mv_p(F,x)),[],[],1);
  e = e/snorm(N,@mv,[],[],1);
  fprintf('mv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)

  % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
  tic
  ntrials = 5;
  for i=1:ntrials
    srskelf_sv_p(F,X);
  end
  t = toc / ntrials;
  [e,niter] = snorm(N,@(x)(x - mv(srskelf_sv_p(F,x))),[],[],1);
  fprintf('sv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)
  
  % run unpreconditioned PCG
  Y = rand(N,1);
  Y = Y/norm(Y);
  X = mv(Y);
[~,~,~,iter] = pcg(@mv,X,1e-12,128);

% run preconditioned PCG
  tic
  [Z,~,~,piter] = pcg(@mv,X,1e-12,32,@(x)(srskelf_sv_p(F,x)));
  t = toc;
  e1 = norm(Z - Y)/norm(Z);
  e2 = norm(X - mv(Z))/norm(X);
  fprintf('cg: %10.4e / %10.4e / %4d (%4d) / %10.4e (s)\n',e1,e2, ...
          piter,iter,t)

  % kernel function
  function K = Kfun(x,y)
    dx = bsxfun(@minus,x(1,:)',y(1,:));
    dy = bsxfun(@minus,x(2,:)',y(2,:));
    K = -1/(2*pi)*log(sqrt(dx.^2 + dy.^2));
  end

  % matrix entries
  function A = Afun(i,j)
    A = Kfun(x(:,i),x(:,j))/N;
    [I,J] = ndgrid(i,j);
    A(I == J) = intgrl;
  end

  % proxy function
  function [Kpxy,nbr] = pxyfun(x,slf,nbr,l,ctr)
    pxy = bsxfun(@plus,proxy*l,ctr');
    Kpxy = Kfun(pxy,x(:,slf))/N;
  end

  % FFT multiplication
  function y = mv(x)
    y = ifft2(G.*fft2(reshape(x,n,n),2*n-1,2*n-1));
    y = reshape(y(1:n,1:n),N,1);
  end
end