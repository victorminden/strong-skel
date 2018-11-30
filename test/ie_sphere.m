function ie_sphere(n,nquad,occ,p,rank_or_tol,store,method)
% IE_SPHERE  An example usage of strong skeletonization, solving a
%  second-kind integral equation (Laplace double-layer potential) on the 
%  unit sphere.  Sane defaults are provided for all parameters.
%
%  IE_SPHERE(N,NQUAD,OCC,P,RANK_OR_TOL,SKIP,METHOD) runs the example on a 
%  sphere discretized using N points (see code), where the near-field 
%  entries of the discretization are corrected with a Gauss-Legendre 
%  quadrature using NQUAD points in each direction and the other parameters
%  are defined as follows:
%  - OCC:         The occupancy parameter, specifying the maximum number of 
%                 points a node in the octree can contain before it is
%                 subdivided.  This therefore gives an upper bound on the 
%                 number of points in a leaf node.
%  - P:           The number of proxy points to use to discretize the proxy 
%                 sphere used during skeletonization.
%  - RANK_OR_TOL: If a natural number, the maximum number of skeletons to
%                 select during a single box of skeletonization.  If a
%                 float between 0 and 1, an approximate relative tolerance
%                 used to automatically select the number of skeletons.
%  - STORE:       Which interactions to store in the fast multipole method 
%                 used to apply the forward operator.  See the IFMM method
%                 in the FLAM library.
%  - METHOD:      The type of strong skeletonization to use.  Options are
%                 'srskelf' (standard) or 'srskelf_hybrid' (alternating 
%                 strong and weak skeletonization).

  % Set sane default parameters
  if nargin < 1 || isempty(n)
    n = 20480;
  end
  if nargin < 2 || isempty(nquad)
    nquad = 4;
  end
  if nargin < 3 || isempty(occ)
    occ = 256;
  end
  if nargin < 4 || isempty(p)
    p = 512;
  end
  if nargin < 5 || isempty(rank_or_tol)
    rank_or_tol = 1e-3;
  end
  if nargin < 6 || isempty(store)
    store = 'a';
  end
  if nargin < 7 || isempty(method)
    method = 'srskelf';
  end


  % Initialize the discretization points on the sphere and the proxy 
  % sphere.  We use random points on the proxy sphere for simplicity.
  ifmmtol = rank_or_tol/1e3;
  [V,F] = trisphere_subdiv(n);
  [x,nu,area] = tri3geom(V,F);
  N = size(x,2);
  proxy = randn(3,p);
  proxy = 1.5*bsxfun(@rdivide,proxy,sqrt(sum(proxy.^2)));

  % Compute the quadrature corrections
  tic
  if nquad > 0
    % Generate reference transformations for each triangle
    [trans,rot,V2,V3] = tri3transrot(V,F);

    % Initialize quadrature on the unit square
    [xq,wq] = glegquad(nquad,0,1);
    [xq,yq] = meshgrid(xq);
    wq = wq*wq';
    xq = xq(:);
    yq = yq(:);
    wq = wq(:);

    % Find neighbors of each triangle
    nlvl = 0;
    l = 2;
    h = sqrt(8*pi/N);
    while l > h
      nlvl = nlvl + 1;
      l = 0.5*l;
    end
    nlvl = max(1,nlvl);
    T = hypoct(x,0,nlvl);

    % Initialize storage
    nz = 0;
    for i = T.lvp(nlvl)+1:T.lvp(nlvl+1)
      node = T.nodes(i);
      nslf = length(node.xi);
      nnbr = length([T.nodes(node.nbor).xi]);
      nz = nz + nslf*(nslf + nnbr - 1);
    end
    I = zeros(nz,1);
    J = zeros(nz,1);
    S = zeros(nz,1);

    % Compute near-field quadratures
    nz = 0;
    for i = T.lvp(nlvl)+1:T.lvp(nlvl+1)
      node = T.nodes(i);
      nbor = [T.nodes(node.nbor).xi];
      for j = node.xi
        [X,Y,W] = qmap_sqtri2(xq,yq,wq,V2(j),V3(:,j));
        for k = [node.xi nbor]
          if j == k
            continue
          end
          trg = rot(:,:,j)*(x(:,k) + trans(:,j));
          q = W'*quadfun(X,Y,trg);
          nz = nz + 1;
          I(nz) = k;
          J(nz) = j;
          S(nz) = q;
        end
      end
    end
  else
    I = [];
    J = [];
    S = [];
  end
  S = sparse(I,J,S,N,N);
  P = zeros(N,1);
  t = toc;
  w = whos('S');
  fprintf('quad: %10.4e (s) / %6.2f (MB)\n',t,w.bytes/1e6)
  clear V F trans rot V2 V3 T I J

  % Factor the matrix using skeletonization (verbose mode)
  opts = struct('verb',1,'symm','n');
  if strcmp(method,'srskelf')
      F = srskelf_asym(@Afun,x,occ,rank_or_tol,@pxyfun,opts);
  elseif strcmp(method,'srskelf_hybrid')
      F = srskelf_asymhybrid(@Afun,x,occ,rank_or_tol,@pxyfun,opts);
  else
     error('No valid method specified for factorization!');
  end
  w = whos('F');
  fprintf([repmat('-',1,80) '\n'])
  fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)

  % Compress the matrix using interpolative FMM for forward-operator apply
  opts = struct('store',store,'verb',1);
  G = ifmm(@Afun,x,x,1024,ifmmtol,@pxyfun_ifmm,opts);
  w = whos('G');
  fprintf([repmat('-',1,80) '\n'])
  fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)

  % Test accuracy of forward and inverse operator using randomized power
  % method
  X = rand(N,1);
  X = X/norm(X);

  % For the forward operator, we approximate NORM(A - F)/NORM(A)
  tic 
  ntrials = 5;
  for i=1:ntrials
    srskelf_mv_nn(F,X);
  end
  t1 = toc / ntrials;
  tic
  ifmm_mv(G,X,@Afun);
  t2 = toc;
  [e,niter] = snorm(N,@(x)(ifmm_mv(G,x,@Afun,'n') - srskelf_mv_nn(F,x)), ...
                      @(x)(ifmm_mv(G,x,@Afun,'c') - srskelf_mv_nc(F,x)));
  e = e/snorm(N,@(x)(ifmm_mv(G,x,@Afun,'n')),@(x)(ifmm_mv(G,x,@Afun,'c')));
  fprintf('mv: %10.4e / %4d / %10.4e (s) / %10.4e (s)\n',e,niter,t1,t2)

  % For the inverse operator, we approximate the upper bound 
  % NORM(I - A*INV(F)) >= NORM(INV(A) - INV(F))/NORM(INV(A))
  tic
  for i =1:ntrials
    srskelf_sv_nn(F,X);
  end
  t = toc/ntrials;
  [e,niter] = snorm(N,@(x)(x - ifmm_mv(G,srskelf_sv_nn(F,x),@Afun,'n')), ...
                      @(x)(x - srskelf_sv_nc(F,ifmm_mv(G,x,@Afun,'c'))));
  fprintf('sv: %10.4e / %4d / %10.4e (s)\n',e,niter,t)

  % To validate the whole discretization, we generate a field from some
  % exterior sources and then solve for the field inside the sphere.
  m = 16;
  src = randn(3,m);
  src = 2*bsxfun(@rdivide,src,sqrt(sum(src.^2)));
  q = rand(m,1);
  B = Kfun(x,src,'s')*q;

  % Solve for surface density
  X = srskelf_sv_nn(F,B);

  % Evaluate field at interior targets
  trg = randn(3,m);
  trg = 0.5*bsxfun(@rdivide,trg,sqrt(sum(trg.^2)));
  Y = bsxfun(@times,Kfun(trg,x,'d',nu),area)*X;

  % Compare against exact field
  Z = Kfun(trg,src,'s')*q;
  e = norm(Z - Y)/norm(Z);
  fprintf('pde: %10.4e\n',e)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Functions used to generate entries of the matrix to be factorized
  function f = quadfun(x,y,trg)
    % QUADFUN(X,Y,TRG)
    dx = trg(1) - x;
    dy = trg(2) - y;
    dz = trg(3);
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    f = 1/(4*pi).*dz./dr.^3;
  end

  function K = Kfun(x,y,lp,nu)
    % KFUN(X,Y,LP,NU) computes the Laplace potential evaluated
    % pairwise between points in X and points in Y (does not handle the
    % singularity).  If LP is 's' then a single-layer potential is
    % computed, but if LP is 'd' then a double-layer potential is computed
    % using the surface normal vectors in NU.
    if nargin < 4
      nu = [];
    end
    dx = bsxfun(@minus,x(1,:)',y(1,:));
    dy = bsxfun(@minus,x(2,:)',y(2,:));
    dz = bsxfun(@minus,x(3,:)',y(3,:));
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    if strcmpi(lp,'s')
      K = 1/(4*pi)./dr;
    elseif strcmpi(lp,'d')
      rdotn = bsxfun(@times,dx,nu(1,:)) + bsxfun(@times,dy,nu(2,:)) + ...
              bsxfun(@times,dz,nu(3,:));
      K = 1/(4*pi).*rdotn./dr.^3;
    end
    K(dr == 0) = 0;
  end

  function A = Afun(i,j)
    % AFUN(I,J) computes entries of the matrix A to be factorized at the
    % index sets I and J.  This handles the near-field correction.
    if isempty(i) || isempty(j)
      A = zeros(length(i),length(j));
      return
    end
    [I,J] = ndgrid(i,j);
    A = bsxfun(@times,Kfun(x(:,i),x(:,j),'d',nu(:,j)),area(j));
    M = spget(i,j);
    idx = M ~= 0;
    A(idx) = M(idx);
    A(I == J) = -0.5;
  end

  % proxy function
  function [Kpxy,nbr] = pxyfun(x,slf,nbr,l,ctr)
    % PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
    % X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
    % appropriately contain a box at level L centered at CTR and then
    % calling KFUN
    pxy = bsxfun(@plus,proxy*l,ctr');
    Kpxy = bsxfun(@times,Kfun(pxy,x(:,slf),'d',nu(:,slf)),area(slf));
    dx = x(1,nbr) - ctr(1);
    dy = x(2,nbr) - ctr(2);
    dz = x(3,nbr) - ctr(3);
    dist = sqrt(dx.^2 + dy.^2 + dz.^2);
    nbr = nbr(dist/l < 1.5);
  end

  function [Kpxy,nbr] = pxyfun_ifmm(rc,rx,cx,slf,nbr,l,ctr)
    % PXYFUN(RC,RX,CX,SLF,NBR,L,CTR) is analogous to PXYFUN but is used for
    % the interpolative fast multipole method.  RC is either 'r' or 'c' to
    % select whether we are proxying row or column interactions, and rx and
    % cx contain the corresponding points for rows and columns.
    pxy = bsxfun(@plus,proxy*l,ctr');
    if strcmpi(rc,'r')
      Kpxy = Kfun(rx(:,slf),pxy,'s')*(4*pi/N);
      dx = cx(1,nbr) - ctr(1);
      dy = cx(2,nbr) - ctr(2);
    elseif strcmpi(rc,'c')
      Kpxy = bsxfun(@times,Kfun(pxy,cx(:,slf),'d',nu(:,slf)),area(slf));
      dx = rx(1,nbr) - ctr(1);
      dy = rx(2,nbr) - ctr(2);
    end
    dist = sqrt(dx.^2 + dy.^2);
    nbr = nbr(dist/l < 1.5);
  end

  function A = spget(I_,J_)
    % SPGET(I_,J_) computes entries of a sparse matrix of near-field
    % corrections that should be added to the kernel matrix, as used in
    % AFUN.
    m_ = length(I_);
    n_ = length(J_);
    [I_sort,E] = sort(I_);
    P(I_sort) = E;
    A = zeros(m_,n_);
    [I_,J_,S_] = find(S(:,J_));
    idx = ismembc(I_,I_sort);
    I_ = I_(idx);
    J_ = J_(idx);
    S_ = S_(idx);
    A(P(I_) + (J_ - 1)*m_) = S_;
  end
end
