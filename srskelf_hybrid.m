function F = srskelf_hybrid(A,x,occ,rank_or_tol,pxyfun,opts)
% SRSKELF_HYBRID   Strong recursive skeletonization factorization with
% hybrid admissibility (symmetric positive definite only).
%
%    F = SRSKELF_HYBRID(A,X,OCC,RANK_OR_TOL,PXYFUN) produces a 
%    factorization F of the interaction matrix A on the points X using tree 
%    occupancy parameter OCC, local precision parameter RANK_OR_TOL, and 
%    proxy function PXYFUN to capture the far field. This is a function of 
%    the form
%
%      [KPXY,NBR] = PXYFUN(X,SLF,NBR,L,CTR)
%
%    that is called for every block, where
%
%      - KPXY: interaction matrix against artificial proxy points
%      - NBR:  block neighbor indices (can be modified)
%      - X:    input points
%      - SLF:  block indices
%      - L:    block size
%      - CTR:  block center
%
%    See the examples for further details.
%
%    F = SRSKELF_HYBRID(A,X,OCC,RANK_OR_TOL,PXYFUN,OPTS) also passes 
%    various options to the algorithm. Valid options include:
%
%      - EXT: set the root node extent to [EXT(I,1) EXT(I,2)] along 
%             dimension I. If EXT is empty (default), then the root extent 
%             is calculated from the data.
%
%      - LVLMAX: maximum tree depth (default: LVLMAX = Inf).
%
%      - VERB: display status of the code if VERB = 1 (default: VERB = 0).

  start = tic;
  % Set sane default parameters
  if nargin < 5
    pxyfun = [];
  end % if
  if nargin < 6
    opts = [];
  end % if
  if ~isfield(opts,'ext')
    opts.ext = [];
  end % if
  if ~isfield(opts,'lvlmax')
    opts.lvlmax = Inf;
  end % if
  if ~isfield(opts,'verb')
    opts.verb = 0;
  end % if
  if ~isfield(opts,'skip')
    opts.skip = 1;
  end % if

  if opts.verb
    disp(['This is symmetric positive definite srskelf with hybrid',  ...
        ' admissibility (RS-WS).']);
    disp('Diagonal blocks will be factorized with Cholesky.');
    
  end % if

  % Build tree to hold the discretization points
  N = size(x,2);
  tic
  t = shypoct(x,occ,opts.lvlmax,opts.ext);

  if opts.verb
    fprintf(['-'*ones(1,80) '\n'])
    fprintf('%5s | %6s | %8s | %8s | %8s | %8s | %10s (s)\n', ...
              'lvl','nblk','nRemIn','nRemOut','inRatio','outRatio','time')
    % Print summary information about tree construction
    fprintf(['-'*ones(1,80) '\n'])
    fprintf('  %3s | %63.2e (s)\n','-',toc)

    % Count the nonempty boxes at each level
    pblk = zeros(t.nlvl+1,1);
    for lvl = 1:t.nlvl
      pblk(lvl+1) = pblk(lvl);
      for i = t.lvp(lvl)+1:t.lvp(lvl+1)
        if ~isempty(t.nodes(i).xi)
          pblk(lvl+1) = pblk(lvl+1) + 1;
        end % if
      end % for
    end % for
  end % if

  % Initialize the data structure holding the factorization
  nbox = t.lvp(end);

  e = cell(nbox,1);
  % Each element of F.factors will contian the following data for one box:
  %   - sk: the skeleton DOF indices
  %   - rd: the redundant DOF indices
  %   - nbr: the neighbor (near-field) DOF indices
  %   - T: the interpolation matrix mapping redundant to skeleton
  %   - E: the left factor of the (symmmetric) Schur complement update to 
  %        sk
  %   - L: the Cholesky factor of the diagonal block
  %   - C: the left factor of the (symmetric) Schur complement update to
  %         nbr
  F = struct('sk',e,'rd',e,'nbr',e,'T',e,'E',e,'L',e,'C',e);
  F = struct('N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F);
  nlvl = 0;
  n = 0;
  % Mark every DOF as "remaining", i.e., not yet eliminated.
  rem = true(N,1);
  lookup_list = zeros(nbox,2);
  

  % Loop over the levels of the tree from bottom to top
  for lvl = t.nlvl:-1:1
    % For each box, pull up information about skeletons from child boxes
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      t.nodes(i).xi = [t.nodes(i).xi [t.nodes(t.nodes(i).chld).xi]];
    end % for
    % We factorize both half-integer levels for every level lower than 2 
    % except the first level otherwise we just have weak skeletonization
    if lvl <= 2
        ub = 1;
    else
        ub = 2;
    end % if
    if lvl == t.nlvl
        lb = 2;
    else
        lb = 1;
    end % if
    % Loop over half-integer levels
    for pass = lb:ub
      time = tic;
      nlvl = nlvl + 1;
      nrem1 = sum(rem);
        
      % Loop over each box in this level
      for i = t.lvp(lvl)+1:t.lvp(lvl+1)
        slf = t.nodes(i).xi;
        nbr = [t.nodes(t.nodes(i).nbor).xi];

        nslf = length(slf);
        slf = sort(slf);

        nnbr = length(nbr);
        nbr = sort(nbr);

        if pass == 1
          lst = nbr;
          l = t.lrt/2^(lvl - 1);
        else
          lst = [t.nodes(t.nodes(i).ilist).xi];
          l = t.lrt/2^(lvl - 1) * 5/3;
        end % if

        % compute proxy interactions and subselect neighbors
        Kpxy = zeros(0,nslf);
        if lvl > 2
          [Kpxy,lst] = pxyfun(x,slf,lst,l,t.nodes(i).ctr);
        end % if

        nlst = length(lst);
        lst = sort(lst);


        %compute interaction matrix
        K1 = full(A(lst,slf));
        K2 = spget('lst','slf');
          
        K = [K1 + K2; Kpxy];
        %skeletonize
        [sk,rd,T] = id(K,rank_or_tol);

        % move on if no compression
        if isempty(rd)
          continue
        end % if

        % compute factors
        K  = full(A(slf,slf)) + spget('slf','slf');
        if pass == 2
          K2 = full(A(nbr,slf)) + spget('nbr','slf');
          K2(:,rd) = K2(:,rd) - K2(:,sk)*T; 
        else
          K2 = [];
        end % if

        K(rd,:) =  K(rd,:) - T'*K(sk,:);
        K(:,rd) = K(:,rd) - K(:,sk)*T;

        L = chol(K(rd,rd),'lower');
        E = K(sk,rd)/L';
        if pass == 1
          C = zeros(0,length(rd));
        else
          C = K2(:,rd)/L';
        end % if


        % store matrix factors
        n = n + 1;
        F.factors(n).sk  = slf(sk);
        F.factors(n).rd  = slf(rd);
        if pass == 2
          F.factors(n).nbr = nbr;
        else
          F.factors(n).nbr = [];
        end % if
            
        F.factors(n).T = T;
        F.factors(n).E = E;
        F.factors(n).L = L;
        F.factors(n).C = C;

        lookup_list(i,pass) = n;

        t.nodes(i).xi = slf(sk);
        rem(slf(rd)) = 0;
      end % for
      F.lvp(nlvl+1) = n;

      % print summary
      if opts.verb
        nrem2 = sum(rem);
        nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);
        fprintf('%3d-%1d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
                lvl,pass,nblk,nrem1,nrem2,nrem1/nblk,nrem2/nblk,toc(time))
      end % if
    end % for
  end % for

  % finish
  F.factors = F.factors(1:n);
  if opts.verb
    fprintf(['-'*ones(1,80) '\n'])
    toc(start)
  end % if
  
  % sparse matrix access function (native MATLAB is slow for large matrices)
  function A = spget(Ityp,Jtyp)
    if strcmpi(Ityp,'slf')
      I_ = slf;
      m_ = nslf;
    elseif strcmpi(Ityp,'nbr')
      I_ = nbr;
      m_ = nnbr;
    elseif strcmpi(Ityp,'lst')
      I_ = lst;
      m_ = nlst;
    end 
    if strcmpi(Jtyp,'slf')
      J_ = slf;
      n_ = nslf;
    elseif strcmpi(Jtyp,'nbr')
      J_ = nbr;
      n_ = nnbr;
    elseif strcmpi(Jtyp,'lst')
      J_ = lst;
      n_ = nlst;
    end
   
    A = zeros(m_,n_);
    update_list = false(nbox,1);
    get_update_list(i);
    update_list = lookup_list(flip(find(update_list)'),:);
    update_list = update_list(update_list ~=0);
    update_list = update_list(:);

    for jj = update_list'
        g = F.factors(jj);
        xj = [g.sk, g.nbr];
        f = length(g.sk);

        if strcmpi(Ityp,Jtyp)
            idxI = ismembc2(xj,I_);

            tmp1 = idxI~=0;

            subI = idxI(tmp1);
            idxI1 = tmp1(1:f);
            idxI2 = tmp1(f+1:end);

            tmp1 = [g.E(idxI1,:); ...
                   g.C(idxI2,:)];

            A(subI, subI) = A(subI,subI) - tmp1*tmp1';
        else
            idxI = ismembc2(xj,I_);
            idxJ = ismembc2(xj,J_);

            tmp1 = idxI~=0;
            tmp2 = idxJ~=0;

            subI = idxI(tmp1);
            subJ = idxJ(tmp2);
            idxI1 = tmp1(1:f);
            idxI2 = tmp1(f+1:end);
            idxJ1 = tmp2(1:f);
            idxJ2 = tmp2(f+1:end);

            tmp1 = [g.E(idxI1,:); ...
                   g.C(idxI2,:)];
            tmp2 = [g.E(idxJ1,:); ...
                   g.C(idxJ2,:)]';
            A(subI, subJ) = A(subI,subJ) - tmp1*tmp2;
        end
                
    end

  

    function get_update_list(node_idx)
        update_list(node_idx) = 1;
        update_list(t.nodes(node_idx).snbor) = 1;
        for k = t.nodes(node_idx).chld
            get_update_list(k);
        end
    end
  end
end