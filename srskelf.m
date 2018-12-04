
function F = srskelf(A,x,occ,rank_or_tol,pxyfun,opts)
% SRSKELF   Strong recursive skeletonization factorization (symmetric 
%           positive definite).
%
%    F = SRSKELF(A,X,OCC,RANK_OR_TOL,PXYFUN) produces a factorization F of 
%    the interaction matrix A on the points X using tree occupancy 
%    parameter OCC, local precision parameter RANK_OR_TOL, and proxy 
%    function PXYFUN to capture the far field. This is a function of the 
%    form
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
%    F = SRSKELF(A,X,OCC,RANK_OR_TOL,PXYFUN,OPTS) also passes various 
%    options to the algorithm. Valid options include:
%
%      - EXT: set the root node extent to [EXT(I,1) EXT(I,2)] along 
%             dimension I. If EXT is empty (default), then the root extent 
%             is calculated from the data.
%
%      - LVLMAX: maximum tree depth (default: LVLMAX = Inf).
%
%      - VERB: display status of the code if VERB = 1 (default: VERB = 0).
%

  start = tic;
  % Set sane default parameters
  if nargin < 5
    pxyfun = [];
  end % if
  if nargin < 6
    opts = [];
  end % if
  if ~isfield(opts,'skip')
    opts.skip = 0;
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

  
  if opts.verb
      disp('This is standard symmetric positive definite srskelf (RS-S).');
      disp('Diagonal blocks will be factorized with Cholesky.');
  end % if

  % Build tree to hold the discretization points
  N = size(x,2);
  tic
  t = shypoct(x,occ,opts.lvlmax,opts.ext);

  % Print summary information
  % TODO(victorminden): include header row
  if opts.verb
    fprintf(['-'*ones(1,80) '\n'])
    fprintf('%3s | %63.2e (s)\n','-',toc)

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
  % TODO(victorminden): document the various things that show up here
  F = struct('sk',e,'rd',e,'nbr',e,'T',e,'E',e,'L',e,'C',e);
  F = struct('N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F, ...
             'symm','p','skip',opts.skip);
  nlvl = 0;
  n = 0;
  % Mark every DOF as "remaining", i.e., not yet eliminated.
  rem = true(N,1);
  lookup_list = zeros(nbox,1);
  
  % Loop over the levels of the tree from bottom to top
  for lvl = t.nlvl:-1:1
    time = tic;
    nlvl = nlvl + 1;
    nrem1 = sum(rem);

    % For each box, pull up information about skeletons from child boxes
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      t.nodes(i).xi = [t.nodes(i).xi [t.nodes(t.nodes(i).chld).xi]];
    end % for

    % Loop over each box in this level
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      slf = t.nodes(i).xi;
      nbr = [t.nodes(t.nodes(i).nbor).xi];

      nslf = length(slf);
      % Sorting not necessary, but makes debugging easier
      slf = sort(slf);
      
      nnbr = length(nbr);
      % Sorting not necessary, but makes debugging easier
      nbr = sort(nbr);

      % TODO(victorminden): what do we want to do with the skip parameter?
      % We can optionally choose to skip levels of skeletonization in which
      % case we can use some dummy values.  If we are at the second level
      % (i.e., the first level we reach in a bottom-to-top loop in which do
      % not exist pairs of non-adjacent boxes) then we can do weak
      % skeletonization, so instead of the interaction list we skeletonize 
      % against the neighbor set.
      if t.nlvl-lvl < opts.skip || lvl == 2
        lst = nbr;
        l = t.lrt/2^(lvl - 1);
      else
        lst = [t.nodes(t.nodes(i).ilist).xi];
        l = t.lrt/2^(lvl - 1) * 5/3;
      end % if

      % Compute proxy interactions and subselect neighbors
      Kpxy = zeros(0,nslf);
      if lvl > 2
        if isempty(pxyfun)
          assert(0);
          lst = setdiff(find(rem),slf);
        else
          [Kpxy,lst] = pxyfun(x,slf,lst,l,t.nodes(i).ctr);
        end % if
      end % if

      nlst = length(lst);
      lst = sort(lst);
      

      % Compute interaction matrix between box and far-field (except level
      % 2).
      K1 = full(A(lst,slf));
      K2 = spget('lst','slf');

      K = [K1 + K2; Kpxy];
      % Compute the skeleton/redundant points and interpolation matrix
      [sk,rd,T] = id(K,rank_or_tol);

      % Move on to next box if no compression for this box
      if isempty(rd)
        continue
      end % if

      % Otherwise, compute the diagonal and off-diagonal blocks for this 
      % box
      K  = full(A(slf,slf)) + spget('slf','slf');
      if t.nlvl-lvl >= opts.skip
        K2 = full(A(nbr,slf)) + spget('nbr','slf');
      end % if
      % Skeletonize
      K(rd,:) =  K(rd,:) - T'*K(sk,:);
      K(:,rd) = K(:,rd) - K(:,sk)*T;
      if t.nlvl-lvl >= opts.skip
        K2(:,rd) = K2(:,rd) - K2(:,sk)*T; 
      end % if
      % Cholesky factor of diagonal block
      L = chol(K(rd,rd),'lower');
      % Throw Cholesky onto intermediate factors
      E = K(sk,rd)/L';
      if t.nlvl-lvl >= opts.skip
        C = K2(:,rd)/L';
      else
          C = zeros(0,length(rd));
      end % if

 
      % Store matrix factors for this box
      n = n + 1;
      F.factors(n).sk  = slf(sk);
      F.factors(n).rd  = slf(rd);
      if t.nlvl-lvl >= opts.skip || lvl == 2
        F.factors(n).nbr = nbr;
      else
          F.factors(n).nbr = [];
      end % if
      F.factors(n).T = T;
      F.factors(n).E = E;
      F.factors(n).L = L;
      F.factors(n).C = C;
      % Box number i is at index n (more sensible for non-uniform case)
      lookup_list(i) = n;
    
      t.nodes(i).xi = slf(sk);
      rem(slf(rd)) = 0;

    end % for
    F.lvp(nlvl+1) = n;
    
    % Print summary for the latest level
    if opts.verb
      nrem2 = sum(rem);
      nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);
      fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
              lvl,nblk,nrem1,nrem2,nrem1/nblk,nrem2/nblk,toc(time))
    end % if
  end % for

  % Truncate extra storage, and we are done.
  F.factors = F.factors(1:n);
  if opts.verb
    fprintf(['-'*ones(1,80) '\n'])
    toc(start)
  end % if
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function A = spget(Ityp,Jtyp)
    % A = SPGET(ITYP,JTYP) Sparse matrix access function (native MATLAB is 
    % slow for large matrices).  We grab the accumulated Schur complement
    % updates to a block of the matrix from previously-skeletonized 
    % levels.  Index sets ITYP and JTYP can be 'slf', 'nbr', or 'lst'.
     
    % Translate input strings to index sets (and their lengths)
    if strcmpi(Ityp,'slf')
      I_ = slf;
      m_ = nslf;
    elseif strcmpi(Ityp,'nbr')
      I_ = nbr;
      m_ = nnbr;
    elseif strcmpi(Ityp,'lst')
      I_ = lst;
      m_ = nlst;
    end % if
    
    if strcmpi(Jtyp,'slf')
      J_ = slf;
      n_ = nslf;
    elseif strcmpi(Jtyp,'nbr')
      J_ = nbr;
      n_ = nnbr;
    elseif strcmpi(Jtyp,'lst')
      J_ = lst;
      n_ = nlst;
    end % if
    
    % Initialize an empty matrix to store updates
    A = zeros(m_,n_);
    % Find the updates, modifying update_list in the function call
    update_list = false(nbox,1);
    get_update_list(i);
    % Translate boxes (indexed relative to tree) to factors (indexed
    % relative to factorization)
    update_list = lookup_list(flip(find(update_list)'));
    % TODO(victorminden): why was this necessary?
    update_list = update_list(update_list~=0)';
        
    for jj = update_list
      g = F.factors(jj);
      xj = [g.sk, g.nbr];
      f = length(g.sk);
            
      if strcmpi(Ityp,Jtyp)
        % If this is a diagonal block, TODO(victorminden)
        idxI = ismembc2(xj,I_);
        tmp1 = idxI~=0;
        subI = idxI(tmp1);
        idxI1 = tmp1(1:f);
        idxI2 = tmp1(f+1:end);
        tmp1 = [g.E(idxI1,:); g.C(idxI2,:)];
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

        tmp1 = [g.E(idxI1,:); g.C(idxI2,:)];
        tmp2 = [g.E(idxJ1,:); g.C(idxJ2,:)]';
        A(subI, subJ) = A(subI,subJ) - tmp1*tmp2;
      end % if         
    end % for

    function get_update_list(node_idx)
      % GET_UPDATE_LIST(NODE_IDX) Recursively get the list of all nodes in
      % the tree that could have generated Schur complement updates to
      % points in node NODE_IDX
      update_list(node_idx) = 1;
      update_list(t.nodes(node_idx).snbor) = 1;
      for k = t.nodes(node_idx).chld
        get_update_list(k);
      end % for
    end % get_update_list
  end % spget
end % srskelf