function F = srskelf_asym(A,x,occ,rank_or_tol,pxyfun,opts)
% SRSKELF_ASYM   Asymmetric strong recursive skeletonization factorization.
%
%    F = SRSKELF_ASYM(A,X,OCC,RANK_OR_TOL,PXYFUN) produces a factorization 
%    F of the interaction matrix A on the points X using tree occupancy 
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
%    F = SRSKELF_ASYM(A,X,OCC,RANK_OR_TOL,PXYFUN,OPTS) also passes various 
%    options to the algorithm. Valid options include:
%
%      - EXT: set the root node extent to [EXT(I,1) EXT(I,2)] along 
%             dimension I.  If EXT is empty (default), then the root extent
%             is calculated from the data.
%
%      - LVLMAX: maximum tree depth (default: LVLMAX = Inf).
%
%      - SYMM: assume that the matrix is asymmetric if SYMM = 'N' and 
%              Hermitian positive definite if SYMM = 'P' (default: SYMM = 
%              'N'). If SYMM = 'N', then local factors are computed using 
%              the LU decomposition; if SYMM = 'P', the Cholesky 
%              decomposition.
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
  if ~isfield(opts,'ext')
    opts.ext = [];
  end % if
  if ~isfield(opts,'lvlmax')
    opts.lvlmax = Inf;
  end % if
  if ~isfield(opts,'symm')
    opts.symm = 'n';
  end % if
  if ~isfield(opts,'verb')
    opts.verb = 0;
  end % if
  
  if opts.verb
    disp('This is standard asymmetric srskelf (RS-S).');
  end

  % Check inputs are sensible
    assert(strcmpi(opts.symm,'p') || strcmpi(opts.symm,'n'), ...
         'RSS:srskelf_asym:invalidSymm', ...
         'Symmetry parameter must be ''p'' or ''n''.');

  % Build tree to hold the discretization points
  N = size(x,2);
  tic
  t = shypoct(x,occ,opts.lvlmax,opts.ext);

  if opts.verb
    fprintf(['-'*ones(1,80) '\n'])
    fprintf('%3s | %6s | %8s | %8s | %8s | %8s | %10s (s)\n', ...
            'lvl','nblk','nRemIn','nRemOut','inRatio','outRatio','time')
    % Print summary information about tree construction
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
  % Each element of F.factors will contain the following data for one box:
  %   - sk: the skeleton DOF indices
  %   - rd: the redundant DOF indices
  %   - nbr: the neighbor (near-field) DOF indices
  %   - T: the interpolation matrix mapping redundant to skeleton
  %   - E: the left factor of the Schur complement update to sk
  %   - F: the right factor of the Schur complement update to sk
  %   - L: the left factor of the diagonal block
  %   - U: the right factor of the diagonal block
  %   - C: the left factor of the Schur complement update to nbr
  %   - D: the right factor of the Schur complement update to nbr
  F = struct('sk',e,'rd',e,'nbr',e,'T',e,'E',e,'F',e,'L',e,'U',e,'C',e,...
             'D',e);
  F = struct('N',N,'nlvl',t.nlvl,'lvp',zeros(1,t.nlvl+1),'factors',F,...
             'symm',opts.symm);
  nlvl = 0;
  n = 0;
  % Mark every DOF as "remaining", i.e., not yet eliminated
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

      % If we are at the second level (i.e., the first level we reach in a 
      % bottom-to-top loop in which there do not exist pairs of 
      % non-adjacent boxes) then we can do weak skeletonization, so instead 
      % of the interaction list we skeletonize against the neighbor set.
      if lvl == 2
        lst = nbr;
        l = t.lrt/2^(lvl - 1);
      else
        lst = [t.nodes(t.nodes(i).ilist).xi];
        l = t.lrt/2^(lvl - 1) * 5/3;
      end % if

      % Compute proxy interactions and subselect neighbors
      Kpxy = zeros(0,nslf);
      if lvl > 2
        [Kpxy,lst] = pxyfun(x,slf,lst,l,t.nodes(i).ctr);
      end % if

      nlst = length(lst);
      % Sorting not necessary, but makes debugging easier
      lst = sort(lst);
      
      % Compute interaction matrix between box and far-field (except level
      % 2, where near-field is included).
      K1 = full(A(lst,slf));
      if strcmpi(opts.symm,'n')
        K1 = [K1; full(A(slf,lst))'];
      end % if

      K2 = spget('lst','slf');
      if strcmpi(opts.symm,'n')
          K2 = [K2; spget('slf','lst')'];
      end % if

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
      K2 = full(A(nbr,slf)) + spget('nbr','slf');
      if strcmpi(opts.symm,'n')
        K3 = full(A(slf,nbr)) + spget('slf','nbr');
      end % if
      
      % Skeletonize
      K(rd,:) =  K(rd,:) - T'*K(sk,:);
      K(:,rd) = K(:,rd) - K(:,sk)*T;
      K2(:,rd) = K2(:,rd) - K2(:,sk)*T; 
      if strcmpi(opts.symm,'n')
        K3(rd,:) = K3(rd,:) - T'*K3(sk,:); 
      end % if
      
      if strcmpi(opts.symm,'p')
        % Cholesky for positive definite input
        L = chol(K(rd,rd),'lower');
        U = [];
        E = K(sk,rd)/L';
        G = [];
        C = K2(:,rd)/L';
        D = [];
      elseif strcmpi(opts.symm,'n')
        % Otherwise, LU
        [L,U] = lu(K(rd,rd));
        E = K(sk,rd)/U;
        G = L\K(rd,sk);
        C = K2(:,rd)/U;
        D = L\K3(rd,:);
      end % if
 
      % Store matrix factors for this box
      n = n + 1;
      F.factors(n).sk  = slf(sk);
      F.factors(n).rd  = slf(rd);
      F.factors(n).nbr = nbr;
      F.factors(n).T = T;
      F.factors(n).E = E;
      F.factors(n).F = G;
      F.factors(n).L = L;
      F.factors(n).U = U;
      F.factors(n).C = C;
      F.factors(n).D = D;
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

   % Truncate extra storage, and we are done
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
    
    A = zeros(m_,n_);
    update_list = false(nbox,1);
    get_update_list(i);
    update_list = lookup_list(flip(find(update_list)'));
    update_list = update_list(update_list~=0)';
    for jj = update_list
      g = F.factors(jj);
      xj = [g.sk, g.nbr];
      f = length(g.sk);
            
      if strcmpi(Ityp,Jtyp)
        % For diagonal block
        idxI = ismembc2(xj,I_);
        tmp1 = idxI~=0;
        subI = idxI(tmp1);
        idxI1 = tmp1(1:f);
        idxI2 = tmp1(f+1:end);
        tmp1 = [g.E(idxI1,:); g.C(idxI2,:)];
        % Different factorization depending on symmetry
        if strcmpi(opts.symm,'p')
          A(subI, subI) = A(subI,subI) - tmp1*tmp1';
        elseif strcmpi(opts.symm,'n')
          tmp2 = [g.F(:,idxI1), g.D(:,idxI2)];
          A(subI, subI) = A(subI,subI) - tmp1*tmp2;
        end % if
      else
        % For off-diagonal block
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
        % Different factorization depending on symmetry
        if strcmpi(opts.symm,'p')
          tmp2 = [g.E(idxJ1,:); g.C(idxJ2,:)]';
        elseif strcmpi(opts.symm,'n')
          tmp2 = [g.F(:,idxJ1), g.D(:,idxJ2)];
        end % if
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
end % srskelf_asym