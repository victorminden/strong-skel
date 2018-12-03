function Y = srskelf_sv_p(F,X)
% Y = SRSKELF_SV_P(F,X) The strong skeletonization factorization stored in
%  F is used to solve FY=X, where the "_P" suffix indicates that the
%  factorization in F exploits the fact that F is positive-definite.  The 
%  input X is a vector or matrix of appropriate size.

  % Initialize
  n = F.lvp(end);
  Y = X;
  % Upward sweep, applying the factors for each box from bottom to top
  for i = 1:n
    sk  = F.factors(i).sk;
    rd  = F.factors(i).rd;
    nbr = F.factors(i).nbr;
    
    Y(rd,:) = Y(rd,:) - F.factors(i).T'*Y(sk,:);
    Y(rd,:) = F.factors(i).L\Y(rd,:);
    Y(sk,:) = Y(sk,:) - F.factors(i).E*Y(rd,:);
    Y(nbr,:) = Y(nbr,:) - F.factors(i).C*Y(rd,:);
  end

  % Downward sweep, applying the factors for each box from top to bottom.
  for i = n:-1:1
    sk  = F.factors(i).sk;
    rd  = F.factors(i).rd;
    nbr = F.factors(i).nbr;

    Y(rd,:) = Y(rd,:) - F.factors(i).E'*Y(sk,:);
    Y(rd,:) = Y(rd,:) - F.factors(i).C'*Y(nbr,:);
    Y(rd,:) = F.factors(i).L'\Y(rd,:);
    Y(sk,:) = Y(sk,:) - F.factors(i).T*Y(rd,:);
  end
end