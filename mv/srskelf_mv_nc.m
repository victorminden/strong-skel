function Y = srskelf_mv_nc(F,X)
% Y = SRSKELF_MV_NC(F,X) The strong skeletonization factorization stored in
%  F is used to compute Y=FX, where the "_NC" suffix indicates that the
%  factorization in F is not symmetric and we are taking the conjugate
%  transpose.  The input X is a vector or matrix of appropriate size.

  % Initialize
  n = F.lvp(end);
  Y = X;
  % Upward sweep, applying the factors for each box from bottom to top
  for i = 1:n
    sk  = F.factors(i).sk;
    rd  = F.factors(i).rd;
    nbr = F.factors(i).nbr;

    Y(sk,:) = Y(sk,:) + F.factors(i).T*Y(rd,:);
    Y(rd,:) = F.factors(i).L'*Y(rd,:);
    Y(rd,:) = Y(rd,:) + F.factors(i).E'*Y(sk,:);
    Y(rd,:) = Y(rd,:) + F.factors(i).C'*Y(nbr,:);
  end

  % Downward sweep, applying the factors for each box from top to bottom
  for i = n:-1:1
    sk = F.factors(i).sk;
    rd = F.factors(i).rd;
    nbr = F.factors(i).nbr;
    
    Y(sk,:)  = Y(sk,:)  + F.factors(i).F'*Y(rd,:);
    Y(nbr,:) = Y(nbr,:) + F.factors(i).D'*Y(rd,:);
    Y(rd,:) = F.factors(i).U'*Y(rd,:);
    Y(rd,:) = Y(rd,:) + F.factors(i).T'*Y(sk,:);
  end
end