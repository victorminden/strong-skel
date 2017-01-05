% RSKELF_MV_N  Dispatch for RSKELF_MV with F.SYMM = 'N'.
%
%    See also RSKELF, RSKELF_MV.

function Y = srskelf_mv_nc(F,X)

  % initialize
  n = F.lvp(end);
  Y = X;
  % upward sweep
  for i = 1:n
    sk  = F.factors(i).sk;
    rd  = F.factors(i).rd;
    nbr = F.factors(i).nbr;
    % came from (rows) rd <- rd - T^T*sk
    Y(sk,:) = Y(sk,:) + F.factors(i).T*Y(rd,:);
    Y(rd,:) = F.factors(i).L'*Y(rd,:);
    % came from (rows) sk <- sk-E*
    Y(rd,:) = Y(rd,:) + F.factors(i).E'*Y(sk,:);
    Y(rd,:) = Y(rd,:) + F.factors(i).C'*Y(nbr,:);
  end

  % downward sweep
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