% SRSKELF_SV_N  Dispatch for RSKELF_SV with F.SYMM = 'N'.
%
%    See also RSKELF, RSKELF_SV.

function Y = srskelf_sv_nc(F,X)

  % initialize
  n = F.lvp(end);
  Y = X;
  % upward sweep
  for i = 1:n
    sk  = F.factors(i).sk;
    rd  = F.factors(i).rd;
    nbr = F.factors(i).nbr;
    Y(rd,:) = Y(rd,:) - F.factors(i).T'*Y(sk,:);
    Y(rd,:) = F.factors(i).U'\Y(rd,:);
    Y(sk,:) = Y(sk,:) - F.factors(i).F'*Y(rd,:);
    Y(nbr,:) = Y(nbr,:) - F.factors(i).D'*Y(rd,:);
  end

  % downward sweep
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