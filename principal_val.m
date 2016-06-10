function [pval,pdir] = principal_val(s)

% function [pval,pdir] = principal_val(s)
%
% Calculates the principal values pval and the principal directions pdir
%

if ( length(s)==3 ) % 2D state
  A=s([1 3;3 2]);
else
  A=s([1 5 4;6 2 4;5 4 3]);
end

[pdir,d]=eig(A);
n=size(d,1);
pval=zeros(3,1);
pval(1:n)=d*ones(n,1);
