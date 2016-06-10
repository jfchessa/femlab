function d=diracr(eta,m,k)

% function  d=diracr(x,m,k)
%
% returns a regularized Dirac delta function of the variable
% x.  The true Dirac function is approximated by a function
% which is nonzero over the open set (-eps,eps), and the 
% integral of the approximation over that open set is unity.
%
% k is the continuity of the regularization and m is the number
% of vanishing moments.  The default is m=1 k=1. Other valid
% options are {m,k}={0,0},{1,0},{1,1},{1,2},{3,0},{3,1},{3,2}

d=0*eta;

if ( nargin == 1 )
  m=1;
  k=1;
end

is=find(abs(eta)<1);

switch 10*m+k
  
case 0
  d(is)=1/2;
  
case 11
  d(is)=15/16*(1-2*eta(is).^2+eta(is).^4);
  
case 12
  d(is)=35/32*(1-3*eta(is).^2+3*eta(is).^4-eta(is).^6);
  
case 30
  d(is)=15/32*(3-10*eta(is).^2+7*eta(is).^4);
  
case 31
  d(is)=105/64*(1-5*eta(is).^2+7*eta(is).^4-3*eta(is).^6);
  
case 32
  d(is)=315/512*(3-20*eta(is).^2+42*eta(is).^4-36*eta(is).^6+11*eta(is).^8);
  
otherwise
  d(is)=3/4*(1-eta(is).^2);
  %disp('invalid m and k for diracr')
  
end
  