function Bmat=bmat_beam1d(L,xi)
%
% function Bmat=bmat_beam1d(L,xi)
%
% Computes the B matrix for a 1d beam element at the point 0<xi<1
%
%
%           1--------------------2
%        xi = 0              xi = 1


s=xi*L;
Bmat = 1/L^3*[ 6*(2*s-L) L*(6*s-4*L) 6*(L-2*s) 2*L*(3*s-L) ];
           
       