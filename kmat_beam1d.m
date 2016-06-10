function [ke,re]=kmat_beam1d(gcoord,EI,q1,q2)
%
% function [ke,re]=kmat_beam1d(gcoord,IE,q1,q2)
%
% Computes the stifness, ke and teh reaction forces due to the  linearly
% variying distributed load defined by q1 and q2 (the values at the nodes)
%
% function ke=kmat_beam1d(gcoord,IE)
%
%

if ( nargin < 3 )   % set defalut values of q1 and q2 if not given
    q1=0;
end
if ( nargin < 4 )
    q2=q1;
end

L = norm( gcoord(1,:) - gcoord(2,:) );  % the length of the beam


ke = EI/L^3*[ 12 6*L -12 6*L;
               6*L 4*L^2 -6*L 2*L^2;
               -12 -6*L 12 -6*L;
               6*L 2*L^2 -6*L 4*L^2 ];
           
re = [ 7*L/20 3*L/20;
       L^2/20 L^2/30;
       3*L/20 7*L/20;
      -L^2/30 -L^2/20 ]*[q1;q2];
       