function [Nv,Nt,dNtdx]=cubic_hermite(xi,L)

% function [Nv,Nt,dNtdx]=cubic_hermit(xi,L)
%
% Returns the 1D cubic Hermite shape functions (C1 continuity, e.g. beams)
% at the interpolation point xi, where 0<xi<1.  The scale length of the 
% element, L, needs to be given as well.
%
%  Nv are the shape functions for the primary interpolation
%  Nt are the shape fuctions for the derivative of the primary interp
%  (i.e. the beam rotation)
%
%  returns the derivative of the of teh Nt shape functions with respect to 
%  x.  This is typically related to the B matrix.  
%  Where x is the scaled length - i.e. dx = L dxi
%
%  Written by Jack Chessa
%  jfchessa@utep.edu
%

s=xi*L;

Nv=[ 1+2*s.^3/L^3-3*s.^2/L^2;  s+s.^3/L^2-2*s.^2/L; -2*s.^3/L^3+3*s.^2/L^2;  s.^3/L^2-s.^2/L;   ];
Nt = [ (6*s.^2)/L^3-(6*s)/L^2;	(3*s.^2)/L^2-(4*s)/L+1;	(6*s)/L^2-(6*s.^2)/L^3;	(3*s.^2)/L^2-(2*s)/L ];

dNtdx = [ -((6*(L-2*s))/L^3); (6*s-4*L)/L^2; (6*(L-2*s))/L^3; -((2*(L-3*s))/L^2) ];
