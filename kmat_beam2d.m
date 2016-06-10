function [ke, rq] = kmat_beam2d(coord, EI, AE, q)
% [ke, rq] = kmat_beam2d(coord, EI, AE, q)
%
%   Generates equations for a beam element, ke, and the 
%   EI = beam stiffness
%   AE = axial stiffness
%   q = distributed load
%   coord = coordinates at the element ends
%
%
% ke = kmat_beam2d(coord, EI, AE)
%
%   Generates equations for a beam element, ke, and the 
%   EI = beam stiffness
%   AE = axial stiffness
%   coord = coordinates at the element ends


L=norm(coord(2,:)-coord(1,:));
ke = [ AE/L, 0,    0,      -AE/L,   0, 0; 
      0, (12*EI)/L^3, (6*EI)/L^2, 0, -((12*EI)/L^3),  (6*EI)/L^2;       
      0,(6*EI)/L^2,   (4*EI)/L, 0,  -((6*EI)/L^2),  (2*EI)/L;
	      -AE/L, 0,    0,      AE/L,   0, 0;
      0, -((12*EI)/L^3),  -((6*EI)/L^2), 0, (12*EI)/L^3, -((6*EI)/L^2);
      0, (6*EI)/L^2,  (2*EI)/L, 0, -((6*EI)/L^2),  (4*EI)/L];

c=(coord(2,1)-coord(1,1))/L;  % element cosine
s=(coord(2,2)-coord(1,2))/L;  % element sine

T=[ c -s 0 0 0 0;
   s c 0 0 0 0;
    0 0 1 0 0 0;
	0 0 0 c -s 0;
    0 0 0 s c 0;
    0 0 0 0 0 1 ];

ke=T'*ke*T;

if ( nargin>3 )
  rq = [(L*q)/2; 0; (L^2*q)/12; (L*q)/2; 0; -((L^2*q)/12)];
	rq=T*rq;
else
  rq = zeros(6,1);
end
