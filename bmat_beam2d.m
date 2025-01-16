
function [B,L]=bmat_beam2d(coord,xi)

% function B=bmat_beam2d(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for a 2D beam
% element.  Xi is the parent coordiante xi=[-1,1]
%
%    coord: the nodal coordinates of the element 
%
% function [B,L]=bmat_beam2d(coord,xi)
%
% Computes the B matrix and the element length
%
% So for these elements   
%        w" = B*d
%        M(xi) = EI w" = EI B*d
%        stress bending = M y/I = E y B*d
%        strain bending = y B*d
%
%
% ********* NOTE THIS DOES NOT COMPUTE AXIAL DEFORMATION *********
% *** to do this add the strain component from a 2d truss bmat ***
%
%     strain = ( yB + Btruss2d )*d 
%
% Written by Jack Chessa, jfchessa@utep.edu

x1=coord(1,1); y1=coord(1,2);  
x2=coord(2,1); y2=coord(2,2);  

L=sqrt((x2-x1)^2+(y2-y1)^2);
invL=1/L;
c=(x2-x1)*invL;
s=(y2-y1)*invL;

a=6*xi*invL^2; 

% rev I think there was an extra 1/L here, JFC 12/5/24 should check this
B=[ -a*s, -a*c, (3*xi-1)*invL, a*s, a*c, (1+3*xi)*invL ];
