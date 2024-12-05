
function [B,L]=bmat_beam3d(coord,xi)

% function B=bmat_beam3d(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for a 3D beam
% element.  Xi is the parent coordiante xi=[-1,1]
%
%    coord: the nodal coordinates of the element 
%
% function [B,L]=bmat_beam3d(coord,xi)
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
% ********* NOTE THIS DOES NOT COMPUTE AXIAL OR TORSION EFFECTS *********
%       *** add the strain component from a 3d truss bmat  ***
%
%     strain normal = ( yB + Btruss3d  )*d 
%     strain torsion = r Btorsion3d*d
%
%     stress normal = E * strain normal
%     stress shear = G * strain shear
% 
% The principal and Mises stresses are then
%
%     s1,s2 = stressnormal/2 +/- sqrt( stressnormal^2 + 4*stressshear^2 )
%     svm = 
%
% Written by Jack Chessa, jfchessa@utep.edu

x1=coord(1,1); y1=coord(1,2); z1=coord(1,3);
x2=coord(2,1); y2=coord(2,2); z2=coord(2,3);

L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
invL=1/L;
nx=(x2-x1)*invL;
ny=(y2-y1)*invL;
nz=(z2-z1)*invL;

% rev I think there was an extra 1/L here, JFC 12/5/24
% B=invL*[ a*nx, a*ny, a*nz, b*nx, b*ny, b*nz,  ..
a=6*xi*invL^2; b=(3*xi-1)*invL
B=[ a*nx, a*ny, a*nz, b*nx, b*ny, b*nz,  ...
    -a*nx, -a*ny, -a*nz, -b*nx, -b*ny, -b*nz,];
