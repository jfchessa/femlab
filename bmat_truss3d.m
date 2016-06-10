
function [B,L]=bmat_truss3d(coord,xi)

% function B=bmat_truss3d(coord)
%
% Computes the strain-displacement matrix (B matrix) for a 3D truss
% element.
%
%    coord: the nodal coordinates of the element 
%
% function [B,L]=bmat_truss3d(coord)
%
% Computes the B matrix and the element length
%
% Written by Jack Chessa, jfchessa@utep.edu

x1=coord(1,1); y1=coord(1,2); z1=coord(1,3);
x2=coord(2,1); y2=coord(2,2); z2=coord(2,3);

L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
invL=1/L;
cx=(x2-x1)*invL;
cy=(y2-y1)*invL;
cz=(z2-z1)*invL;
B=invL*[-cx -cy -cz cx cy cz];
