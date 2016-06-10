
function [B,L]=bmat_truss2d(coord,xi)

% function B=bmat_truss2d(coord)
%
% Computes the strain-displacement matrix (B matrix) for a 2D truss
% element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,L]=bmat_tria3(coord)
%
% Computes the B matrix and the element length
%
% Written by Jack Chessa, jfchessa@utep.edu

x1=coord(1,1); y1=coord(1,2);
x2=coord(2,1); y2=coord(2,2);

L=sqrt((x2-x1)^2+(y2-y1)^2);
invL=1/L;
c=(x2-x1)*invL;
s=(y2-y1)*invL;
B=invL*[-c -s c s];
