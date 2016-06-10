function [B,A]=bmat_tria3(coord,xi)

% function B=bmat_tria3(coord)
%
% Computes the strain-displacement matrix (B matrix) for a three node
% triangular element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,A]=bmat_tria3(coord)
%
% Computes the B matrix and the element area
%
% Written by Jack Chessa, jfchessa@utep.edu

x1=coord(1,1); y1=coord(1,2);
x2=coord(2,1); y2=coord(2,2);
x3=coord(3,1); y3=coord(3,2);

b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2;
c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1;
f1 = x2*y3 - x3*y2; f2 = x3*y1 - x1*y3; f3 = x1*y2 - x2*y1;
jac = (f1 + f2 + f3);
A=jac/2;
B = [b1, 0, c1; 0, c1, b1; b2, 0, c2; 0, c2, b2; 
    b3, 0, c3; 0, c3, b3]'/jac;