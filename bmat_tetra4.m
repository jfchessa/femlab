function [B,A]=bmat_tetra4(coord,xi)

% function B=bmat_tetra4(coord)
%
% Computes the strain-displacement matrix (B matrix) for a four node
% tetrahedral element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,V]=bmat_tetra4(coord)
%
% Computes the B matrix and the element volume
%
% Written by Jack Chessa, jfchessa@utep.edu


x1=coord(1,1); y1=coord(1,2); z1=coord(1,3);
x2=coord(2,1); y2=coord(2,2); z2=coord(2,3);
x3=coord(3,1); y3=coord(3,2); z3=coord(3,3);
x4=coord(4,1); y4=coord(4,2); z4=coord(4,3);


jac=det([x2-x1 x3-x1 x4-x1;
         y2-y1 y3-y1 y4-y1;
         z2-z1 z3-z1 z4-z1 ]);

V=0.16666666666667*jac;

dN=inv(jac)*[-1 -1 -1;1 0 0;0 1 0; 0 0 1];      
       
B = [dN(1,1) 0 0 dN(2,1) 0 0 dN(3,1) 0 0 dN(4,1) 0 0 ;
     0 dN(1,2) 0 0 dN(2,2) 0 0 dN(3,2) 0 0 dN(4,2) 0;
     0 0 dN(1,3) 0 0 dN(2,3) 0 0 dN(3,3) 0 0 dN(4,3);
     0  dN(1,3) dN(1,2) 0  dN(2,3) dN(2,2) 0 dN(3,3) dN(3,2) 0 dN(4,3) dN(4,2);
     dN(1,3) 0 dN(1,1) dN(2,3) 0 dN(2,1) dN(3,3) 0 dN(3,1) dN(4,3) 0 dN(4,1);
     dN(1,2) dN(1,1) 0  dN(2,2) dN(2,1) 0 dN(3,2) dN(3,1) 0 dN(4,2) dN(4,1) 0 ];