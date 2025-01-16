function [B,jac]=bmat_quad8(coord,xi)

% function B=bmat_quad8(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for a eight node
% quadrilateral element.
%
%    coord: the nodal coordinates of the element (4x2 matrix)
%
% function [B,jac]=bmat_quad4(coord,xi)
%
% Computes the B matrix and the Jacobian
%
% Written by Jack Chessa, jfchessa@utep.edu

dNxi =dshape_quad8(xi);

J=coord'*dNxi;
jac=det(J);
dN=dNxi/J;

B=zeros(3,16);
for i=1:8
  B(1,2*i-1) = dN(i,1);
  B(2,2*i)   = dN(i,2);
  B(3,2*i-1) = dN(i,2);
  B(3,2*i)   = dN(i,1);
end