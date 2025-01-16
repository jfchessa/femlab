function [B,jac]=bmat_hexa20(coord,xi)

% function B=bmat_hexa20(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for an twenty node
% hexahedral element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,jac]=bmat_hexa20(coord,xi)
%
% Computes the B matrix and the element Jacobian at xi
%
% Written by Jack Chessa, jfchessa@utep.edu

dNdxi = dshape_hexa20(xi);
J=coord'*dNdxi;
jac=det(J);
dN=dNdxi/J;

B=zeros(6,60);
for i=1:20
  B(1,3*i-2)=dN(i,1);
  B(2,3*i-1)=dN(i,2);
  B(3,3*i)  =dN(i,3);
  B(4,3*i-1)=dN(i,3); B(4,3*i)=dN(i,2);
  B(5,3*i-2)=dN(i,3); B(5,3*i)=dN(i,1);
  B(6,3*i-2)=dN(i,2); B(6,3*i-1)=dN(i,1);
end
