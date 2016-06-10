function [B,A]=bmat_tetra10(coord,xi)

% function B=bmat_tetra10(coord)
%
% Computes the strain-displacement matrix (B matrix) for a ten node
% tetrahedral element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,jac]=bmat_tetra10(coord)
%
% Computes the B matrix and the element Jacobian at xi
%
% Written by Jack Chessa, jfchessa@utep.edu

l2 = xi(1)
l3 = xi(2)
l4 = xi(3)
l1=1-l2-l3-l4

dNdxi = [ -(4*l1-1) -(4*l1-1) -(4*l1-1);
  4*l2-1 0.0 0.0 ;
  0.0 4*l3-1 0.0 ;
  0.0 0.0 4*l4-1 ;
  4*(l1-l2) -4*l2 -4*l2;
  4*l3 4*l2 0.0;
 -4*l3 4*(l1-l3) -4*l3;
 -4*l4 -4*l4 4*(l1-l4);
  4*l4 0.0 4*l2;
  0.0 4*l4 4*l3 ];

J=coord'*dNxi;
jac=det(J);
dN=dNxi*inv(J);

B=zeros(6,30);
for i=1:10
  B(1,3*i-2)=dN(i,1);
  B(2,3*i-1)=dN(i,2);
  B(3,3*i)  =dN(i,3);
  B(4,3*i-1)=dN(i,3); B(4,3*i)=dN(i,2);
  B(5,3*i-2)=dN(i,3); B(5,3*i)=dN(i,1);
  B(6,3*i-2)=dN(i,2); B(6,3*i-1)=dN(i,1);
end