function [B,jac]=bmat_hexa8(coord,xi)

% function B=bmat_hexa8(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for an eight node
% hexahedral element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,jac]=bmat_hexa8(coord,xi)
%
% Computes the B matrix and the element Jacobian at xi
%
% Written by Jack Chessa, jfchessa@utep.edu

eta=xi(2);
zeta=xi(3);
xi=xi(1);

dNxi=0.125*[-1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
  1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
  1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
  -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;
  -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
  1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
  1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
  -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ];

J=coord'*dNxi;
jac=det(J);
dN=dNxi/J;

B=zeros(6,24);
for i=1:8
  B(1,3*i-2)=dN(i,1);
  B(2,3*i-1)=dN(i,2);
  B(3,3*i)  =dN(i,3);
  B(4,3*i-1)=dN(i,3); B(4,3*i)=dN(i,2);
  B(5,3*i-2)=dN(i,3); B(5,3*i)=dN(i,1);
  B(6,3*i-2)=dN(i,2); B(6,3*i-1)=dN(i,1);
end