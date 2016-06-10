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

eta=xi(2);
zeta=xi(3);
xi=xi(1);

dNdxi=zeros(20,3);
dNdxi(1,1)=-0.125*(1-eta)*(1-zeta)*(-zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(1-zeta)
dNdxi(1,2)=-0.125*(1-xi)*(1-zeta)*(-zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(1-zeta)
dNdxi(1,3)=-0.125*(1-eta)*(1-xi)*(-zeta-xi-eta-2) -0.125*(1-eta)*(1-xi)*(1-zeta)
dNdxi(2,1)=-0.25*(1-eta)*(1-zeta^2)
dNdxi(2,2)=-0.25*(1-xi)*(1-zeta^2)
dNdxi(2,3)=-0.5*(1-eta)*(1-xi)*zeta
dNdxi(3,1)=-0.125*(1-eta)*(zeta+1)*(zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(zeta+1)
dNdxi(3,2)=-0.125*(1-xi)*(zeta+1)*(zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(zeta+1)
dNdxi(3,3)= 0.125*(1-eta)*(1-xi)*(zeta-xi-eta-2)+0.125*(1-eta)*(1-xi)*(zeta+1)
dNdxi(4,1)=-0.5*(1-eta)*xi*(zeta+1)
dNdxi(4,2)=-0.25*(1-xi^2 )*(zeta+1)
dNdxi(4,3)= 0.25*(1-eta)*(1-xi^2)
dNdxi(5,1)= 0.125*(1-eta)*(zeta+1)*(zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(zeta+1)
dNdxi(5,2)= -0.125*(xi+1)*(zeta+1)*(zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(zeta+1)
dNdxi(5,3)= 0.125*(1-eta)*(xi+1)*(zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(zeta+1)
dNdxi(6,1)= 0.25*(1-eta)*(1-zeta^2)
dNdxi(6,2)=-0.25*(xi+1)*(1-zeta^2)
dNdxi(6,3)=-0.5*(1-eta)*(xi+1)*zeta
dNdxi(7,1)= 0.125*(1-eta)*(1-zeta)*(-zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(1-zeta)
dNdxi(7,2)=-0.125*(xi+1)*(1-zeta)*(-zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(1-zeta)
dNdxi(7,3)=-0.125*(1-eta)*(xi+1)*(-zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(1-zeta)
dNdxi(8,1)=-0.5*(1-eta)*xi*(1-zeta)
dNdxi(8,2)=-0.25*(1-xi^2)*(1-zeta)
dNdxi(8,3)=-0.25*(1-eta)*(1-xi^2)
dNdxi(9,1)=-0.25*(1-eta )*(1-zeta)
dNdxi(9,2)=-0.5*eta*(1-xi)*(1-zeta)
dNdxi(9,3)=-0.25*(1-eta^2 )*(1-xi)
dNdxi(10,1)=-0.25*(1-eta^2)*(zeta+1)
dNdxi(10,2)=-0.5*eta*(1-xi)*(zeta+1)
dNdxi(10,3)= 0.25*(1-eta^2)*(1-xi)
dNdxi(11,1)= 0.25*(1-eta^2)*(zeta+1)
dNdxi(11,2)=-0.5*eta*(xi+1)*(zeta+1)
dNdxi(11,3)= 0.25*(1-eta^2)*(xi+1)
dNdxi(12,1)= 0.25*(1-eta^2)*(1-zeta)
dNdxi(12,2)=-0.5*eta*(xi+1)*(1-zeta)
dNdxi(12,3)=-0.25*(1-eta^2)*(xi+1)
dNdxi(13,1)=-0.125*(eta+1)*(1-zeta)*(-zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(1-zeta)
dNdxi(13,2)= 0.125*(1-xi)*(1-zeta)*(-zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(1-zeta)
dNdxi(13,3)=-0.125*(eta+1)*(1-xi)*(-zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(1-zeta)
dNdxi(14,1)=-0.25*(eta+1)*(1-zeta^2)
dNdxi(14,2)= 0.25*(1-xi)*(1-zeta^2)
dNdxi(14,3)=-0.5*(eta+1)*(1-xi)*zeta
dNdxi(15,1)=-0.125*(eta+1)*(zeta+1)*(zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(zeta+1)
dNdxi(15,2)= 0.125*(1-xi)*(zeta+1)*(zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(zeta+1)
dNdxi(15,3)= 0.125*(eta+1)*(1-xi)*(zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(zeta+1)
dNdxi(16,1)=-0.5*(eta+1)*xi*(zeta+1)
dNdxi(16,2)= 0.25*(1-xi^2)*(zeta+1)
dNdxi(16,3)= 0.25*(eta+1)*(1-xi^2)
dNdxi(17,1)= 0.125*(eta+1)*(zeta+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1)
dNdxi(17,2)= 0.125*(xi+1)*(zeta+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1)
dNdxi(17,3)= 0.125*(eta+1)*(xi+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1)
dNdxi(18,1)= 0.25*(eta+1)*(1-zeta)
dNdxi(18,2)= 0.25*(xi+1)*(1-zeta^2)
dNdxi(18,3)=-0.5*(eta+1)*(xi+1)*zeta
dNdxi(19,1)= 0.125*(eta+1)*(1-zeta)*(-zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(1-zeta)
dNdxi(19,2)= 0.125*(xi+1)*(1-zeta)*(-zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(1-zeta)
dNdxi(19,3)=-0.125*(eta+1)*(xi+1)*(-zeta+xi+eta-2)-0.125*(eta+1)*(xi+1)*(1-zeta)
dNdxi(20,1)=-0.5*(eta+1)*xi*(1-zeta)
dNdxi(20,2)= 0.25*(1-xi^2)*(1-zeta)
dNdxi(20,3)=-0.25*(eta+1)*(1-xi^2)

J=coord'*dNxi;
jac=det(J);
dN=dNdxi*inv(J);

B=zeros(6,60);
for i=1:20
  B(1,3*i-2)=dN(i,1);
  B(2,3*i-1)=dN(i,2);
  B(3,3*i)  =dN(i,3);
  B(4,3*i-1)=dN(i,3); B(4,3*i)=dN(i,2);
  B(5,3*i-2)=dN(i,3); B(5,3*i)=dN(i,1);
  B(6,3*i-2)=dN(i,2); B(6,3*i-1)=dN(i,1);
end