function [B,jac]=bmat_tria6(coord,xi)

% function B=bmat_tria6(coord)
%
% Computes the strain-displacement matrix (B matrix) for a six node
% triangular element.
%
%    coord: the nodal coordinates of the element (3x2 matrix)
%
% function [B,jac]=bmat_tria3(coord)
%
% Computes the B matrix and the jacobian (1/2 element area)
%
% Written by Jack Chessa, jfchessa@utep.edu


eta=xi(2);
xi=xi(1);

dNxi=[4*(xi+eta)-3   4*(xi+eta)-3;
             4*xi-1              0; 
                  0        4*eta-1;
     4*(1-eta-2*xi)          -4*xi;
              4*eta           4*xi;
             -4*eta  4*(1-xi-2*eta)];

J=coord'*dNxi;
jac=det(J);
dN=dNxi*inv(J);

B=zeros(3,12);
for i=1:6
  B(1,2*i-1)=dN(i,1);
  B(2,2*i)=dN(i,2);
  B(3,2*i-1)=dN(i,2);
  B(3,2*i)=dN(i,1);
end