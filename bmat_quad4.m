function [B,jac]=bmat_quad4(coord,xi)

% function B=bmat_quad4(coord,xi)
%
% Computes the strain-displacement matrix (B matrix) for a four node
% quadrilateral element.
%
%    coord: the nodal coordinates of the element (4x2 matrix)
%
% function [B,jac]=bmat_quad4(coord,xi)
%
% Computes the B matrix and the Jacobian
%
% Written by Jack Chessa, jfchessa@utep.edu

eta=xi(2);
xi=xi(1);

dNxi=0.25*[ -(1-eta) -(1-xi);
             (1-eta) -(1+xi);
			 (1+eta)  (1+xi);
			-(1+eta)  (1-xi) ];

J=coord'*dNxi;
jac=det(J);
dN=dNxi*inv(J);

B = [ dN(1,1) 0 dN(2,1) 0 dN(3,1) 0 dN(4,1) 0;
	  0 dN(1,2) 0 dN(2,2) 0 dN(3,2) 0 dN(4,2);
	dN(1,2) dN(1,1) dN(2,2) dN(2,1) dN(3,2) dN(3,1) dN(4,2) dN(4,1) ];